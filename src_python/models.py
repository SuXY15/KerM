from utils import *

""" Mixing Model base class
"""
class MixingModel():
    def __init__(self, particles, weights):
        self.name = "MixingModel"
        self.weights = deepcopy(weights)
        self.phis = deepcopy(particles)
        self.var = var(self.phis, weights)
        self.dphidt = np.zeros_like(particles)
        self.N = len(particles)

    def update(self, Omega_phi, dt):
        self.phis += self.dphidt * dt
        self.var = var(self.phis, self.weights)
        return self.var

    def PDF(self):
        hist, bins = np.histogram(self.phis, weights=self.weights, bins=50, density=True)
        bin_centers = (bins[1:]+bins[:-1]) * 0.5
        return bin_centers, hist


"""  IEM model
""" 
class IEM(MixingModel):
    def __init__(self, particles, weights):
        super(IEM, self).__init__(particles, weights)
        self.name = "IEM"

    def update(self, Omega_phi, dt):
        E_phi = mean(self.phis, self.weights) # Expectations
        self.dphidt = - 1./2. * Omega_phi * (self.phis - E_phi)
        return super(IEM, self).update(Omega_phi, dt)


"""  Modified Curl model
""" 
class MCurl(MixingModel):
    def __init__(self, particles, weights):
        super(MCurl, self).__init__(particles, weights)
        self.name = "MCurl"
        self.N = len(particles)

    def update(self, Omega_phi, dt):
        n = int(3/2 * Omega_phi * self.N * dt + 1)
        p = 3/2 * Omega_phi * self.N * dt / n
        self.dphidt = np.zeros_like(self.phis)
        
        if n>self.N:
            print("Error: n(%d) > N(%d), please use smaller `dt` or `Omega_phi`"%(n,self.N))

        W = self.weights
        Wavr = np.mean(self.weights)
        p_sets = np.random.choice(self.N, n, p=(W+Wavr)/np.sum(W+Wavr))
        q_sets = np.array([np.random.choice(self.N, 1, p=(W[ip]+W)/np.sum(W[ip]+W))
                    for ip in p_sets]).flatten()
        idx = np.where(np.random.rand(n) < p, True, False)
        p_idxs = p_sets[idx]
        q_idxs = q_sets[idx]

        d_pq = self.phis[p_idxs] - self.phis[q_idxs]
        alpha = np.random.rand(*d_pq.shape)
        self.dphidt[p_idxs] = - alpha / 2 * d_pq / dt
        self.dphidt[q_idxs] =   alpha / 2 * d_pq / dt
        return super(MCurl, self).update(Omega_phi, dt)


""" Euclidean Minimum Spanning Tree model
""" 
class EMST(MixingModel):
    def __init__(self, particles, weights):
        super(EMST, self).__init__(particles, weights)
        self.name = "EMST"
        self.N = len(particles)

    def update(self, Omega_phi, dt):
        dt_in = deepcopy(dt)

        sorted_phis = sorted(self.phis.flatten())
        sorted_idxs = np.argsort(self.phis.flatten())

        # get weights
        w = self.weights / np.sum(self.weights[sorted_idxs])
        W = np.cumsum(w)[:-1]
        wv= np.array([min(Wi, 1-Wi) for Wi in W])
        B = (2*wv) # size(Np-1)

        # mixing on the edges
        dphi = np.zeros_like(self.phis)
        for v in range(self.N-1):
            mv = sorted_idxs[v]
            nv = sorted_idxs[v+1]
            dphi[mv] += - B[v] * (self.phis[mv] - self.phis[nv]) / w[mv]
            dphi[nv] += - B[v] * (self.phis[nv] - self.phis[mv]) / w[nv]

        AA = mean(dphi**2, self.weights)
        BB = 2*mean(dphi*self.phis, self.weights)
        CC = Omega_phi*self.var
        
        dt = 1.0 * BB**2/(4*AA*CC)
        dt = min(dt, dt_in)

        alpha = -BB/(2*AA*dt)
        # alpha = -CC/BB
        # alpha = (-BB+np.sqrt(abs(BB**2-4*AA*CC*dt)))/(2*AA*dt)

        # # root finding process, unnecessary for 1D
        # for i in range(2):
        #     new_var = var(self.phis+dphi*alpha*dt, self.weights)
        #     var_decay = 1-np.mean(new_var / self.var)
        #     var_ratio = var_decay / (1-np.exp(-Omega_phi*dt))
        #     alpha = alpha / var_ratio
        #     print("   ", var_decay, var_ratio)
            
        self.dphidt = dphi * alpha

        # print("EMST: dt_in=%6.1e, dt=%6.1e, alpha=%e"%(dt_in, dt, alpha))
        if dt < dt_in:
            super(EMST, self).update(Omega_phi, dt)
            return self.update(Omega_phi, dt_in-dt)
        else:
            # print("EMST: finsih inner loop.")
            return super(EMST, self).update(Omega_phi, dt)


""" Kernel Mixing model, with Localness in X space
"""
class KerM(MixingModel):
    def __init__(self, particles, weights, sigma_x=np.inf):
        super(KerM, self).__init__(particles, weights)
        self.name = "KerM"
        self.N = len(particles)
        self.sigma_x = sigma_x # mixing length scale for gaussian kernel

    def kernel_func(self, d):
        return np.exp(-d**2 / self.sigma_x**2 / 4).reshape(len(d),)

    def update(self, Omega_phi, dt):
        # # quick sort for CDF
        # pos = np.arange(self.N)
        # idx = np.argsort(self.phis)
        # x = np.zeros(self.N)
        # x[idx] = np.cumsum(self.weights[idx])

        # bucket sort for CDF
        Nbin = 50
        minz, maxz = np.min(self.phis) - 1e-8, np.max(self.phis) + 1e-8
        minB = np.arange(Nbin) / Nbin * (maxz - minz) + minz
        maxB = minB + 1 / Nbin * (maxz - minz)
        xbins = np.floor(((self.phis - minz) / (maxz - minz) * Nbin)).astype(int)
        countB = np.bincount(xbins)
        accumB = np.cumsum(countB) - countB
        x =  ((self.phis - minB[xbins]) / (maxB[xbins] - minB[xbins]) * countB[xbins] + accumB[xbins]) / self.N

        W = self.weights
        Wavr = np.mean(self.weights)

        # in python, the weighted sampling is inefficient
        p_idxs = np.random.choice(self.N, self.N, p=(W+Wavr)/np.sum(W+Wavr))
        q_idxs = np.array([np.random.choice(self.N, 1, p=(W[ip]+W)/np.sum(W[ip]+W))
                    for ip in p_idxs]).flatten()
        x_pq = distance(x[p_idxs] - x[q_idxs])
        f_pq = self.kernel_func(x_pq)
        d_pq = distance(self.phis[p_idxs] - self.phis[q_idxs])
        coeff = self.var / np.sum( 0.5 * f_pq * d_pq**2 / len(d_pq))

        # print("Coeff=%.3f"%(coeff))

        n = int(3/2 * Omega_phi * self.N * dt * coeff + 1)
        p = 3/2 * Omega_phi * self.N * dt * coeff / n
        
        if n>self.N:
            print("Error: n(%d) > N(%d), please use smaller `dt` or `Omega_phi`"%(n,self.N))

        self.dphidt = np.zeros_like(self.phis)
        p_sets = np.random.choice(self.N, n, p=(W+Wavr)/np.sum(W+Wavr))
        q_sets = np.array([np.random.choice(self.N, 1, p=(W[ip]+W)/np.sum(W[ip]+W))
                    for ip in p_sets]).flatten()
        x_pq = distance(x[p_sets] - x[q_sets])
        f = self.kernel_func(x_pq)
        f_idx = np.where(np.random.rand(*x_pq.shape) < f, True, False)
        p_idx = np.where(np.random.rand(*x_pq.shape) < p, True, False)
        p_idxs = p_sets[f_idx & p_idx]
        q_idxs = q_sets[f_idx & p_idx]

        d_pq = self.phis[p_idxs] - self.phis[q_idxs]
        x_pq = distance(x[p_idxs] - x[q_idxs])
        f_pq = self.kernel_func(x_pq)
        f_pq = f_pq.reshape(*d_pq.shape)
        alpha = np.random.rand(*d_pq.shape)
        self.dphidt[p_idxs] = - alpha / 2 * d_pq / dt
        self.dphidt[q_idxs] =   alpha / 2 * d_pq / dt
        return super(KerM, self).update(Omega_phi, dt)

