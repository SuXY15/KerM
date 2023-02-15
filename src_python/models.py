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
        Nbin = int(max(20, min(np.floor(np.sqrt(self.N)), 100)))
        hist, bins = np.histogram(self.phis, weights=self.weights, bins=Nbin, density=True)
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


""" Modified Curl model
    To acceralte the sampling process only uniform weights are supported
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
        
        # # in python, the weighted sampling is inefficient
        # p_sets = np.random.choice(self.N, n, p=(W+Wavr)/np.sum(W+Wavr))
        # q_sets = np.array([np.random.choice(self.N, 1, p=(W[ip]+W)/np.sum(W[ip]+W))
        #             for ip in p_sets]).flatten()
        # idx = np.where(np.random.rand(n) < p, True, False)
        # p_idxs = p_sets[idx]
        # q_idxs = q_sets[idx]

        idx = np.where(np.random.rand(n) < p, True, False)
        p_idxs = np.random.choice(self.N, n)[idx]
        q_idxs = np.random.choice(self.N, n)[idx]

        d_pq = self.phis[p_idxs] - self.phis[q_idxs]
        alpha = np.random.rand(*d_pq.shape)
        self.phis[p_idxs] -= alpha / 2 * d_pq 
        self.phis[q_idxs] += alpha / 2 * d_pq 
        self.var = var(self.phis, self.weights)
        return self.var


"""  Mapping Closure model
""" 
class MCMG(MixingModel):
    def __init__(self, particles, weights):
        super(MCMG, self).__init__(particles, weights)
        self.name = "MCMG"
        self.N = len(particles)

    def update(self, Omega_phi, dt):
        rtpi = 0.39894228
        ldab = 2
        N = self.N

        homdt = 0.5 * Omega_phi * dt
        wtsum = np.sum(self.weights)

        cdfs = np.zeros(N+1)
        bphs = np.zeros(N-1)
        Ab = np.zeros((ldab,N))

        favg = np.sum(self.weights * self.phis) / wtsum
        fvar = np.sum(self.weights * self.phis**2) / wtsum - favg**2

        phis = np.sort(self.phis)
        idxs = np.argsort(self.phis)
        cdfs[:N-1] = np.cumsum(self.weights / wtsum)[:N-1]
        cdfs[N-1] = 0.5 * cdfs[0]
        cdfs[N] = 0.5 * (cdfs[N-2] + 1)

        etas = norm.ppf(cdfs)

        #--- B_1+1/2
        etai   = etas[N-1]
        etaph  = etas[0]
        etap   = 0.5 * (etas[0] + etas[1])
        gph    = rtpi * np.exp(-0.5 * etaph**2)
        bphs[0] = float(N) * gph / (etap - etai)

        #--- B_i+1/2
        etai = 0.5 * (etas[0:N-2] + etas[1:N-1])
        etaph = etas[1:N-1]
        etap = 0.5 * (etas[1:N-1] + etas[2:N])
        gph = rtpi * np.exp(-0.5 * etaph**2)
        bphs[1:N-1] = float(N) * gph / (etap - etai)
        
        #--- B_npt-1/2
        etai      = 0.5 * (etas[N-3] + etas[N-2])
        etaph     = etas[N-2]
        etap      = etas[N]
        gph       = rtpi * np.exp(-0.5 * etaph**2)
        bphs[N-2] = float(N) * gph / (etap - etai)

        #--- construct matrix (I-dt*A)
        Ab = np.zeros((ldab,N))
        Ab[1,0] = 1. + bphs[0] * homdt
        Ab[1,1:N-1] = 1. + (bphs[0:N-2] + bphs[1:N-1]) * homdt
        Ab[1,N-1] = 1. + bphs[N-2] * homdt
        
        Ab[0,1:N] = -bphs[0:N-1] * homdt

        c, x, info = lapack.dpbsv(Ab, phis, lower=0)
        new_phis = np.zeros(N)
        new_phis[idxs] = x
        new_favg = np.sum(self.weights * new_phis) / wtsum
        new_fvar = np.sum(self.weights * new_phis**2) / wtsum - new_favg**2

        fac = np.sqrt(fvar/(new_fvar+1e-30)) * np.exp(-homdt)
        self.phis = favg + (new_phis - new_favg) * fac
        self.var = var(self.phis, self.weights)
        return self.var


""" Euclidean Minimum Spanning Tree model
""" 
class EMST(MixingModel):
    def __init__(self, particles, weights):
        super(EMST, self).__init__(particles, weights)
        self.name = "EMST"
        self.N = len(particles)

    def update(self, Omega_phi, dt):
        dt_in = deepcopy(dt)
        N = self.N

        sorted_phis = np.sort(self.phis.flatten())
        sorted_idxs = np.argsort(self.phis.flatten())

        # get weights
        w = self.weights / np.sum(self.weights[sorted_idxs])
        W = np.cumsum(w)[:-1]
        wv= np.array([min(Wi, 1-Wi) for Wi in W])
        B = (2*wv) # size(Np-1)

        # mixing on the edges
        dphi = np.zeros_like(self.phis)
        v = np.arange(0,N-1)
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

        # # root finding process, unnecessary for 1D EMST
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
    To acceralte the sampling process only uniform weights are supported
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
        Nbin = int(max(20, min(np.floor(np.sqrt(self.N)), 100)))
        minz, maxz = np.min(self.phis) - 1e-8, np.max(self.phis) + 1e-8
        minB = np.arange(Nbin) / Nbin * (maxz - minz) + minz
        maxB = minB + 1 / Nbin * (maxz - minz)
        xbins = np.floor(((self.phis - minz) / (maxz - minz) * Nbin)).astype(int)
        countB = np.bincount(xbins)
        accumB = np.cumsum(countB) - countB
        x =  ((self.phis - minB[xbins]) / (maxB[xbins] - minB[xbins]) * countB[xbins] + accumB[xbins]) / self.N

        W = self.weights
        Wavr = np.mean(self.weights)

        # this sampling process is only suitable for uniform weights 
        p_idxs = np.random.choice(self.N, self.N, replace=True)
        q_idxs = np.random.choice(self.N, self.N, replace=True)

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
        
        # this sampling process is only suitable for uniform weights 
        p_sets = np.random.choice(self.N, n, replace=True)
        q_sets = np.random.choice(self.N, n, replace=True)

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

