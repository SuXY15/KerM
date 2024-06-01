import os, sys, time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import normaltest
from copy import deepcopy

matplotlib.rc('font', **{'size':14})
np.random.seed(0x7777777)

def mean(data, weights=None):
    if weights is None:
        return np.mean(data, axis=0)
    w = weights / np.sum(weights)
    return np.sum(data*w, axis=0)

def var(data, weights=None):
    if weights is None:
        return np.var(data, axis=0)
    w = weights / np.sum(weights)
    m = mean(data, w)
    return np.sum(w*(data-m)**2, axis=0)

""" Get Norm2 distance between phi_p and phi_q
"""
def distance(d_pq):
    if len(d_pq.shape)==1:
        return np.abs(d_pq)
    if len(d_pq.shape)>1:
        return np.linalg.norm(d_pq, axis=len(d_pq.shape)-1).reshape(len(d_pq),)

""" Mixing Model base class, for 2D particles
"""
class MixingModel():
    def __init__(self, particles, weights=None):
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


""" Kernel Mixing model, with Localness in X space
    To acceralte the sampling process only uniform weights are supported
"""
class KerM(MixingModel):
    def __init__(self, particles, weights=None, sigma_x=np.inf):
        super(KerM, self).__init__(particles, weights)
        self.name = "KerM"
        self.N = len(particles)
        self.sigma_x = sigma_x # mixing length scale for gaussian kernel

    def kernel_func(self, d):
        return np.exp(-d**2 / self.sigma_x**2 / 4).reshape(len(d),)

    def update(self, Omega_phi, dt):
        Nc = self.N*2

        # this sampling process is only suitable for uniform weights 
        p_idxs = np.random.choice(N, Nc, replace=True)
        q_idxs = np.random.choice(N, Nc, replace=True)

        d_pq = distance(self.phis[p_idxs] - self.phis[q_idxs])
        f_pq = self.kernel_func(d_pq)
        coeff = self.var / np.sum( 0.5 * f_pq * d_pq**2 / len(d_pq))

        # print("coeff:", coeff)
        # print("Coeff=%.3f"%(coeff))

        n = int(3/2 * Omega_phi * self.N * dt * np.mean(coeff) + 1)
        p = 3/2 * Omega_phi * self.N * dt * np.mean(coeff) / n
        
        if n>self.N:
            print("Error: n(%d) > N(%d), please use smaller `dt` or `Omega_phi`"%(n,self.N))

        self.dphidt = np.zeros_like(self.phis)
        
        # this sampling process is only suitable for uniform weights 
        p_sets = np.random.choice(self.N, n, replace=True)
        q_sets = np.random.choice(self.N, n, replace=True)

        d_pq = distance(self.phis[p_sets] - self.phis[q_sets])
        f = self.kernel_func(d_pq)
        f_idx = np.where(np.random.rand(*d_pq.shape) < f, True, False)
        p_idx = np.where(np.random.rand(*d_pq.shape) < p, True, False)
        p_idxs = p_sets[f_idx & p_idx]
        q_idxs = q_sets[f_idx & p_idx]

        d_pq = self.phis[p_idxs] - self.phis[q_idxs]
        alpha = np.random.rand(d_pq.shape[0])
        self.dphidt[p_idxs] = - np.multiply(d_pq.T, alpha / 2).T / dt
        self.dphidt[q_idxs] =   np.multiply(d_pq.T, alpha / 2).T / dt
        return super(KerM, self).update(Omega_phi, dt)

N = 5000  # number of particles
Omega_phi = 2
dt = 1.e-3
max_steps = 2000

# generate initial PDF
v = np.random.normal(0, 1, size=(N*9))
v = v[np.logical_and(v>-3, v<3)]
z1 = v[:N*2].reshape(N,2)
z2 = v[N*2:N*4].reshape(N,2)
z3 = v[N*4:N*6].reshape(N,2)
z1[:,0] = z1[:,0]*0.1-1.0
z1[:,1] = z1[:,1]*0.1-0.5773
z2[:,0] = z2[:,0]*0.1+1.0
z2[:,1] = z2[:,1]*0.1-0.5773
z3[:,0] = z3[:,0]*0.1
z3[:,1] = z3[:,1]*0.1+ 1.1547
particles = np.vstack([z1, z2, z3])

m = KerM ( particles, weights=None, sigma_x=0.25*2)
print(m.var)
m_var0 = np.mean(m.var)

plt.ion()
fig, axs = plt.subplots(2,2,figsize=(8,6))
fig.subplots_adjust(left=0.07, bottom=0.1, top=0.92, right=0.95, wspace=0.3, hspace=0.3)

k = 0
for i in range(max_steps):
    PHI = np.sqrt(np.mean(m.var)/m_var0)
    if (1-PHI)*10 >= k:
        row = k//2
        col = k%2
        k = k+1
        axs[row,col].plot(m.phis[:,0], m.phis[:,1], 'C0.', ms=2, alpha=0.2)
        axs[row,col].set_xlim([-1.4, 1.4])
        axs[row,col].set_ylim([-1, 1.6])
        axs[row,col].set_title(r"$\Phi={:.2f}$".format(PHI))
        plt.draw()
        plt.pause(1e-8)
    if k==4:
        break

    var_i = m.update(Omega_phi, dt)
        
plt.ioff()
plt.show()