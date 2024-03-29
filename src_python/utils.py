import os, sys, time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.linalg import lapack
from copy import deepcopy
from itertools import accumulate
from scipy.interpolate import interp1d

color_arr = ('k','r','b','m','y','g','k','r','b','m','y')
symbol_arr = ('s','o','v','^','p','*','s','o','v','^','p')

def mean(data, weights=None):
    # return np.mean(data)
    w = weights / np.sum(weights)
    return np.sum(data*w)

def var(data, weights=None):
    # return np.var(data)
    w = weights / np.sum(weights)
    m = mean(data, w)
    return np.sum(w*(data-m)**2)

def PDF(data, weights):
    NBins = min(20, max(int(np.sqrt(len(data))), 100))
    hist, bins = np.histogram(data, weights=weights, bins=NBins, density=True)
    bin_centers = (bins[1:]+bins[:-1]) * 0.5
    return bin_centers, hist

""" Central differential of data
"""
def cdiff(data):
    if len(data)==1: return np.array([0])
    d = np.diff(data)
    return np.array([d[0]] + list((d[1:]+d[:-1])/2) + [d[-1]])

""" Central mean of data
"""
def cmean(data):
    if len(data)==1: return np.array([0])
    d = data
    return np.array([d[0]] + list((d[1:]+d[:-1])/2) + [d[-1]])

""" Get Norm2 distance between phi_p and phi_q
"""
def distance(d_pq):
    if len(d_pq.shape)==1:
        return np.abs(d_pq)
    if len(d_pq.shape)>1:
        return np.linalg.norm(d_pq, axis=len(d_pq.shape)-1).reshape(len(d_pq),)

""" Accept-Reject Sampling method
    f(x) is the PDF function
        0 <= f(x) <= 1 should be satisfied
        max(f(x))=1 would be better
"""
def acceptSampling_i(f, xr):
    xmin, xmax = xr
    while True:
        x = np.random.rand()*(xmax-xmin) + xmin
        if np.random.rand() < f(x):
            return x

""" Sampling for given size
"""
def acceptSampling(f, xr, size=(1,)):
    data = np.zeros(size)
    if len(size)>1:
        data = [acceptSampling(f, xr, size=size[1:]) for d in data]
    else:
        data = [acceptSampling_i(f, xr) for d in data]
    return np.array(data)

""" PDF statistics: mean and variance
"""
def PDFstat(xi, yi):
    yi = (cdiff(xi)*yi)
    yi = yi/np.sum(yi)
    mu = np.sum(yi*xi)
    var = np.sum(yi*(xi-mu)**2)
    return mu, var

def smooth(data, window=5): 
    # data: NumPy 1-D array containing the data to be smoothed 
    # window: smoothing window size needs, which must be odd number
    out0 = np.convolve(data,np.ones(window,dtype=int),'valid')/window  
    r = np.arange(1,window-1,2) 
    start = np.cumsum(data[:window-1])[::2]/r 
    stop = (np.cumsum(data[:-window:-1])[::2]/r)[::-1] 
    return np.concatenate(( start, out0, stop )) 

def genSamples(PDFxy, N, method="uniform"):
    xr = np.min(PDFxy[:,0]), np.max(PDFxy[:,0])
    f_PDF = interp1d(PDFxy[:,0], PDFxy[:,1]/np.max(PDFxy[:,1]))
    particles = np.sort(acceptSampling(f_PDF, xr, size=(N,)))

    # uniform weights
    if method == "uniform":
        weights = np.ones(N)/N
    else:
        dx = cdiff(particles)
        weights = f_PDF(particles)*smooth(dx,3)
        weights = weights / np.sum(weights)
    idx = np.arange(N)
    np.random.shuffle(idx)
    return particles[idx], weights[idx]

if __name__ == "__main__":
    PDFxy = np.array([[0,0],[1,1]])
    xi, pi = PDF(*genSamples(PDFxy, N=4000, method="uniform"))
    plt.plot(xi, pi, '-', label='uniform')
    xi, pi = PDF(*genSamples(PDFxy, N=4000, method="weighted"))
    plt.plot(xi, pi, '--', label='weighted')
    plt.legend()
    plt.show()