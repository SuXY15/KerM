from models import *
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import normaltest

matplotlib.rc('font', **{'size':14})
np.random.seed(0x7777777)

Omega_phi = 2
dt = 4.e-3
max_steps = 10000

# ==================================
# DNS datasets
# # 1D inert: F2c in DNS_PoF_1988, sigma_x = 0.3
# casename = "PoF_1988_F2c"
# N = 10000    # number of particles
# var_arr = [0.92, 0.80, 0.54, 0.35, 0.28]
# ddt_arr = [0.22, 0.42, 0.83, 1.28, 1.49]
# namerule = "./data/inert/PoF_DNS_1988_F2c_%.2f.txt"
# ylim_0 = [0,3]

# # 1D inert F2e in DNS_PoF_1988
# casename = "PoF_1988_F2e"
# N = 10000  # number of particles
# var_arr = [0.94, 0.76, 0.54, 0.38, 0.27]
# ddt_arr = [0,    0.2,  0.61, 1.06, 1.27]
# namerule = "./data/inert/PoF_DNS_1988_F2e_%.2f.txt"
# ylim_0 = [0,3]

# 1D inert DNS PoF Juneja 1996 Phi2, sigma_x = 0.3
casename = "PoF_1996_Fig9b"
N = 50000  # number of particles
var_arr = [1.0, 0.8, 0.6, 0.5, 0.4, 0.3]
ddt_arr = [0,   0.2, 0.4, 0.5, 0.6, 0.7]
namerule = "./data/inert/Juneja_1996_PoF_Phi2_Var_%.1f.txt"
ylim_0 = [0,4]

# ==================================
# Get DNS data
dns_data = {}
dns_vars = []
for v in var_arr:
    dns_data[v] = deepcopy(np.loadtxt(namerule%v))
    mu, var = PDFstat(dns_data[v][:,0], dns_data[v][:,1])
    dns_vars.append(deepcopy(var))
    print("mu, var, normstd = %.4f %.4f %.4f"%(mu, var, np.sqrt(var) / 0.90))
print()

# ==================================
# Generating particles with given PDF
init_PDF = dns_data[var_arr[0]]
xr = np.min(init_PDF[:,0]), np.max(init_PDF[:,0])
f_init_PDF = interp1d(init_PDF[:,0], init_PDF[:,1]/np.max(init_PDF[:,1]))
particles = np.sort(acceptSampling(f_init_PDF, xr, size=(N,)))

# uniform weights
particles, weights = genSamples(init_PDF, N, method="uniform")

# # weighted PDF, something wrong in python code, not used now
# particles, weights = genSamples(init_PDF, N, method="weighted")

# saving samples to file
samples = np.vstack([particles, weights]).T
np.savetxt("./data/%s_samples.txt"%casename, samples)
np.savetxt("./data/%s_variances.txt"%casename, [len(var_arr)] + var_arr)
sys.exit()

# ==================================
# Select Mixing Models
# TEST1: for different mixing models
models = [ 
    IEM  ( particles, weights ),
    MCurl( particles, weights ),
    EMST ( *genSamples(init_PDF, 400) ),
    KerM ( particles, weights, sigma_x=0.25 ),
]

# # TEST2: for different mixing parameter sigma
# models = [ # mixing models
#     KerMX(particles, sigma_x=0.05),
#     KerMX(particles, sigma_x=0.10),
#     KerMX(particles, sigma_x=0.25),
#     KerMX(particles, sigma_x=1.00),
# ]
# for model in models:
#     model.name = model.name + r"($\sigma_k=%.2f$)"%model.sigma_x

# ==================================
# Simulation the pariticle mixing
y0lim = [1,0]
y1lim = [1,0]
fig, axs = plt.subplots(2, len(models), figsize=(4*len(models),7))
for j, m in enumerate(models):
    print("Start mixing in python: %s (N=%d)"%(m.name, m.phis.shape[0]))
    
    m_vars = [m.var]
    
    j_data = 0 # has reached the jth dns data
    doplot = False
    for i in range(max_steps):

        if j_data<len(var_arr) and m.var < dns_vars[j_data] + 1e-8 or i%10000 == 0:
            doplot = True

        if doplot:
            print("Plot j = %d; i =%-5d, i/ddt=%.1f"%(j_data, i, 0 if j_data==0 else i/ddt_arr[j_data]))
            dns_PDF = dns_data[var_arr[j_data]]
            
            xi, pi = dns_PDF[:,0], dns_PDF[:,1]
            axs[0, j].plot(xi, pi, '--', c=color_arr[j_data], alpha=0.8, label="")

            xi, pi = m.PDF()
            axs[0, j].plot(xi, pi, '-', c=color_arr[j_data], alpha=0.8, label="t=%.2f"%ddt_arr[j_data])

        var_i = m.update(Omega_phi, dt)
        m_vars.append(var_i)

        if doplot:
            j_data += 1
            doplot = False
            if j_data == len(var_arr):
                break

    xi, pi = m.PDF()
    
    axs[0, j].set_title(m.name)
    axs[0, j].set_xlabel(r"$\phi$")

    y0lim[0] = min(y0lim[0], np.min(pi)*0.9)
    y0lim[1] = max(y0lim[1], np.max(pi)*1.1)

    var_i = m.update(Omega_phi, dt)
    m_vars.append(var_i)

    axs[1, j].plot(np.arange(len(m_vars))/2*Omega_phi, m_vars)
    axs[1, j].set_xlabel("Timesteps")
    axs[1, j].set_yscale("log")
    y1lim[0] = min(y1lim[0], np.min(m_vars)*0.9)
    y1lim[1] = max(y1lim[1], np.max(m_vars)*1.1)

    if j!=0:
        axs[0, j].get_yaxis().set_visible(False)
        axs[1, j].get_yaxis().set_visible(False)
    
    print()

axs[0, -1].legend(ncol=2, loc="upper center", frameon=False, 
                    handletextpad=0.2, columnspacing=0.4)
axs[0, 0].set_ylabel("PDF evolution")
axs[1, 0].set_ylabel("Scalar Variance")
for j in range(len(models)):
    axs[0, j].set_ylim(ylim_0)
    axs[1, j].set_ylim(y1lim)

fig.subplots_adjust(left=0.06,bottom=0.10,top=0.9,right=0.98, wspace=0.0, hspace=0.35)
plt.savefig("figs/comparison_python_%s.png"%casename)

plt.show()
