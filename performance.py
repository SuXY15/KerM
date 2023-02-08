from src_python.utils import *
matplotlib.rc('font', **{'size':14})

PATH = "./data/"

# # 1D inert: F2c in DNS_PoF_1988
# casename = "PoF_1988_F2c"
# var_arr = [0.92, 0.80, 0.54, 0.35, 0.28]
# ddt_arr = [0.22, 0.42, 0.83, 1.28, 1.49]
# namerule = "./data/inert/PoF_DNS_1988_F2c_%.2f.txt"
# ylim_0 = [0,3]

# # 1D inert F2e in DNS_PoF_1988
# casename = "PoF_1988_F2e"
# var_arr = [0.94, 0.76, 0.54, 0.38, 0.27]
# ddt_arr = [0,    0.2,  0.61, 1.06, 1.27]
# namerule = "./data/inert/PoF_DNS_1988_F2e_%.2f.txt"
# ylim_0 = [0,3]

# sigma_x = 0.25
# 1D inert DNS PoF Juneja 1996 Phi2
casename = "PoF_1996_Fig9b"
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

def fit(N,T,k,log=False):
    i = len(N)>>1
    y = N**k / N[i]**k*T[i]
    return y*np.log(N) / np.log(N[i]) if log else y

def compare_cost(lang="cpp"):
    if lang=="cpp":
        MMs = ["IEM", "MC",  "EMST-1D", "KerM"]
        Slp = [1.,    1.,    3.,        1.    ]  # slope in log scale
        Log = [False, False, True,      False ]
    elif lang=="fortran":
        MMs = ["IEM", "MC",  "MCMG", "EMST", "KerM",]
        Slp = [1,     1,     1,      2,      1      ]  # slope in log scale
        Log = [False, False, False,   False,  False  ]
    else:
        print("ERROR::only cpp and fortran are supported")

    # ==================================
    # Get and show COST info
    print("Hanlding the costs and plots ...")
    plt.figure(figsize=(6,4))
    Tmin, Tmax = 1e10, -1e10
    for i,mm in enumerate(MMs):
        data = np.loadtxt(PATH + casename + "_%s_"%lang + mm + "_costs.txt")
        N, T, k = data[:,0], data[:,1], Slp[i]
        Nuniq = np.array(sorted(set(N)))
        Tmean = [np.mean(T[N==n]) for n in Nuniq]
        Tstd = [np.std(T[N==n]) for n in Nuniq]
        T_min = [np.min(T[N==n]) for n in Nuniq]
        N,T = Nuniq, T_min
        # plt.errorbar(N, T, yerr=Tstd, fmt=color_arr[i]+symbol_arr[i], fillstyle="none", label=mm)
        plt.plot(N, T, color_arr[i]+symbol_arr[i], fillstyle="none", label=mm)
        plt.plot(N, fit(N,T,k,Log[i]), color_arr[i]+'--', label="$N"+("^%d"%k if k!=1 else "") + ("log(N)$" if Log[i] else "$"))
        Tmin, Tmax = min(Tmin, np.min(T)), max(Tmax, np.max(T))

    plt.legend(frameon=False, ncol=len(MMs), loc="upper center",
                handletextpad=0.15, columnspacing=0.3)
    plt.ylim([Tmin / 2, Tmax*200])
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("N")
    plt.ylabel("Time Cost [s]")
    plt.subplots_adjust(left=0.15,bottom=0.15,top=0.95,right=0.98, wspace=0.0)
    plt.savefig("figs/performance_%s_%s.png"%(lang,casename))

    # ==================================
    # Get and show PDF info
    fig, axs = plt.subplots(1, len(MMs), figsize=(4*len(MMs),4))
    for i,mm in enumerate(MMs):
        print("Handling data and plots for %s ..."%mm)
        for j,var in enumerate(var_arr):
            dns_PDF = dns_data[var]
            xi, pi = dns_PDF[:,0], dns_PDF[:,1]
            axs[i].plot(xi, pi, '--', c=color_arr[j], alpha=0.8, label="")

            data = np.loadtxt(PATH+"%s_%s_%s_%.6f.txt"%(casename,lang,mm,var))
            xi, pi = PDF(data[:,0], data[:,1])
            axs[i].plot(xi, pi, '-', c=color_arr[j], alpha=0.8, label="t=%.2f"%ddt_arr[j])

        axs[i].set_title(mm)
        axs[i].set_xlabel(r"$\phi$")
        if i==0:
            axs[0].set_ylabel("PDF evolution")
        else:
            axs[i].get_yaxis().set_visible(False)
        axs[i].set_ylim(ylim_0)
        axs[i].legend(ncol=2, loc="upper right", frameon=False, 
                        handletextpad=0.2, columnspacing=0.4)
    fig.subplots_adjust(left=0.05,bottom=0.15,top=0.9,right=0.98, wspace=0.0)
    plt.savefig("figs/comparison_%s_%s.png"%(lang,casename))

if len(sys.argv)>1:
    compare_cost(sys.argv[1])
else:
    compare_cost(lang="cpp")

plt.show()
