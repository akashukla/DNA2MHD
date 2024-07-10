import numpy as np
import matplotlib.pyplot as plt
import dna2mhd_utils as dn
from matplotlib.ticker import LogFormatter
lpath_base = "/scratch/08929/echansen/dna2mhdrunhc"
Ns = [16,32,64]
dts = [0.05,0.1,0.2]
dts_i = ["1","2","3"]
orders = ["2","3","4","5"]
fmts = ["rs","bs","rv","bv"]
labels = ["RK2","RK3","RK4","DP5"]
zs = ['','z']

for z in zs:
    fig,ax = plt.subplots(1,3)
    for i in range(3):
        ax[i].set_title("N = "+str(Ns[i]))
        for j in range(3):
            for k in range(4):
                lpath = lpath_base + str(Ns[i])+str(orders[k])+str(dts_i[j])+z
                print(lpath)
                de = dn.enheldev(lpath)[0]
                print(de)
                if i == 0 and j == 0:
                    ax[i].plot(dts[j],np.abs(de),fmts[k],label=labels[k])
                else:
                    ax[i].plot(dts[j],np.abs(de),fmts[k])
        ax[i].set_xscale("log")
        ax[i].set_yscale("log")
        ax[i].set_xlim(0.01,1)
        ax[i].set_ylim(10.0**(-16.0),10.0**(-6.0))
        ax[i].set_xlabel("dt ($\omega_c^{-1}$)")
        ax[i].xaxis.set_major_formatter(LogFormatter(labelOnlyBase=True))
        ax[i].xaxis.set_minor_formatter("")
    ax[0].set_ylabel("Energy Deviation")
    fig.legend(loc=4)
    fig.suptitle("Energy Deviations")
    fig.tight_layout()
    plt.savefig("convgdt0627"+z)
    
