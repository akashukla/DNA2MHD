import numpy as np
import matplotlib.pyplot as plt
import dna2mhd_utils as dn

tempstr = "/scratch/08929/echansen/dna2mhdrun"

Ns = []
for i in range(4,8):
    Ns.append(2**i)

ss = ["17","18","19"]
ts = ["0.005","0.01","0.1"]

plt.figure()
for i in Ns:
    for j in range(3):
        lpath = tempstr+str(i)+"s"+ss[j]
        de,dmh,dch = dn.enheldev(lpath)
        if i == 16:
            plt.plot(i,np.abs(de),label=ts[j])
        else:
            plt.plot(i,np.abs(de))
        plt.plot()

plt.xlabel("N")
plt.ylabel("Abs Energy Deviation")
plt.title("SphLike Coordinate Energy Deviation for Varied Time Step")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc=1)
plt.savefig("endev1")


