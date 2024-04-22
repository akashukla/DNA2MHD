import numpy as np
import matplotlib.pyplot as plt
import dna2mhd_utils as dn

tempstr = "/scratch/08929/echansen/dna2mhdrun"

start4 = 986
start3 = 1036
start2 = 1136

start4s = range(start4,start4+50,10)
start3s = range(start3,start3+50,10)
start2s = [1146,1136,1156,1166,1176]

intes = [2,3,4]
intmap = {2:start2s,3:start3s,4:start4s}
fmtmap = {2:"m+",3:"ro",4:"b*"}
legmap = {2:"Ralston RK2",3:"Ralston RK3",4:"Classic RK4"}
Ns = [128,64,32,16,8]

plt.figure()
for i in intes:
    means = []
    for j in range(len(intmap[i])):
        des = []
        for k in range(10):
            lpath = tempstr+str(intmap[i][j]+k)
            print(lpath)
            de,dmh,dch = dn.enheldev(lpath)
            des.append(de)
        debar = np.mean(np.abs(des))
        means.append(debar)
        destdm = np.std(np.abs(des))/np.sqrt(10)
    plt.plot(Ns,means,fmtmap[i],label=legmap[i])

plt.xlabel("N")
plt.ylabel("Abs Energy Deviation")
plt.title("Energy Deviation for Varied Integrators at T = 1, dt = 0.1")
plt.xscale("log")
plt.yscale("log")
plt.ylim(10**(-9),10**(-5))
plt.legend(loc=1)
plt.savefig("endev1")


