import dna2mhd_utils as dn
from output2 import maxinds
import numpy as np

lpaths = ["/scratch/08929/echansen/dna2mhdrun442"]

fname = {lpaths[0]:lpaths[0]+"/DNAHD.out11725125"}
a = np.random.randint(0,31,size=4)

for lpath in lpaths:
    dn.plot_energy(lpath,3)
    dn.plot_enspec(lpath,4,zz=-1,log=True,newload=True,show=False)
    for zz in a:
        dn.plot_enspec(lpath,5,zz=zz,log=True,newload=True,show=False,tmaxfac=2)

par = {}
opts = ['b','v','bdv','vdb','cbdb','bdcb','vdv','bdb']

for lpath in lpaths:

    inds = maxinds(fname[lpath],32)
    dn.getb(lpath)
    dn.getv(lpath)
    dn.read_parameters(lpath)

    for ind in inds:
        for i in range(3):
            dn.plot_bv(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,ask=False)
            dn.plot_nls(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,ask=False)
            for opt in opts:
                dn.plot_bspectrum(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,opt=opt)
        print("Plotted ",ind)
