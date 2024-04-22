import dna2mhd_utils as dn
from output2 import maxinds
import numpy as np

lpaths = ["/scratch/08929/echansen/dna2mhdrun745"]
#lpaths = ["/scratch/08929/echansen/dna2mhdrun452","/scratch/08929/echansen/dna2mhdrun453"]

a = np.random.randint(0,63,size=4)

for lpath in lpaths:
    dn.plot_energy(lpath,3,show=False)
    dn.plot_enspec(lpath,4,zz=-1,log=True,newload=True,show=False)
    dn.plot_xispec(lpath,4,zz=-1,log=True,show=False)
    dn.mode_break(lpath,show=False)
    for zz in a:
        dn.plot_enspec(lpath,5,zz=zz,log=True,newload=True,show=False,tmaxfac=2)
        dn.plot_xispec(lpath,5,zz=zz,log=True,show=False,tmaxfac=2)
    dn.planeplotter(lpath,0,show=False)
    dn.planeplotter(lpath,-1,show=False)

x = """
fname = {lpaths[0]:lpaths[0]+"/DNAHD.out1467543",lpaths[1]:lpaths[1]+"/DNAHD.out1472208"}
opts = ['b','v','bdv','vdb','cbdb','bdcb','vdv','bdb']

for lpath in lpaths:

    inds = maxinds(fname[lpath],64)
    dn.getb(lpath)
    dn.getv(lpath)
    dn.read_parameters(lpath)

    for ind in inds:
        dn.plot_profile_xi(lpath,ind[0][1],ind[0][2],ind[0][3],show=False)
        for t in ['b','t','v']:
            dn.plot_profile_enspec(lpath,ind[0][1],ind[0][2],ind[0][3],tbv=t,show=False)
        for i in range(3):
            dn.plot_bv(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,ask=False)
            dn.plot_nls(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,ask=False)
            for opt in opts:
                dn.plot_bspectrum(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,opt=opt)
        print("Plotted ",ind)
"""
