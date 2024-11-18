import dna2mhd_utils as dn
from output2 import maxinds
import numpy as np

lpaths = ["/pscratch/sd/e/echansen/DNA2MHDruns/fullwhisti"]

#lpaths = ["/scratch/08929/echansen/dna2mhdrunDPP128LC",
#          "/scratch/08929/echansen/dna2mhdrunDPP128LCi"]

for lpath in lpaths:
    dn.plot_energy(lpath,3,show=False)
    print("Through energy")
    dn.plot_enspec(lpath,npt=4,zz=-1,show=False,log=True,linplot=False,newload=True,fullspec=False,old=True,tmaxfac=1,tmax=2000000)
    dn.plot_enspec(lpath,npt=4,zz=-1,show=False,log=True,linplot=False,newload=True,fullspec=False,old=True,tmaxfac=1,tmax=2000000,version=0)
    print("Through Energy Spec")
    dn.mode_break(lpath,show=False)
    print("Through Mode Breakdown")
    dn.structurefunction(lpath,tmax=2*10**10)
    print("Through Structure Functions")
    dn.mode_nlparam(lpath,0,1)
    dn.mode_nlparam(lpath,0,3)
    dn.mode_nlparam(lpath,-1,1)
    dn.mode_nlparam(lpath,-1,3)
    print("Through Nonlinearity Parameter")

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
