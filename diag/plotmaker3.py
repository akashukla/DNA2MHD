import dna2mhd_utils as dn
from output2 import maxinds
import numpy as np

lpaths = ["/pscratch/sd/e/echansen/DNA2MHDruns/siam_whist",
          "/pscratch/sd/e/echansen/DNA2MHDruns/siam_whist_ideal"]

for lpath in lpaths:
    dn.plot_energy(lpath)
    print("\nThrough energy\n")
     
    dn.plot_enspec(lpath,zz=-1,version=3)
    dn.plot_enspec(lpath,zz=2,version=0)
    print("\nThrough Energy Spec\n")
    x = """
    dn.nlparam(lpath)
    print("\nThrough NL Param\n")
    """
    
    dn.mode_break(lpath,show=False)
    print("\nThrough Mode Breakdown\n")
    dn.structurefunction(lpath,tmax=2*10**10)
    print("\nThrough Structure Functions\n")
    x = """
    dn.mode_nlparam(lpath,0,1)
    dn.mode_nlparam(lpath,0,3)
    dn.mode_nlparam(lpath,-1,1)
    dn.mode_nlparam(lpath,-1,3)
    print("Through Nonlinearity Parameter")
    """

x = """inds = []
for i in range(10):
    inds.append(np.random.randint(0,16,3))

for lpath in lpaths:

    inds = []
    for i in range(10):
        inds.append(np.random.randint(0,16,3))

    for ind in inds:
        for i in range(3):
            dn.plot_bv(lpath,ind[0],ind[1],ind[2],i,show=False)
            dn.plot_bvspectrum(lpath,"b",ind[0],ind[1],ind[2],i,show=False)
            dn.plot_bvspectrum(lpath,"v",ind[0],ind[1],ind[2],i,show=False)
                
        print("Plotted ",ind)
"""
