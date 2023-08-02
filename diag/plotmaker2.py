import dna2mhd_utils as dn
import output
import numpy as np

lpaths = ["/scratch/08929/echansen/dna2mhdrun419","/scratch/08929/echansen/dna2mhdrun420"]
#fname = {lpaths[0]:"../input_files/DNAHD.out11499571"}
fname = {lpaths[0]:lpaths[0]+"/DNAHD.out11523448",lpaths[1]:lpaths[1]+"/DNAHD.out11523449"}
a = np.random.randint(0,62,size=4)

for lpath in lpaths:
    dn.plot_energy(lpath,3)
    dn.plot_enspec(lpath,4,zz=-1,log=True,newload=True,show=False)
    for zz in a:
        dn.plot_enspec(lpath,5,zz=zz,log=True,newload=True,show=False,tmaxfac=2)

for lpath in lpaths:

    inds = output.maxinds(fname[lpath])
    dn.getb(lpath)
    dn.getv(lpath)
    
    for ind in inds:
        for i in range(3):
            dn.plot_bv(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,ask=False)
            dn.plot_bspectrum(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False)
            dn.plot_vspectrum(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False)
        print("Plotted ",ind)
