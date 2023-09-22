import dna2mhd_utils as dn
import output
import numpy as np

lpaths = ["/scratch/08929/echansen/dna2mhdrun433","/scratch/08929/echansen/dna2mhdrun434","/scratch/08929/echansen/dna2mhdrun435","/scratch/08929/echansen/dna2mhdrun436","/scratch/08929/echansen/dna2mhdrun437"]
#fname = {lpaths[0]:"../input_files/DNAHD.out11499571"}
fname = {lpaths[0]:lpaths[0]+"/DNAHD.out11715556",lpaths[1]:lpaths[1]+"/DNAHD.out11715737",lpaths[2]:lpaths[2]+"/DNAHD.out11716570",lpaths[3]:lpaths[3]+"/DNAHD.out11717610",lpaths[4]:lpaths[4]+'/DNAHD.out11723153'}
a = np.random.randint(0,62,size=4)
for lpath in lpaths:
    dn.plot_energy(lpath,3,tmax=50)
    dn.plot_enspec(lpath,4,zz=-1,log=True,newload=True,show=False,tmax=50)
    for zz in a:
        dn.plot_enspec(lpath,5,zz=zz,log=True,newload=True,show=False,tmaxfac=2)

par = {}
opts = ['b','v','bdv','vdb','cbdb','bdcb','vdv','bdb']

for lpath in lpaths:

    inds = output.maxinds(fname[lpath])
    if lpath == lpaths[4]:
        for j in [[['V',26,10,63],'V'],[['V',3,27,50],'V'],[['V',7,36,35],'V'],[['V',1,37,38],'V'],[['V',10,19,1],'V'],[['V',30,38,2],'V']]:
            inds.append(j)

    dn.getb(lpath,tmax=50)
    dn.getv(lpath,tmax=50)
    dn.read_parameters(lpath)

    for ind in inds:
        for i in range(3):
            dn.plot_bv(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,ask=False)
            dn.plot_nls(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,ask=False,tmax=50)
            for opt in opts:
                dn.plot_bspectrum(lpath,ind[0][1],ind[0][2],ind[0][3],i,show=False,opt=opt)
        print("Plotted ",ind)
