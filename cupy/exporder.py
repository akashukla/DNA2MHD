from dna2mhd_exp import DNA2MHD
import dna2mhd_utils as dn
import numpy as np
import os
import sys

start = np.int32(sys.argv[1])

factor = 1
N = 128
nkx0 = N*factor
nky0 = N*factor
nkz0 = N*factor

kxmin = 0.025
kymin = 0.025
kzmin = 0.005

nu = 0.1
eta = 1.0
hyper = 2
nu *= (factor)**(2*hyper)

dt = 0.05

a = np.load("norm35sample.npy")
n = np.size(a)//12

for i in range(start,start+1):

    lpath = "/pscratch/sd/e/echansen/64interactionangle"+str(int(np.floor(a[i,9])))
    print(lpath)
    if not os.path.exists(lpath):
        os.mkdir(lpath)
    iterations = int(np.floor(a[i,-1]/dt * 0.001))
    triplet = a[i,:9].tolist()
    
    solver = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                     dt,iterations,lpath,linear=False,explicitrk4=False,
                     initcond="threewave",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                     forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],hyper=hyper,hallparam=1.0,
                     solveprec=16,maxwallclock=86200,triplet=triplet,records=50,bittype=64)

    #solver.ralstonrk2()
    #solver.etdrk2()
    solver.dp547s()
    #solver.gauss2()

    dn.plot_energy(lpath,checkenergyonly=True)
    #dn.enheldev(lpath)
    #dn.plot_energy(lpath)
    #dn.plot_enspec(lpath,zz=-1,version=0)
    #dn.mode_break(lpath,show=False)
    #dn.threewaveenergy(lpath)
    #dn.nonlinearities(lpath)
    
