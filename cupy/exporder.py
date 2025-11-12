from dna2mhd_exp import DNA2MHD
import numpy as np
import os

N = 128
nkx0 = N
nky0 = N
nkz0 = 128

kxmin = 0.025
kymin = 0.025
kzmin = 0.005

nu = 0.1
eta = 1.0
hyper = 2
dt = 0.1

a = np.load("perpsample.npy")
n = np.size(a)//11

for i in range(0,n):

    lpath = "/pscratch/sd/e/echansen/interactionangle"+str(int(np.floor(a[i,9])))
    if not os.path.exists(lpath):
        os.mkdir(lpath)
    iterations = int(np.floor(a[i,10]/dt * 0.01))
    triplet = a[i,:9].tolist()
    
    solver = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                 dt,iterations,lpath,linear=False,explicitrk4=False,
                 initcond="threewave",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                 forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],hyper=hyper,hallparam=1.0,
                 solveprec=16,maxwallclock=86200,triplet=triplet,records=20)

    solver.ralstonrk2()

#solver.etdrk2()

#solver.dp547s()
