"""
Script to order one or more three wave interaction simulations from submitruns.py
"""


from dna2mhd_exp import DNA2MHD
import numpy as np
import os
import sys

start = np.int32(sys.argv[1])

factor = 1
N = 256
nkx0 = N*factor
nky0 = N*factor
nkz0 = N*factor

for i in range(start,start+1):

    
    # Use h5py file for initial conditions
    # This file must contain - wavenumber perp and parallel grid scales
    # Viscosity determined with microscale in simulation domain
    # Resistivity determined from microscale viscosity and "realistic" Prandtl number (cross field viscosity to Spitzer resistivity)
    # Number of iterations to simulate
    # Time step for plasma domain - determined from fastest Alfvén whistler scale
    # Resonant wave interaction wavenumbers in ascending wavenumber order
    # If using other waves than positive whistlers, add three mode types
    
    with h5py.File("initconds032426.hdf5","r") as f:

        # Setting number - specific to file with 5 interactions, 5 settings
        # 0 - LAPD, 1 - HSX, 2 - DIII-D, 3 - Coronal Loop, 4 - Solar Wind
        seti = i // 5
        # Interaction number - specific to file with 5 interactions, 5 settings
        intnum = i % 5
        
        triplets = f["waves"][seti,intnum,:]

        kzmin = f["kzmins"][seti]
        kxmin = f["kpmins"][seti]
        iterations = f["iterations"][seti]
        dt = f["timesteps"][seti]
        nu = f["viscs"][seti]
        eta = f["etas"][seti]
        fname = f["fnames"][seti]

    hyper = 1
    lpath = "/pscratch/sd/e/echansen/threewaves032426/"+fname+str(intnum)+"/"
    print(lpath,nu,eta,kzmin,kxmin)
    if not os.path.exists(lpath):
        os.mkdir(lpath)
    test_iterations = iterations//1000
    triplet = triplets.tolist()
    
    solver = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                     dt,iterations,lpath,
                     linear=False,
                     initcond="threewave",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                     forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],hyper=hyper,hallparam=1.0,
                     solveprec=16,maxwallclock=86200,triplet=triplet,records=50,bittype=64,exactnueta=True)

    solver.gauss2()
