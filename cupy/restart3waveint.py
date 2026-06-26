"""
Script to order one or more three wave interaction simulations from submitruns.py
"""


from dna2mhd import DNA2MHD
import numpy as np
import os
import sys
import h5py

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

     # Setting number - specific to file with 5 interactions, 5 settings                                                                                                                                                                       
        # 0 - LAPD, 1 - HSX, 2 - DIII-D, 3 - Coronal Loop, 4 - Solar Wind
    

    seti = i // 5

    # Interaction number - specific to file with 5 interactions, 5 settings 
    intnum = i % 5
    
    with h5py.File("initconds060426.hdf5","r") as f:

        triplets = f["waves"][seti,intnum,:]

        kzmin = f["kzmins"][seti]
        kxmin = f["kpmins"][seti]
        kymin = kxmin
        iterations = f["iterations"][seti]
        fname = f["fnames"][seti].astype("T") # Convert from numpy bytes object to string

    lpath = "/pscratch/sd/e/echansen/threewaves060426/"+fname+str(N)+"set"+str(intnum)+"/"

    with h5py.File(lpath+"output.hdf5","r") as f:

        dt = f["dt"][0]
        nu = f["nu"][0]
        eta = f["eta"][0]
        hyper = f["hyper"][0]
        
        
    print(lpath,nu,eta,kzmin,kxmin)
    triplet = triplets.tolist()
    
    solver = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                     dt,iterations,lpath,
                     linear=False,
                     initcond="checkpoint",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                     forcetype="threewave",hyper=hyper,hallparam=1.0,
                     solveprec=16,maxwallclock=86200,triplet=triplet,records=50,bittype=64,exactnueta=True)

    solver.gauss2split()
