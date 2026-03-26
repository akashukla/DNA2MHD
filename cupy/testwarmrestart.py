"""
Script to order one or more three wave interaction simulations from submitruns.py
Testing for the 
"""


from dna2mhd import DNA2MHD
import numpy as np
import os
import sys
import h5py

start = np.int32(sys.argv[1])

factor = 1
N = 64
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
        kymin = kxmin
        
        iterations = f["iterations"][seti]
        dt = f["timesteps"][seti]
        nu = f["viscs"][seti]
        eta = f["etas"][seti]
        fname = f["fnames"][seti]

    hyper = 1
    lpath1 = "/pscratch/sd/e/echansen/warmrestart1/"
    print(lpath1,nu,eta,kzmin,kxmin)
    if not os.path.exists(lpath1):
        os.mkdir(lpath1)
    test_iterations = iterations//1000
    triplet = triplets.tolist()
    
    
    solver1 = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                     dt,test_iterations,lpath1,
                     linear=False,
                     initcond="threewave",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                     forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],hyper=hyper,hallparam=1.0,
                     solveprec=16,maxwallclock=86200,triplet=triplet,records=50,bittype=64,exactnueta=True)

    solver1.gauss2()

    # Now define second simulation with half time steps which will be done twice
    lpath2 = "/pscratch/sd/e/echansen/warmrestart2/"
    if not os.path.exists(lpath2):
        os.mkdir(lpath2)

    solver2a = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                     dt,test_iterations//2,lpath2,
                     linear=False,
                     initcond="threewave",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                     forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],hyper=hyper,hallparam=1.0,
                     solveprec=16,maxwallclock=86200,triplet=triplet,records=25,bittype=64,exactnueta=True)

    solver2a.gauss2()

    # Use checkpoint initial conditions to finish up to test_iterations
    solver2b = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                     dt,test_iterations//2,lpath2,
                     linear=False,
                     initcond="checkpoint",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                     forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],hyper=hyper,hallparam=1.0,
                     solveprec=16,maxwallclock=86200,triplet=triplet,records=25,bittype=64,exactnueta=True)

    solver2b.gauss2()

    def readoutput(f):

        written = f["written"][0]
        
        b1 = f["magneticfields"][written,:,:,:,:]
        v1 = f["velocityfields"][written,:,:,:,:]
        print(f["time"][written])

        return(b1,v1)
    
    with h5py.File(lpath1+"output.hdf5","r") as f:

        b1,v1 = readoutput(f)
        
    with h5py.File(lpath2+"output.hdf5","r") as f:

        b2,v2 = readoutput(f)

    print("Magnetic Field Diffs",np.amax(np.abs(b1-b2)))
    print("Velocity Field Diffs",np.amax(np.abs(v1-v2)))
