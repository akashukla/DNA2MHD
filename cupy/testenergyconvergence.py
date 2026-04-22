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
        fname = f["fnames"][seti].astype("T") # Convert from numpy bytes object to string

    hyper = 1
    
    test_iterations = iterations//10000 * 8
    triplet = triplets.tolist()

    # Adjust viscosities from N = 256 value if needed - place microscale at N/sqrt(2)
    nu *= 1/((N/256)**(4/3))
    eta *= 1/((N/256)**(4/3))

    # Adjust time step - the file time steps are small by 2.5 because of 3/2 padding
    # But adjust if using stiff integrator
    # dt *= 1/3
    dt0 = dt * (256/N)**2 * 3/4 * 1/5

    #Set nu and eta to zero for energy convergence test
    nu = 0.0
    eta = 0.0

    #Loop to test time step convergence
    dtlist = []
    energylist = []
    blist = []
    vlist = []
    enfieldlist = []
    finenlist = []
    hlist = []
    def enfromarray(a):

        zeros = np.sum(np.abs(a[:,:,:,:,0])**2,axis=(1,2,3))
        nonzeros = 2*np.sum(np.abs(a[:,:,:,:,1:])**2,axis=(1,2,3,4))
        return(4*np.pi**3 *(zeros+nonzeros))

    for dti in range(0,4):

        test_iterations_dti = test_iterations * 2**dti
        dt = dt0 / (2**dti)
        dtlist.append(dt)
        lpath = "/pscratch/sd/e/echansen/threewaves032426/"+fname+str(intnum)+"test"+str(dti)+"/"
        print(lpath,nu,eta,kzmin,kxmin)
        if not os.path.exists(lpath):
            os.makedirs(lpath)


        solver = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                     dt,test_iterations_dti,lpath,
                     linear=False,
                     initcond="threewave",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                     forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],hyper=hyper,hallparam=1.0,
                     solveprec=16,maxwallclock=86200,triplet=triplet,records=20,bittype=64,exactnueta=True)
        
        #solver.ralstonrk2()
        solver.explicit(order=4)
        #solver.exptime_ideal(order=2,kt=0)
        #solver.gauss2split()
        with h5py.File(lpath+"/output.hdf5","a") as f:
            energy = f["enval"][:,0]
            b = f["magneticfields"][:,:,:,:,:]
            v = f["velocityfields"][:,:,:,:,:]
            helicity = f["enval"][:,1]
            mhelcorr = f["enval"][:,-1]

        energylist.append(np.amax(np.abs(energy-energy[0])))
        hlist.append(np.amax(np.abs(helicity+mhelcorr-helicity[0])))
        blist.append(b)
        vlist.append(v)
        finenlist.append(energy[-1])

        energyfield = enfromarray(b)+enfromarray(v)
        print(energyfield[1:]-energyfield[:-1])
        enfieldlist.append(np.amax(energyfield-energyfield[0]))
        if np.amax(energyfield-energy) > 1e-16:
            print("Energy diag broken")
            print("diag ",energy)
            print("field ",energyfield)
            print("diff ",np.abs(energyfield-energy))

    diffblist = []
    diffvlist = []
    diffenlist = []

    for i in range(0,3):
        diffblist.append(np.amax(np.abs(blist[i]-blist[-1])))
        diffvlist.append(np.amax(np.abs(vlist[i]-vlist[-1])))
        diffenlist.append(np.amax(np.abs(finenlist[i]-finenlist[-1])))

    print(dtlist)
    print(energylist)
    print(diffblist)
    print(diffvlist)
    print(enfieldlist)
    print(diffenlist)
    print(hlist)