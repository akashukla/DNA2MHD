"""
Script to turn off and compare nonlinearites and compare invariants
"""

from dna2mhd import DNA2MHD
import numpy as np
import os
import sys
import h5py

start = np.int32(sys.argv[1])

N = 256
nkx0 = N
nky0 = N
nkz0 = N//2

wave = 27
    
with h5py.File("initconds061726.hdf5","r") as f:

    # Setting number - specific to file with 5 interactions, 5 settings
    # 0 - LAPD, 1 - HSX, 2 - DIII-D, 3 - Coronal Loop, 4 - Solar Wind, 5 - MHD-like range
    seti = wave // 5
    # Interaction number - specific to file with 5 interactions, 5 settings
    intnum = wave % 5
    triplets = f["waves"][seti,intnum,:]

    kzmin = f["kzmins"][seti]
    kxmin = f["kpmins"][seti]
    kymin = kxmin
    dt = f["timesteps"][seti]
    nu = f["viscs"][seti]
    eta = f["etas"][seti]
    fname = f["fnames"][seti].astype("T") # Convert from numpy bytes object to string

iterations = 100000

lpath = "/pscratch/sd/e/echansen/threewaves061726/"+fname+str(N)+"set"+str(intnum)+"/"
print(lpath,nu,eta,kzmin,kxmin)
if not os.path.exists(lpath):
    os.makedirs(lpath)
test_iterations = iterations//10000
triplet = triplets.tolist()

#    # Adjust viscosities from N = 256 value if needed - place microscale at N/3
#    nu *= 1/((N/256 * np.sqrt(2)/3)**(4/3))
#    eta *= 1/((N/256 * np.sqrt(2)/3)**(4/3))

# There was a problem in the spreadsheet producing time steps since the length/k scales were set wrong

def alpha(k):
    return((k+np.sqrt(4+k**2))/2)

dt = (2/3)/(nkz0 * kzmin * alpha(2*nkx0*kxmin))

runnames = ["vdvhalloff","vdvoffhallon","mhd","hmhdp1","hmhdp10","hmhdp01"]

lpath = "/pscratch/sd/e/echansen/threewaves061726/"+runnames[start]+str(N)+"/"
print(lpath,nu,eta,kzmin,kxmin)
if not os.path.exists(lpath):
    os.makedirs(lpath)

if start == 0:
    vdv = 0.0
    hp = 0.0
    ene = 1
if start == 1:
    vdv = 0.0
    hp = 1.0
    ene = 1
if start == 2: # MHD
    vdv = 1.0
    hp = 0.0
    ene = 1
else:
    vdv = 1.0
    hp = 1.0
    ene = 2

if start <= 3:
    eta = 1.0
if start == 4:
    eta = 10.0
if start == 5:
    eta = 0.1
    
solver = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                     dt,iterations,lpath,linear=False,
                     initcond="threewave",energystart=0.1,init_kolm=0,hmhdwave=[1,0,0,0],
                     forcetype="threewave",hallparam=hp,vdv=vdv,
                     solveprec=16,maxwallclock=170000,triplet=triplet,records=20,bittype=64,exactnueta=ene)

solver.gauss2split()
