"""
Script to order one or more three wave interaction simulations from submitruns.py
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

runnames = ["vdvhalloff","vdvoffhallon","mhd","hmhdp1","hmhdp10","hmhdp01"]
lpath = "/pscratch/sd/e/echansen/threewaves061726/"+runnames[start]+str(N)+"/"

with h5py.File(lpath+"output.hdf5","r") as f:

    dt = f["dt"][0]
    nu = f["nu"][0]
    eta = f["eta"][0]
    hyper = f["hyper"][0]
    hp = f["hall"][0]
    vdv = f["vdv"][0]
    triplets = f["threewaves"][:]
    kxmin = f["kx"][1]
    kymin = f["ky"][1]
    kzmin = f["kz"][1]

iterations = 100000

triplet = triplets.tolist()

solver = DNA2MHD(nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                     dt,iterations,lpath,linear=False,
                     initcond="checkpoint",energystart=0.1,init_kolm=0,hmhdwave=[1,0,0,0],
                     forcetype="threewave",hallparam=hp,vdv=vdv,
                     solveprec=16,maxwallclock=170000,triplet=triplet,records=20,bittype=64,exactnueta=0)

solver.gauss2split()
