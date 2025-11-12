import numpy as np
import h5py
import matplotlib.pyplot as plt
lpath = "/pscratch/sd/e/echansen/interactionangle51"

with h5py.File(lpath+"/output.hdf5","r") as f:

    ['crosshel', 'hamiltonian', 'helicitycorr', 'kinenergy', 'leftcyclo', 'leftwhistler', 'magenergy', 'maghel', 'magneticfields', 'rightcyclo', 'rightwhistler', 'time', 'velocityfields', 'written']
    
    print(f.keys())
    
    written = f["written"][0]
    time = f["time"][:]
    ham = f["hamiltonian"][:]
    mhel = f["maghel"][:]
    chel = f["crosshel"][:]
    ke = f["kinenergy"][:]
    me = f["magenergy"][:]
    lw = f["leftwhistler"][:]
    lc = f["leftcyclo"][:]
    rw = f["rightwhistler"][:]
    rc = f["rightcyclo"][:]
    mhc = f["helicitycorr"][:]

    b1 = f["magneticfields"][written,:,:,:,:]
    v1 = f["velocityfields"][written,:,:,:,:]

print(time)
print(ham)
print(mhel)
print(chel)
