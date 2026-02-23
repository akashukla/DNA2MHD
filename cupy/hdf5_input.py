import numpy as np
import h5py
import matplotlib.pyplot as plt
import os


lpath = "/pscratch/sd/e/echansen/interactionangle51"
with h5py.File(lpath+"/output.hdf5","r") as f:

    print(f.keys())
    
    written = f["written"][0]
    time = f["time"][:]
    enval = f["enval"][:,:]

    b1 = f["magneticfields"][written,:,:,:,:]
    v1 = f["velocityfields"][written,:,:,:,:]


if not os.path.exists(lpath+"/eplots/"):
    os.mkdir(lpath+"/eplots/")

x = """
# Energy
fig,ax = plt.subplots(1)
ax.plot(time,ham/(4*np.pi**3),color="black",label="Total")
ax.plot(time,ke/(4*np.pi**3),color="tomato",label="Kinetic")
ax.plot(time,me/(4*np.pi**3),color="turquoise",label="Magnetic")
ax.set_ylabel("Energy (Guide Field Energy)",size="large")
ax.set_xlabel("Time ($\omega_c^{-1}$)",size="large")
ax.legend()
fig.suptitle("Energy Evolution")
plt.savefig(lpath+"eplots/energy.png",bbox_inches="tight")
"""

print(time)
print(enval[:,0])
print(enval[:,1]+enval[:,-1])
print(enval[:,2]+enval[:,-1])
