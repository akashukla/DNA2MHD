import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

modes = np.load("paraparsample.npy")
#modes = np.loadtxt("modes128k2.csv",delimiter=",")
print(modes.shape)

Nmodes = modes.shape[0]

mode_values = np.zeros([Nmodes,9])
# perp scale length
mode_values[:,0:2] = modes[:,0:2] * modes[:,6,None]
mode_values[:,6:8] = modes[:,6:8] * modes[:,6,None]
mode_values[:,3:5] = mode_values[:,6:8]-mode_values[:,0:2]

# parallel scale length
mode_values[:,2] = modes[:,2] * modes[:,7]
mode_values[:,8] = modes[:,8] * modes[:,7]
mode_values[:,5] = mode_values[:,8] - mode_values[:,2]

# norms
norms = np.zeros([Nmodes,4])
norms[:,0] =np.sqrt(np.sum(mode_values[:,0:3]**2.0,axis=1))
norms[:,1] =np.sqrt(np.sum(mode_values[:,3:6]**2.0,axis=1))
norms[:,2] =np.sqrt(np.sum(mode_values[:,6:9]**2.0,axis=1))
norms[:,3] = np.sum(norms[:,0:3],axis=1)

# mus
# note hybrid mode is in the middle

mus = np.zeros([Nmodes,3])
mus[:,0] = - (norms[:,1]-norms[:,2])/(2*norms[:,0]) * norms[:,3]
mus[:,2] = - (norms[:,1]-norms[:,0])/(2*norms[:,2]) * norms[:,3]
mus[:,1] = - (norms[:,0]-norms[:,2])/(2*norms[:,1]) * norms[:,3]

# positive curl eigenstates

# shear waves first
shear_wave = np.zeros([Nmodes,9])
shear_wave[:,0::3] = mode_values[:,1::3]/np.sqrt(mode_values[:,0::3]**2 + mode_values[:,1::3]**2)
shear_wave[:,1::3] = -mode_values[:,0::3]/np.sqrt(mode_values[:,0::3]**2 + mode_values[:,1::3]**2)

# pseudo waves from shear waves
pseudo_wave = np.zeros([Nmodes,9])
pseudo_wave[:,0:3] = np.cross(mode_values[:,0:3],shear_wave[:,0:3],axis=1)/norms[:,0,None]
pseudo_wave[:,3:6] = np.cross(mode_values[:,3:6],shear_wave[:,3:6],axis=1)/norms[:,1,None]
pseudo_wave[:,6:9] = np.cross(mode_values[:,6:9],shear_wave[:,6:9],axis=1)/norms[:,2,None]

pcurleig = np.zeros([Nmodes,9],dtype="complex128")
pcurleig = shear_wave + 1j * pseudo_wave

# Interaction Kernels for Wave Times

# J21 = |(k2 dot e(1-2))(e2 dot e1*) = |(k2 dot e(2-1)*)|| (e2 dot e1*)|

J21 = np.abs(np.sum(mode_values[:,6:9]*np.conj(pcurleig[:,3:6]),axis=1))
J21 *= np.abs(np.sum(pcurleig[:,6:9]*np.conj(pcurleig[:,0:3]),axis=1))

# J2K = |k2 dot e(2-K)||e2 dot eK*| = |k2 dot e1| |e2 dot eK*|

J2K = np.abs(np.sum(mode_values[:,6:9]*pcurleig[:,0:3],axis=1))
J2K *= np.abs(np.sum(pcurleig[:,6:9]*np.conj(pcurleig[:,3:6]),axis=1))

# Jn1K = |-k1 dot e(-1-K)| |e-1 dot eK*| = |k1 dot e2*| |e1* dot eK*|

Jn1K = np.abs(np.sum(mode_values[:,0:3]*np.conj(pcurleig[:,6:9]),axis=1))
Jn1K *= np.abs(np.sum(np.conj(pcurleig[:,0:3])*np.conj(pcurleig[:,3:6]),axis=1))

# us
# middle mode is hybrid

us = np.zeros([Nmodes,3])
us[:,0] = mus[:,0]/2 * J21
us[:,2] = mus[:,2]/2 * J21
us[:,1] = mus[:,1]/2 * (J2K-Jn1K)

psi1_0 = 4/3.5 * 0.01 # initial energy of mode times 2 from sim
psiK_0 = 2/3.5 * 0.01 # initial energy of hybrid mode from sim times 2

# normalized initial energies
Psi1_0 = psi1_0/np.abs(us[:,0])  
PsiK_0 = psiK_0/np.abs(us[:,1]) 

# Mahajan's constant E
E = Psi1_0 + PsiK_0

print(E,us[:,0],us[:,1],us[:,2])

# Wave Time
T = np.sqrt(E*np.abs(us[:,0]*us[:,1]*us[:,2]))

# Perp angle calculation between waves
dp = np.sum(mode_values[:,0:2]*mode_values[:,6:8],axis=1) / \
    np.sqrt(np.sum(mode_values[:,0:2]**2,axis=1)*np.sum(mode_values[:,6:8]**2,axis=1))
angle = np.arccos(dp)


plt.scatter(angle,1/T,s=50,marker="s",cmap="tab20b",c=norms[:,3])
plt.colorbar()
plt.xlabel("Planar Angle Between Initial and Final Wave")
plt.ylabel("Mahajan 2021 Energy Depletion Time")
plt.yscale("log")
plt.show()

