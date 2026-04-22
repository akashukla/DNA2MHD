import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py

if True:    

    with h5py.File("initconds032426.hdf5","r") as f:
        triplets = f["waves"][:,:,:]
        mode_values = np.reshape(triplets,shape=[25,9])
else:
    modes = np.load("aps2025/norm35sample.npy")
    ii = [0,2,4,5,6,7,9,10,11,12,13,14,15,16,17,19,20]

    modes = modes[ii,:]
    mode_values[:,0:9] = modes[:,0:9]
    #modes = np.loadtxt("modes128k2.csv",delimiter=",")

# norms
Nmodes = mode_values.shape[0]
norms = np.zeros([Nmodes,3])
norms[:,0] =np.sqrt(np.sum(mode_values[:,0:3]**2.0,axis=1))
norms[:,1] =np.sqrt(np.sum(mode_values[:,3:6]**2.0,axis=1))
norms[:,2] =np.sqrt(np.sum(mode_values[:,6:9]**2.0,axis=1))

alphas = norms/2 + np.sqrt((norms/2)**2 + 1)

# Determine which mode is the sum of the other two

s0 = np.sum(mode_values[:,0:3]-mode_values[:,3:6]-mode_values[:,6:9],axis=-1)
s1 = np.sum(mode_values[:,3:6]-mode_values[:,6:9]-mode_values[:,0:3],axis=-1)
s2 = np.sum(mode_values[:,6:9]-mode_values[:,0:3]-mode_values[:,3:6],axis=-1)

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

mus = np.zeros([Nmodes,3])
Js = np.zeros([Nmodes,3])
us = np.zeros([Nmodes,3])
# alpha * g
mus += alphas * (1/np.roll(alphas,1,axis=-1)-np.roll(norms,1,axis=-1) - (1/np.roll(alphas,2,axis=-1)-np.roll(norms,2,axis=-1)))
# add h and divide by sqrt(1+m^2 / 4)
mus += (-np.roll(norms,2,axis=-1)+np.roll(norms,1,axis=-1))/norms * (1-1/(np.roll(alphas,1,axis=-1)*np.roll(alphas,2,axis=-1)))
mus /= np.sqrt(1 + norms**2 / 4)


def magintkernel(u,ediff,eu,ev):
    """Compute magnitude of interaction kernel for numpy three vectors u, v """
    """ abs(u dot e_(v-u)) abs(e_u dot e_v^*)"""

    return(np.abs(np.sum(u*ediff)*np.sum(eu*np.conj(ev))))
    
# Have to use order that the modes are sums to calculate J
for i in range(Nmodes):
    if s0[i] == 0:
        sumvar1 = 0
        sumvar2 = 1
        sumvarh = 2
    elif s1[i] == 0:
        sumvar1 = 1
        sumvar2 = 2
        sumvarh = 0
    else:
        sumvar1 = 2
        sumvar2 = 0
        sumvarh = 1

    m = mode_values[i,sumvar1*3:sumvar1*3+3]
    em = pcurleig[i,sumvar1*3:sumvar1*3+3]
    n = mode_values[i,sumvar2*3:sumvar2*3+3]
    en = pcurleig[i,sumvar2*3:sumvar2*3+3]
    h = mode_values[i,sumvarh*3:sumvarh*3+3]
    eh = pcurleig[i,sumvarh*3:sumvarh*3+3]

    Js[:,sumvar1] = magintkernel(n,eh,en,em)
    Js[:,sumvar2] = magintkernel(-h,em,-np.conj(eh),en)
    Js[:,sumvarh] = magintkernel(m,-np.conj(en),em,eh)

Js -= np.roll(Js,1,axis=-1)

us = mus * Js/2

psi1_0 = (2/3.5)**2 * 0.01 # initial energy of mode times 2 from sim
psiK_0 = (1/3.5)**2 * 0.01 # initial energy of hybrid mode from sim times 2

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
angle = np.acos(dp)*180/np.pi

vtime = [1200,1800,1400,1600,1600,1800,1400,1400,1600,1600,1200,1200,1800,1200,3600,1400,1200]
tmtotime = [1400,3000,2600,2400,1600,2400,2000,1600,3800,2000,2200,2000,2000,3400,5600,2600,6200]

vtime2 = [600,600,600,600,600,900,600,900,600,900,600,600,600,600,600,600,600]
tmto2 = [1200,4800,2100,4800,1200,1200,6000,1500,6900,1200,1200,900,1200,8600,4500,2000,8900]

plt.plot(angle,1/T,"k8",label="Jacobi")
plt.plot(angle,vtime2,"r8",label="Vorticity")
plt.plot(angle,tmto2,"b8",label="Two Mode Drop")
#plt.scatter(angle,1/T,s=50,marker="s",cmap="tab20b",c=norms[:,3],label="Jacobi")
#plt.scatter(angle,vtime,s=50,marker="^",cmap="tab20b",c=norms[:,3],label="Vorticity")
#plt.scatter(angle,tmtotime,s=50,marker="8",cmap="tab20b",c=norms[:,3],label="Two Mode Drop")
#plt.colorbar()
plt.xlabel("Planar Angle Between Initial and Final Wave (Degrees)",size="large")
plt.ylabel("Estimated Time for Wave Interaction",size="large")
plt.title("Turbulent Time Scale Estimates",size="large")
plt.yscale("log")
plt.legend()
plt.savefig("timescale")
plt.show()

vtime = np.array(vtime)
tmtotime = np.array(tmtotime)

print(np.corrcoef(-np.log10(T),np.log10(vtime2)))
print(np.corrcoef(-np.log10(T),np.log10(tmto2)))
print(np.corrcoef(-np.log10(tmtotime),np.log10(tmto2)))


