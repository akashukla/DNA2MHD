import numpy as np
from numpy import format_float_positional as ff
import matplotlib.pyplot as plt

perpscale = 0.05
parscale = 0.01

def resonancecondition(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2):

    k = np.sqrt(perpscale**2.0 * (kx**2 + ky**2) + parscale**2.0 * kz**2 )
    wk = h1 * parscale * kz * (k / 2 + b1 * np.sqrt(1+(k/2)**2.0))
    l = np.sqrt(perpscale**2.0 * (lx**2 + ly**2) + parscale**2.0 * lz**2 )
    wl = h1 * parscale * lz * (l / 2 + b1 * np.sqrt(1+(l/2)**2.0))

    m = np.sqrt(perpscale**2.0 * ((kx+lx)**2 + (ky+ly)**2) + parscale**2.0 * (kz+lz)**2)
    wm = h2 * parscale * (kz+lz) * (m / 2 + b2 * np.sqrt(1+(m/2)**2.0))

    if min(np.abs(wk),np.abs(wl),np.abs(wm)) == 0:
        print(kx,ky,kz,lx,ly,lz,b1,b2,h1,h2)

    #if np.abs(wm - wk -wl) < 0.005:
    #    print(wk,wl,wm,wk+wl-wm,lx,ly,lz)
    return(np.abs(wm-wl-wk)/min(np.abs(wk),np.abs(wl),np.abs(wm)),1.0/np.abs(wm))

test_ks = [[0,35,9],[0,15,9],[0,10,15],[0,47,3],[0,40,17],[0,15,19]]
ind = 0
kx = perpscale*test_ks[ind][0]
ky = perpscale*test_ks[ind][1]
kz = parscale*test_ks[ind][2]

Nx = 64
Ny = 64
Nz = 64

lx = np.arange(0,Nx)*perpscale
ly = np.hstack((np.arange(0,Ny),np.arange(1-Ny,0)))*perpscale
lz = np.hstack((np.arange(0,Nz),np.arange(1-Nz,0)))*parscale

LX,LY = np.meshgrid(lx,ly,indexing="ij")

k = np.sqrt(kx**2.0+ky**2.0+kz**2.0)
l = np.sqrt(lz[None,None,:]**2.0 + lx[:,None,None]**2.0 + ly[None,:,None]**2.0)
m = np.sqrt((lz[None,None,:]+kz)**2.0 + (lx[:,None,None]+kx)**2.0 + (ly[None,:,None]+ky)**2.0)

h1 = 1
h2 = 1
b1 = 1
b2 = 1

wk = h1*(k/2 + b1 * np.sqrt(1+(k/2)**2.0))*kz
wl = h1*(l/2 + b1 * np.sqrt(1+(l/2)**2.0))*lz[None,None,:]
wm = h2*(m/2 + b2 * np.sqrt(1+(m/2)**2.0))*(kz+lz[None,None,:])

freqdiff = np.abs(wm-wl-wk)
scaledfreqdiff = np.abs(wm-wl-wk)/np.minimum(np.abs(wk),np.abs(wl),np.abs(wm))

plt.contourf(LX,LY,np.log10(np.amin(freqdiff,axis=-1)))
plt.colorbar()
plt.show()
plt.close()
