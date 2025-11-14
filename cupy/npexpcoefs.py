import numpy as np
#import cupy as cp
#from cupy import float32
#from numba import njit,jit
#from mpi4py import MPI
#import nvmath.distributed
from numpy import pi
import matplotlib.pyplot as plt

kxmin = 0.025
kymin = 0.025
kzmin = 0.005

nx0_big = 192
ny0_big = 192
nz0_big = 192

nu = 0.05
eta = 1
hyper = 2

dt = 0.1

kxgrid = np.hstack((np.arange(0,nx0_big//2,dtype="float32"),np.arange(-nx0_big//2,0,dtype="float32")))*kxmin
kygrid = np.hstack((np.arange(0,ny0_big//2,dtype="float32"),np.arange(-ny0_big//2,0,dtype="float32")))*kymin
kzgrid = np.arange(0,nz0_big//2+1,dtype="float32")*kzmin

kmags = np.sqrt(kxgrid[:,None,None]**2.0 + kygrid[None,:,None]**2.0+kzgrid[None,None,:]**2.0).astype("float32")
kmax = np.amax(kmags)

vnu = nu/(kmax**(2*hyper))
etab = eta*vnu
hallparam = 1

print("Dissipation Factors",vnu,etab)

eig = np.zeros([nx0_big,ny0_big,nz0_big//2+1,4],dtype="complex64")
eig[:,:,:,0] = - ((vnu+etab)*kmags**(2*hyper) + 1j * kmags * kzgrid[None,None,:]*hallparam)/2
eig[:,:,:,2] = - ((vnu+etab)*kmags**(2*hyper) - 1j * kmags * kzgrid[None,None,:]*hallparam)/2
eig[:,:,:,1] = np.sqrt(eig[:,:,:,0]**2 - vnu * etab * kmags**(4*hyper) - kzgrid[None,None,:]**2 - 1j * kmags * kzgrid[None,None,:] * vnu * kmags**(2*hyper))
eig[:,:,:,3] = np.sqrt(eig[:,:,:,2]**2 - vnu * etab * kmags**(4*hyper) - kzgrid[None,None,:]**2 + 1j * kmags * kzgrid[None,None,:] * vnu * kmags**(2*hyper))

eig[:,:,:,0] += eig[:,:,:,1]
eig[:,:,:,2] += eig[:,:,:,3]
eig[:,:,:,1] *= -2.0
eig[:,:,:,3] *= -2.0
eig[:,:,:,1] += eig[:,:,:,0]
eig[:,:,:,3] += eig[:,:,:,2]


r = np.abs(eig).flatten()
ii = np.argsort(r)

x = np.abs(np.real(eig)).flatten()
y = np.abs(np.imag(eig)).flatten()

plt.scatter(dt*x[ii[::20]],dt*y[ii[::20]],c=dt*r[ii[::20]],cmap="berlin",s=0.01)
plt.xscale("log")
plt.yscale("log")
plt.show()

def exponentialcoefficients(kxgrid,kygrid,kzgrid,eta,nu,hyp,dt,Nquad,kmags,eig,hall,taylor=False,pade=False):

    padNx = np.size(kxgrid)//2
    padNy = np.size(kygrid)//2
    padNz = np.size(kzgrid)

    shape4 = [2*padNx,2*padNy,padNz,4]

    # Time integration variables

    # Construct diagonal time integration variables from contour quadrature
    expLdiag = np.exp(eig*dt)
    coef1diag = np.zeros(shape4,dtype="complex64")
    coef2diag = np.zeros(shape4,dtype="complex64")

    # Contour integral for the needed exponential time integrator functions
    if not taylor and not pade:
        for q in range(Nquad):
        
            r = np.exp((q+0.5)/Nquad * 1j * 2*pi)
        
            # Initialize each nontrivial eigenvalue
            coef1diag += (expLdiag*np.exp(r*dt) - 1)/((eig+r))
            coef2diag += (expLdiag*np.exp(r*dt)-1-(eig+r)*dt)/(dt*(eig+r)**2)

        coef1diag /= Nquad
        coef2diag /= Nquad

    if taylor:
        # Ten term Taylor
        coef1diag += dt
        coef1diag += dt**2 * eig / 2
        coef1diag += dt**3 * eig**2 / 6
        coef1diag += dt**4 * eig**3 / 24
        coef1diag += dt**5 * eig**4 / 120
        coef1diag += dt**6 * eig**5 / 720
        coef1diag += dt**7 * eig**6 / 5040
        coef1diag += dt**8 * eig**7 / 40320
        coef1diag += dt**9 * eig**8 / 362880
        coef1diag += dt**10 * eig**9 / 3628800
        coef1diag += dt**11 * eig**10 / 39916800

        coef2diag += dt/2
        coef2diag += dt**2 * eig / 6
        coef2diag += dt**3 * eig**2 / 24
        coef2diag += dt**4 * eig**3 / 120
        coef2diag += dt**5 * eig**4 / 720
        coef2diag += dt**6 * eig**5 / 5040
        coef2diag += dt**7 * eig**6 / 40320
        coef2diag += dt**8 * eig**7 / 362880
        coef2diag += dt**9 * eig**8 / 3628800
        coef2diag += dt**10 * eig**9 / 39916800
        coef2diag += dt**11 * eig**10 / 479001600

    # Get needed values from normal mode decomposition
    # Only using incompressible Hall MHD so not including the other modes

    if pade:
        coef1diag += (1 + (dt*eig)/14 + (dt*eig)**2 / 42 + (dt*eig)**3 / 840)
        coef1diag /= (1 - 3/7 * (dt*eig) + (dt*eig)**2 / 14 - (dt*eig)**3 / 210)
        coef1diag *= dt

        coef2diag += (1/2 - (dt*eig)/48 + (dt*eig)**2 / 168 + (dt*eig)**3 / 6720)
        coef2diag /= (1 - 3/8 * (dt * eig) + 3/56 * (dt*eig)**2 - 1/336 * (dt*eig)**3)
        coef2diag *= dt
        

    # Matrix is block 2x2 diagonal in k and curl eigenstates

    shape4full = [2*padNx,2*padNy,padNz,8]

    def changeofbasis(diagcalc):

        # Perform PDPinv calculation for adjusted Hall MHD proportionality factors between b and v from diagonal basis
        
        result = np.zeros(shape4full,dtype="complex64")

        alpha = hall*kmags/2 + np.sqrt(1 + (hall*kmags)**2  / 4)
        for it in [0,4]:
        
            if nu != eta:

                result[:,:,:,it] += (nu * kmags**(2*hyp) + eig[:,:,:,it//2]) * 1j * kzgrid[None,None,:] * diagcalc[:,:,:,it//2]
                result[:,:,:,it] -= (nu * kmags**(2*hyp) + eig[:,:,:,it//2 + 1]) * 1j * kzgrid[None,None,:] * diagcalc[:,:,:,it//2 + 1]
                
                result[:,:,:,it+1] -= (nu * kmags**(2*hyp) + eig[:,:,:,it//2])*(nu * kmags**(2*hyp) + eig[:,:,:,it//2 + 1])*diagcalc[:,:,:,it//2]
                result[:,:,:,it+1] += (nu * kmags**(2*hyp) + eig[:,:,:,it//2])*(nu * kmags**(2*hyp) + eig[:,:,:,it//2 + 1])*diagcalc[:,:,:,it//2 + 1]

                result[:,:,:,it+2] -= kzgrid[None,None,:]**2.0 * diagcalc[:,:,:,it//2]
                result[:,:,:,it+2] += kzgrid[None,None,:]**2.0 * diagcalc[:,:,:,it//2 + 1]

                result[:,:,:,it+3] -= (nu * kmags**(2*hyp) + eig[:,:,:,it//2+1]) * 1j * kzgrid[None,None,:] * diagcalc[:,:,:,it//2]
                result[:,:,:,it+3] += (nu * kmags**(2*hyp) + eig[:,:,:,it//2]) * 1j * kzgrid[None,None,:] * diagcalc[:,:,:,it//2 + 1]

                result[:,:,:,it:it+4] /= (eig[:,:,:,it//2] - eig[:,:,:,it//2+1])[:,:,:,None] * 1j * kzgrid[None,None,:,None]

            else:

                alpha **= 1-it/2 
                result[:,:,:,it] += diagcalc[:,:,:,it//2] + alpha**2 * diagcalc[:,:,:,it//2+1]
                result[:,:,:,it+1] += alpha * (diagcalc[:,:,:,it//2]-diagcalc[:,:,:,it//2+1])
                result[:,:,:,it+2] += alpha * (diagcalc[:,:,:,it//2]-diagcalc[:,:,:,it//2+1])
                result[:,:,:,it+3] += alpha**2 * diagcalc[:,:,:,it//2] + diagcalc[:,:,:,it//2 + 1]
                
                result[:,:,:,it:it+4] /= alpha[:,:,:,None]**2 + 1

        return(result)

    # Structure results: plus eigenstate bb, bv, vb, vv, minus bb, bv, vb, vv

    expL = changeofbasis(expLdiag)
    coef1 = changeofbasis(coef1diag)
    coef2 = changeofbasis(coef2diag)

    return(expL,coef1,coef2)
        
    # Can I vectorize this loop? The main thing that comes to mind is einsum
"""
    for ikx in range(0,2*padNx):
        for iky in range(0,2*padNy):
            print(iky)
            for ikz in range(0,padNz):
                
                # Positive curl eigenstate modes
                evmat_plus = cp.array([[nu*kmags[ikx,iky,ikz]**(2*hyp)+eig[ikx,iky,ikz,:2],[1j*kzgrid[ikz],1j*kzgrid[ikz]]]],dtype="complex64")
                expk = evmat_plus @ cp.diag(expLdiag[ikx,iky,ikz,:2]) @ cp.linalg.inv(evmat_plus)
                coef1k = evmat_plus @ cp.diag(coef1diag[ikx,iky,ikz,:2]) @ cp.linalg.inv(evmat_plus)
                coef2k = evmat_plus @ cp.diag(coef2diag[ikx,iky,ikz,:2]) @ cp.linalg.inv(evmat_plus)
                expL[ikx,iky,ikz,:4] = expk.flatten().copy()
                coef1[ikx,iky,ikz,:4] = coef1k.flatten().copy()
                coef2[ikx,iky,ikz,:4] = coef2k.flatten().copy()     
                
                # Negative curl eigenstate modes 
                evmat_minus = cp.array([[nu*kmags[ikx,iky,ikz]**(2*hyp)+eig[ikx,iky,ikz,2:],[1j*kzgrid[ikz],1j*kzgrid[ikz]]]],dtype="complex64")
                expk = evmat_minus @ cp.diag(expLdiag[ikx,iky,ikz,2:]) @ cp.linalg.inv(evmat_minus)
                coef1k = evmat_minus @ cp.diag(coef1diag[ikx,iky,ikz,2:]) @ cp.linalg.inv(evmat_minus)
                coef2k = evmat_minus @ cp.diag(coef2diag[ikx,iky,ikz,2:]) @ cp.linalg.inv(evmat_minus)
                expL[ikx,iky,ikz,4:] = expk.flatten().copy()
                coef1[ikx,iky,ikz,4:] = coef1k.flatten().copy()
                coef2[ikx,iky,ikz,4:] = coef2k.flatten().copy()
    return(expL,coef1,coef2)
"""
a1,b1,c1 = exponentialcoefficients(kxgrid,kygrid,kzgrid,etab,vnu,hyper,dt,128,kmags,eig,hallparam)
a2,b2,c2 = exponentialcoefficients(kxgrid,kygrid,kzgrid,etab,vnu,hyper,dt,128,kmags,eig,hallparam,taylor=True)
a3,b3,c3 = exponentialcoefficients(kxgrid,kygrid,kzgrid,etab,vnu,hyper,dt,128,kmags,eig,hallparam,pade=True)

print("Diff Coef1 ",np.amax(b2-b1))
print("Diff Coef2 ",np.amax(c2-c1))

iii = np.argmax(c2-c1)

print(c2.flatten()[iii],c1.flatten()[iii],c3.flatten()[iii])
print(x[iii],y[iii])
