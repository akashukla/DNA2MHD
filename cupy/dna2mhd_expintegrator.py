import numpy as np
import cupy as cp
from cupy import float32
from numba import njit,jit
from mpi4py import MPI
import nvmath.distributed
from numpy import pi

"""
dt = float32(0.1)
Nquad = 32

Nx = 128
Ny = 128
Nz = 128
padNx = Nx*3//2
padNy = Ny*3//2
padNz = Nz*3//2

Lx = 100
Ly = 100
Lz = 1000
nu = float32(0.001)

kxgrid = cp.arange(0,padNx,dtype="float32")*2*pi/Lx
kxgrid = cp.hstack((kxgrid,-kxgrid[:0:-1]))
kygrid = cp.arange(0,padNy,dtype="float32")*2*pi/Ly
kygrid = cp.hstack((kygrid,-kygrid[:0:-1]))
kzgrid = cp.arange(0,padNz,dtype="float32")*2*pi/Lz

nu = 0.0001
eta = nu
hyp = 2
"""

def exponentialcoefficients(kxgrid,kygrid,kzgrid,eta,nu,hyp,dt,Nquad,kmags,eig):

    padNx = cp.size(kxgrid)//2
    padNy = cp.size(kygrid)//2
    padNz = cp.size(kzgrid)

    shape4 = [2*padNx,2*padNy,padNz,4]

    # Time integration variables

    # Construct diagonal time integration variables from contour quadrature
    expLdiag = cp.exp(eig*dt)
    coef1diag = cp.zeros(shape4,dtype="complex64")
    coef2diag = cp.zeros(shape4,dtype="complex64")

    # Contour integral for the needed exponential time integrator functions
    for q in range(Nquad):
        
        r = cp.exp((q+0.5)/Nquad * 1j * 2*pi)
        
        # Initialize each nontrivial eigenvalue
        coef1diag += (expLdiag*cp.exp(r*dt) - 1)/((eig+r))
        coef2diag += (expLdiag*cp.exp(r*dt)-1-(eig+r)*dt)/(dt*(eig+r)**2)

    coef1diag /= Nquad
    coef2diag /= Nquad

    # Get needed values from normal mode decomposition
    # Only using incompressible Hall MHD so not including the other modes

    # Matrix is block 2x2 diagonal in k and curl eigenstates

    shape4full = [2*padNx,2*padNy,padNz,8]

    # Structure results: plus eigenstate bb, bv, vb, vv, minus bb, bv, vb, vv

    def changeofbasis(diagcalc):

        # Perform PDPinv calculation for adjusted Hall MHD proportionality factors between b and v from diagonal basis
        
        result = cp.zeros(shape4full,dtype="complex64")

        alpha = kmags/2 + cp.sqrt(1 + (kmags **2)  / 4)
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
                
                alpha **= 1-it/2 # for negative curl eigenstates send alpha to 1/alpha
                
                result[:,:,:,it] += diagcalc[:,:,:,it//2] + alpha**2 * diagcalc[:,:,:,it//2+1]
                result[:,:,:,it+1] += alpha * (diagcalc[:,:,:,it//2]-diagcalc[:,:,:,it//2+1])
                result[:,:,:,it+2] += alpha * (diagcalc[:,:,:,it//2]-diagcalc[:,:,:,it//2+1])
                result[:,:,:,it+3] += alpha**2 * diagcalc[:,:,:,it//2] + diagcalc[:,:,:,it//2 + 1]

                result[:,:,:,it:it+4] /= alpha[:,:,:,None]**2 + 1
                
        return(result)

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
