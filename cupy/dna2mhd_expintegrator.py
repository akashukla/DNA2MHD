import numpy as np
import cupy as cp
from cupy import float32
from numba import njit,jit
from mpi4py import MPI
import nvmath.distributed
from numpy import pi


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

def exponentialcoefficients(kxgrid,kygrid,kzgrid,eta,nu,hyp,dt,Nquad):

    padNx = np.size(kxgrid)//2
    padNy = np.size(kygrid)//2
    padNz = np.size(kzgrid)

    kmags = cp.sqrt(kxgrid[:,None,None]**2 + kygrid[:,None,:]**2 + kzgrid[:,:,None]**2,dtype="float32")
    shape4 = [2*padNx,2*padNy,padNz,4]

    # Construct array of dissipative eigenvalues
    eig = cp.zeros(shape4,dtype="complex64")
    eig[:,:,:,0] = - ((nu+eta)*kmags**(2*hyp) + 1j * kmags * kzgrid[None,None,:])/2
    eig[:,:,:,2] = - ((nu+eta)*kmags**(2*hyp) - 1j * kmags * kzgrid[None,None,:])/2
    eig[:,:,:,1] = cp.sqrt(eig[:,:,:,0]**2 - nu * eta * kmags**(4*hyp) - kzgrid[None,None,:]**2 - 1j * kmags * kzgrid[None,None,:] * nu * kmags**(2*hyp))
    eig[:,:,:,3] = cp.sqrt(eig[:,:,:,2]**2 - nu * eta * kmags**(4*hyp) - kzgrid[None,None,:]**2 + 1j * kmags * kzgrid[None,None,:] * nu * kmags**(2*hyp))

    eig[:,:,:,0] += eig[:,:,:,1]
    eig[:,:,:,2] += eig[:,:,:,3]
    eig[:,:,:,1] *= -2
    eig[:,:,:,3] *= -2
    eig[:,:,:,1] += eig[:,:,:,0]
    eig[:,:,:,3] += eig[:,:,:,2]
    
    # Time integration variables

    # Construct diagonal time integration variables from contour quadrature
    expLdiag = cp.exp(eig*dt)
    coef1diag = cp.zeros(shape4,dtype="complex64")
    coef2diag = cp.zeros(shape4,dtype="complex64")

    for q in range(Nquad):
        
        r = np.exp(q/Nquad * 1j * 2*pi)
        
        # Initialize each nontrivial eigenvalue
        coef1diag += (expL*cp.exp(r*dt) - 1)/((eig+r))
        coef2diag += (expL*cp.exp(r+dt)-1-(eig+r)*dt)/((eig+r)**2)

    coef1diag /= Nquad
    coef2diag /= Nquad

    # Get needed values from normal mode decomposition
    # Only using incompressible Hall MHD so not including the other modes

    # Matrix is block 2x2 diagonal in k and curl eigenstates

    shape4full = [2*padNx,2*padNy,padNz,8]

    # Structure results: plus eigenstate bb, bv, vb, vv, minus bb, bv, vb, vv
    
    expL = cp.zeros(shape4full,dtype="complex64")
    coef1 = cp.zeros(shape4full,dtype="complex64")
    coef2 = cp.zeros(shape4full,dtype="complex64")
    
    # Can I vectorize this loop? The main thing that comes to mind is einsum

    for ikx in range(0,2*padNx):
        for iky in range(0,2*padNy):
            for ikz in range(0,2*padNz):
                
                # Positive curl eigenstate modes
                evmat_plus = cp.array([[nu*kmags**(2*hyp)+eig[ikx,iky,ikz,:2],[1j*kzgrid[ikz],1j*kzgrid[ikz]]]],dtype="complex64")
                expk = evmat_plus @ cp.diag(expLdiag[ikx,iky,ikz,:2]) @ cp.linalg.inv(evmat_plus)
                coef1k = evmat_plus @ cp.diag(coef1diag[ikx,iky,ikz,:2]) @ cp.linalg.inv(evmat_plus)
                coef2k = evmat_plus @ cp.diag(coef2diag[ikx,iky,ikz,:2]) @ cp.linalg.inv(evmat_plus)
                expL[ikx,iky,ikz,:4] = expk.flatten().copy()
                coef1[ikx,iky,ikz,:4] = coef1k.flatten().copy()
                coef2[ikx,iky,ikz,:4] = coef2k.flatten().copy()     
                
                # Negative curl eigenstate modes 
                evmat_minus = cp.array([[nu*kmags**(2*hyp)+eig[ikx,iky,ikz,2:],[1j*kzgrid[ikz],1j*kzgrid[ikz]]]],dtype="complex64")
                expk = evmat_minus @ cp.diag(expLdiag[ikx,iky,ikz,2:]) @ cp.linalg.inv(evmat_minus)
                coef1k = evmat_minus @ cp.diag(coef1diag[ikx,iky,ikz,2:]) @ cp.linalg.inv(evmat_minus)
                coef2k = evmat_minus @ cp.diag(coef2diag[ikx,iky,ikz,2:]) @ cp.linalg.inv(evmat_minus)
                expL[ikx,iky,ikz,4:] = expk.flatten().copy()
                coef1[ikx,iky,ikz,4:] = coef1k.flatten().copy()
                coef2[ikx,iky,ikz,4:] = coef2k.flatten().copy()

    return(expL,coef1,coef2)
