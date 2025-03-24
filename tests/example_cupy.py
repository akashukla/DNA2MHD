import cupy as cp
from cupy.fft import rfftn,irfftn
import time

nkx0 = 128

nx0_big = 3*nkx0//2
ny0_big = 3*nkx0//2
nz0_big = 3*nkx0//2

realarray = cp.zeros([nx0_big,ny0_big,nz0_big],dtype="float64")
specarray = cp.zeros([nx0_big,ny0_big,nz0_big//2 + 1],dtype="complex128")

t1 = time.perf_counter()

for repeat in range(0,5000):

    ii = cp.arange(0,nx0_big)
    jj = cp.arange(0,ny0_big)
    kk = cp.arange(0,nz0_big//2 + 1)

    specarray = ii[:,None,None]/cp.sqrt((2*ii[:,None,None]+2)**2.0 + (3+3*jj[None,:,None])**3.0 + (4+4*kk[None,None,:])**4.0)

    specarray = 1j * kk[None,None,:] * specarray - jj[None,:,None]*kk[None,None,:]*specarray
                
    realarray = irfftn(specarray)
    realarray = realarray**2.0

    specarray = rfftn(realarray)
    specarray = 1j * jj[None,:,None] * specarray/(nx0_big*ny0_big*nz0_big) - jj[None,:,None]*ii[:,None,None]*specarray

    if repeat % 300 == 0:
        print("Iteration ",repeat)

t2 = time.perf_counter()

print("Time for Operations ",t2-t1)





                               
