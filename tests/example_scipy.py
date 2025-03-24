import numpy as np
from scipy.fft import rfftn,irfftn
import time

nx0_big = 192
ny0_big = 192
nz0_big = 192

realarray = np.zeros([nx0_big,ny0_big,nz0_big],dtype="float64")
specarray = np.zeros([nx0_big,ny0_big,nz0_big//2 + 1],dtype="complex128")

t1 = time.perf_counter()

for repeat in range(0,500):

    ii = np.arange(0,nx0_big)
    jj = np.arange(0,ny0_big)
    kk = np.arange(0,nz0_big//2 + 1)

    specarray = ii[:,None,None]/np.sqrt((2*ii[:,None,None]+2)**2.0 + (3+3*jj[None,:,None])**3.0 + (4+4*kk[None,None,:])**4.0)

    specarray = 1j * kk[None,None,:] * specarray - jj[None,:,None]*kk[None,None,:]*specarray
                
    realarray = irfftn(specarray)
    realarray = realarray**2.0

    specarray = rfftn(realarray)
    specarray = 1j * jj[None,:,None] * specarray/(nx0_big*ny0_big*nz0_big) - jj[None,:,None]*ii[:,None,None]*specarray

    if repeat % 30 == 0:
        print("Iteration ",repeat)

t2 = time.perf_counter()

print("Time for Operations ",t2-t1)





                               
