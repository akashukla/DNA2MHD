import cupy as cp
from cupy.fft import rfftn,irfftn
import time
#from dna2mhd_coefs import exponentialcoefficients
import os
from shutil import rmtree
#from cupy import float32,complex64,int32,float64,complex128
from cupy import int32
import h5py
import numpy as np
from numpy import float32 as float32c
from numpy import complex64 as complex64c
from numpy import int32 as int32c

class RHS:

    def idealhallrhs(self,binput,vinput):

        # Evolution of incompressible ideal Hall MHD system, allowing for a splitting to handle dissipation for current implicit schemes

        bout,vout = self.convolutions(binput,vinput)
        bout,vout = self.removediv(bout,vout)

        return(bout,vout)
    
    def hallrhs(self,binput,vinput):

        # Evolution of incompressible Hall MHD system including dissipation
        bout,vout = self.rhslinear(binput,vinput)
        bout,vout = self.convolutions(binput,vinput,bout,vout)

        bout -= self.etab * self.kmags**(2*self.hyper) * binput
        vout -= self.vnu * self.kmags**(2*self.hyper) * vinput

        mhc = self.helicitycorrection(binput,vinput)

        self.removediv(bout,vout)
        
        return(bout,vout,mhc)

    def rhslinear(self,binput,vinput):

        # For reference, linear evolution of the Hall MHD system prescribed from Mahajan 2021, etc.

        bout = cp.zeros_like(self.b1)
        vout = cp.zeros_like(self.v1)
        
        bout[0,:,:,:] = (1j * vinput[0,:,:,:] + self.kygrid[None,:,None]*binput[2,:,:,:] - self.kzgrid[None,None,:]*binput[1,:,:,:])
        bout[1,:,:,:] = (1j * vinput[1,:,:,:] + self.kzgrid[None,None,:]*binput[0,:,:,:] - self.kxgrid[:,None,None]*binput[2,:,:,:])
        bout[2,:,:,:] = (1j * vinput[2,:,:,:] + self.kxgrid[:,None,None]*binput[1,:,:,:] - self.kygrid[None,:,None]*binput[0,:,:,:])
        bout *= self.kzgrid[None,None,None,:]

        vout[0,:,:,:] = 1j * self.kzgrid[None,None,:] * binput[0,:,:,:] - 1j * binput[2,:,:,:] * self.kxgrid[:,None,None]
        vout[1,:,:,:] = 1j * self.kzgrid[None,None,:] * binput[1,:,:,:] - 1j * binput[2,:,:,:] * self.kygrid[None,:,None]
        vout[2,:,:,:] = 1j * self.kzgrid[None,None,:] * binput[2,:,:,:] - 1j * binput[2,:,:,:] * self.kzgrid[None,None,:]

        return(bout,vout)

    def convolutions(self,binput,vinput,bout,vout):

        # Carry out convolution sums for nonlinear evolution using FFTs

        bx = irfftn(binput[0,:,:,:])
        by = irfftn(binput[1,:,:,:])
        bz = irfftn(binput[2,:,:,:])

        vx = irfftn(vinput[0,:,:,:])
        vy = irfftn(vinput[1,:,:,:])
        vz = irfftn(vinput[2,:,:,:])

        # Compute curl operators spectrally before FFTs

        dum = 1j*self.kygrid[None,:,None]*binput[2,:,:,:]-1j*self.kzgrid[None,None,:]*binput[1,:,:,:]
        curlbx = irfftn(dum)

        dum = 1j*self.kygrid[None,:,None]*vinput[2,:,:,:]-1j*self.kzgrid[None,None,:]*vinput[1,:,:,:]
        curlvx = irfftn(dum)

        dum = 1j*self.kzgrid[None,None,:]*binput[0,:,:,:]-1j*self.kxgrid[:,None,None]*binput[2,:,:,:]
        curlby = irfftn(dum)

        dum = 1j*self.kzgrid[None,None,:]*vinput[0,:,:,:]-1j*self.kxgrid[:,None,None]*vinput[2,:,:,:]
        curlvy = irfftn(dum)

        dum = 1j*self.kxgrid[:,None,None]*binput[1,:,:,:]-1j*self.kygrid[None,:,None]*binput[0,:,:,:]
        curlbz = irfftn(dum)

        dum = 1j*self.kxgrid[:,None,None]*vinput[1,:,:,:]-1j*self.kygrid[None,:,None]*vinput[0,:,:,:]
        curlvz = irfftn(dum)

        # Compute the nonlinear RHS and transform into real space

        dum = vy * bz - vz * by - self.hallparam * ( curlby * bz - curlbz * by )
        dum2 = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big
        bout[1,:,:,:] += 1j * self.kzgrid[None,None,:] * dum2
        bout[2,:,:,:] += -1j * self.kygrid[None,:,None] * dum2

        dum = vy * curlvz - vz * curlvy + (curlby * bz - curlbz * by)
        vout[0,:,:,:] = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big

        dum = vz * bx - vx * bz - self.hallparam * ( curlbz * bx - curlbx * bz )
        dum2 = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big
        bout[0,:,:,:] += -1j * self.kzgrid[None,None,:] * dum2
        bout[2,:,:,:] += 1j * self.kxgrid[:,None,None] * dum2

        dum = vz * curlvx - vx * curlvz + (curlbz * bx - curlbx * bz)
        vout[1,:,:,:] = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big

        dum = vx * by	- vy * bx - self.hallparam * ( curlbx * by - curlby * bx )
        dum2 = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big
        bout[0,:,:,:] += 1j * self.kygrid[None,:,None] * dum2
        bout[1,:,:,:] += -1j * self.kxgrid[:,None,None] * dum2
        
        dum = vx * curlvy - vy * curlvx + (curlbx * by - curlby * bx)
        vout[2,:,:,:] = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big

        # Convolution sums are dealiased using the 3/2 rule
        
        bout *= self.padding[None,:,:,:]
        vout *= self.padding[None,:,:,:]

        #print("Through NL")

        return(bout,vout)

    def reset_phase(self):

        # Randomly change forcing phase

        return(cp.exp(1j*cp.pi*2*cp.random.rand(self.nx0_big,self.ny0_big,self.nz0_big//2 + 1)).astype(self.ctype))
    
    def force(self):

        if self.forcetype == "hallwave":

            # Inject energy using prescribed type of wave at selected wavenumbers from forcemask, with random phases
            
            phase = self.reset_phase()
            self.b1 += self.dt/2 * self.forcemask[None,:,:,:] * self.forcewave[0] \
                * self.pcurleig * self.alphaleftwhist[None,:,:,:] / cp.sqrt(self.alphaleftwhist[None,:,:,:]**2.0 + 1) * phase[None,:,:,:]
            self.v1 += self.dt/2 * self.forcemask[None,:,:,:] * self.forcewave[0] \
                * self.pcurleig / cp.sqrt(self.alphaleftwhist[None,:,:,:]**2.0 + 1) * phase[None,:,:,:]

            phase = self.reset_phase()
            self.b1 += self.dt/2 * self.forcemask[None,:,:,:] * self.forcewave[2] \
                * self.pcurleig * -self.alphaleftwhist[None,:,:,:] / cp.sqrt(self.alphaleftwhist[None,:,:,:]**2.0 + 1) * phase[None,:,:,:]
            self.v1 += self.dt/2 * self.forcemask[None,:,:,:] * self.forcewave[2] \
                * self.pcurleig / cp.sqrt(self.alphaleftwhist[None,:,:,:]**2.0 + 1) * phase[None,:,:,:]

            phase = self.reset_phase()
            self.b1 += self.dt/2 * self.forcemask[None,:,:,:] * self.forcewave[1] \
		* self.pcurleig / cp.sqrt(self.alphaleftwhist[None,:,:,:]**2.0 + 1) * phase[None,:,:,:]
            self.v1 += self.dt/2 * self.forcemask[None,:,:,:] * self.forcewave[1] \
                * self.pcurleig * -self.alphaleftwhist[None,:,:,:]/ cp.sqrt(self.alphaleftwhist[None,:,:,:]**2.0 + 1) * phase[None,:,:,:]
            
            phase = self.reset_phase()
            self.b1 += self.dt/2 * self.forcemask[None,:,:,:] * self.forcewave[3] \
                * self.pcurleig / cp.sqrt(self.alphaleftwhist[None,:,:,:]**2.0 + 1) * phase[None,:,:,:]
            self.v1 += self.dt/2 * self.forcemask[None,:,:,:] * self.forcewave[3] \
                * self.pcurleig * self.alphaleftwhist[None,:,:,:]/ cp.sqrt(self.alphaleftwhist[None,:,:,:]**2.0 + 1) * phase[None,:,:,:]

        return(None)
        
    def removediv(self,binput,vinput):

        # Maintain incompressibility by solving the Poisson equation for pressure at each time step
        
        divergence_v = self.kxgrid[:,None,None]*vinput[0,:,:,:]+self.kygrid[None,:,None]*vinput[1,:,:,:]+ \
            self.kzgrid[None,None,:]*vinput[2,:,:,:]
        
        avgvelocity = cp.copy(vinput[:,0,0,0])
        vinput[0,:,:,:] -= divergence_v/(self.kmags**2.0) * self.kxgrid[:,None,None]
        vinput[1,:,:,:] -= divergence_v/(self.kmags**2.0) * self.kygrid[None,:,None]
        vinput[2,:,:,:] -= divergence_v/(self.kmags**2.0) * self.kzgrid[None,None,:]
        vinput[:,0,0,0] = avgvelocity

        if cp.any(cp.isnan(vinput)):
            print("Remove Div shows nans ")
            print(cp.nonzero(cp.isnan(vinput)))
            quit()

        return(binput,vinput)

class DIAGS:

    def hamiltonian(self,option):

        ham = 0.0

        # Sums to compute the magnetic and kinetic energies
        # Notice kz = 0 terms only contribute once to the sum since negative kz terms are found from the complex conjugate of nonzero kz
        
        if option == 0 or option == 2 :
            ham += cp.sum(cp.abs(self.b1[:,:,:,0])**2.0,axis=(0,1,2))
            ham += 2*cp.sum(cp.abs(self.b1[:,:,:,1:])**2.0,axis=(0,1,2,3))

        if option == 0 or option == 1 :
            ham += cp.sum(cp.abs(self.v1[:,:,:,0])**2.0,axis=(0,1,2))
            ham += 2*cp.sum(cp.abs(self.v1[:,:,:,1:])**2.0,axis=(0,1,2,3))

        ham *= 4*cp.pi**3

        return(ham)
        
    def helicity(self):

        # Compute vector potential from laplacian inverse curl b, for Coulomb gauge

        AVP = cp.stack((1j*self.kygrid[None,:,None]*self.b1[2,:,:,:]-1j*self.kzgrid[None,None,:]*self.b1[1,:,:,:],
                        1j*self.kzgrid[None,None,:]*self.b1[0,:,:,:]-1j*self.kxgrid[:,None,None]*self.b1[2,:,:,:],
                        1j*self.kxgrid[:,None,None]*self.b1[1,:,:,:]-1j*self.kygrid[None,:,None]*self.b1[0,:,:,:]),\
                       axis=0)
        cb2 = 2*cp.sum(cp.abs(AVP[:,:,:,1:])**2.0)+cp.sum(cp.abs(AVP[:,:,:,0])**2.0)
        AVP /= self.kmags[None,:,:,:]**2.0
        AVP[:,0,0,0] = 0.0

        # Compute vorticity
        WVORT = cp.stack((1j*self.kygrid[None,:,None]*self.v1[2,:,:,:]-1j*self.kzgrid[None,None,:]*self.v1[1,:,:,:],
                          1j*self.kzgrid[None,None,:]*self.v1[0,:,:,:]-1j*self.kxgrid[:,None,None]*self.v1[2,:,:,:],
                          1j*self.kxgrid[:,None,None]*self.v1[1,:,:,:]-1j*self.kygrid[None,:,None]*self.v1[0,:,:,:]),\
                         axis=0)
        vort = 2*cp.sum(cp.abs(WVORT[:,:,:,1:])**2.0)+cp.sum(cp.abs(WVORT[:,:,:,0])**2.0)

        # Sums to compute the magnetic helicity and canonical helicity from the vector potential and vorticity
        # Notice that kz = 0 terms only contribute once to the sum

        mhel = cp.real(2*cp.sum((AVP[:,:,:,1:])*cp.conj(self.b1[:,:,:,1:])))
        mhel += cp.real(cp.sum(AVP[:,:,:,0]*cp.conj(self.b1[:,:,:,0])))

        canhel = cp.real(2*cp.sum((AVP[:,:,:,1:]+self.hallparam*self.v1[:,:,:,1:])*cp.conj(self.b1[:,:,:,1:]+self.hallparam*WVORT[:,:,:,1:])))
        canhel += cp.real(cp.sum((AVP[:,:,:,0]+self.hallparam*self.v1[:,:,:,0])*cp.conj(self.b1[:,:,:,0]+self.hallparam*WVORT[:,:,:,0]))) + cp.real(self.hallparam*self.v1[2,0,0,0])

        mhel *= 8 * cp.pi**3
        canhel *= 8 * cp.pi**3
        
        return(mhel,canhel,vort,cb2)

    def tripletamps(self,lwk):

        # Find three wave amplitudes from individual wave energies found in hallmodeenergies()

        lw1 = 8.0*cp.pi**3 * cp.abs(lwk[self.itriplet[0],self.itriplet[1],self.itriplet[2]])**2.0
        lw2 = 8.0*cp.pi**3 * cp.abs(lwk[self.itriplet[3],self.itriplet[4],self.itriplet[5]])**2.0
        lw3 = 8.0*cp.pi**3 * cp.abs(lwk[self.itriplet[6],self.itriplet[7],self.itriplet[8]])**2.0

        return(lw1,lw2,lw3)

    def hallmodeenergies(self):

        # Energies in each type of inviscid normal mode for kz != 0
        
        self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
        self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
        self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
        self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

        # Compute normal mode coefficients by a dot product with the orthogonal normal mode basis
        # This computation must change for non-ideal Hall MHD
        
        lwk = (self.alphaleftwhist * self.bplus + self.vplus)/cp.sqrt(self.alphaleftwhist**2 + 1)
        lw = 2.0*cp.sum(cp.abs(lwk[:,:,1:])**2.0)
        lw *= (4.0*cp.pi**3)
        
        if self.initcond == "threewave":
            lw1,lw2,lw3 = self.tripletamps(lwk)

        lwk = (-1/self.alphaleftwhist * self.bplus + self.vplus)/cp.sqrt((1/self.alphaleftwhist)**2 + 1)
        lc = 2.0 * cp.sum(cp.abs(lwk[:,:,1:])**2.0)
        lc *= 4.0 * cp.pi**3
        if self.initcond == "threewave":
            lc1,lc2,lc3 = self.tripletamps(lwk)

        lwk = (-self.alphaleftwhist * self.bminus + self.vminus)/cp.sqrt(self.alphaleftwhist**2 + 1)
        rw = 2.0 * cp.sum(cp.abs(lwk[:,:,1:])**2.0)
        rw *= 4.0 * cp.pi**3
        if self.initcond == "threewave":
            rw1,rw2,rw3 = self.tripletamps(lwk)

        lwk = (1/self.alphaleftwhist * self.bminus + self.vminus)/cp.sqrt((1/self.alphaleftwhist)**2 + 1)
        rc = 2.0 * cp.sum(cp.abs(lwk[:,:,1:])**2.0)
        rc *= 4.0 * cp.pi**3
        if self.initcond == "threewave":
            rc1,rc2,rc3 = self.tripletamps(lwk)

        # Prepare triplet wave energies if a three wave simulation
        if self.initcond == "threewave":
            
            self.tripletenergies = cp.array([lw1,lc1,rw1,rc1,lw2,lc2,rw2,rc2,lw3,lc3,rw3,rc3])

        return(lw,lc,rw,rc)

    def helicitycorrection(self,binput,vinput):

        # Compute the correction to the magnetic helicity needed for computation, as in Stribling 1994
        correction = cp.real(binput[0,:,:,:]*cp.conj(vinput[1,:,:,:])-binput[1,:,:,:]*cp.conj(vinput[0,:,:,:]))

        mhc = 2.0*cp.sum(correction[:,:,1:])+cp.sum(correction[:,:,0])
        mhc *= -16*cp.pi**3
        
        return(mhc)

    def createhdf5file(self):

        # Define a new HDF5 simulation data file

        # Remove past data if needed 
        if os.path.exists(self.lpath+"/output.hdf5") and (self.itime == 0):

            rmtree(self.lpath)
            os.mkdir(self.lpath)

        # Send grid and hall parameter data to CPU
        
        kxgridcpu = cp.asnumpy(self.kxgrid)
        kygridcpu = cp.asnumpy(self.kygrid)
        kzgridcpu = cp.asnumpy(self.kzgrid)
        hcpu = cp.asnumpy(self.hallparam)
        dtcpu = cp.asnumpy(self.dt)
        
        if not os.path.exists(self.lpath+"/output.hdf5"):

            # Define datasets for new HDF5 file, with identifying simulation parameters and grids
            with h5py.File(self.lpath+"/output.hdf5","a") as f:

                kx = f.create_dataset("kx",(np.size(kxgridcpu),),data=kxgridcpu)
                ky = f.create_dataset("ky",(np.size(kygridcpu),),data=kygridcpu)
                kz = f.create_dataset("kz",(np.size(kzgridcpu),),data=kzgridcpu)
                written = f.create_dataset("written",(1,),dtype="int32",data=int32c(0))
                dt = f.create_dataset("dt",(1,),dtype="float64",data=dtcpu)
                hall = f.create_dataset("hall",(1,),data=hcpu)

                itime = f.create_dataset("itime",(1,),dtype="int32",data=0)
                times = f.create_dataset("time",(self.record+1,),dtype=self.rtype,maxshape=(None))
                enval = f.create_dataset("enval",(self.record+1,12),dtype=self.rtype,maxshape=(None,12))

                if self.initcond == "threewave":

                    tripletcpu = cp.asnumpy(self.triplet)
                    threewaves = f.create_dataset("threewaves",(12,),dtype=self.rtype,data=tripletcpu)
                    threewaveenergies = f.create_dataset("threewaveenergies",(self.record+1,12),dtype=self.rtype,maxshape=(None,12))

                magneticfields = f.create_dataset("magneticfields",(self.record+1,3,self.nx0_big,self.ny0_big,self.nz0_big//2+1),dtype=self.ctype,maxshape=(None,3,self.nx0_big,self.ny0_big,self.nz0_big//2+1))
                velocityfields = f.create_dataset("velocityfields",(self.record+1,3,self.nx0_big,self.ny0_big,self.nz0_big//2+1),dtype=self.ctype,maxshape=(None,3,self.nx0_big,self.ny0_big,self.nz0_big//2+1))

        return(None)
    
    def writefile(self):

        # Compute diagnostics
        ham = self.hamiltonian(0)
        kinham = self.hamiltonian(1)
        magham = self.hamiltonian(2)
        mh,ch,vort,cb2 = self.helicity()
        lw,lc,rw,rc = self.hallmodeenergies()


        # Send data to CPU
        self.writtencpu = cp.asnumpy(self.written)
        self.timecpu = cp.asnumpy(self.time)
        hamcpu = cp.asnumpy(ham)
        kinhamcpu = cp.asnumpy(kinham)
        maghamcpu = cp.asnumpy(magham)
        mhcpu = cp.asnumpy(mh)
        chcpu = cp.asnumpy(ch)
        lwcpu = cp.asnumpy(lw)
        lccpu = cp.asnumpy(lc)
        rwcpu = cp.asnumpy(rw)
        rccpu = cp.asnumpy(rc)
        vortcpu = cp.asnumpy(vort)
        cb2cpu = cp.asnumpy(cb2)
        
        b1cpu = cp.asnumpy(self.b1)
        v1cpu = cp.asnumpy(self.v1)
        mhccpu = cp.asnumpy(self.mhelcorr)
        itimecpu = cp.asnumpy(self.itime)

        if self.initcond == "threewave":
            waveenergiescpu = cp.asnumpy(self.tripletenergies)        

        # Use CPU to write HDF5 file 
        with h5py.File(self.lpath+"/output.hdf5","a") as f:

            dset = f["written"]
            print("Written shape ",self.timecpu)

            f["itime"][0] = itimecpu
            f["written"][0] = self.writtencpu
            f["time"][self.writtencpu] = self.timecpu
            f["enval"][self.writtencpu,0] = hamcpu
            f["enval"][self.writtencpu,1] = mhcpu
            f["enval"][self.writtencpu,2] = chcpu
            f["enval"][self.writtencpu,3] = kinhamcpu
            f["enval"][self.writtencpu,4] = maghamcpu
            f["enval"][self.writtencpu,5] = lwcpu
            f["enval"][self.writtencpu,6] = lccpu
            f["enval"][self.writtencpu,7] = rwcpu
            f["enval"][self.writtencpu,8] = rccpu
            f["enval"][self.writtencpu,9] = vortcpu
            f["enval"][self.writtencpu,10] = cb2cpu
            f["enval"][self.writtencpu,11] = mhccpu
            
            f["magneticfields"][self.writtencpu,:,:,:,:] = b1cpu
            f["velocityfields"][self.writtencpu,:,:,:,:] = v1cpu

            if self.initcond == "threewave":
                f["threewaveenergies"][self.writtencpu,:] = waveenergiescpu
                
        return(None)

class DNA2MHD(RHS,DIAGS):
    def __init__(self,nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                 dt,iterations,lpath,linear=False,
                 initcond="hallwave",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                 forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],hyper=1,hallparam=1.0,
                 solveprec=16,maxwallclock=86200,triplet=None,records=20,bittype=32,exactnueta=False):

        # Define working real/complex data type
        if bittype == 32:
            floattype = cp.float32
            complextype = cp.complex64
        else:
            floattype = cp.float64
            complextype = cp.complex128

        if lpath[-1] == "/":
            lpath = lpath[:-1]
        self.lpath = lpath

        self.ctype = "complex"+str(2*bittype)
        self.rtype = "float"+str(bittype)

        # Intended runtime parameters
        self.maxwallclock = maxwallclock
        self.dt = floattype(dt)
        self.iterations = iterations
        self.itime = int32(0)
        self.itime_start = int32(0)
        self.time = floattype(0)
        self.solveprec = solveprec
        
        # Define grid parameters

        # Grid sizes in each direction
        self.nkx0 = int32(nkx0)
        self.nky0 = int32(nky0)
        self.nkz0 = int32(nkz0)

        self.nx0_big = 3 * nkx0//2
        self.ny0_big = 3 * nky0//2
        self.nz0_big = 3 * nkz0//2

        # Smallest k scales resolved        
        self.kxmin = floattype(kxmin)
        self.kymin = floattype(kymin)
        self.kzmin = floattype(kzmin)

        self.kxgrid = cp.concatenate((cp.arange(0,self.nx0_big//2,dtype=self.rtype),cp.arange(-self.nx0_big//2,0,dtype=self.rtype)))*kxmin
        self.kygrid = cp.concatenate((cp.arange(0,self.ny0_big//2,dtype=self.rtype),cp.arange(-self.ny0_big//2,0,dtype=self.rtype)))*kymin
        self.kzgrid = cp.arange(0,self.nz0_big//2+1,dtype=self.rtype)*kzmin

        KX,KY,KZ = cp.meshgrid(self.kxgrid,self.kygrid,self.kzgrid,indexing="ij")

        self.kmags = cp.sqrt(self.kxgrid[:,None,None]**2.0 + self.kygrid[None,:,None]**2.0+self.kzgrid[None,None,:]**2.0).astype(self.rtype)
        self.kmax = cp.amax(self.kmags)

        # Positive curl eigenstates defined
        zvec = cp.array([0,0,1],dtype=self.ctype)
        ks = cp.stack((KX,KY,KZ)).astype(self.ctype)

        self.pcurleig = cp.cross(ks,zvec[:,None,None,None],axis=0)
        self.pcurleig += 1j* cp.cross(ks,self.pcurleig,axis=0)/self.kmags[None,:,:,:]
        self.pcurleig *= cp.sqrt(2)/(2*self.kmags[:,:,0][None,:,:,None])
        self.pcurleig[:,0,0,:] = cp.array([1,1j,0],dtype=self.ctype)[:,None]/cp.sqrt(2.0)
        self.pcurleig[:,0,0,0] = 0.0

        # Initial condition parameters 

        self.initcond = initcond
        self.energystart = energystart
        
        self.init_kolm = floattype(init_kolm) # Power law decay of initial modes for full mode initialization
        
        # Initial and forced normal modes: 0 whistler, 1 left cyclotron, 2 right whistler, 3 right cyclotron  
        self.hmhdwave = hmhdwave 
        self.forcetype = forcetype
        
        self.forceamp = floattype(forceamp)
        self.nforce = nforce
        self.forcewave = forcewave
        
        self.nu = floattype(nu) # Dissipation at largest scale
        self.eta = floattype(eta) # Magnetic Prandtl number
        self.hyper = hyper

        # Three wavenumbers for resonant interaction in threewave simulation
        # Either positive whistlers with 9 wavenumbers (wave x3, wave2 x3, wave3 x3); or 12 with above mode types specified after wavenumbers
        self.triplet = triplet
        
        self.exactnueta = exactnueta        
        if exactnueta == False:
            # Normalized viscosity and resistivity from largest scale dissipation and Prandtl number
            self.vnu = nu/(self.kmax**(2*hyper))
            self.etab = eta*self.vnu
        else:
            # Specify viscosity, resistivity in di^2 wci, mu0 di^2 wci; e.g.
            self.vnu = nu
            self.etab = eta
        
        print("Dissipation Factors",self.vnu,self.etab)

        # Ratio of frequency to kz for ideal positive helicity cyclotron modes; the other frequency ratios can be determined from this
        self.alphaleftwhist = - (hallparam*self.kmags/2 + cp.sqrt(1+ (hallparam*self.kmags/2)**2.0))

        # Calculation of non-ideal Hall MHD eigenvalues
        
        eigb2 = - ((self.vnu+self.etab)*self.kmags**(2*hyper) + 1j * self.kmags * self.kzgrid[None,None,:]*hallparam)/2
        eigsqrt = eigb2**2 - self.kzgrid[None,None,:]**2.0 - self.vnu * self.etab * self.kmags**(4*hyper) \
            - 1j * self.vnu * self.kmags**(2*hyper) * self.kzgrid[None,None,:] * self.kmags
        
        self.eig = cp.zeros([self.nx0_big,self.ny0_big,self.nz0_big//2+1,4],dtype=self.ctype)
        self.eig[:,:,:,0] = eigb2 - cp.sqrt(eigsqrt)
        self.eig[:,:,:,1] = eigb2 + cp.sqrt(eigsqrt)
        self.eig[:,:,:,2] = cp.conj(eigb2) + cp.sqrt(cp.conj(eigsqrt))
        self.eig[:,:,:,3] = cp.conj(eigb2) - cp.sqrt(cp.conj(eigsqrt))
        

        # Padding array for dealiasing
        self.padding = cp.ones_like(self.kmags,dtype="int32")
        self.padding[self.nkx0//2:1-self.nkx0//2,:,:] = 0
        self.padding[:,self.nky0//2:1-self.nky0//2,:] = 0
        self.padding[:,:,self.nkz0:] = 0

        # Mask of forcing waves beyond certain mode amplitudes
        self.forcemask = cp.zeros_like(self.padding,dtype=self.rtype)
        self.forcemask[1:nforce+1,1:nforce+1,1:nforce+1] = 1
        self.forcemask[-nforce:,-nforce:,1:nforce+1] = 1
        self.forcemask[1:nforce+1,-nforce:,1:nforce+1] = 1
        self.forcemask[-nforce:,1:nforce+1,1:nforce+1] = 1
        self.forcemask *= forceamp
        
        # Initialize fields
        self.b1 = cp.zeros_like(self.pcurleig)
        self.v1 = cp.zeros_like(self.pcurleig)
        self.mhelcorr = cp.float32(0.0)
        self.fieldsetup()

        # Testing linear simulation
        self.linear = linear
        # Ratio of ion skin depth to length scale; this should be 1 with normalization used, and length scale enters through k
        self.hallparam = hallparam

        # Output parameters - filename, how many times to store fields, and set up HDF5 output file
        self.record = records
        self.irecord = self.iterations // records
        self.iterations += self.itime_start
        print("Recording Time ",self.irecord,"dt ",self.dt)
        self.written = 0
        self.createhdf5file()
        
        return(None)

    def fieldsetup(self):

        # Specify fields as a combination of all ideal HMHD normal modes
        if self.initcond == "hallwave":
            
            phase = self.reset_phase()[1:,1:,1:]
            self.b1[:,1:,1:,1:] = self.pcurleig[:,1:,1:,1:] * self.alphaleftwhist[None,1:,1:,1:] \
                / cp.sqrt(self.alphaleftwhist[None,1:,1:,1:]**2.0 + 1) * phase * self.hmhdwave[0]
            self.v1[:,1:,1:,1:] = self.pcurleig[:,1:,1:,1:] / cp.sqrt(self.alphaleftwhist[None,1:,1:,1:]**2.0 + 1) * phase * self.hmhdwave[0]
            
            phase = self.reset_phase()[1:,1:,1:]
            self.b1[:,1:,1:,1:] += self.pcurleig[:,1:,1:,1:] * -self.alphaleftwhist[None,1:,1:,1:] \
                / cp.sqrt(self.alphaleftwhist[None,1:,1:,1:]**2.0 + 1) * phase * self.hmhdwave[2]
            self.v1[:,1:,1:,1:] += self.pcurleig[:,1:,1:,1:] / cp.sqrt(self.alphaleftwhist[None,1:,1:,1:]**2.0 + 1) * phase * self.hmhdwave[2]
            
            phase = self.reset_phase()[1:,1:,1:]
            self.b1[:,1:,1:,1:] += self.pcurleig[:,1:,1:,1:] / cp.sqrt(self.alphaleftwhist[None,1:,1:,1:]**2.0 + 1) * phase * self.hmhdwave[1]
            self.v1[:,1:,1:,1:] += self.pcurleig[:,1:,1:,1:] * -self.alphaleftwhist[None,1:,1:,1:] \
                / cp.sqrt(self.alphaleftwhist[None,1:,1:,1:]**2.0 + 1) * phase * self.hmhdwave[1]
            
            phase = self.reset_phase()[1:,1:,1:]
            self.b1[:,1:,1:,1:] += self.pcurleig[:,1:,1:,1:] / cp.sqrt(self.alphaleftwhist[None,1:,1:,1:]**2.0 + 1) * phase * self.hmhdwave[3]
            self.v1[:,1:,1:,1:] += self.pcurleig[:,1:,1:,1:] * self.alphaleftwhist[None,1:,1:,1:] \
                / cp.sqrt(self.alphaleftwhist[None,1:,1:,1:]**2.0 + 1) * phase * self.hmhdwave[3]
            
        self.b1[:,1:,1:,1:] *= self.padding[None,1:,1:,1:] * self.kmags[None,1:,1:,1:]**(-self.init_kolm)
        self.v1[:,1:,1:,1:] *= self.padding[None,1:,1:,1:] * self.kmags[None,1:,1:,1:]**(-self.init_kolm)

        # Specify fields as a three wave interaction
        
        if self.initcond == "threewave" or self.initcond == "onewave":

            self.itriplet = []

            if self.triplet == None or  (len(self.triplet) != 9 and len(self.triplet) != 12):
                
                raise ValueError("LW Triplet or triplet and normal mode types must be specified for three wave initial condition")

            # Relative amplitudes
            if self.initcond == "threewave":
                amp0 = 2.0
                amp1 = 1.0
                amp2 = 0.5
            else:
                amp0 = 1.0
                amp1 = 0.0
                amp2 = 0.0

            # Specify triplet as either list of three wavevectors or list of three wave vectors and normal mode types
            # Wave triplet - 1 + whistler, 0 + cyclotron, 2 - whistler, 3 - cyclotron

            if len(self.triplet) == 9:
                # Either initialize wave triplets as positive helicity whistlers from wavenumbers
                
                for i in range(3):
                    self.triplet.append(0)
                self.initializewave(self.triplet[0:3],amp0,0)
                if amp1 > 0:
                    self.initializewave(self.triplet[3:6],amp1,0)
                    self.initializewave(self.triplet[6:9],amp2,0)

            else:

                # Or initialize triplets from wavenumbers and the last three digits of self.triplet specify the different waves
                self.initializewave(self.triplet[0:3],amp0,self.triplet[9])
                if amp1 > 0:
                    self.initializewave(self.triplet[3:6],amp1,self.triplet[10])
                    self.initializewave(self.triplet[6:9],amp2,self.triplet[11])

            self.triplet = cp.array(self.triplet)

        if self.energystart != None:

            # Normalize data to a given initial energy
            
            energy = cp.sum(cp.abs(self.b1[:,:,:,0])**2.0+cp.abs(self.v1[:,:,:,0])**2.0)
            energy += 2* cp.sum(cp.abs(self.b1[:,:,:,1:])**2.0+cp.abs(self.v1[:,:,:,1:])**2.0)

            self.b1 *= cp.sqrt(self.energystart/energy)
            self.v1 *= cp.sqrt(self.energystart/energy)
                
        if self.initcond == "checkpoint":

            # Use the last step of a previously written simulation to start data using HDF5

            # CPU read data
            with h5py.File(self.lpath+"/output.hdf5","r") as f:

                written = f["written"][0]
                b1cpu = f["magneticfields"][written,:,:,:,:]
                v1cpu = f["velocityfields"][written,:,:,:,:]
                timecpu = f["time"][written]
                dtcpu = f["dt"][0]
                itimecpu = f["itime"][0]
                mhc = f["enval"][written,-1]

            # Send to GPU storage
            self.b1 = cp.asarray(b1cpu)
            self.v1 = cp.asarray(v1cpu)
            self.time = cp.asarray(timecpu)
            self.itime = cp.asarray(itimecpu)
            self.itime_start = cp.copy(self.itime)
            self.dt = cp.asarray(dtcpu)
            self.mhelcorr = cp.asarray(mhc)

        return(None)

    def initializewave(self,wave,amp,modeno=0):

        # Set up normal mode type modeno for indices given by rescaling the wavenumber wave

        # Make kz index of wave positive for positive kz wavenumber storage
        if wave[2] < 0:
            wave[0] *= -1
            wave[1] *= -1
            wave[2] *= -1

        ix = np.round(wave[0],decimals=3)//self.kxmin
        iy = np.round(wave[1],decimals=3)//self.kymin
        iz = np.round(wave[2],decimals=3)//self.kzmin

        self.itriplet.append(ix)
        self.itriplet.append(iy)
        self.itriplet.append(iz)
        
        if modeno < 2:

            self.b1[:,ix,iy,iz] = self.pcurleig[:,ix,iy,iz]
            self.v1[:,ix,iy,iz] = self.pcurleig[:,ix,iy,iz]

        else:

            self.b1[:,ix,iy,iz] = cp.conj(self.pcurleig[:,ix,iy,iz])
            self.v1[:,ix,iy,iz] = cp.conj(self.pcurleig[:,ix,iy,iz])

        self.b1[:,ix,iy,iz] *= (self.alphaleftwhist[ix,iy,iz] ** (1 - 2 * (modeno % 2))) * ((-1)**(modeno//2))
        self.v1[:,ix,iy,iz] *= 1

        # Set wave amplitude to amp
        
        self.b1[:,ix,iy,iz] *= amp/cp.sqrt((self.alphaleftwhist[ix,iy,iz] ** (1 - 2 * (modeno % 2)))**2 + 1)
        self.v1[:,ix,iy,iz] *= amp/cp.sqrt((self.alphaleftwhist[ix,iy,iz] ** (1 - 2 * (modeno % 2)))**2 + 1)

        return(None)

    def fixedpointidealgauss2(self):

        a11 = 0.5

        self.b2 = self.b1.copy()
        self.v2 = self.v1.copy()
        self.brhs2,self.vrhs2 = self.idealhallrhs(self.b2,self.v2)
            
        self.brhs1 = self.brhs2.copy()
        self.vrhs1 = self.vrhs2.copy()

        self.b2 += a11 * self.dt * self.brhs1
        self.v2 += a11 * self.dt * self.vrhs1

        self.brhs2,self.vrhs2 = self.idealhallrhs(self.b2,self.v2)

        linferror = max(cp.amax(cp.abs(self.brhs2-self.brhs1)),cp.amax(cp.abs(self.vrhs2-self.vrhs1)))
            
        solveiteration = 0
        # Typical convergent mins are 1e-9, so set a bit above that
        while solveiteration < 40 and linferror > 10.0**(-8.0):

            self.brhs1 = cp.copy(self.brhs2)
            self.vrhs1 = cp.copy(self.vrhs2)

            self.b2 = self.b1 + a11 * self.dt * self.brhs1
            self.v2 = self.v1 + a11 * self.dt * self.vrhs1

            self.brhs2,self.vrhs2 = self.idealhallrhs(self.b2,self.v2)
                
            linferror = max(cp.amax(cp.abs(self.brhs2-self.brhs1)),cp.amax(cp.abs(self.vrhs2-self.vrhs1)))
            solveiteration += 1

        if self.itime % 10 == 0:
            print("Iterations itime ",self.itime," = ",solveiteration)

        if solveiteration == 40:
            print("Fixed point solver fail error ",linferror," itime ",self.itime)

        self.b1 += self.dt * self.brhs2
        self.v1 += self.dt * self.vrhs2
        self.mhelcorr += self.helicitycorrection(self.b2,self.v2) * self.dt

        return(None)
            
    def gauss2split(self): 
        # Splitting method for dissipation, implicit midpoint nonlinear

        self.startruntime = time.time()
        self.runtime = 0

        # Add guide field to do linear and nonlinear calculation at once                                                                                                                   
        self.b1[2,0,0,0] = 1.0

        a11 = 0.5

        while self.itime < self.iterations and self.runtime < self.maxwallclock: 

            if self.irecord > 0 and (self.itime - self.itime_start) % self.irecord == 0 and (self.itime > self.itime_start or self.itime_start == 0):
                self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
                self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
                self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
                self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

                self.writefile()
                self.written += 1

            self.b1 *= cp.exp(-self.etab * self.kmags[None,:,:,:]**(2*self.hyper) * self.dt/2)
            self.v1 *= cp.exp(-self.vnu * self.kmags[None,:,:,:]**(2*self.hyper) * self.dt/2)
            
            self.fixedpointidealgauss2()

            self.b1 *= cp.exp(-self.etab * self.kmags**(2*self.hyper) * self.dt/2)
            self.v1 *= cp.exp(-self.vnu * self.kmags**(2*self.hyper) * self.dt/2)

            self.itime += 1
            self.time += self.dt
            self.runtime = time.time()-self.startruntime

        self.writefile()
        print(self.runtime)

        return(None)

    def linearupdate(self,dt):
        """Exact linear update for splitting or exponential method"""

        if self.itime == 0:
            # Compute exponential for positive helicity branches; negative branch is inverse of positive branch
            self.expwhistler = cp.exp(1j * dt * self.kzgrid[None,None,:] * self.alphaleftwhist)
            self.expcyclotron = cp.exp(-1j * dt * self.kzgrid[None,None,:] / self.alphaleftwhist)

        self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
        self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
        self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
        self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

        # Positive helicity whistler
        # Compute normal mode coefficients by a dot product with the orthogonal normal mode basis
        lwk = (self.alphaleftwhist * self.bplus + self.vplus)/cp.sqrt(self.alphaleftwhist**2 + 1)
        # Advance linear modes through multiplication with exponential of frequency
        lwk *= self.expwhistler
        # Reconstruct b and v from the evolved normal modes
        self.b1 = self.pcurleig * lwk * self.alphaleftwhist /cp.sqrt(self.alphaleftwhist**2 + 1)
        self.v1 = self.pcurleig * lwk /cp.sqrt(self.alphaleftwhist**2 + 1)

        # Positive helicity cyclotron
        lwk = (-1/self.alphaleftwhist * self.bplus + self.vplus)/cp.sqrt((1/self.alphaleftwhist)**2 + 1)
        lwk *= self.expcyclotron
        self.b1 -= self.pcurleig * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)
        self.v1 += self.pcurleig * lwk * self.alphaleftwhist / cp.sqrt(self.alphaleftwhist**2 + 1)

        # Negative helicity whistler
        lwk = (-self.alphaleftwhist * self.bminus + self.vminus)/cp.sqrt(self.alphaleftwhist**2 + 1)
        lwk *= 1/self.expwhistler
        self.b1 -= self.alphaleftwhist * cp.conj(self.pcurleig) * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)
        self.v1 += cp.conj(self.pcurleig) * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)

        # Negative helicity cyclotron
        lwk = (1/self.alphaleftwhist * self.bminus + self.vminus)/cp.sqrt((1/self.alphaleftwhist)**2 + 1)
        lwk *= 1/self.expcyclotron
        self.b1 += cp.conj(self.pcurleig) * lwk /cp.sqrt(self.alphaleftwhist**2 + 1)
        self.v1 += cp.conj(self.pcurleig) * self.alphaleftwhist * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)

        return(None)

    def idealRalston2(self):

        self.b2 = self.b1.copy()
        self.v2 = self.v1.copy()
        self.brhs2,self.vrhs2 = self.idealhallrhs(self.b2,self.v2)
        mhc1 = self.helicitycorrection(self.b2,self.v2)
            
        self.brhs1 = self.brhs2.copy()
        self.vrhs1 = self.vrhs2.copy()

        self.b2 += 2/3 * self.dt * self.brhs1
        self.v2 += 2/3 * self.dt * self.vrhs1

        self.brhs2,self.vrhs2 = self.idealhallrhs(self.b2,self.v2)
        mhc2 = self.helicitycorrection(self.b2,self.v2)

        self.b1 += self.dt / 4 * (self.brhs1 + 3 * self.brhs2) 
        self.v1 += self.dt / 4 * (self.vrhs1 + 3 * self.vrhs2)
        self.mhelcorr += self.dt / 4 * (mhc1 + 3 * mhc2)
        
        return(None)

    def LinearSplitNLR2(self):

        """Splitting dissipation and linear terms, nonlinear terms Ralston RK2"""

        self.startruntime = time.time()
        self.runtime = 0

        while self.itime < self.iterations and self.runtime < self.maxwallclock: 

            if self.irecord > 0 and (self.itime - self.itime_start) % self.irecord == 0 and (self.itime > self.itime_start or self.itime_start == 0):
                self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
                self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
                self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
                self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

                self.writefile()
                self.written += 1

            # First step - half dissipation
            self.b1 *= cp.exp(-self.etab * self.kmags[None,:,:,:]**(2*self.hyper) * self.dt/2)
            self.v1 *= cp.exp(-self.vnu * self.kmags[None,:,:,:]**(2*self.hyper) * self.dt/2)

            # Second stage - half step ideal linear
            self.linearupdate(self.dt/2)

            # Third stage - nonlinear terms Ralston RK2
            self.fixedpointidealgauss2()

            # Fourth stage - half step ideal linear
            self.linearupdate(self.dt/2)

            # Final stage - half dissipation
            self.b1 *= cp.exp(-self.etab * self.kmags**(2*self.hyper) * self.dt/2)
            self.v1 *= cp.exp(-self.vnu * self.kmags**(2*self.hyper) * self.dt/2)

            self.itime += 1
            self.time += self.dt
            self.runtime = time.time()-self.startruntime

        self.writefile()
        print(self.runtime)

        return(None)

    def exptime_ideal(self,records=200,order=1,cutoff=0.01,kt=0):

        """Use exponential first order integrator for ideal Hall MHD"""
    
        self.startruntime = time.time()
        self.runtime = 0

        # Write eigenvalues for the system - positive whistler, positive cyclo, negative whistler, negative cyclo
        self.eig[:,:,:,0] = 1j * self.kzgrid[None,None,:] * self.alphaleftwhist
        self.eig[:,:,:,1] = -1j * self.kzgrid[None,None,:] / self.alphaleftwhist
        self.eig[:,:,:,2] = -1j * self.kzgrid[None,None,:] * self.alphaleftwhist
        self.eig[:,:,:,3] = 1j * self.kzgrid[None,None,:] / self.alphaleftwhist

        h = self.eig * self.dt

        print("Range eig dt ",cp.amin(cp.abs(h[h>0])),cp.amax(cp.abs(h)))

        # Exponential coefficient (exp(L dt) - I)/L
        coef1 = cp.zeros([self.nx0_big,self.ny0_big,self.nz0_big//2+1,4],dtype=self.ctype)
        ii = cp.abs(h) < cutoff
        coef1[~ii] = self.dt * (cp.exp(h[~ii]) - 1)/h[~ii]
        # Find coefficients with cutoff for small eig dt with (4,4) Pade approximant h[ii]
        coef1[ii] = self.dt*(1 + 1/26 * h[ii] +5/156 * h[ii]**2 +1/858 *h[ii]**3 + 1/5720 *h[ii]**4 + 1/205920 *h[ii]**5 +1/8648640 *h[ii]**6)
        coef1[ii] /= (1 - 6/13*h[ii] +5/52 *h[ii]**2 -5/429 * h[ii]**3 +1/1144 *h[ii]**4 - 1/25740 *h[ii]**5+1/1235520 *h[ii]**6)

        coef1[ii] = self.dt * (1 +1/10 * h[ii] +1/60 * h[ii]**2 )
        coef1[ii] /= (1 - 2/5 * h[ii] +1/20 * h[ii]**2 )

        #coef1 = self.dt * cp.ones_like(self.eig)
        print(cp.amax(cp.abs(coef1-self.dt)))
        print(cp.any(cp.isnan(coef1)))

        if order == 2 : # Second order coefficient (exp(L dt) - I - L dt)/(dt L**2)
            coef2 = cp.zeros_like(coef1)
            coef2[~ii] = self.dt * (cp.exp(h[~ii]) - 1 - h[~ii])/(h[~ii]**2)

            coef2[ii] = self.dt * (1/2 - 1/21 * h[ii] + 25/2184 *h[ii]**2 - 1/2730 *h[ii]**3 +3/80080 *h[ii]**4 +1/121080960 *h[ii]**6)
            coef2[ii] /= (1 - 3/7 * h[ii] +15/182 *h[ii]**2 - 5/546 *h[ii]**3 + 5/8008 *h[ii]**4 - 1/40040 *h[ii]**5+1/2162160 *h[ii]**6)

            coef2[ii] = (1/2 + 1/360 * h[ii]**2 )*self.dt
            coef2[ii] /= (1 - 1/3 * h[ii] +1/30 * h[ii]**2 )

        if kt > 0: # Get the coefficients using Kassam Trefethen contour integration - kt is number of quad points
            
            coef1 *= 0.0
            r = cp.exp(1j * (cp.arange(kt)/kt + 1/cp.pi ) * 2 * cp.pi)

            for ir in r:
                coef1 += (cp.exp(ir + h)-1)/(ir + h)

            coef1 *= self.dt / kt
            if order == 2:
                coef2 *= 0.0

                for ir in r:
                    coef2 += (cp.exp(ir + h)-1-(ir + h))/((ir + h)**2)
                
                coef2 *= self.dt/kt
    
        # Time loop
        while self.itime < self.iterations and self.runtime < self.maxwallclock:

            # Obtain positive and negative curl eigenstate projections
            self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
            self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
            self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
            self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

            #dot1 -= time.time()
            #print("dot time ", np.abs(dot1))

            # File output
            if self.irecord > 0 and (self.itime - self.itime_start) % self.irecord== 0:

                self.writefile()

                self.written += 1
                
            # Get nonlinear terms and helicity correction
            self.brhs1,self.vrhs1 = self.idealhallrhs(self.b1,self.v1)
            mhc = self.helicitycorrection(self.b1,self.v1)
            self.mhelcorr += self.dt * mhc

            # Linear advance of fields with matrix exp - done in linear update
            self.linearupdate(self.dt)

            # Decompose the nonlinear terms into curl eigenstates
            self.bplus = cp.sum(cp.conj(self.pcurleig)*self.brhs1,axis=0)
            self.vplus = cp.sum(cp.conj(self.pcurleig)*self.vrhs1,axis=0)
            self.bminus = cp.sum(self.pcurleig*self.brhs1,axis=0)
            self.vminus = cp.sum(self.pcurleig*self.vrhs1,axis=0)

            # Increment b and v according to eigenvalue coefficients
            # Positive helicity whistler
            lwk = (self.alphaleftwhist * self.bplus + self.vplus)/cp.sqrt(self.alphaleftwhist**2 + 1)
            lwk *= coef1[:,:,:,0]
            self.b1 += self.pcurleig * lwk * self.alphaleftwhist /cp.sqrt(self.alphaleftwhist**2 + 1)
            self.v1 += self.pcurleig * lwk /cp.sqrt(self.alphaleftwhist**2 + 1)

            # Positive helicity cyclotron
            lwk = (-1/self.alphaleftwhist * self.bplus + self.vplus)/cp.sqrt((1/self.alphaleftwhist)**2 + 1)
            lwk *= coef1[:,:,:,1]
            self.b1 -= self.pcurleig * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)
            self.v1 += self.pcurleig * lwk * self.alphaleftwhist / cp.sqrt(self.alphaleftwhist**2 + 1)

            # Negative helicity whistler
            lwk = (-self.alphaleftwhist * self.bminus + self.vminus)/cp.sqrt(self.alphaleftwhist**2 + 1)
            lwk *= coef1[:,:,:,2]
            self.b1 -= self.alphaleftwhist * cp.conj(self.pcurleig) * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)
            self.v1 += cp.conj(self.pcurleig) * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)

            # Negative helicity cyclotron
            lwk = (1/self.alphaleftwhist * self.bminus + self.vminus)/cp.sqrt((1/self.alphaleftwhist)**2 + 1)
            lwk *= coef1[:,:,:,3]
            self.b1 += cp.conj(self.pcurleig) * lwk /cp.sqrt(self.alphaleftwhist**2 + 1)
            self.v1 += cp.conj(self.pcurleig) * self.alphaleftwhist * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)
            
            if order == 2:
                # ETD2RK time stepping
                self.brhs2,self.vrhs2 = self.idealhallrhs(self.b1,self.v1)
                mhc2 = self.helicitycorrection(self.b1,self.v1)

                self.mhelcorr += (mhc2 - mhc) * self.dt/2

                self.brhs1 = self.brhs2 - self.brhs1
                self.vrhs1 = self.vrhs2 - self.vrhs1

                # Decompose again
                self.bplus = cp.sum(cp.conj(self.pcurleig)*self.brhs1,axis=0)
                self.vplus = cp.sum(cp.conj(self.pcurleig)*self.vrhs1,axis=0)
                self.bminus = cp.sum(self.pcurleig*self.brhs1,axis=0)
                self.vminus = cp.sum(self.pcurleig*self.vrhs1,axis=0)

                lwk = (self.alphaleftwhist * self.bplus + self.vplus)/cp.sqrt(self.alphaleftwhist**2 + 1)
                lwk *= coef2[:,:,:,0]
                self.b1 += self.pcurleig * lwk * self.alphaleftwhist /cp.sqrt(self.alphaleftwhist**2 + 1)
                self.v1 += self.pcurleig * lwk /cp.sqrt(self.alphaleftwhist**2 + 1)

                lwk = (-1/self.alphaleftwhist * self.bplus + self.vplus)/cp.sqrt((1/self.alphaleftwhist)**2 + 1)
                lwk *= coef2[:,:,:,1]
                self.b1 -= self.pcurleig * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)
                self.v1 += self.pcurleig * lwk * self.alphaleftwhist / cp.sqrt(self.alphaleftwhist**2 + 1)

                lwk = (-self.alphaleftwhist * self.bminus + self.vminus)/cp.sqrt(self.alphaleftwhist**2 + 1)
                lwk *= coef2[:,:,:,2]
                self.b1 -= self.alphaleftwhist * cp.conj(self.pcurleig) * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)
                self.v1 += cp.conj(self.pcurleig) * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)

                lwk = (1/self.alphaleftwhist * self.bminus + self.vminus)/cp.sqrt((1/self.alphaleftwhist)**2 + 1)
                lwk *= coef2[:,:,:,3]
                self.b1 += cp.conj(self.pcurleig) * lwk /cp.sqrt(self.alphaleftwhist**2 + 1)
                self.v1 += cp.conj(self.pcurleig) * self.alphaleftwhist * lwk / cp.sqrt(self.alphaleftwhist**2 + 1)

            self.itime += 1
            self.time += self.dt
            self.runtime = time.time()-self.startruntime

            #if self.itime % max(self.iterations//50,10) == 0 and self.itime > 0:
            #    self.steadystate()

        self.writefile()
        print(self.runtime)

        return(None)

    def explicit(self,records=200,order=1):

        """Library of explicit schemes I've used for benchmarking"""

        self.startruntime = time.time()
        self.runtime = 0

        # Add guide field to do linear and nonlinear calculation at once

        while self.itime < self.iterations and self.runtime < self.maxwallclock:
            
            if self.irecord > 0 and (self.itime - self.itime_start) % self.irecord== 0 :

                self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
                self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
                self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
                self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

                self.writefile()
                self.written += 1

            if order == 2: # Ralston RK2 

                # First RHS evaluation
                k1b,k1v,mhc1 = self.hallrhs(self.b1,self.v1)

                # Store original state
                b1_old = self.b1.copy()
                v1_old = self.v1.copy()
                
                # Intermediate state
                perturbation = 2/3 * self.dt * k1b
                b2_temp = b1_old + perturbation
                v2_temp = v1_old + 2/3 * self.dt * k1v

                # Debug: check the perturbation size
                if self.itime == self.itime_start:
                    pert_norm = float(cp.linalg.norm(perturbation))
                    diff_norm = float(cp.linalg.norm(b2_temp - b1_old))
                    print(f"DEBUG Ralston: perturbation norm={pert_norm:.6e}, diff norm={diff_norm:.6e}")
                    print(f"DEBUG Ralston: dt={self.dt:.6e}, k1_norm={float(cp.linalg.norm(k1b)):.6e}")
                    print(f"DEBUG Ralston: 2/3*dt*k1_expected = {2/3 * self.dt * float(cp.linalg.norm(k1b)):.6e}")

                # Second RHS evaluation at intermediate state
                k2b,k2v,mhc2 = self.hallrhs(b2_temp,v2_temp)
                
                # Debug: check if k1 and k2 are different
                if self.itime == self.itime_start:
                    k1_norm = float(cp.linalg.norm(k1b))
                    k2_norm = float(cp.linalg.norm(k2b))
                    diff_k = float(cp.linalg.norm(k2b - k1b))
                    print(f"DEBUG Ralston: k1_norm={k1_norm:.6e}, k2_norm={k2_norm:.6e}, ||k2-k1||={diff_k:.6e}")
                
                # Final update: y_{n+1} = y_n + dt*(1/4*k1 + 3/4*k2)
                self.b1 = b1_old + 1/4 * self.dt * k1b + 3/4 * self.dt * k2b
                self.v1 = v1_old + 1/4 * self.dt * k1v + 3/4 * self.dt * k2v
                self.mhelcorr += (1/4 * self.dt * mhc1 + 3/4 * self.dt * mhc2)

            if order == 1: # Euler's method

                self.brhs1,self.vrhs1,mhc = self.hallrhs(self.b1,self.v1)

                self.b1 += self.dt * self.brhs1
                self.v1 += self.dt * self.vrhs1
                self.mhelcorr += self.dt * mhc

            if order == 4: #Classic RK4

                self.b2 = self.b1.copy()
                self.v2 = self.v1.copy()

                self.brhs1,self.vrhs1,mhc1 = self.hallrhs(self.b2,self.v2)
                self.b2 = self.b1 + self.dt/2 * self.brhs1
                self.v2 = self.v1 + self.dt/2 * self.vrhs1

                self.brhs2,self.vrhs2,mhc2 = self.hallrhs(self.b2,self.v2)
                self.b2 = self.b1 + self.dt/2 * self.brhs2
                self.v2 = self.v1 + self.dt/2 * self.vrhs2

                self.brhs3,self.vrhs3,mhc3 = self.hallrhs(self.b2,self.v2)
                self.b2 = self.b1 + self.dt * self.brhs3
                self.v2 = self.v1 + self.dt * self.vrhs3

                self.brhs4,self.vrhs4,mhc4 = self.hallrhs(self.b2,self.v2)

                self.b1 += self.dt/6 * (self.brhs1 + self.brhs4 + 2*(self.brhs2 + self.brhs3))
                self.v1 += self.dt/6 * (self.vrhs1 + self.vrhs4 + 2*(self.vrhs2 + self.vrhs3))
                self.mhelcorr += self.dt/6 * (mhc1 + mhc4 + 2 * (mhc2 + mhc3)) 

            if order == 5:
                """5th order stability region Dormand Prince scheme https://doi.org/10.1016/0771-050X(80)90013-3"""

                self.b2 = self.b1.copy()
                self.v2 = self.v1.copy()

                self.brhs1,self.vrhs1,mhc1 = self.hallrhs(self.b1,self.v1)

                self.b2 = self.b1 + self.dt * 2/9 * self.brhs1
                self.v2 = self.v1 + self.dt * 2/9 * self.vrhs1
                
                self.brhs2,self.vrhs2,mhc2 = self.hallrhs(self.b2,self.v2)

                self.b2 = self.b1 + self.dt * (1/12 * self.brhs1 + 1/4 * self.brhs2)
                self.v2 = self.v1 + self.dt * (1/12 * self.vrhs1 + 1/4 * self.vrhs2)

                self.brhs3,self.vrhs3,mhc3 = self.hallrhs(self.b2,self.v2)

                self.b2 = self.b1 + self.dt * (55/324 * self.brhs1 - 25/108 * self.brhs2 + 50/81 * self.brhs3)
                self.v2 = self.v1 + self.dt * (55/324 * self.vrhs1 - 25/108 * self.vrhs2 + 50/81 * self.vrhs3)

                self.brhs4,self.vrhs4,mhc4 = self.hallrhs(self.b2,self.v2)

                self.b2 = self.b1 + self.dt * (83/330 * self.brhs1 - 13/22 * self.brhs2 + 61/66 * self.brhs3 + 9/110 * self.brhs4)
                self.v2 = self.v1 + self.dt * (83/330 * self.vrhs1 - 13/22 * self.vrhs2 + 61/66 * self.vrhs3 + 9/110 * self.vrhs4)

                self.brhs5,self.vrhs5,mhc5 = self.hallrhs(self.b2,self.v2)

                self.b2 = self.b1 + self.dt * (-19/28 * self.brhs1 + 9/4 * self.brhs2 + 1/7 * self.brhs3 + -27/7 * self.brhs4 + 22/7 * self.brhs5)
                self.v2 = self.v1 + self.dt * (-19/28 * self.vrhs1 + 9/4 * self.vrhs2 + 1/7 * self.vrhs3 + -27/7 * self.vrhs4 + 22/7 * self.vrhs5)

                self.brhs6,self.vrhs6,mhc6 = self.hallrhs(self.b2,self.v2)

                self.b1 += (19/200 * self.brhs1 + 3/5 * self.brhs3 - 243/400 * self.brhs4 + 33/40 * self.brhs5 + 7/80  * self.brhs6) * self.dt
                self.v1 += (19/200 * self.vrhs1 + 3/5 * self.vrhs3 - 243/400 * self.vrhs4 + 33/40 * self.vrhs5 + 7/80  * self.vrhs6) * self.dt 
                self.mhelcorr += (19/200 * mhc1 + 3/5 * mhc3 -243/400 * mhc4 + 33/40 * mhc5 + 7/80 * mhc6) * self.dt
                
            self.itime += 1
            self.time += self.dt
            self.runtime = time.time()-self.startruntime

        self.writefile()
        print(self.runtime)

        return(None)
