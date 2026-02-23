import cupy as cp
from cupy.fft import rfftn,irfftn
import time
from dna2mhd_expintegrator import exponentialcoefficients
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

        bout,vout = self.convolutions(binput,vinput)
        bout,vout = self.removediv(bout,vout)

        return(bout,vout)
    
    def hallrhs(self,binput,vinput):
        
        bout,vout = self.convolutions(binput,vinput)

        bout -= self.etab * self.kmags**(2*self.hyper) * binput
        vout -= self.vnu * self.kmags**(2*self.hyper) * vinput

        mhc = self.helicitycorrection(binput,vinput)

        self.removediv(bout,vout)
        
        return(bout,vout,mhc)

    def rhslinear(self,binput,vinput):

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

    def convolutions(self,binput,vinput):

        # Save some memory by doing for h non0
        #  dA/dt    = (v - h curl b) x b  dv/dt = v x curl v + curl b x b
        # d(A+h v)/dt = (v x (h curl v + b))

        bout = cp.zeros_like(binput)
        vout = cp.zeros_like(vinput)

        bx = irfftn(binput[0,:,:,:])
        by = irfftn(binput[1,:,:,:])
        bz = irfftn(binput[2,:,:,:])

        vx = irfftn(vinput[0,:,:,:])
        vy = irfftn(vinput[1,:,:,:])
        vz = irfftn(vinput[2,:,:,:])

        dum = 1j*self.kygrid[None,:,None]*binput[2,:,:,:]-1j*self.kzgrid[None,None,:]*binput[1,:,:,:]
        curlbx = irfftn(dum)

        dum = 1j*self.kzgrid[None,None,:]*binput[0,:,:,:]-1j*self.kxgrid[:,None,None]*binput[2,:,:,:]
        curlby = irfftn(dum)

        dum = 1j*self.kxgrid[:,None,None]*binput[1,:,:,:]-1j*self.kygrid[None,:,None]*binput[0,:,:,:]
        curlbz = irfftn(dum)

        dum = 1j*self.kygrid[None,:,None]*vinput[2,:,:,:]-1j*self.kzgrid[None,None,:]*vinput[1,:,:,:]
        curlvx = irfftn(dum)

        dum = 1j*self.kzgrid[None,None,:]*vinput[0,:,:,:]-1j*self.kxgrid[:,None,None]*vinput[2,:,:,:]
        curlvy = irfftn(dum)

        dum = 1j*self.kxgrid[:,None,None]*vinput[1,:,:,:]-1j*self.kygrid[None,:,None]*vinput[0,:,:,:]
        curlvz = irfftn(dum)

        dum = vy * bz - vz * by + self.hallparam * ( curlby * bz - curlbz * by )
        dum2 = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big
        bout[1,:,:,:] += 1j * self.kzgrid[None,None,:] * dum2
        bout[2,:,:,:] += -1j * self.kygrid[None,:,None] * dum2        

        dum = vz * bx	- vx * bz + self.hallparam * ( curlbz * bx - curlbx * bz )
        dum2 = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big
        bout[0,:,:,:] += -1j * self.kzgrid[None,None,:] * dum2
        bout[2,:,:,:] += 1j * self.kxgrid[:,None,None] * dum2

        dum = vx * by	- vy * bx + self.hallparam * ( curlbx * by - curlby * bx )
        dum2 = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big
        bout[0,:,:,:] += 1j * self.kygrid[None,:,None] * dum2
        bout[1,:,:,:] += -1j * self.kxgrid[:,None,None] * dum2
        
        dum = vy * curlvz - vz * curlvy + (curlby * bz - curlbz * by)
        vout[0,:,:,:] = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big
        
        dum = vz * curlvx - vx * curlvz + (curlbz * bx - curlbx * bz)
        vout[1,:,:,:] = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big
        
        dum = vx * curlvy - vy * curlvx + (curlbx * by - curlby * bx)
        vout[2,:,:,:] = rfftn(dum) * self.nx0_big * self.ny0_big * self.nz0_big
        
        bout *= self.padding[None,:,:,:]
        vout *= self.padding[None,:,:,:]

        #print("Through NL")

        return(bout,vout)

    def reset_phase(self):

        return(cp.exp(1j*cp.pi*2*cp.random.rand(self.nx0_big,self.ny0_big,self.nz0_big//2 + 1)).astype(self.ctype))
    
    def force(self):

        if self.forcetype == "hallwave":
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

        if option == 0 or option == 2 :
            ham += cp.sum(cp.abs(self.b1[:,:,:,0])**2.0,axis=(0,1,2))
            ham += 2*cp.sum(cp.abs(self.b1[:,:,:,1:])**2.0,axis=(0,1,2,3))

        if option == 0 or option == 1 :
            ham += cp.sum(cp.abs(self.v1[:,:,:,0])**2.0,axis=(0,1,2))
            ham += 2*cp.sum(cp.abs(self.v1[:,:,:,1:])**2.0,axis=(0,1,2,3))

        ham *= 4*cp.pi**3

        return(ham)
        
    def helicity(self):

        AVP = cp.stack((1j*self.kygrid[None,:,None]*self.b1[2,:,:,:]-1j*self.kzgrid[None,None,:]*self.b1[1,:,:,:],
                        1j*self.kzgrid[None,None,:]*self.b1[0,:,:,:]-1j*self.kxgrid[:,None,None]*self.b1[2,:,:,:],
                        1j*self.kxgrid[:,None,None]*self.b1[1,:,:,:]-1j*self.kygrid[None,:,None]*self.b1[0,:,:,:]),\
                       axis=0)
        cb2 = 2*cp.sum(cp.abs(AVP[:,:,:,1:])**2.0)+cp.sum(cp.abs(AVP[:,:,:,0])**2.0)
        AVP /= self.kmags[None,:,:,:]**2.0
        AVP[:,0,0,0] = 0.0
        
        WVORT = cp.stack((1j*self.kygrid[None,:,None]*self.v1[2,:,:,:]-1j*self.kzgrid[None,None,:]*self.v1[1,:,:,:],
                          1j*self.kzgrid[None,None,:]*self.v1[0,:,:,:]-1j*self.kxgrid[:,None,None]*self.v1[2,:,:,:],
                          1j*self.kxgrid[:,None,None]*self.v1[1,:,:,:]-1j*self.kygrid[None,:,None]*self.v1[0,:,:,:]),\
                         axis=0)
        vort = 2*cp.sum(cp.abs(WVORT[:,:,:,1:])**2.0)+cp.sum(cp.abs(WVORT[:,:,:,0])**2.0)

        mhel = cp.real(2*cp.sum((AVP[:,:,:,1:])*cp.conj(self.b1[:,:,:,1:])))
        mhel += cp.real(cp.sum(AVP[:,:,:,0]*cp.conj(self.b1[:,:,:,0])))

        canhel = cp.real(2*cp.sum((AVP[:,:,:,1:]+self.hallparam*self.v1[:,:,:,1:])*cp.conj(self.b1[:,:,:,1:]+self.hallparam*WVORT[:,:,:,1:])))
        canhel += cp.real(cp.sum((AVP[:,:,:,0]+self.hallparam*self.v1[:,:,:,0])*cp.conj(self.b1[:,:,:,0]+self.hallparam*WVORT[:,:,:,0]))) + cp.real(self.hallparam*self.v1[2,0,0,0])

        mhel *= 8 * cp.pi**3
        canhel *= 8 * cp.pi**3
        
        return(mhel,canhel,vort,cb2)

    def tripletamps(self,lwk):

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

        # Write three wave energy file if three wave simulation
        if self.initcond == "threewave":
            
            self.tripletenergies = cp.array([lw1,lc1,rw1,rc1,lw2,lc2,rw2,rc2,lw3,lc3,rw3,rc3])

        return(lw,lc,rw,rc)

    def helicitycorrection(self,binput,vinput):

        correction = cp.real(binput[0,:,:,:]*cp.conj(vinput[1,:,:,:])-binput[1,:,:,:]*cp.conj(vinput[0,:,:,:]))

        mhc = 2.0*cp.sum(correction[:,:,1:])+cp.sum(correction[:,:,0])
        mhc *= -16*cp.pi**3
        
        return(mhc)

    def checkpointfile(self):

        f = open(self.lpath+"/s_checkpoint.dat","wb")

        f.write(self.itime)
        f.write(self.nkx0)
        f.write(self.nky0)
        f.write(self.nkz0)
        f.write(self.time)
        f.write(self.dt)
        f.write(self.b1)
        f.write(self.v1)
        f.write(self.mhelcorr)

        f.close()
        
        return(None)

    def createhdf5file(self):

        if os.path.exists(self.lpath+"/output.hdf5") and (self.itime == 0):

            rmtree(self.lpath)
            os.mkdir(self.lpath)

        kxgridcpu = cp.asnumpy(self.kxgrid)
        kygridcpu = cp.asnumpy(self.kygrid)
        kzgridcpu = cp.asnumpy(self.kzgrid)
        hcpu = cp.asnumpy(self.hallparam)
        
        if not os.path.exists(self.lpath+"/output.hdf5"):

            with h5py.File(self.lpath+"/output.hdf5","a") as f:

                kx = f.create_dataset("kx",(np.size(kxgridcpu),),data=kxgridcpu)
                ky = f.create_dataset("ky",(np.size(kygridcpu),),data=kygridcpu)
                kz = f.create_dataset("kz",(np.size(kzgridcpu),),data=kzgridcpu)
                written = f.create_dataset("written",(1,),dtype="int32",data=int32c(0))
                times = f.create_dataset("time",(self.record+1,),dtype=self.rtype)
                enval = f.create_dataset("enval",(self.record+1,12),dtype=self.rtype)
                #totalenergy = f.create_dataset("hamiltonian",(self.record+1,),dtype=self.rtype)
                #magnetichelicity = f.create_dataset("maghel",(self.record+1,),dtype=self.rtype)
                #crosshel = f.create_dataset("crosshel",(self.record+1,),dtype=self.rtype)
                #kinetic = f.create_dataset("kinenergy",(self.record+1,),dtype=self.rtype)
                #magnetic = f.create_dataset("magenergy",(self.record+1,),dtype=self.rtype)
                #leftwhistler = f.create_dataset("leftwhistler",(self.record+1,),dtype=self.rtype)
                #leftcyclo = f.create_dataset("leftcyclo",(self.record+1,),dtype=self.rtype)
                #rightwhistler = f.create_dataset("rightwhistler",(self.record+1,),dtype=self.rtype)
                #rightcyclo = f.create_dataset("rightcyclo",(self.record+1,),dtype=self.rtype)
                #helicitycorr = f.create_dataset("helicitycorr",(self.record+1,),dtype=self.rtype)
                hall = f.create_dataset("hall",(1,),data=hcpu)

                if self.initcond == "threewave":

                    tripletcpu = cp.asnumpy(self.triplet)
                    threewaves = f.create_dataset("threewaves",(12,),dtype=self.rtype,data=tripletcpu)
                    threewaveenergies = f.create_dataset("threewaveenergies",(self.record+1,12),dtype=self.rtype)

                magneticfields = f.create_dataset("magneticfields",(self.record+1,3,self.nx0_big,self.ny0_big,self.nz0_big//2+1),dtype=self.ctype)
                velocityfields = f.create_dataset("velocityfields",(self.record+1,3,self.nx0_big,self.ny0_big,self.nz0_big//2+1),dtype=self.ctype)

        return(None)
    
    def writefile(self):

        ham = self.hamiltonian(0)
        kinham = self.hamiltonian(1)
        magham = self.hamiltonian(2)
        mh,ch,vort,cb2 = self.helicity()
        lw,lc,rw,rc = self.hallmodeenergies()

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
        if self.initcond == "threewave":
            waveenergiescpu = cp.asnumpy(self.tripletenergies)        
        
        with h5py.File(self.lpath+"/output.hdf5","a") as f:

            dset = f["written"]
            print("Written shape ",self.timecpu)

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
            
            #f["hamiltonian"][self.writtencpu] = hamcpu
            #f["maghel"][self.writtencpu] = mhcpu
            #f["crosshel"][self.writtencpu] = chcpu
            #f["kinenergy"][self.writtencpu] = kinhamcpu
            #f["magenergy"][self.writtencpu] = maghamcpu
            #f["leftwhistler"][self.writtencpu] = lwcpu
            #f["leftcyclo"][self.writtencpu] = lccpu 
            #f["rightwhistler"][self.writtencpu] = rwcpu
            #f["rightcyclo"][self.writtencpu] = rccpu
            #f["helicitycorr"][self.writtencpu] = mhccpu

            f["magneticfields"][self.writtencpu,:,:,:,:] = b1cpu
            f["velocityfields"][self.writtencpu,:,:,:,:] = v1cpu

            if self.initcond == "threewave":
                f["threewaveenergies"][self.writtencpu,:] = waveenergiescpu
                
            #times[self.written] = self.timecpu
            ##totalenergy[self.written] = hamcpu
            
        return(None)

    def energyfile(self):        
                
        if self.itime == 0:

            if os.path.exists(self.lpath):
                rmtree(self.lpath)
                os.mkdir(self.lpath)
            f = open(self.lpath+"/energy.dat","wb")
        else:
            f = open(self.lpath+"energy.dat","ab")

        mh,ch = self.helicity()
        lw,lc,rw,rc = self.hallmodeenergies()

        f.write(self.time)
        f.write(self.hamiltonian(0))
        f.write(mh)
        f.write(ch)
        f.write(self.hamiltonian(1))
        f.write(self.hamiltonian(2))
        f.write(lw)
        f.write(lc)
        f.write(rw)
        f.write(rc)
        f.write(self.mhelcorr)
            
        f.close()
        
        return(None)

    def fulloutputfile(self):

        if self.itime == 0:
            f = open(self.lpath+"/allfields","wb")
        else:
            f = open(self.lpath+"/allfields","ab")

        f.write(self.time)
        f.write(self.b1)
        f.write(self.v1)

        f.close()
        
        return(None)
    
    def steadystate(self):

        self.newspec = cp.sum(cp.abs(self.b1)**2.0 + cp.abs(self.v1)**2.0,axis=0)/2

        ii = cp.nonzero(self.oldspec)

        convmost = cp.amax(cp.abs(cp.log10(self.newspec[ii]/self.oldspec[ii])))
        convleast = cp.amin(cp.abs(cp.log10(self.newspec[ii]/self.oldspec[ii])))

        print("Convergence of Energy Spectrum Itime ",self.itime," Best Log Ratio ",convleast," Worst Log Ratio ",convmost)
        
        self.oldspec = cp.copy(self.newspec)

        return(None)
    
        
class DNA2MHD(RHS,DIAGS):
    def __init__(self,nkx0,nky0,nkz0,kxmin,kymin,kzmin,nu,eta,
                 dt,iterations,lpath,linear=False,explicitrk4=False,
                 initcond="hallwave",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                 forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],hyper=1,hallparam=1.0,
                 solveprec=16,maxwallclock=86200,triplet=None,records=20,bittype=32):

        self.maxwallclock = maxwallclock
        self.linear = linear
        
        self.nkx0 = int32(nkx0)
        self.nky0 = int32(nky0)
        self.nkz0 = int32(nkz0)

        self.nx0_big = 3 * nkx0//2
        self.ny0_big = 3 * nky0//2
        self.nz0_big = 3 * nkz0//2

        if bittype == 32:
            floattype = cp.float32
            complextype = cp.complex64
        else:
            floattype = cp.float64
            complextype = cp.complex128

        self.ctype = "complex"+str(2*bittype)
        self.rtype = "float"+str(bittype)
                    
        self.kxmin = floattype(kxmin)
        self.kymin = floattype(kymin)
        self.kzmin = floattype(kzmin)

        self.explicitrk4 = explicitrk4

        self.initcond = initcond
        self.energystart = energystart
        self.init_kolm = floattype(init_kolm)
        self.hmhdwave = hmhdwave
        
        self.forcetype = forcetype
        self.forceamp = floattype(forceamp)
        self.nforce = nforce
        self.forcewave = forcewave
        
        self.hallparam = hallparam

        if lpath[-1] == "/":
            lpath = lpath[:-1]
        self.lpath = lpath
        
        self.nu = floattype(nu) # Dissipation at largest scale
        self.eta = floattype(eta) # Magnetic Prandtl number
        self.hyper = hyper

        self.dt = floattype(dt)
        self.iterations = iterations
        self.itime = int32(0)
        self.itime_start = int32(0)
        self.time = floattype(0)
        self.solveprec = solveprec

        self.kxgrid = cp.concatenate((cp.arange(0,self.nx0_big//2,dtype=self.rtype),cp.arange(-self.nx0_big//2,0,dtype=self.rtype)))*kxmin
        self.kygrid = cp.concatenate((cp.arange(0,self.ny0_big//2,dtype=self.rtype),cp.arange(-self.ny0_big//2,0,dtype=self.rtype)))*kymin
        self.kzgrid = cp.arange(0,self.nz0_big//2+1,dtype=self.rtype)*kzmin 

        KX,KY,KZ = cp.meshgrid(self.kxgrid,self.kygrid,self.kzgrid,indexing="ij")
        
        self.kmags = cp.sqrt(self.kxgrid[:,None,None]**2.0 + self.kygrid[None,:,None]**2.0+self.kzgrid[None,None,:]**2.0).astype(self.rtype)
        self.kmax = cp.amax(self.kmags)

        self.vnu = nu/(self.kmax**(2*hyper))
        self.etab = eta*self.vnu
        print("Dissipation Factors",self.vnu,self.etab)

        self.alphaleftwhist = - (hallparam*self.kmags/2 + cp.sqrt(1+ (hallparam*self.kmags/2)**2.0))
        self.triplet = triplet

        eigb2 = - ((self.vnu+self.etab)*self.kmags**(2*hyper) + 1j * self.kmags * self.kzgrid[None,None,:]*hallparam)/2
        eigsqrt = eigb2**2 - self.kzgrid[None,None,:]**2.0 - self.vnu * self.etab * self.kmags**(4*hyper) \
            - 1j * self.vnu * self.kmags**(2*hyper) * self.kzgrid[None,None,:] * self.kmags
        
        self.eig = cp.zeros([self.nx0_big,self.ny0_big,self.nz0_big//2+1,4],dtype=self.ctype)
        self.eig[:,:,:,0] = eigb2 - cp.sqrt(eigsqrt)
        self.eig[:,:,:,1] = eigb2 + cp.sqrt(eigsqrt)
        self.eig[:,:,:,2] = cp.conj(eigb2) + cp.sqrt(cp.conj(eigsqrt))
        self.eig[:,:,:,3] = cp.conj(eigb2) - cp.sqrt(cp.conj(eigsqrt))
        
        zvec = cp.array([0,0,1],dtype=self.ctype)
        ks = cp.stack((KX,KY,KZ)).astype(self.ctype)

        self.pcurleig = cp.cross(ks,zvec[:,None,None,None],axis=0)
        self.pcurleig += 1j* cp.cross(ks,self.pcurleig,axis=0)/self.kmags[None,:,:,:]
        self.pcurleig *= cp.sqrt(2)/(2*self.kmags[:,:,0][None,:,:,None])
        self.pcurleig[:,0,0,:] = cp.array([1,1j,0],dtype=self.ctype)[:,None]/cp.sqrt(2.0)
        self.pcurleig[:,0,0,0] = 0.0

        self.padding = cp.ones_like(self.kmags,dtype="int32")
        self.padding[self.nkx0//2:1-self.nkx0//2,:,:] = 0
        self.padding[:,self.nky0//2:1-self.nky0//2,:] = 0
        self.padding[:,:,self.nkz0:] = 0

        self.forcemask = cp.zeros_like(self.padding,dtype=self.rtype)
        self.forcemask[1:nforce+1,1:nforce+1,1:nforce+1] = 1
        self.forcemask[-nforce:,-nforce:,1:nforce+1] = 1
        self.forcemask[1:nforce+1,-nforce:,1:nforce+1] = 1
        self.forcemask[-nforce:,1:nforce+1,1:nforce+1] = 1
        self.forcemask *= forceamp
        
        self.b1 = cp.zeros_like(self.pcurleig)
        self.v1 = cp.zeros_like(self.pcurleig)

        self.mhelcorr = cp.float32(0.0)
        
        self.fieldsetup()

        self.oldspec = cp.sum(cp.abs(self.b1)**2.0+cp.abs(self.v1)**2.0,axis=0)/2
        self.newspec = 2*cp.ones_like(self.kmags)

        self.record = records
        self.irecord = self.iterations // records
        self.written = 0
        self.createhdf5file()
        
        return(None)

    def fieldsetup(self):

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

        if self.initcond == "threewave":

            self.itriplet = []

            if self.triplet == None or  (len(self.triplet) != 9 and len(self.triplet) != 12):
                
                raise ValueError("LW Triplet or triplet and normal mode types must be specified for three wave initial condition")

            # Specify triplet as either list of three wavevectors or list of three wave vectors and normal mode types
            # Wave triplet - 1 + whistler, 0 + cyclotron, 2 - whistler, 3 - cyclotron
            
            if len(self.triplet) == 9:
                for i in range(3):
                    self.triplet.append(0)
                self.initializewave(self.triplet[0:3],2.0,0)
                self.initializewave(self.triplet[3:6],1.0,0)
                self.initializewave(self.triplet[6:9],0.5,0)

            else:
                self.initializewave(self.triplet[0:3],2.0,self.triplet[9])
                self.initializewave(self.triplet[3:6],1.0,self.triplet[10])
                self.initializewave(self.triplet[6:9],0.5,self.triplet[11])

            self.triplet = cp.array(self.triplet)

        if self.energystart != None:
            energy = cp.sum(cp.abs(self.b1[:,:,:,0])**2.0+cp.abs(self.v1[:,:,:,0])**2.0)
            energy += 2* cp.sum(cp.abs(self.b1[:,:,:,1:])**2.0+cp.abs(self.v1[:,:,:,1:])**2.0)

            self.b1 *= cp.sqrt(self.energystart/energy)
            self.v1 *= cp.sqrt(self.energystart/energy)
                
        if self.initcond == "checkpoint":

            with h5py.File(self.lpath+"/output.hdf5","r") as f:

                dset = f["written"]
                print("Written shape ",self.timecpu)

                written = f["written"][0]
                b1cpu = f["magneticfields"][written,:,:,:,:]
                v1cpu = f["velocityfields"][written,:,:,:,:]
                timecpu = f["time"][written]
                dt = f["time"][written]-f["time"][written-1]
                mhc = f["enval"][written,-1]

            self.b1 = cp.asarray(b1cpu)
            self.v1 = cp.asarray(v1cpu)
            self.time = cp.asarray(timecpu)
            self.itime = cp.asarray(timecpu/dt)
            self.itime_start = cp.copy(self.itime)
            self.dt = cp.asarray(dt)
            self.mhelcorr = cp.asarray(mhc)

        return(None)

    def initializewave(self,wave,amp,modeno=0):

        # Set up normal mode type modeno for index ix, iy, iz with relative amplitude amp

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
        self.b1[:,ix,iy,iz] *= amp/cp.sqrt((self.alphaleftwhist[ix,iy,iz] ** (1 - 2 * (modeno % 2)))**2 + 1)
        self.v1[:,ix,iy,iz] *= amp/cp.sqrt((self.alphaleftwhist[ix,iy,iz] ** (1 - 2 * (modeno % 2)))**2 + 1)

        return(None)

    def etdrk2(self,records=200):
    
        self.startruntime = time.time()
        self.runtime = 0

        debugnote = """"Zero eigenvalue coefficients for reference"""
        expL = cp.ones([self.nx0_big,self.ny0_big,self.nz0_big//2+1,8])
        expL[:,:,:,1:3] = 0.0
        expL[:,:,:,5:7] = 0.0
        coef1 = expL * self.dt
        coef2 = coef1/2

        expL,coef1,coef2 = exponentialcoefficients(self.kxgrid,self.kygrid,self.kzgrid,self.etab,self.vnu,self.hyper,self.dt,128,self.kmags,self.eig,self.hallparam,
                                                   pade=True,ctype=self.ctype)
        #blurb = Have to adjust the zero mode separately because the curl eigenstates are undefined - set coefficients there to zero
        expL[0,0,0,:] = cp.array([1.0,0.0,0.0,1.0,1.0,0.0,0.0,1.0])
        coef1[0,0,0,:] = expL[0,0,0,:]*self.dt
        coef2[0,0,0,:] = coef1[0,0,0,:]/2.0

        print("Coef1 ",cp.amax(cp.abs(coef1/self.dt)))
        print("Coef2 ",cp.amax(cp.abs(coef2/self.dt)))
        
        while self.itime < self.iterations and self.runtime < self.maxwallclock:

            #dot1 = time.time()
            
            self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
            self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
            self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
            self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

            #dot1 -= time.time()
            #print("dot time ", np.abs(dot1))

            if self.irecord > 0 and self.itime % self.irecord == 0:

                self.writefile()
                #self.energyfile()
                #self.checkpointfile()
                #self.fulloutputfile()

                self.written += 1
                
            avgb = self.b1[:,0,0,0].copy()
            avgv = self.v1[:,0,0,0].copy()

            #rhs1time = time.time()
            self.brhs1,self.vrhs1 = self.convolutions(self.b1,self.v1)
            #rhs1time -= time.time()
            #print("rhs1 time ", np.abs(rhs1time))
            mhc = self.helicitycorrection(self.b1,self.v1)
            
            #expmult = time.time()
            self.b1 = self.pcurleig * (expL[:,:,:,0] * self.bplus + expL[:,:,:,1] * self.vplus)
            self.v1 = self.pcurleig * (expL[:,:,:,2] * self.bplus + expL[:,:,:,3] * self.vplus)
            self.b1 += cp.conj(self.pcurleig) * (expL[:,:,:,4] * self.bminus + expL[:,:,:,5] * self.vminus)
            self.v1 += cp.conj(self.pcurleig) * (expL[:,:,:,6] * self.bminus + expL[:,:,:,7] * self.vminus)
            self.b1[:,0,0,0] = avgb
            self.v1[:,0,0,0] = avgv
            #expmult -= time.time()
            #print("expmult time", np.abs(expmult))

            #coef1time = time.time()
            self.bplus = cp.sum(cp.conj(self.pcurleig)*self.brhs1,axis=0)
            self.vplus = cp.sum(cp.conj(self.pcurleig)*self.vrhs1,axis=0)
            self.bminus = cp.sum(self.pcurleig*self.brhs1,axis=0)
            self.vminus = cp.sum(self.pcurleig*self.vrhs1,axis=0)

            self.b1 += self.pcurleig * (coef1[:,:,:,0] * self.bplus + coef1[:,:,:,1] * self.vplus)
            self.v1 += self.pcurleig * (coef1[:,:,:,2] * self.bplus + coef1[:,:,:,3] * self.vplus)
            self.b1 += cp.conj(self.pcurleig) * (coef1[:,:,:,4] * self.bminus + coef1[:,:,:,5] * self.vminus)
            self.v1 += cp.conj(self.pcurleig) * (coef1[:,:,:,6] * self.bminus + coef1[:,:,:,7] * self.vminus)

            #coef1time -= time.time()
            #print("coef1 time", np.abs(coef1time))

            self.b1 += self.dt * self.brhs1
            self.v1 += self.dt * self.vrhs1

            b1zeroderiv = self.brhs1[:,0,0,0].copy()
            v1zeroderiv = self.vrhs1[:,0,0,0].copy()
            self.mhelcorr += self.dt * mhc
            mhc1 = mhc

            self.b1[:,0,0,0] += self.dt*self.brhs1[:,0,0,0]
            self.v1[:,0,0,0] += self.dt*self.vrhs1[:,0,0,0]

            self.brhs1,self.vrhs1 = self.convolutions(self.b1,self.v1)
            mhc = self.helicitycorrection(self.b1,self.v1)

            self.bplus -= cp.sum(cp.conj(self.pcurleig)*(self.brhs1),axis=0) 
            self.vplus -= cp.sum(cp.conj(self.pcurleig)*(self.vrhs1),axis=0)
            self.bminus -= cp.sum(self.pcurleig*self.brhs1,axis=0)
            self.vminus -= cp.sum(self.pcurleig*self.vrhs1,axis=0)
                            
            # Now self.bplus etc equal to rhs1- rhs2 in basis, reverse sign of sum below

            #coef2time = time.time()
            self.b1 -= self.pcurleig * (coef2[:,:,:,0] * self.bplus + coef2[:,:,:,1] * self.vplus)
            self.v1 -= self.pcurleig * (coef2[:,:,:,2] * self.bplus + coef2[:,:,:,3] * self.vplus)
            self.b1 -= cp.conj(self.pcurleig) * (coef2[:,:,:,4] * self.bminus + coef2[:,:,:,5] * self.vminus)
            self.v1 -= cp.conj(self.pcurleig) * (coef2[:,:,:,6] * self.bminus + coef2[:,:,:,7] * self.vminus)
            self.b1[:,0,0,0] += self.dt/2 * (self.brhs1[:,0,0,0]-b1zeroderiv)
            self.v1[:,0,0,0] += self.dt/2 * (self.vrhs1[:,0,0,0]-v1zeroderiv)
            self.mhelcorr += self.dt/2 * (mhc-mhc1)
            #coef2time -= time.time()
            #print("coef2 time",np.abs(coef2time))
            
            self.itime += 1
            self.time += self.dt
            self.runtime = time.time()-self.startruntime

            #if self.itime % max(self.iterations//50,10) == 0 and self.itime > 0:
            #    self.steadystate()

        self.writefile()
        #self.energyfile()
        #self.checkpointfile()
        #self.fulloutputfile()
        print(self.runtime)

        return(None)

    def ralstonrk2(self,records=200):

        self.startruntime = time.time()
        self.runtime = 0

        # Add guide field to do linear and nonlinear calculation at once
        self.b1[2,0,0,0] = 1.0

        while self.itime < self.iterations and self.runtime < self.maxwallclock:

            if self.irecord > 0 and self.itime % self.irecord == 0:

                self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
                self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
                self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
                self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

                self.writefile()
                self.written += 1

            self.brhs1,self.vrhs1,mhc = self.hallrhs(self.b1,self.v1)
            #mhc = self.helicitycorrection(self.b1,self.v1)

            self.b2 = self.b1 + 2/3 * self.dt * self.brhs1
            self.v2 = self.v1 + 2/3 * self.dt * self.vrhs1

            self.b1 += 1/4 * self.dt * self.brhs1
            self.v1 += 1/4 * self.dt * self.vrhs1
            self.mhelcorr += 1/4 * self.dt * mhc
            
            self.brhs1,self.vrhs1,mhc = self.hallrhs(self.b2,self.v2)
            #mhc = self.helicitycorrection(self.b1,self.v2)
            
            self.b1 += 3/4 * self.dt * self.brhs1
            self.v1 += 3/4 * self.dt * self.vrhs1
            self.mhelcorr += 3/4 * self.dt * mhc

            self.itime += 1
            self.time += self.dt
            self.runtime = time.time()-self.startruntime

        self.writefile()
        print(self.runtime)

    def dp547s(self):

        """5th order stability region Dormand Prince scheme"""
        """https://doi.org/10.1016/0771-050X(80)90013-3"""

        self.startruntime = time.time()
        self.runtime = 0

        # Add guide field to do linear and nonlinear calculation at once
        self.b1[2,0,0,0] = 1.0

        while self.itime < self.iterations and self.runtime < self.maxwallclock:

            if self.irecord > 0 and self.itime % self.irecord == 0:
                self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
                self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
                self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
                self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

                self.writefile()
                self.written += 1

            self.brhs1,self.vrhs1,mhc1 = self.hallrhs(self.b1,self.v1)
            mhc1 = self.helicitycorrection(self.b1,self.v1)

            dumb = self.b1 + self.dt * 2/9 * self.brhs1
            dumv = self.v1 + self.dt * 2/9 * self.vrhs1
            
            self.brhs2,self.vrhs2,mhc2 = self.hallrhs(dumb,dumv)

            dumb = self.b1 + self.dt * (1/12 * self.brhs1 + 1/4 * self.brhs2)
            dumv = self.v1 + self.dt * (1/12 * self.vrhs1 + 1/4 * self.vrhs2)

            self.brhs3,self.vrhs3,mhc3 = self.hallrhs(dumb,dumv)

            dumb = self.b1 + self.dt * (55/324 * self.brhs1 - 25/108 * self.brhs2 + 50/81 * self.brhs3)
            dumv = self.v1 + self.dt * (55/324 * self.vrhs1 - 25/108 * self.vrhs2 + 50/81 * self.vrhs3)

            self.brhs4,self.vrhs4,mhc4 = self.hallrhs(dumb,dumv)

            dumb = self.b1 + self.dt * (83/330 * self.brhs1 - 13/22 * self.brhs2 + 61/66 * self.brhs3 + 9/110 * self.brhs4)
            dumv = self.v1 + self.dt * (83/330 * self.vrhs1 - 13/22 * self.vrhs2 + 61/66 * self.vrhs3 + 9/110 * self.vrhs4)

            self.brhs5,self.vrhs5,mhc5 = self.hallrhs(dumb,dumv)

            dumb = self.b1 + self.dt * (-19/28 * self.brhs1 + 9/4 * self.brhs2 + 1/7 * self.brhs3 + -27/7 * self.brhs4 + 22/7 * self.brhs5)
            dumv = self.v1 + self.dt * (-19/28 * self.vrhs1 + 9/4 * self.vrhs2 + 1/7 * self.vrhs3 + -27/7 * self.vrhs4 + 22/7 * self.vrhs5)

            self.brhs6,self.vrhs6,mhc6 = self.hallrhs(dumb,dumv)

            self.b1 += (19/200 * self.brhs1 + 3/5 * self.brhs3 - 243/400 * self.brhs4 + 33/40 * self.brhs5 + 7/80  * self.brhs6) * self.dt
            self.v1 += (19/200 * self.vrhs1 + 3/5 * self.vrhs3 - 243/400 * self.vrhs4 + 33/40 * self.vrhs5 + 7/80  * self.vrhs6) * self.dt 
            self.mhelcorr += (19/200 * mhc1 + 3/5 * mhc3 -243/400 * mhc4 + 33/40 * mhc5 + 7/80 * mhc6) * self.dt

            self.itime += 1
            self.time += self.dt
            self.runtime = time.time()-self.startruntime

        self.writefile()
        print(self.runtime)

        return(None)
            
    def gauss2(self): # What we started with - splitting method for dissipation, implicit midpoint nonlinear

        self.startruntime = time.time()
        self.runtime = 0

        # Add guide field to do linear and nonlinear calculation at once                                                                                                                   
        self.b1[2,0,0,0] = 1.0

        a11 = 0.5

        while self.itime < self.iterations and self.runtime < self.maxwallclock: 

            if self.irecord > 0 and self.itime % self.irecord == 0:
                self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
                self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
                self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
                self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

                self.writefile()
                self.written += 1

            self.b1 *= cp.exp(-self.etab * self.kmags[None,:,:,:]**(2*self.hallparam) * self.dt/2)
            self.v1 *= cp.exp(-self.vnu * self.kmags[None,:,:,:]**(2*self.hallparam) * self.dt/2)
            
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

            self.b1 *= cp.exp(-self.etab * self.kmags**(2*self.hallparam) * self.dt/2)
            self.v1 *= cp.exp(-self.vnu * self.kmags**(2*self.hallparam) * self.dt/2)

            self.itime += 1
            self.time += self.dt
            self.runtime = time.time()-self.startruntime

        self.writefile()
        print(self.runtime)

        return(None)
