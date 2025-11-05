import cupy as cp
from cupy.fft import rfftn,irfftn
import time
from dna2mhd_expintegrator import exponentialcoefficients
import os
from shutil import rmtree

class RHS:

    def hallrhs(self,binput,vinput):

        if self.linear:
            bout,vout = self.rhslinear(binput,vinput)
        else:
            bout,vout = self.convolutions(binput,vinput)
        bout,vout = self.removediv(bout,vout)

        return(bout,vout)

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

        bout = cp.zeros_like(self.b1)
        vout = cp.zeros_like(self.v1)

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

        return(cp.exp(1j*cp.pi*2*cp.random.rand(self.nx0_big,self.ny0_big,self.nz0_big//2 + 1)))
    
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
                       axis=0)/(self.kmags[None,:,:,:]**2.0)
        AVP[:,0,0,0] = 0.0

        WVORT = cp.stack((1j*self.kygrid[None,:,None]*self.v1[2,:,:,:]-1j*self.kzgrid[None,None,:]*self.v1[1,:,:,:],
                          1j*self.kzgrid[None,None,:]*self.v1[0,:,:,:]-1j*self.kxgrid[:,None,None]*self.v1[2,:,:,:],
                          1j*self.kxgrid[:,None,None]*self.v1[1,:,:,:]-1j*self.kygrid[None,:,None]*self.v1[0,:,:,:]),\
                         axis=0)

        mhel = cp.real(2*cp.sum((AVP[:,:,:,1:])*cp.conj(self.b1[:,:,:,1:])))
        mhel += cp.real(cp.sum(AVP[:,:,:,0]*cp.conj(self.b1[:,:,:,0])))

        canhel = cp.real(2*cp.sum((AVP[:,:,:,1:]+self.hallparam*self.v1[:,:,:,1:])*cp.conj(self.b1[:,:,:,1:]+self.hallparam*WVORT[:,:,:,1:])))
        canhel += cp.real(cp.sum((AVP[:,:,:,0]+self.hallparam*self.v1[:,:,:,0])*cp.conj(self.b1[:,:,:,0]+self.hallparam*WVORT[:,:,:,0]))) + cp.real(self.hallparam*self.v1[2,0,0,0])

        mhel *= 8 * cp.pi**3
        canhel *= 8 * cp.pi**3
        
        return(mhel,canhel)

    def tripletamps(self,lwk):

        lw1 = 8.0*cp.pi**3 * cp.abs(lwk[self.itriplet[0],self.itriplet[1],self.itriplet[2]])**2.0
        lw2 = 8.0*cp.pi**3 * cp.abs(lwk[self.itriplet[1],self.itriplet[4],self.itriplet[5]])**2.0
        lw3 = 8.0*cp.pi**3 * cp.abs(lwk[self.itriplet[2],self.itriplet[7],self.itriplet[8]])**2.0

        return(lw1,lw2,lw3)

    def hallmodeenergies(self):

        # Dot product

        det = ( (self.vnu * self.kmags**(2*self.hyper)+self.eig[:,:,:,0])* (self.vnu * self.kmags**(2*self.hyper) + self.eig[:,:,:,1])+ kzgrid[None,None,:]**2.0 )
        
        lwk = (1j * self.kzgrid[None,None,:] * self.bplus - (self.eig[:,:,:,1] + self.vnu * self.kmags**(2*self.hyper)) * self.vplus)/det 
        lwk *= cp.sqrt(cp.abs(self.eig[:,:,:,0] + self.vnu * self.kmags**(2*self.hyper))**2.0 + self.kzgrid[None,None,:]**2.0)
        lw = 2.0*cp.sum(cp.abs(lwk[:,:,1:])**2.0) + cp.sum(cp.abs(lwk[:,:,0])**2.0)        
        lw *= (4.0*cp.pi**3)

        if self.initcond == "threewave":
            lw1,lw2,lw3 = self.tripletamps(lwk)
            
        lwk = (-1j * self.kzgrid[None,None,:] * self.bplus + (self.eig[:,:,:,0] + self.vnu * self.kmags**(2*self.hyper)) * self.vplus)/det
        lwk *= cp.sqrt(cp.abs(self.eig[:,:,:,1] + self.vnu * self.kmags**(2*self.hyper))**2.0 + self.kzgrid[None,None,:]**2.0)
        lc = 2.0 * cp.sum(cp.abs(lwk[:,:,1:])**2.0) + cp.sum(cp.abs(lwk[:,:,0])**2.0)
        lc *= 4.0 * cp.pi**3
        if self.initcond == "threewave":
            lc1,lc2,lc3 = self.tripletamps(lwk)

        det = ( (self.vnu * self.kmags**(2*self.hyper)+self.eig[:,:,:,2])* (self.vnu * self.kmags**(2*self.hyper) + self.eig[:,:,:,3])+ kzgrid[None,None,:]**2.0 )
        
        lwk = (1j * self.kzgrid[None,None,:] * self.bminus - (self.eig[:,:,:,3] + self.vnu * self.kmags**(2*self.hyper)) * self.vminus)/det
        lwk *= cp.sqrt(cp.abs(self.eig[:,:,:,2] + self.vnu * self.kmags**(2*self.hyper))**2.0 + self.kzgrid[None,None,:]**2.0)
        rw = 2.0 * cp.sum(cp.abs(lwk[:,:,1:])**2.0) + cp.sum(cp.abs(lwk[:,:,0])**2.0)
        rw *= 4.0 * cp.pi**3
        if self.initcond == "threewave":
            rw1,rw2,rw3 = self.tripletamps(lwk)

        lwk = (-1j * self.kzgrid[None,None,:] * self.bminus + (self.eig[:,:,:,2] + self.vnu * self.kmags**(2*self.hyper)) * self.vminus)/det
        lwk *= cp.sqrt(cp.abs(self.eig[:,:,:,3] + self.vnu * self.kmags**(2*self.hyper))**2.0 + self.kzgrid[None,None,:]**2.0)
        rc = 2.0 * cp.sum(cp.abs(lwk[:,:,1:])**2.0) + cp.sum(cp.abs(lwk[:,:,0])**2.0)
        rc *= 4.0 * cp.pi**3
        if self.initcond == "threewave":
            rc1,rc2,rc3 = self.tripletamps(lwk)

        # Write three wave energy file if three wave simulation
        if self.initcond == "threewave":

            if self.itime == 0:
                f = open(self.lpath+"/threewave.dat","wb")
            else:
                f = open(self.lpath+"/threewave.dat","ab")

            f.write(self.itime)
            f.write(lw1)
            f.write(lc1)
            f.write(rw1)
            f.write(rc1)
            f.write(lw2)
            f.write(lc2)
            f.write(rw2)
            f.write(rc2)
            f.write(lw3)
            f.write(lc3)
            f.write(rw3)
            f.write(rc3)
            f.close()

        return(lw,lc,rw,rc)

    def helicitycorrection(self):

        correction = cp.real(self.b2[0,:,:,:]*cp.conj(self.v2[1,:,:,:])-self.b2[1,:,:,:]*cp.conj(self.v2[0,:,:,:]))

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
    def __init__(self,nkx0,nky0,nkz0,kxmin,kymin,kzmin,
                 nu,eta,
                 dt,iterations,
                 lpath,linear=False,explicitrk4=False,
                 initialcondition="hallwave",energystart=0.01,init_kolm=0,hmhdwave=[1,0,0,0],
                 forcetype="hallwave",forceamp=0.0,nforce=4,forcewave=[1,0,0,0],
                 hyper=1,hallparam=1.0,
                 solveprec=16,maxwallclock=86200,
                 triplet=None):

        self.maxwallclock = maxwallclock
        self.linear = linear
        
        self.nkx0 = nkx0
        self.nky0 = nky0
        self.nkz0 = nkz0

        self.nx0_big = 3 * nkx0//2
        self.ny0_big = 3 * nky0//2
        self.nz0_big = 3 * nkz0
        
        self.kxmin = kxmin
        self.kymin = kymin
        self.kzmin = kzmin

        self.explicitrk4 = explicitrk4

        self.initialcondition = initialcondition
        self.energystart = energystart
        self.init_kolm = init_kolm
        self.hmhdwave = hmhdwave
        
        self.forcetype = forcetype
        self.forceamp = forceamp
        self.nforce = nforce
        self.forcewave = forcewave
        
        self.hallparam = hallparam

        if lpath[-1] == "/":
            lpath = lpath[:-1]
        self.lpath = lpath
        
        self.nu = nu # Dissipation at largest scale
        self.eta = eta # Magnetic Prandtl number
        self.hyper = hyper

        self.dt = dt
        self.iterations = iterations
        self.itime = 0
        self.itime_start = 0
        self.time = 0
        self.solveprec = solveprec

        self.kxgrid = cp.hstack((cp.arange(0,self.nx0_big//2,dtype="float64"),cp.arange(-self.nx0_big//2,0,dtype="float64")))*kxmin
        self.kygrid = cp.hstack((cp.arange(0,self.ny0_big//2,dtype="float64"),cp.arange(-self.ny0_big//2,0,dtype="float64")))*kymin
        self.kzgrid = cp.arange(0,self.nz0_big//2+1,dtype="float64")*kzmin        

        KX,KY,KZ = cp.meshgrid(self.kxgrid,self.kygrid,self.kzgrid,indexing="ij")
        
        self.kmags = cp.sqrt(self.kxgrid[:,None,None]**2.0 + self.kygrid[None,:,None]**2.0+self.kzgrid[None,None,:]**2.0)
        self.kmax = cp.amax(self.kmags)

        self.vnu = nu/(self.kmax**(2*hyper))
        self.etab = eta*self.vnu
        print("Dissipation Factors",self.vnu,self.etab)

        self.alphaleftwhist = - (hallparam*self.kmags/2 + cp.sqrt(1+ (hallparam*self.kmags/2)**2.0))
        self.triplet = triplet

        self.eig[:,:,:,0] = - ((self.vnu+self.etab)*self.kmags**(2*hyper) + 1j * self.kmags * self.kzgrid[None,None,:])/2
        self.eig[:,:,:,2] = - ((self.vnu+self.etab)*self.kmags**(2*hyper) - 1j * self.kmags * self.kzgrid[None,None,:])/2
        self.eig[:,:,:,1] = cp.sqrt(self.eig[:,:,:,0]**2 - self.vnu * self.etab * self.kmags**(4*hyper) - self.kzgrid[None,None,:]**2 - 1j * self.kmags * self.kzgrid[None,None,:] * self.vnu * self.kmags**(2*hyper))
        self.eig[:,:,:,3] = cp.sqrt(self.eig[:,:,:,2]**2 - self.vnu * self.etab * self.kmags**(4*hyper) - self.kzgrid[None,None,:]**2 + 1j * self.kmags * self.kzgrid[None,None,:] * self.vnu * self.kmags**(2*hyper))
        
        self.eig[:,:,:,0] += self.eig[:,:,:,1]
        self.eig[:,:,:,2] += self.eig[:,:,:,3]
        self.eig[:,:,:,1] *= -2.0
        self.eig[:,:,:,3] *= -2.0
        self.eig[:,:,:,1] += self.eig[:,:,:,0]
        self.eig[:,:,:,3] += self.eig[:,:,:,2]


        zvec = cp.array([0,0,1],dtype="complex128")
        ks = cp.stack((KX,KY,KZ)).astype("complex128")

        self.pcurleig = cp.cross(ks,zvec[:,None,None,None],axis=0)
        self.pcurleig += 1j* cp.cross(ks,self.pcurleig,axis=0)/self.kmags[None,:,:,:]
        self.pcurleig *= cp.sqrt(2)/(2*self.kmags[:,:,0][None,:,:,None])
        self.pcurleig[:,0,0,:] = cp.array([1,1j,0],dtype="complex128")[:,None]

        self.padding = cp.ones_like(self.kmags,dtype="int32")
        self.padding[self.nkx0//2:1-self.nkx0//2,:,:] = 0
        self.padding[:,self.nky0//2:1-self.nky0//2,:] = 0
        self.padding[:,:,self.nkz0:] = 0

        self.forcemask = cp.zeros_like(self.padding,dtype="float64")
        self.forcemask[1:nforce+1,1:nforce+1,1:nforce+1] = 1
        self.forcemask[-nforce:,-nforce:,1:nforce+1] = 1
        self.forcemask[1:nforce+1,-nforce:,1:nforce+1] = 1
        self.forcemask[-nforce:,1:nforce+1,1:nforce+1] = 1
        self.forcemask *= forceamp
        
        self.b1 = cp.zeros_like(self.pcurleig)
        self.v1 = cp.zeros_like(self.pcurleig)

        self.mhelcorr = cp.float64(0.0)
        
        self.fieldsetup()

        self.oldspec = cp.sum(cp.abs(self.b1)**2.0+cp.abs(self.v1)**2.0,axis=0)/2
        self.newspec = 2*cp.ones_like(self.kmags)

        return(None)

    def fieldsetup(self):

        if self.initialcondition == "hallwave":
            
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

        if self.initialcondition == "threewave":

            self.itriplet = []

            if self.triplet == None or  (len(self.triplet) != 9 and len(self.triplet) != 12):
                
                raise ValueError("Triplet must be specified for three wave initial condition")

            # Specify triplet as either list of three wavevectors or list of three wave vectors and normal mode types

            if len(self.triplet) == 9:

                self.initializewave(self.triplet[0:3],2.0,0)
                self.initializewave(self.triplet[3:6],1.0,0)
                self.initializewave(self.triplet[6:9],0.5,0)

            else:
                self.initializewave(self.triplet[0:3],2.0,self.triplet[9])
                self.initializewave(self.triplet[3:6],1.0,self.triplet[10])
                self.initializewave(self.triplet[6:9],0.5,self.triplet[11])

        if self.energystart != None:
            energy = cp.sum(cp.abs(self.b1[:,:,:,0])**2.0+cp.abs(self.v1[:,:,:,0])**2.0)
            energy += 2* cp.sum(cp.abs(self.b1[:,:,:,1:])**2.0+cp.abs(self.v1[:,:,:,1:])**2.0)

            self.b1 *= cp.sqrt(self.energystart/energy)
            self.v1 *= cp.sqrt(self.energystart/energy)
                
        if self.initialcondition == "checkpoint":

            pastresults = cp.load(self.lpath+"/s_checkpoint.npz")
            if self.nkx0 != pastresults["nkx0"] or self.nky0 != pastresults["nky0"] or self.nkz0 != pastresults["nkz0"]:
                print("Loaded run has different size than desired simulation")
                quit()
            self.itime = pastresults["itime"]
            self.itime_start = pastresults["itime"]
            self.time = pastresults["time"]
            self.dt = pastresults["dt"]
            self.b1 = pastresults["b1"]
            self.v1 = pastresults["v1"]
            self.mhelcorr = pastresults["mhelcorr"]

        return(None)

    def initializewave(self,wave,amp,modeno=0):

        # Set up normal mode type modeno for index ix, iy, iz with relative amplitude amp

        ix = wave[0]//self.kxmin
        iy = wave[1]//self.kymin
        iz = wave[2]//self.kzmin

        self.itriplet.append(ix)
        self.itriplet.append(iy)
        self.itriplet.append(iz)
        
        if modeno < 2:

            self.b1[:,ix,iy,iz] = self.pcurleig[:,ix,iy,iz]
            self.v1[:,ix,iy,iz] = self.pcurleig[:,ix,iy,iz]

        else:

            self.b1[:,ix,iy,iz] = cp.conj(self.pcurleig[:,ix,iy,iz])
            self.v1[:,ix,iy,iz] = cp.conj(self.pcurleig[:,ix,iy,iz])

        self.b1[:,ix,iy,iz] *= (self.vnu*self.kmags[ix,iy,iz]**(2*self.hyper)+self.eig[ix,iy,iz,modeno])
        self.v1[:,ix,iy,iz] *= 1j * self.kzgrid[iz]

        # Normalize amplitude of wave to amp
        length = cp.sum(cp.abs(self.b1[:,ix,iy,iz])**2.0 + cp.abs(self.v1[:,ix,iy,iz])**2.0)
        self.b1 *= amp/length
        self.v1 *= amp/length

        # Include complex conjugate

        self.b1[:,-ix,-iy,iz] = cp.conj(self.b1[:,ix,iy,iz])
        self.v1[:,-ix,-iy,iz] = cp.conj(self.v1[:,ix,iy,iz])

        return(None)

    def etdrk2(self,records=200):
    
        self.startruntime = time.time()
        self.runtime = 0

        irecord = self.iteration // records

        expL,coef1,coef2 = exponentialcoefficients(self.kxgrid,self.kygrid,self.kzgrid,self.eta,self.nu,self.hyp,self.dt,32)

        # Have to adjust the zero mode separately because the curl eigenstates are undefined
        # So set the exponential method impacts to zero mode to be zero
        expL[0,0,0,:] = 0.0
        coef1[0,0,0,:] = 0.0
        coef2[0,0,0,:] = 0.0

        while self.itime < self.iterations and self.runtime < self.maxwallclock:

            self.bplus = cp.sum(cp.conj(self.pcurleig)*self.b1,axis=0)
            self.vplus = cp.sum(cp.conj(self.pcurleig)*self.v1,axis=0)
            self.bminus = cp.sum(self.pcurleig*self.b1,axis=0)
            self.vminus = cp.sum(self.pcurleig*self.v1,axis=0)

            if self.itime % irecord == 0:
                self.energyfile()
                self.checkpointfile()
                self.fulloutputfile()

            avgb = self.b1[:,0,0,0]
            avgv = self.v1[:,0,0,0]
            
            self.brhs1,self.vrhs1 = self.hallrhs(self.b1,self.v1)
            mhc = self.helicitycorrection()
            
            self.b1 = self.pcurleig * (expL[:,:,:,0] * self.bplus + expL[:,:,:,1] * self.vplus)
            self.v1 = self.pcurleig * (expL[:,:,:,2] * self.bplus + expL[:,:,:,3] * self.vplus)
            self.b1 += cp.conjg(self.pcurleig) * (expL[:,:,:,4] * self.bminus + expL[:,:,:,5] * self.vminus)
            self.v1 += cp.conjg(self.pcurleig) * (expL[:,:,:,6] * self.bminus + expL[:,:,:,7] * self.vminus)
            self.b1[:,0,0,0] = avgb
            self.v1[:,0,0,0] = avgv

            self.bplus = cp.sum(cp.conj(self.pcurleig)*self.brhs1,axis=0)
            self.vplus = cp.sum(cp.conj(self.pcurleig)*self.vrhs1,axis=0)
            self.bminus = cp.sum(self.pcurleig*self.brhs1,axis=0) 
            self.vminus = cp.sum(self.pcurleig*self.vrhs1,axis=0)

            self.b1 += self.pcurleig[:,:,:,:] * (coef1[:,:,:,0] * self.bplus + coef1[:,:,:,1] * self.vplus)
            self.v1 += self.pcurleig[:,:,:,:] * (coef1[:,:,:,2] * self.bplus + coef1[:,:,:,3] * self.vplus)
            self.b1 += cp.conjg(self.pcurleig) * (coef1[:,:,:,4] * self.bminus + coef1[:,:,:,5] * self.vminus)
            self.v1 += cp.conjg(self.pcurleig) * (coef1[:,:,:,6] * self.bminus + coef1[:,:,:,7] * self.vminus)

            # What do we do about zero? Eigenvalue zero so just integrate as normal

            self.b1[:,0,0,0] += self.dt/2*self.brhs1[:,0,0,0]
            self.v1[:,0,0,0] += self.dt/2*self.vrhs1[:,0,0,0]
            self.mhelcorr += self.dt/2 * mhc
       
            self.brhs1,self.vrhs1 = self.hallrhs(self.b1,self.v1)
            mhc = self.helicitycorrection()

            self.bplus -= cp.sum(cp.conj(self.pcurleig)*(self.brhs1,axis=0) 
            self.vplus -= cp.sum(cp.conj(self.pcurleig)*(self.vrhs1,axis=0)
            self.bminus -= cp.sum(self.pcurleig*self.brhs1,axis=0) 
            self.vminus -= cp.sum(self.pcurleig*self.vrhs1,axis=0)
                            
            # Now self.bplus etc equal to rhs1- rhs2 in basis, reverse sign of sum below

            self.b1 = b2 - self.pcurleig[:,:,:,:] * (coef2[:,:,:,0] * self.bplus + coef2[:,:,:,1] * self.vplus)
            self.v1 = v2 - self.pcurleig[:,:,:,:] * (coef2[:,:,:,2] * self.bplus + coef2[:,:,:,3] * self.vplus)
            self.b1 -= cp.conjg(self.pcurleig) * (coef2[:,:,:,4] * self.bminus + coef2[:,:,:,5] * self.vminus)
            self.v1 -= cp.conjg(self.pcurleig) * (coef2[:,:,:,6] * self.bminus + coef2[:,:,:,7] * self.vminus)
            self.b1[:,0,0,0] += self.dt/2 * (self.brhs1[:,0,0,0])
            self.v1[:,0,0,0] += self.dt/2 * (self.vrhs1[:,0,0,0])
            self.mhelcorr += self.dt/2 * mhc                            
                            
            self.itime += 1
            self.time += self.dt
            self.runtime = time.time()-self.startruntime

            if self.itime % max(self.iterations//50,10) == 0 and self.itime > 0:
                self.steadystate()

	self.energyfile()
        self.checkpointfile()
        self.fulloutputfile()

        print(self.runtime)

        return(None)
