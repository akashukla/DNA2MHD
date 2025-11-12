import cupy as cp
from cupy.fft import rfftn,irfftn
import time

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

        print(ham)

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

    def hallmodeenergies(self):

        # Dot product
        lwk = cp.abs(cp.sum(cp.conj(self.pcurleig)*(self.alphaleftwhist[None,:,:,:]*self.b1 + self.v1),axis=0) \
                     /cp.sqrt(self.alphaleftwhist**2.0+1))**2.0
        lw = 2.0*cp.sum(lwk[:,:,1:]) + cp.sum(lwk[:,:,0])
        lw *= (4.0*cp.pi**3)
        
        lwk = cp.abs(cp.sum(self.pcurleig*(-self.alphaleftwhist[None,:,:,:]*self.b1 + self.v1),axis=0) \
                     /cp.sqrt(self.alphaleftwhist**2.0+1))**2.0
        rw = 2.0 * cp.sum(lwk[:,:,1:]) + cp.sum(lwk[:,:,0])
        rw *= 4.0 * cp.pi**3
        
        lwk = cp.abs(cp.sum(cp.conj(self.pcurleig)*(self.b1 - self.alphaleftwhist[None,:,:,:]*self.v1),axis=0) \
                     /cp.sqrt(self.alphaleftwhist**2.0+1))**2.0
        lc = 2.0 * cp.sum(lwk[:,:,1:]) + cp.sum(lwk[:,:,0])
        lc *= 4.0 * cp.pi**3

        lwk = cp.abs(cp.sum(self.pcurleig*(self.b1 + self.alphaleftwhist[None,:,:,:]*self.v1),axis=0) \
                     /cp.sqrt(self.alphaleftwhist**2.0+1))**2.0
        rc = 2.0 * cp.sum(lwk[:,:,1:]) + cp.sum(lwk[:,:,0])
        rc *= 4.0 * cp.pi**3

        return(lw,lc,rw,rc)

    def helicitycorrection(self):

        correction = cp.real(self.b2[0,:,:,:]*cp.conj(self.v2[1,:,:,:])-self.b2[1,:,:,:]*cp.conj(self.v2[0,:,:,:]))

        mhc = 2.0*cp.sum(correction[:,:,1:])+cp.sum(correction[:,:,0])
        mhc *= -16*cp.pi**3
        
        return(mhc)

    def checkpointfile(self):

        cp.savez(self.lpath+"/s_checkpoint",itime=self.itime,nkx0=self.nkx0,\
                 nky0=self.nky0,nkz0=self.nkz0,time=self.time,dt=self.dt,\
                 b1=self.b1,v1=self.v1,mhelcorr=self.mhelcorr)

        return(None)

    def energyfile(self):        

        
        if self.itime == 0:
            self.energyarray = cp.array([])
        elif self.itime == self.itime_start:
            self.energyarray = cp.load(self.lpath+"/energy.npy")

        mh,ch = self.helicity()
        lw,lc,rw,rc = self.hallmodeenergies()

        energyarraytime = cp.zeros(11)
        energyarraytime[0] = self.time
        energyarraytime[1] = self.hamiltonian(0)
        energyarraytime[2] = mh
        energyarraytime[3] = ch

        energyarraytime[4] = self.hamiltonian(1)
        energyarraytime[5] = self.hamiltonian(2)

        energyarraytime[6] = lw
        energyarraytime[7] = lc
        energyarraytime[8] = rw
        energyarraytime[9] = rc

        energyarraytime[10] = self.mhelcorr
        
        self.energyarray = cp.append(self.energyarray,energyarraytime)

        if self.itime == self.iterations:
            cp.save(self.lpath+"/energy.npy",self.energyarray)
        
        return(None)

    def fulloutputfile(self):

        cp.savez(self.lpath+"/checkpoint"+str(self.itime),time=self.time,b1=self.b1,v1=self.v1)

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
                 solveprec=16,maxwallclock=86200):

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

        """
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
        """

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

        # Put the guide field in b1 - this picks up linear terms with the convolution
        self.b1[2,0,0,0] = 1.0

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

    def rk4(self):

        
        self.brhs1,self.vrhs1 = self.hallrhs(self.b1,self.v1)
        self.b2 = self.b1 + self.dt * 1/6 * self.brhs1
        self.v2 = self.v1 + self.dt * 1/6 * self.vrhs1
        
        self.brhs2,self.vrhs2 = self.hallrhs(self.b1+1/2*self.dt*self.brhs1,self.v1+1/2*self.dt*self.vrhs1)
        self.b2 += self.dt * 1/3 * self.brhs2
        self.v2 += self.dt * 1/3 * self.vrhs2

        self.brhs1,self.vrhs1 = self.hallrhs(self.b1+1/2*self.dt*self.brhs2,self.v1+1/2*self.dt*self.vrhs2)
        self.b2 += self.dt * 1/3 * self.brhs1
        self.v2 += self.dt * 1/3 * self.vrhs1

        self.brhs2,self.vrhs2 = self.hallrhs(self.b1+self.dt*self.brhs1,self.v1+self.dt*self.vrhs1)
        self.b1 = self.b2 + self.dt * 1/6 * self.brhs2
        self.v1 = self.v2 + self.dt * 1/6 * self.vrhs2

        return(None)
        
    def gauss2(self):

        a11 = 1/2
        self.b2 = self.b1
        self.v2 = self.v1 
        self.brhs2,self.vrhs2 = self.hallrhs(self.b2,self.v2)
        
        self.brhs1 = self.brhs2
        self.vrhs1 = self.vrhs2 

        self.b2 = self.b1 + a11 * self.dt * self.brhs1
        self.v2 = self.v1 + a11 * self.dt * self.vrhs1

        self.brhs2,self.vrhs2 = self.hallrhs(self.b2,self.v2)

        linferror = max(cp.amax(cp.abs(self.brhs2-self.brhs1)),cp.amax(cp.abs(self.vrhs2-self.vrhs1)))

        solveiteration = 0
        while solveiteration < 40 and linferror > 10.0**(-16.0):

            self.brhs1 = cp.copy(self.brhs2)
            self.vrhs1 = cp.copy(self.vrhs2)

            self.b2 = self.b1 + a11 * self.dt * self.brhs1
            self.v2 = self.v1 + a11 * self.dt * self.vrhs1

            self.brhs2,self.vrhs2 = self.hallrhs(self.b2,self.v2)
            linferror = max(cp.amax(cp.abs(self.brhs2-self.brhs1)),cp.amax(cp.abs(self.vrhs2-self.vrhs1)))
            solveiteration += 1

        if self.itime % 10 == 0:
            print("Iterations itime ",self.itime," = ",solveiteration)

        self.b1 += self.dt * self.brhs2
        self.v1 += self.dt * self.vrhs2

        self.mhelcorr += self.helicitycorrection() * self.dt

        return(None)

    def splitsimulation(self):

        self.startruntime = time.time()
        self.runtime = 0

        while self.itime < self.iterations and self.runtime < self.maxwallclock:

            # print("Bz ",self.b1[2,0,0,0])
            
            if self.itime % max(self.iterations//500,1) == 0:
                self.energyfile()
            if self.itime % 200 == 0: #max(self.iterations//50,10) == 0:
                self.checkpointfile()
            if self.itime % 200 == 0: #max(self.iterations//5,100) == 0:
                self.fulloutputfile()

            #print("Through diagnostics ",self.itime)
            
            self.b1 *= cp.exp(-self.vnu * self.kmags**2.0 * self.dt/2)
            self.v1 *= cp.exp(-self.etab * self.kmags**2.0 * self.dt/2)

            #print("Through dissipation ",self.itime)

            self.force()

            #print("Through force ",self.itime)

            if self.explicitrk4:
                self.rk4()
            else:
                self.gauss2()

            #print("Through Gauss2 ",self.itime)

            self.force()

            #print("Through force ",self.itime)            
            
            self.b1 *= cp.exp(-self.vnu * self.kmags**2.0 * self.dt/2)
            self.v1 *= cp.exp(-self.etab * self.kmags**2.0 * self.dt/2)

            #print("Through dissipation ",self.itime)
            
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

    

        

        
