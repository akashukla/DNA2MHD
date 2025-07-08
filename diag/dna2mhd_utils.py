 #/usr/bin/env python
# File: steps_updated.py

#from start_nl import *
import numpy as np
import matplotlib.pyplot as plt
#from dna_diags import read_parameters, get_time_from_gout,read_time_step_g,get_grids
import os
import re
import multiprocessing as mp
import sys
import scipy.fft
import scipy.signal
from scipy.signal import find_peaks
from scipy.fft import fft,fftfreq,fftshift,irfftn
import scipy.optimize as spo
import matplotlib.animation as anim
from matplotlib.ticker import ScalarFormatter,MultipleLocator,MaxNLocator,StrMethodFormatter,LogLocator,LogFormatterExponent,LinearLocator
from numpy import format_float_positional as ff

par={}       #Global Variable to hold parameters once read_parameters is called
namelists={}

def read_parameters(lpath):
    """Reads parameters from parameters.dat \n
    The parameters are in a dictionary call par \n
    and can be accessed via par['parameter_name']"""
    if lpath==None:
        parfile=open('./parameters.dat','r')
    else:
        parfile=open(lpath+'/parameters.dat', 'r')
    parameters_in=parfile.read()
    lines=parameters_in.split('\n')
    #    parameters={}
    #note: par comes from config.py
    num_lines=len(lines)
    print( "Number of lines", num_lines)
    print(lines[0])
    for i in range(num_lines):
         temp=lines[i].split()
         if temp:
              str_check_namelist=re.match("&",temp[0])
         if str_check_namelist:
              current_namelist=temp[0]
              print(current_namelist)
              namelists[current_namelist]=" "
         if len(temp)>2:
              #if (re.match(\d):
              str_check_sn=re.match("\d*\.?\d*[eE]-?\+?\d*",temp[2])
              str_check_int=re.match("\d*",temp[2])
              str_check_float=re.match("\d*\.\d*",temp[2])
              if (str_check_sn and str_check_sn.end()==len(temp[2])):
                   par[temp[0]]=float(temp[2])
                   namelists[current_namelist]=namelists[current_namelist]+" "+temp[0]
              elif (str_check_float and str_check_float.end()==len(temp[2])):
                   par[temp[0]]=float(temp[2])
                   namelists[current_namelist]=namelists[current_namelist]+" "+temp[0]
              elif (str_check_int and str_check_int.end()==len(temp[2])):
                   float_temp=float(temp[2])
                   par[temp[0]]=int(float_temp)
                   namelists[current_namelist]=namelists[current_namelist]+" "+temp[0]
              else:
                   par[temp[0]]=temp[2]
                   namelists[current_namelist]=namelists[current_namelist]+" "+temp[0]

    #par['kxmax']=(par['nkx0']-1)*par['kxmin']
    #par['kymax']=(par['nky0']/2-1)*par['kymin']
    #par['kzmax']=(par['nkz0']/2-1)*par['kzmin']
    par['ky_nyq']=(par['nky0']//2)*par['kymin']
    par['kz_nyq']=(par['nkz0']//2)*par['kzmin']
    par['nx0_big'] = 1 + 3 * par['nkx0'] // 2
    par['ny0_big'] = 3 * par['nky0'] // 2
    par['nz0_big'] = 3 * par['nkz0'] // 2
    if par['etg_factor'] != 0.0:
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print( "Warning! field solver in dna diags not implement for ky=0 and etg_factor != 0.")
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print( "!!!!!!!!!!!!!!!!!!!!!!!!!!")


def get_grids():
    """Returns kx,ky,kz grids in the same form as used in the code \n
    kxgrid = 0, kxmin, . . . kxmax \n
    kygrid = 0, kymin, . . . kymax, kymax+kymin, -kymax, . . . -kymin """

    kxgrid=np.arange((par['nx0_big']))
    kxgrid=kxgrid*par['kxmin']
    kygrid=np.empty(par['ny0_big'])
    kzgrid=np.empty(par['nz0_big'])
    herm_grid=np.arange(2)
    herm_grid=1.0*herm_grid

    for i in range(par['ny0_big']//2):
        kygrid[par['ny0_big']-1-i]=-float(i+1)*par['kymin']
        kygrid[i]=float(i)*par['kymin']
    kygrid[par['ny0_big']//2]=par['ny0_big']/2*par['kymin']
    
    for i in range(par['nz0_big']//2):
        kzgrid[par['nz0_big']-1-i]=-float(i+1)*par['kzmin']
        kzgrid[i]=float(i)*par['kzmin']
    kzgrid[par['nz0_big']//2]=par['nz0_big']//2*par['kzmin']
    return kxgrid,kygrid,kzgrid

def read_time_step_bv(bv,which_itime,swap_endian=False,):
   """Reads a time step from bv_out.dat.  Time step determined by \'which_itime\'"""
   file_name = par['diagdir'][1:-1]+'/'+bv+'_out.dat'
   f = open(file_name,'rb')
   ntot=par['nx0_big']*par['ny0_big']*par['nz0_big']*3#par['nv0']
   mem_tot=ntot*16
   gt0=np.empty((3,par['nz0_big'],par['ny0_big'],par['nx0_big']))
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='complex128',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   #print sum(gt0)
   f.close()
   return gt0

def read_time_step_energy(which_itime,swap_endian=False):
   """Reads a time step from energy_out.dat.  Time step determined by \'which_itime\'"""
   file_name = par['diagdir'][1:-1]+'/energy_out.dat'
   f = open(file_name,'rb')
   gt0=np.empty((1))
   ntot = 12
   mem_tot = (ntot)*8
   gt0 = np.empty(ntot)
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='float64',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   #print sum(gt0)                                                                                      \
   f.close()
   return gt0

def get_time_from_bvout(bv,swap_endian=False,tmax=2000000):
   """Returns time array taken from b_out.dat"""
   file_name = par['diagdir'][1:-1]+ '/'+bv+'_out.dat'
   f = open(file_name,'rb')
   ntot=par['nx0_big']*par['ny0_big']*par['nz0_big']*3#par['nv0']
   mem_tot=ntot*16
   time=np.empty(0)
   continue_read=1
   i=0
   while (continue_read):
     f.seek(i*(mem_tot+8))
     i=i+1
     inp=np.fromfile(f,dtype='float64',count=1)
     if swap_endian:
         inp=inp.newbyteorder()
     #print inp
     
     if inp==0 or inp:
         time = np.append(time,inp)
     else:
         continue_read=0
     if inp >= tmax:
         continue_read=0
     
   f.close()
   return time

def get_time_from_energyout(swap_endian=False,tmax=2000000):
   """Returns time array taken from v_out.dat"""
   file_name = par['diagdir'][1:-1]+ '/energy_out.dat'
   f = open(file_name,'rb')
   ntot = 12
   mem_tot=ntot*8
   time=np.empty(0)
   continue_read=1
   i=0
   while (continue_read):
     f.seek(i*(mem_tot+8))
     i=i+1
     inp=np.fromfile(f,dtype='float64',count=1)
     if swap_endian:
         inp=inp.newbyteorder()
     #print inp                                                                                  
     if inp==0 or inp:
         time = np.append(time,inp)
     else:
         continue_read=0
     if inp >= tmax:
         continue_read=0

   print(time)
   # work = input('Proceed? Y/N ')
   # if work == 'N':
   #    quit('Wrong times')
   f.close()
   return time

def getbv(lpath,bv,tmax=2000000):
    """Saves b_out.dat (located in the directory specified by lpath) into a python-readable format b_xyz.dat
    which will also be located in the lpath directory.
    """
    read_parameters(lpath)
    time = get_time_from_bvout(bv,tmax=tmax)
    #time=time[:1000]
    kx,ky,kz=get_grids()
    i_n=[0,1,2]
    savepath = lpath + '/'+bv+'_xyz.dat'
    #g=np.zeros((len(time)-1,len(kx),len(ky,),len(kz),len(i_n)), dtype='complex64')
    #print('allocating array')
    g=np.memmap(savepath,dtype='complex64',mode='w+', shape=(len(time),len(kx),len(ky,),len(kz),len(i_n)) )
     #g=np.zeros((len(time),len(kx),len(ky,),len(kz),len(i_n)), dtype='complex64')
    #print('starting loop')
    print(par)
    print('time length = ', len(time))
    for t in range(len(time)):
        if(t%1000==0):
            print(str(t))
        gt = read_time_step_bv(bv,t)
        gt = np.reshape(gt,(par['nx0_big'],par['ny0_big'],par['nz0_big'],3),order='F')
        g[t] = gt
    #np.save(lpath+'/g_allk_g04',g)
    #print('finished loop')
    np.save(lpath+'/time'+bv+'.npy',time)
    np.save(lpath+'/'+bv+'shape.npy',g.shape)
    return time, g

def getenergy(lpath,tmax=2000000):
    time = get_time_from_energyout(tmax=tmax)
    #time=time[:1000] 

    kx,ky,kz=get_grids()
    i_n=[0,1,2]
    ntot = 12
    savepath = lpath+'/energy_xyz.dat'
    #g=np.zeros((len(time)-1,len(kx),len(ky,),len(kz),len(i_n)), dtype='complex64') 
   #print('allocating array') 
    g=np.memmap(savepath,dtype='float64',mode='w+', shape=(len(time),ntot))
    np.save(lpath+'/energyshape.npy',g.shape)
    np.save(lpath+'/timeenergy.npy',time)
    #g=np.zeros((len(time),len(kx),len(ky,),len(kz),len(i_n)), dtype='complex64')
    #print('starting loop') 

    print(par)
    print('time length = ', len(time))
    for t in range(len(time)):
        if(t%20==0):
            print(str(t))
        gt = read_time_step_energy(t)
        gt = np.reshape(gt,ntot,order='F')
        g[t] = gt
    #np.save(lpath+'/g_allk_g04',g)

    #print('finished loop')

    np.save(lpath+'/energyshape.npy',g.shape)
    np.save(lpath+'/timeenergy.npy',time)
    
    return time,g

def load_bv(lpath,bv):
    """
    This method can only be run after getb has been called at least once to save the b_xyz.dat file_name
    This quickly loads the b array which will have indices [time,kx,ky,kz, x/y/z]
    """
    read_parameters(lpath)
    time = np.load(lpath+'/time'+bv+'.npy')
    bload=np.memmap(lpath+'/'+bv+'_xyz.dat',dtype='complex64',mode='r',shape=tuple(np.load(lpath+'/bshape.npy')))
    return time, bload

def load_energy(lpath):
    """  This method can only be run after getv has been called at least once to save the v_xyz.dat file_name    
    This quickly loads the v array which will have indices [time,kx,ky,kz, x/y/z]"""
    read_parameters(lpath)
    time = np.load(lpath+'/timeenergy.npy')
    enload=np.memmap(lpath+'/energy_xyz.dat',dtype='float64',mode='r',shape=tuple(np.load(lpath+'/energyshape.npy')))
    return time, enload

def read_checkpoint(lpath):

   """Reads a time step from bv_out.dat.  Time step determined by \'which_itime\'"""
   file_name = par['diagdir'][1:-1]+'/s_checkpoint.dat'
   f = open(file_name,'rb')
   
   ntot=par['nx0_big']*par['ny0_big']*par['nz0_big']*3
   mem_tot=ntot*16

   itime = np.fromfile(f,dtype='int32',count=1)   
   dt = np.fromfile(f,dtype='float64',count=1)   
   nkx0 = np.fromfile(f,dtype='int32',count=1)
   nky0 = np.fromfile(f,dtype='int32',count=1)
   nkz0 = np.fromfile(f,dtype='int32',count=1)
   time = np.fromfile(f,dtype='float64',count=1)

   gt0=np.empty((3,par['nz0_big'],par['ny0_big'],par['nx0_big']))
   gt0=np.fromfile(f,dtype='complex128',count=ntot)
   b1= np.reshape(gt0,(par['nx0_big'],par['ny0_big'],par['nz0_big'],3),order='F')

   gt0=np.empty((3,par['nz0_big'],par['ny0_big'],par['nx0_big']))
   gt0=np.fromfile(f,dtype='complex128',count=ntot)
   v1= np.reshape(gt0,(par['nx0_big'],par['ny0_big'],par['nz0_big'],3),order='F')

   mhc = np.fromfile(f,dtype='float64',count=1)
      
   return itime,dt,nkx0,nky0,nkz0,time,b1,v1,mhc

def index_shift(ix,iy,iz):
    if iy > par['nky0']/2:
        iy = iy + par['nky0']/2
    if iz > par['nkz0']/2:
        iz = iz + par['nkz0']/2
    return (ix,iy,iz)

def plot_bv(lpath,ix,iy,iz,ind,show=True):
    """
    This is an example method that plots the timetraces of b and v at the specified wavevector (kx[ix],ky[iy],kz[iz]).
    ind specifies whether you want the x(0),y(1), or z(2) component.
    """
    ind_strings= ['x','y','z']
    ind_string=ind_strings[ind]
    ix,iy,iz = index_shift(ix,iy,iz)

    if lpath[-1] == '/':
        lpath = lpath[:-1]
    read_parameters(lpath)
    
    if os.path.isfile(lpath+"/timeb.npy") and os.path.isfile(lpath+"/timev.npy"):
        timeb,b=load_bv(lpath,"b")
        timev,v=load_bv(lpath,"v")
    else:
        timeb,b = getbv(lpath,"b")
        timev,v = getbv(lpath,"v")
    
    fig,ax=plt.subplots(2)
    ax[0].plot(timeb,b[:,ix,iy,iz,ind].real,label='Re')
    ax[0].plot(timeb,b[:,ix,iy,iz,ind].imag,label='Im')
    ax[0].set_ylabel('b_%s'%ind_string,size="large")
    ax[1].plot(timev,v[:,ix,iy,iz,ind].real,label='Re')
    ax[1].plot(timev,v[:,ix,iy,iz,ind].imag,label='Im')
    ax[1].set_ylabel('v_%s'%ind_string,size="large")
    ax[0].set_ylim(-3*np.median(np.abs(b[:,ix,iy,iz,ind])),3*np.median(np.abs(3*b[:,ix,iy,iz,ind])))
    ax[1].set_ylim(-3*np.median(np.abs(v[:,ix,iy,iz,ind])),3*np.median(np.abs(3*v[:,ix,iy,iz,ind])))
    ax[0].legend()
    ax[1].legend()
    kx,ky,kz=get_grids()
    fig.suptitle('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]),size="large")
    fig.supxlabel('Time ($\omega_c^{-1}$)',size="large")
    if lpath[-1] != '/':
        lpath = lpath + '/'
    if not os.path.exists(lpath + 'bvs/'):
        os.mkdir(lpath + 'bvs/')
    plt.savefig(lpath+'bvs/bv_%s_%d_%d_%d'%(ind_string,ix,iy,iz),bbox_inches='tight')
    if show == True:
        plt.show()

    return timeb,b,timev,v

def center_width(data_width,data_min):
    # Get preferred center and half-width for LinearLocator data

    plot_width = 10**(np.ceil(np.log10(data_width)))
    if plot_width > 2 * data_width:
        plot_width = plot_width/2
    if plot_width > 2 * data_width:
        plot_width = plot_width/2
    if plot_width > 2 * data_width:
        plot_width = 10**(np.ceil(np.log10(data_width))) * 0.15

    plot_center = (np.ceil(10*data_min/plot_width)+4)*plot_width/10

    print(plot_width,plot_center)
    
    return(plot_width,plot_center)

def plot_energy(lpath,xb=1,tmax=2000000):
    """ Plots Scalars Written in energy_out.dat 
    Written in Order : 
    Energy, Magnetic Helicity, Canonical Helicity,
    Kinetic Energy, Magnetic Energy, 
    LW Energy, LC Energy, RW Energy, RC Energy,
    Magnetic Helicity Bound, Canonical Helicity Bound, Magnetic Helicity Correction"""
    
    read_parameters(lpath)
    if lpath[-1] == "/":
        lpath = lpath[:-1]
        
    if False : # os.path.isfile(lpath+'/timeenergy.npy'):
        timeen,enval = load_energy(lpath)
    else:
        timeen,enval = getenergy(lpath,tmax=tmax)

    #shapes = {1:(1,1),2:(2,1),3:(2,2),4:(2,2),5:(2,3),6:(2,3),7:(3,3),8:(3,3),9:(3,3)}
    #s = shapes[ntp+1]
    
    if np.any(np.isnan(enval)):
        return(timeen,-99.0*np.ones(np.shape(enval)))
    if not os.path.exists(lpath + '/eplots/'):
        os.mkdir(lpath + '/eplots/')
    xx = len(timeen)
    timeen = timeen[range(0,xx,xb)]
    enval = enval[range(0,xx,xb),:]

    # Total Energy Plot
    fig,ax = plt.subplots(1)
    ev = enval[:,0]/(4* np.pi**3)
    ax.plot(timeen,ev,"ks",markersize=1)
    r = np.max(ev) - np.min(ev)
    ax.set_ylabel("Total Energy / Guide Field Energy",size="large")
    ax.set_xlabel("Time ($\omega_c^{-1}$)",size="large")
    
    ax.set_ylim(bottom=min(0.5*(ev[0]),0.8*(np.amin(ev))),top=max(1.2*np.amax(ev),1.5*ev[0]))
    
    ax.yaxis.set_major_locator(LinearLocator())
    form = ScalarFormatter()
    ax.yaxis.set_major_formatter(form)

    fig.suptitle("Total Energy Evolution")
    plt.savefig(lpath+"/eplots/energy",bbox_inches='tight')

    ax.set_yscale("log")
    plt.savefig(lpath+"/eplots/energylog",bbox_inches='tight')
    plt.close()

    # Energy Discrepancy Plot
    fig,ax = plt.subplots(1)
    ev = np.abs((enval[:,0]-enval[0,0])/(4*np.pi**3))
    ax.plot(timeen,ev,color="tomato",marker="s",markersize=1)
    ax.set_ylabel("Energy - Initial Energy / Guide Field Energy",size="large")
    ax.set_xlabel("Time ($\omega_c^{-1}$)",size="large")
    ax.set_ylim(10**(-20.0),10**0.0)
    plt.tight_layout()
    
    ax.set_yscale("log")
    fig.suptitle("Error in Energy Conservation")
    plt.savefig(lpath+"/eplots/energyerr",bbox_inches='tight')
    plt.close()    

    # Magnetic Helicity Plot
    if not par["init_null"]:
        ev = (enval[:,1]+enval[:,-1])/enval[0,-3]
        ev0 = ev-enval[:,-1]/enval[0,-3]
    else:
        ev = enval[:,1]+enval[:,-1]
        ev0 = enval[:,1]
    fig,ax = plt.subplots(1)
    ax.plot(timeen,ev,"bs",markersize=1,label="Transformed")
    ax.plot(timeen,ev0,"rs",markersize=1,label="Original")
    r = np.max(ev) - np.min(ev)
    ax.set_ylabel("Helicity / Initial Helicity Bound",size="large")
    ax.set_xlabel("Time ($\omega_c^{-1}$)",size="large")
    ax.legend(loc = "upper right")
    
    ax.set_ylim(bottom=min(0.5*(ev[0]),(np.amin(ev0))),top=max(np.amax(ev0),1.5*ev[0]))
    ax.yaxis.set_major_locator(LinearLocator())
    form = ScalarFormatter()
    ax.yaxis.set_major_formatter(form)

    fig.suptitle("Magnetic Helicity Evolution")
    plt.savefig(lpath+"/eplots/maghcty",bbox_inches='tight')
    plt.close()

    if not par["init_null"]:
        ev = np.abs((enval[:,1]+enval[:,-1])/enval[0,-3])
        ev0 = np.abs(ev-enval[:,-1]/enval[0,-3])
    else:
        ev = np.abs(enval[:,1]+enval[:,-1])
        ev0 = np.abs(enval[:,1])
    fig,ax = plt.subplots(1)
    ax.plot(timeen,ev,"bs",markersize=1,label="Transformed")
    ax.plot(timeen,ev0,"rs",markersize=1,label="Original")
    r = np.max(ev) - np.min(ev)
    ax.set_ylabel("Helicity / Initial Helicity Bound",size="large")
    ax.set_xlabel("Time ($\omega_c^{-1}$)",size="large")
    ax.legend(loc = "upper right")


    ax.set_yscale("log")
    plt.savefig(lpath+"/eplots/maghctylog",bbox_inches='tight')
    plt.close()

    # Canonical Helicity Plot
    fig,ax = plt.subplots(1)
    if not par["init_null"]:
        ev = (enval[:,2]+enval[:,-1])/enval[0,-2]
        ev0 = enval[:,2]/enval[0,-2]
    else:
        ev = enval[:,2]+enval[:,-1]
        ev0 = enval[:,2]
    ax.plot(timeen,ev,"bs",markersize=1,label="Transformed")
    ax.plot(timeen,ev0,"rs",markersize=1,label="Original")
    r = np.max(ev) - np.min(ev)
    ax.set_ylabel("Helicity / Initial Helicity Bound",size="large")
    ax.set_xlabel("Time ($\omega_c^{-1}$)",size="large")
    ax.legend(loc="upper right")

    ax.set_ylim(bottom=min(0.5*(ev[0]),(np.amin(ev0))),top=max(np.amax(ev0),1.5*ev[0]))
    ax.yaxis.set_major_locator(LinearLocator())
    form = ScalarFormatter()
    ax.yaxis.set_major_formatter(form)
    
    fig.suptitle("Canonical Helicity Evolution")
    plt.savefig(lpath+"/eplots/canhcty",bbox_inches='tight')
    plt.close()
    
    if not par["init_null"]:
        ev = np.abs((enval[:,2]+enval[:,-1])/enval[0,-2])
        ev0 = np.abs(enval[:,2]/enval[0,-2])
    else:
        ev = np.abs(enval[:,2]+enval[:,-1])
        ev0 = np.abs(enval[:,2])
    fig,ax = plt.subplots(1)
    ax.plot(timeen,ev,"bs",markersize=1,label="Transformed")
    ax.plot(timeen,ev0,"rs",markersize=1,label="Original")
    r = np.max(ev) - np.min(ev)
    ax.set_ylabel("Helicity / Initial Helicity Bound",size="large")
    ax.set_xlabel("Time ($\omega_c^{-1}$)",size="large")
    ax.legend(loc = "upper right")

    ax.set_yscale("log")
    plt.savefig(lpath+"/eplots/canhctylog",bbox_inches='tight')
    plt.close()

    # Helicity Discrepancies
    fig,ax = plt.subplots(1)
    evm = np.abs(enval[:,1]+enval[:,-1]-enval[0,1])
    evc = np.abs(enval[:,2]+enval[:,-1]-enval[0,2])
    ax.plot(timeen,evm,color="tomato",marker="s",markersize=1,label="Magnetic")
    ax.plot(timeen,evc,color="mediumturquoise",marker="s",markersize=1,label="Canonical")
    ax.set_ylabel("Helicity - Initial Helicity",size="large")
    ax.set_xlabel("Time ($\omega_c^{-1}$)",size="large")
    ax.set_ylim(10**(-20.0),10**0.0)
    ax.set_yscale("log")
    ax.legend(loc="upper right")
    fig.suptitle("Error in Helicity Conservation")
    plt.tight_layout()
    plt.savefig(lpath+"/eplots/helicityerr",bbox_inches='tight')
    plt.close()    
    
    # Magnetic vs Kinetic Energy Plot
    fig,ax = plt.subplots(1)
    ax.plot(timeen,enval[:,3]/enval[:,0],"rs",markersize=1,label="Kinetic Energy")
    ax.plot(timeen,enval[:,4]/enval[:,0],"bs",markersize=1,label="Magnetic Energy")
    ax.set_ylabel('Energy Component / Total Energy',size="large")
    ax.set_xlabel('Time ($\omega_c^{-1}$)',size="large")
    ax.set_ylim(bottom=0,top=1)
    ax.legend()
    fig.suptitle("Kinetic and Magnetic Energies")
    plt.savefig(lpath+'/eplots/spliten',bbox_inches='tight')
    plt.close()

    # Mode Energies Plot
    fig,ax = plt.subplots(1)
    fmts = ['b--','b:','r--','r:']
    labels = ['+ Helicity Whistler','+ Helicity Cyclotron','- Helicity Whistler','- Helicity Cyclotron']
    for i in range(4):
        ax.plot(timeen,enval[:,5+i]/(4*np.pi**3),fmts[i],label=labels[i])
    ax.set_ylabel("Mode Energy / Guide Field Energy",size="large")
    ax.set_xlabel("Time ($\omega_c^{-1}$)",size="large")
    ax.set_ylim(10**(-8),np.amax(np.abs(enval[:,0])))
    ax.set_yscale("log")
    ax.legend()
    fig.suptitle("Hall MHD Normal Mode Energy Distribution")
    plt.savefig(lpath+"/eplots/modeen",bbox_inches='tight')
    plt.close()

    return timeen,enval

def numpy_enspec(lpath):
    # Save magnetic and kinetic energy spectrum data for future access

    if lpath[-1] == "/":
        lpath = lpath[-1]
        
    if not os.path.exists(lpath + '/eplots/'):
        os.mkdir(lpath + '/eplots/')

    itime,dt,nkx0,nky0,nkz0,time,b1,v1,mhc = read_checkpoint(lpath)
    if not os.path.isfile(lpath+"/eplots/enspecs.npz"):
        
        ekb = np.sum(np.abs(b1)**2.0,3)
        ekv = np.sum(np.abs(v1)**2.0,3)

        np.savez(lpath+"/eplots/enspecs",time=time,itime=itime,ekb=ekb,ekv=ekv)
        
    else:
        data = np.load(lpath+"/eplots/enspecs.npz")
        time = data["time"]
        itime = data["itime"]
        ekb = data["ekb"]
        ekv = data["ekv"]
            
    return(time,itime,ekb,ekv)
        
def plot_enspec(lpath,zz=-1,version=3,show=False):
    # Plot energy spectrum at checkpoint time
    read_parameters(lpath)
    kx,ky,kz = get_grids()
    Ky,Kx,Kz = np.meshgrid(ky,kx,kz)

    print("Maximum k Values\n")
    print(np.amax(ky),np.amax(kx),np.amax(kz)) 
    print(np.amax(Ky),np.amax(Kx),np.amax(Kz))

    kmag = np.sqrt(np.abs(Kx)**2 + np.abs(Ky)**2 + np.abs(Kz)**2)
    
    if lpath[-1] == "/":
        lpath = lpath[:-1]
        
    if not os.path.exists(lpath + '/eplots/'):
        os.mkdir(lpath + '/eplots/')        

    time,itime,ekb,ekv = numpy_enspec(lpath)

    ekb = ekb[1:,1:,1:]
    ekv = ekv[1:,1:,1:]
    kmag = kmag[1:,1:,1:]
    xmax = 2 * np.amax(kmag)
        
    def enspec_format(fig,ax,prefix1,prefix2,xmin = 0.1,xmax = 10,ymin=10**(-6),ymax=10):
        ax.set_ylim(ymin/2,ymax*2)
        ax.set_xlim(xmin/2,xmax*2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        if zz == -1:
            fig.suptitle(prefix1+" Energy Spectrum "+ff(time,2)+"($\omega_c^{-1}$)")
            ax.set_ylabel(prefix2+" Energy Spectrum",size="large")
            fig.supxlabel("|k| ($d_i^{-1}$)",size="large")
        else:
            fig.suptitle(prefix1+" Perpendicular Energy Spectrum "+ff(time,2)+"($\omega_c^{-1}$)")
            ax.set_ylabel(prefix2+" Perpendicular Energy Spectrum",size="large")
            fig.supxlabel("|$k_\perp$| ($d_i^{-1}$)",size="large")
        plt.tight_layout()
        return(fig,ax)

    ek = ekb + ekv

    # Get x values
    if (zz == -1):
        x = np.reshape(kmag,np.size(kmag))
        a = np.argsort(x)
        yb = np.reshape(ekb,np.size(kmag))
        yv = np.reshape(ekv,np.size(kmag))
        yt = np.reshape(ek,np.size(kmag))
    else:
        x = np.reshape(kmag[:,:,zz],np.size(kmag[:,:,zz]))
        a = np.argsort(x)
        yb = np.reshape(ekb[:,:,zz],np.size(kmag[:,:,0]))
        yv = np.reshape(ekv[:,:,zz],np.size(kmag[:,:,0]))
        yt = np.reshape(ek[:,:,zz],np.size(kmag[:,:,0]))

    xmin = x[a[1]]/2

    ymax = 10**0.0
    ymin = 10**(-20.0)

    print("X Limits",xmin,xmax)
    
    # Obtain needed phase information to correct initial spectrum
    
    fig,ax = plt.subplots(1)
    ax.plot(x[a[::101]],yb[a[::101]],"ks",markersize=1)
    fig,ax = enspec_format(fig,ax,"Magnetic 3D","Magnetic",xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
    fig.savefig(lpath+'/eplots/'+str(itime[0])+'benspec'+str(zz+1)+'.png',bbox_inches='tight')
    plt.close()
    
    fig,ax = plt.subplots(1)
    ax.plot(x[a[::101]],yv[a[::101]],"ks",markersize=1)
    fig,ax = enspec_format(fig,ax,"Kinetic 3D","Kinetic",xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
    fig.savefig(lpath+'/eplots/'+str(itime[0])+'venspec'+str(zz+1)+'.png',bbox_inches='tight')
    plt.close()    
    
    fig,ax = plt.subplots(1)
    ax.plot(x[a[::101]],yt[a[::101]],"ks",markersize=1)
    fig,ax = enspec_format(fig,ax,"3D","",xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
    fig.savefig(lpath+'/eplots/'+str(itime[0])+'enspec'+str(zz+1)+'.png',bbox_inches='tight')
    plt.close()
    
    fig,ax = plt.subplots(1)
    kperps,ekti = integrated_spectrum_1d(ek,lpath,v=version)
    ax.plot(kperps,ekti,"ks",markersize=1,label="Total")
    #kperps,ekbi = integrated_spectrum_1d(ekb,lpath,v=version)
    #ax.plot(kperps,ekbi,"b",label="Magnetic")
    #kperps,ekvi = integrated_spectrum_1d(ekv,lpath,v=version)
    #ax.plot(kperps,ekbi,"v",label="Kinetic")
    fig.suptitle("Integrated Total Energy Spectra")
    plt.tight_layout()
    ax.set_ylabel("Total Energy Spectra",size="large")

    ymax = 10**1.0
    ymin = 10**(-16.0)
    
    ax.set_xlim(xmin/2,2*xmax)
    ax.set_ylim(ymin/2,2*ymax)
    ax.yaxis.set_major_locator(LogLocator())
    ax.yaxis.set_minor_locator(LogLocator(subs=[2,3,4,5,6,7,8,9]))
    form = LogFormatterExponent(minor_thresholds="all")
    form2 = LogFormatterExponent(minor_thresholds="auto")
    ax.yaxis.set_major_formatter(form)
    ax.yaxis.set_minor_formatter(form2)
            
    ax.set_xlabel("$k_{\perp} (d_i^{-1})$",size="large")
    ax.set_yscale("log")
    ax.set_xscale("log")
    fig.savefig(lpath+'/eplots/t1denspec'+str(version)+'.png',bbox_inches="tight")
    if show == True:
        plt.show()
    else:
        plt.close()

    return(time,ek)

#if __name__ == '__main__':
#    #count = mp.cpu_count()
#    #start = 1
#    #stop = start+count
#    params = [(12,0.05), (15,0.20), (15,0.50), (5,0.00), (6, 0.00), (6,0.01), (6,0.05), (6,0.10), (6,0.50), (7,0.00), (7,0.01), (7,0.05), (7,0.50), (8,0.00), (8,0.01), (8,0.10), (8,0.50), (9,0.00), (9,0.01), (9,0.20), (9,0.50)]
#    count = len(params)
#    print('params = ', params)
#    print('count = %d'%count)
#    p = mp.Pool(count)
#    p.starmap(saveg,params)
#    #scores = p.map(gbmerror, range(start, stop))
#    #scores = np.array(scores)
#    #np.save('scores', scores)
#    p.close()
#    p.join()
#    print('all done')


#iif __name__ == '__main__':
#    omt = int(sys.argv[1])
#    nu = float(sys.argv[2])
#    style = str(sys.argv[3])
#    print(omt,nu, style)
#    print(type(omt), type(nu))
#    saveg(omt,nu,style)
#

def analytical_omega(lpath,ix,iy,iz):
    read_parameters(lpath)
    kx,ky,kz = get_grids()
    k = np.sqrt(kx[ix]**2 + ky[iy]**2 + kz[iz]**2)
    wp = kz[iz] * np.sqrt((1+0.5*(par["hall"]*k)**2) + np.sqrt((1+0.5*(par["hall"]*k)**2)**2 - 1))
    wm = kz[iz] * np.sqrt((1+0.5*(par["hall"]*k)**2) - np.sqrt((1+0.5*(par["hall"]*k)**2)**2 - 1))
    #wp = kz[iz]*(-np.sqrt(kx[ix]**2+ky[iy]**2+kz[iz]**2)/2 + np.sqrt(1+ (kx[ix]**2+ky[iy]**2+kz[iz]**2)/4))
    #wm = kz[iz]*(-np.sqrt(kx[ix]**2+ky[iy]**2+kz[iz]**2)/2 - np.sqrt(1+ (kx[ix]**2+ky[iy]**2+kz[iz]**2)/4))
    return wp,wm

def enheldev(lpath,local=0):
    te,e = plot_energy(lpath)
    dt = te[-1]-te[0]
    de = e[:,0]-e[0,0]
    dmh = (e[:,1]-e[0,1]) + (e[:,-1] - e[0,-1])
    dch = (e[:,2]-e[0,2]) + (e[:,-1] - e[0,-1])

    print(dt)

    # print(de[0])                                                                                                                                                                                                    
    # print(de[1]+de[-3])                                                                                                                                                                                             
    # print(de[4]+de[-3])                                                                                                                                                                                             
    if local > 0:
        dei = np.max(np.abs(de[local:]-de[:-local]))
        dmhi = np.max(np.abs(dmh[local:]-dmh[:-local]))
        dchi = np.max(np.abs(dch[local:]-dch[:-local]))
        return(dei,dmhi,dchi)
    else:
        return np.max(np.abs(de)),np.max(np.abs(dmh)),np.max(np.abs(dch))

def convert_spec_to_real(lpath,spectra):
    """Performs IFFTs needed to get real space values on one time b or mode"""

    read_parameters(lpath)
    print(len(np.shape(spectra)))
    if len(np.shape(spectra)) == 3:
        if par['splitx']:
            ts = np.transpose(spectra,axes=(2,1,0))
            ftrs = irfftn(ts,axes=(0,1,2))
            res = np.transpose(ftrs,(2,1,0))
        else:
            res = irfftn(spectra,axes=(0,1,2))
    else:
        if par['splitx']:
            ts = np.transpose(spectra,axes=(2,1,0,3))
            ftrs = irfftn(ts,axes=(0,1,2))
            res = np.transpose(ftrs,(2,1,0,3))
        else:
            res = irfftn(spectra,axes=(0,1,2))
    return(res)

def plot_bvspectrum(lpath,bv,ix,iy,iz,ind,show=False):
    """                                                                                                                                                                  
    ix,iy,iz specifies the wavevector                                                                                                                                    
    ind specifies x/y/z (0/1/2) component                                                                                                                                
    This is an example method that performs the fft on the real part of b and plots the result.                                                                          
    It will return an array of the frequencies found.                                                                                                                    
    *** Right now it seems like freqs need to multiplied by 2*pi to get the right dispersion relation.                                                                   
        I think this makes sense because w = 2*pi*f                                                                                                                      
    """
    ind_strings= ['x','y','z']
    ind_string=ind_strings[ind]
    ix,iy,iz = index_shift(ix,iy,iz)
    
    read_parameters(lpath)
    kx,ky,kz=get_grids()
    print(par['dt_max'])
    if bv=='b':
        time,b=load_bv(lpath,"b")
        dt = time[1]-time[0]
        b_k = b[:,ix,iy,iz,ind]
        sp=fftshift(fft(b_k-np.mean(b_k)))    
        freq = fftshift(fftfreq(time.shape[-1],d=dt))
    elif bv=='v':
        time,b=load_bv(lpath,"v")
        dt = time[1]-time[0]
        b_k = b[:,ix,iy,iz,ind]
        sp=fftshift(fft(b_k-np.mean(b_k)))
        freq = fftshift(fftfreq(time.shape[-1],d=dt))
    omega = 2*np.pi*freq
    peaks,_ = find_peaks(np.abs(sp),threshold=1)
    print(freq[peaks])
    print(freq[peaks]*2*np.pi)
    omega_plot = omega[(omega>-2*np.pi)&(omega<2*np.pi)]
    sp_plot= sp[(omega>-2*np.pi)&(omega<2*np.pi)]
    fig,ax= plt.subplots(1)
    ax.plot(omega_plot, np.abs(sp_plot))
    w1,w2 = np.abs(analytical_omega(lpath,ix,iy,iz))
    print(w1,w2)
    MM = 1.2*np.max(np.abs(sp_plot))
    ax.plot([w1,w1],[-MM,MM])
    ax.plot([w2,w2],[-MM,MM])
    ax.plot([-w1,-w1],[-MM,MM])
    ax.plot([-w2,-w2],[-MM,MM])
    ax.set_ylim(0.0,MM)
    ax.set_xlim(max(-6.3,-10*w1),min(6.3,10*w1))
    
    ax.set_ylabel('|FFT('+bv+'_%s)|'%ind_string ,size="large")
    ax.set_xlabel('frequency',size="large")
    fig.suptitle('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
    if lpath[-1] != '/':
        lpath =lpath +'/'
    if not os.path.exists(lpath + bv+'spectra/'):
        os.mkdir(lpath + bv+'spectra/')
    plt.savefig(lpath+bv+'spectra/'+bv+'spectrum_%s_%d_%d_%d'%(ind_string,ix,iy,iz))
    if show == True:
        plt.show()
    plt.close()
    return freq[peaks]*2*np.pi

def integrated_spectrum_1d(spec,lpath,v=3):
    read_parameters(lpath)
    kx,ky,kz = get_grids()
    kygrid,kxgrid = np.meshgrid(ky[1:],kx[1:])
    kpgrid = np.sqrt(kxgrid**2 + kygrid**2)
    kps = kpgrid.flatten()
    kperps = np.unique(kps)

    if v == 3: 
        specperp = np.sum(spec,axis=2)
        a = np.argsort(kps)
        kps_sorted = kps[a]
        spiral = specperp.flatten()[a]

        averages = []
        N = int(np.amax(kperps)//(3*kx[1]))
        kperps = np.linspace(0,np.amax(kperps),num=N)

        for i in range(0,N-1):
            test = np.nonzero((kps_sorted >= kperps[i])*(kps_sorted < kperps[i+1]))
            averages.append(np.average(spiral[test]) * 2 * np.pi * kperps[i+1])
            print("Ring i # Elements: ",np.count_nonzero((kps_sorted >= kperps[i])*(kps_sorted < kperps[i+1])))

        spec1d = np.array(averages)

    else:
        spec1d = 2*np.pi*kpgrid[:-1,0]*np.sum(spec,axis=-1)[:-1,0].flatten()
        kperps = kpgrid[:,0]

    return(kperps[:-1],spec1d)

def nlparam(lpath):

    read_parameters(lpath)
    kx,ky,kz = get_grids()
    t,itime,ekbf,ekvf = numpy_enspec(lpath)

    ekbf = ekbf[1:,1:,1:]
    ekvf = ekvf[1:,1:,1:]
    ekm = ekbf+ekvf
    
    kyy,kxx,kzz = np.meshgrid(ky[1:],kx[1:],kz[1:])
    kmags = np.sqrt(kxx**2 +kyy**2 + kzz**2)
    kperps = np.sqrt(kxx**2 + kyy**2)
    
    sq = np.sqrt(1+kmags**2 / 4)
    whist = kzz * (kmags/2 + sq)
    cyclo = kzz * (sq - kmags/2)
    mhd = kzz

    rmse = np.sqrt(2*ekm)
    xi_mhd = rmse*kperps/mhd
    xi_whist = rmse*kperps/whist
    xi_cyclo = rmse*kperps/cyclo

    plt.plot(kmags.flatten(),xi_cyclo.flatten(),'bs',label="Cyclotron",markersize=1)
    plt.plot(kmags.flatten(),xi_mhd.flatten(),'ks',label="MHD",markersize=1)
    plt.plot(kmags.flatten(),xi_whist.flatten(),'rs',label="Whister",markersize=1)
    plt.ylabel("Nonlinearity Parameter",size="large")
    plt.xlabel("|k| ($d_i^{-1}$)",size="large")
    plt.title("3D Nonlinearity Parameter Spectrum t = "+np.format_float_positional(t,1)+" ($\omega_c^{-1}$)")
    plt.ylim(10**(-5),10**1)
    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    plt.show()

    return(kperps,xi_cyclo,xi_mhd,xi_whist)

def modes_from_check(lpath):
    """Post Process b and v into Normal Modes"""

    read_parameters(lpath)

    if lpath[-1] == "/":
        lpath = lpath[:-1]

    itime,dt,nkx0,nky0,nkz0,time,b1,v1,mhc = read_checkpoint(lpath)

    if os.path.isfile(lpath+"/modes"+str(itime[0])+".npz"):
        data = np.load(lpath+"/modes"+str(itime[0])+".npz")
        time = data["time"]
        itime = data["itime"]
        lwk = data["lwk"]
        lck = data["lck"]
        rwk = data["rwk"]
        rck = data["rck"]
    else:
        kx,ky,kz = get_grids()
        
        Ky,Kx,Kz = np.meshgrid(ky,kx,kz)
        kmags = np.sqrt(Kx**2 + Ky**2 + Kz**2)
        alpha_lw = -(par["hall"]*kmags)/2 - np.sqrt(1+(par["hall"]*kmags)**2 /4)
        alpha_lc = -(par["hall"]*kmags)/2 + np.sqrt(1+(par["hall"]*kmags)**2 /4)
    
        Kvec = np.zeros([par['nx0_big'],par['ny0_big'],par['nz0_big'],3],dtype='complex64')
        Zvec = np.zeros([par['nx0_big'],par['ny0_big'],par['nz0_big'],3],dtype='complex64')
        Kvec[:,:,:,0] = Kx
        Kvec[:,:,:,1] = Ky
        Kvec[:,:,:,2] = Kz
        Zvec[:,:,:,2] = 1
        
        pceig = np.zeros([par['nx0_big'],par['ny0_big'],par['nz0_big'],3],dtype='complex64')
        pceig = np.cross(Kvec,Zvec,axis=-1)
        pceig += 1/kmags[:,:,:,None] * 1.0j * np.cross(Kvec,pceig,axis=-1)
        pceig[0,0,0:par["nz0_big"]//2,:] = np.array([-1j/np.sqrt(2),1/np.sqrt(2),0.0],dtype='complex64')
        pceig[0,0,par["nz0_big"]//2:par["nz0_big"],:] = np.array([1j/np.sqrt(2),1/np.sqrt(2),0.0],dtype='complex64')
                
        #for i in range(par['nx0_big']):
        #    for j in range(par['ny0_big']):
        #        for k in range(par['nz0_big']):
        #            pceig[i,j,k,:] = np.cross(Kvec[i,j,k,:],Zvec[i,j,k,:])+1/kmags[i,j,k] * 1.0j * np.cross(Kvec[i,j,k,:],np.cross(Kvec[i,j,k,:],Zvec[i,j,k,:]))

        pceig = pceig / (np.sqrt(2) * np.sqrt(Kx[:,:,:,None]**2 + Ky[:,:,:,None]**2))
        pceig[0,0,:,:] = 0

        lwk = np.zeros([par['nx0_big'],par['ny0_big'],par['nz0_big']],dtype='complex64')
        lck = np.zeros([par['nx0_big'],par['ny0_big'],par['nz0_big']],dtype='complex64')
        rwk = np.zeros([par['nx0_big'],par['ny0_big'],par['nz0_big']],dtype='complex64')
        rck = np.zeros([par['nx0_big'],par['ny0_big'],par['nz0_big']],dtype='complex64')

        lwk = np.sum(np.conj(pceig[:,:,:,:])*(alpha_lw[:,:,:,None]*b1[:,:,:,:]+v1[:,:,:,:]),axis=-1)/np.sqrt(alpha_lw**2.0 + 1)
        lck = np.sum(np.conj(pceig[:,:,:,:])*(alpha_lc[:,:,:,None]*b1[:,:,:,:]+v1[:,:,:,:]),axis=-1)/np.sqrt(alpha_lc**2.0 + 1)
        rwk = np.sum(pceig[:,:,:,:]*(-alpha_lw[:,:,:,None]*b1[:,:,:,:]+v1[:,:,:,:]),axis=-1)/np.sqrt(alpha_lw**2.0 + 1)
        rck = np.sum(pceig[:,:,:,:]*(-alpha_lc[:,:,:,None]*b1[:,:,:,:]+v1[:,:,:,:]),axis=-1)/np.sqrt(alpha_lc**2.0 + 1)
        
        #for i in range(par['nx0_big']):
        #    for j in range(par['ny0_big']):
        #        for k in range(par['nz0_big']):
        #            lwk[i,j,k] = np.dot(np.conj(pceig[i,j,k,:]),alpha_lw[i,j,k]*b1[i,j,k,:]+v1[i,j,k,:])/np.sqrt(alpha_lw[i,j,k]**2+1)
        #            lck[i,j,k] = np.dot(np.conj(pceig[i,j,k,:]),alpha_lc[i,j,k]*b1[i,j,k,:]+v1[i,j,k,:])/np.sqrt(alpha_lc[i,j,k]**2+1)
        #            rwk[i,j,k] = np.dot((pceig[i,j,k,:]),-alpha_lw[i,j,k]*b1[i,j,k,:]+v1[i,j,k,:])/np.sqrt(alpha_lw[i,j,k]**2+1)
        #            rck[i,j,k] = np.dot((pceig[i,j,k,:]),-alpha_lc[i,j,k]*b1[i,j,k,:]+v1[i,j,k,:])/np.sqrt(alpha_lc[i,j,k]**2+1)
        
        np.savez(lpath+"/modes"+str(itime[0])+".npz",time=time,itime=itime,lwk=lwk,lck=lck,rwk=rwk,rck=rck)

    return(time,itime,lwk,lck,rwk,rck)

def mode_break(lpath,show=False,tmax=200000):

    read_parameters(lpath)
    kx,ky,kz = get_grids()

    time,itime,lwk,lck,rwk,rck = modes_from_check(lpath)

    lwk1 = lwk.flatten()
    lck1 = lck.flatten()
    rwk1 = rwk.flatten()
    rck1 = rck.flatten()

    lwk_inds = np.argsort(np.abs(lwk1)**2.0)
    lck_inds = np.argsort(np.abs(lck1)**2.0)
    rwk_inds = np.argsort(np.abs(rwk1)**2.0)
    rck_inds = np.argsort(np.abs(rck1)**2.0)

    Kx,Ky,Kz = np.meshgrid(kx,ky,kz,indexing="ij")
    
    for i in range(10):
        print("LWK ",i,np.abs(lwk1[lwk_inds[-i]])**2.0,Kx.flatten()[lwk_inds[-i]],Ky.flatten()[lwk_inds[-i]],Kz.flatten()[lwk_inds[-i]])
        print("LCK ",i,np.abs(lck1[lck_inds[-i]])**2.0,Kx.flatten()[lck_inds[-i]],Ky.flatten()[lck_inds[-i]],Kz.flatten()[lck_inds[-i]])
        print("RWK ",i,np.abs(rwk1[rwk_inds[-i]])**2.0,Kx.flatten()[rwk_inds[-i]],Ky.flatten()[rwk_inds[-i]],Kz.flatten()[rwk_inds[-i]])
        print("RCK ",i,np.abs(rck1[rck_inds[-i]])**2.0,Kx.flatten()[rck_inds[-i]],Ky.flatten()[rck_inds[-i]],Kz.flatten()[rck_inds[-i]])

    mode_ks = np.stack((lwk,lck,rwk,rck))
    fmts = ['b--','b:','r--','r:']
    labels = ['+ Helicity Whistler','+ Helicity Cyclotron','- Helicity Whistler','- Helicity Cyclotron']
    
    # Other plots: 1D spectra
    comment = """

    fig,ax = plt.subplots(2)

    for i in range(4):
        kperps,spec1di = integrated_spectrum_1d(mode_ks[i,0,:,:,:],lpath)
        kperps,spec1df = integrated_spectrum_1d(mode_ks[i,-1,:,:,:],lpath)
        ax[0].plot(kperps,spec1di,fmts[i],label=labels[i],markersize=1)
        ax[1].plot(kperps,spec1df,fmts[i],label=labels[i],markersize=1)
    ax[0].set_ylabel("Initial")
    ax[1].set_ylabel("Final")
    ax[1].set_xlabel("$k_\perp$ ($d_i^{-1}$)")
    ax[0].set_yscale("log")
    ax[1].set_yscale("log")
    ax[0].set_xscale("log")
    ax[1].set_xscale("log")
    ax[0].set_ylim(10**(-10),10**1)
    ax[1].set_ylim(10**(-10),10**1)
    ax[0].legend()
    ax[1].legend()
    fig.suptitle("Mode Energy Spectra at t = %.2f and t = %.2f " % (t[0],t[-1]))
    fig.savefig(lpath+"/eplots/modespec.png")
    if show == True:
        plt.show()
    plt.close()"""

    fig,ax = plt.subplots(1)

    cmin = 0
    cmax = 0

    kperps,spec1df1 = integrated_spectrum_1d(0.5 * np.abs(mode_ks[0,:,:,:])**2.0,lpath,v=3)
    kperps,spec1df2 = integrated_spectrum_1d(0.5 * np.abs(mode_ks[1,:,:,:])**2.0,lpath,v=3)
    kperps,spec1df3 = integrated_spectrum_1d(0.5 * np.abs(mode_ks[2,:,:,:])**2.0,lpath,v=3)
    kperps,spec1df4 = integrated_spectrum_1d(0.5 * np.abs(mode_ks[3,:,:,:])**2.0,lpath,v=3)

    ax.plot(kperps,spec1df1,fmts[0],label=labels[0])
    ax.plot(kperps,spec1df2,fmts[1],label=labels[1])
    ax.plot(kperps,spec1df3,fmts[2],label=labels[2])
    ax.plot(kperps,spec1df4,fmts[3],label=labels[3])

    if int(par["init_cond"]) >= 31:
        for i in range(3):
            kxi = par["wave"+str(i+1)+"x"]
            kyi = par["wave"+str(i+1)+"y"]
            kzi = par["wave"+str(i+1)+"z"]
            ki = np.sqrt(kx[kxi-1]**2 + ky[kyi-1]**2)
            kpind = np.argwhere(kperps < ki)[-1]
            if i == 1:
                ax.plot(kperps[kpind],spec1df1[kpind],marker="x",color=fmts[0][0],label="Excited Wave",markersize=10)
            else:
                ax.plot(kperps[kpind],spec1df1[kpind],marker="x",color=fmts[0][0],markersize=10)
            ax.plot(kperps[kpind],spec1df2[kpind],marker="x",color=fmts[1][0],markersize=10)
            ax.plot(kperps[kpind],spec1df3[kpind],marker="x",color=fmts[2][0],markersize=10)
            ax.plot(kperps[kpind],spec1df4[kpind],marker="x",color=fmts[3][0],markersize=10)
    
    m = min(np.amin(spec1df1[np.nonzero(spec1df1)]),np.amin(spec1df2[np.nonzero(spec1df2)]),
            np.amin(spec1df3[np.nonzero(spec1df3)]),np.amin(spec1df4[np.nonzero(spec1df4)]))
    M = max(np.amax(spec1df1[np.nonzero(spec1df1)]),np.amax(spec1df2[np.nonzero(spec1df2)]),
            np.amax(spec1df3[np.nonzero(spec1df3)]),np.amax(spec1df4[np.nonzero(spec1df4)]))

    ax.set_ylabel("Final Mode Energy Spectrum",size="large")
    ax.set_xlabel("$k_\perp$ ($d_i^{-1}$)",size="large")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylim(max(m/10,10**(-16)),M*10)

    if int(par["init_cond"]) >= 31:
        ax.legend(loc="lower right")
    else:
        ax.legend(loc="upper right")
    fig.suptitle("Mode Energy Spectra at t = %.2f $(\omega_c^{-1})$ " % (time))
    fig.savefig(lpath+"/eplots/modespec.png",bbox_inches="tight")
    if show == True:
        plt.show()
    plt.close()

    
    
    return (time,lwk,lck,rwk,rck)

def enheldev(lpath,local=0):
    te,e = plot_energy(lpath)
    dt = te[-1]-te[0]
    de = e[:,0]-e[0,0]
    dmh = (e[:,1]-e[0,1]) + (e[:,-1] - e[0,-1])
    dch = (e[:,2]-e[0,2]) + (e[:,-1] - e[0,-1])
    
    print(dt)
    
    # print(de[0])
    # print(de[1]+de[-3])
    # print(de[4]+de[-3])
    if local > 0:
        dei = np.max(np.abs(de[local:]-de[:-local]))
        dmhi = np.max(np.abs(dmh[local:]-dmh[:-local]))
        dchi = np.max(np.abs(dch[local:]-dch[:-local]))
        return(dei,dmhi,dchi)
    else:
        return np.max(np.abs(de)),np.max(np.abs(dmh)),np.max(np.abs(dch))

def structurefunction(lpath,tmax=2*10**10):

    read_parameters(lpath)
    kx,ky,kz = get_grids()

    x=""" if os.path.isfile(lpath+'/dumlasts.txt'):
        bk = np.load(lpath+'/b_fin.npy')
        vk = np.load(lpath+'v_fin.npy')
    else:
        t,bk,vk = lastbv(lpath)
        bk = np.load(lpath+'/b_fin.npy')
        vk = np.load(lpath+'v_fin.npy')"""

    time,itime,lwk,lck,rwk,rck = modes_from_check(lpath)
    
    labels = ['+ Helicity Whistler','+ Helicity Cyclotron','- Helicity Whistler','- Helicity Cyclotron']
    pmsymbols = ["++","+-","-+","--"]
    stname = ["phw","phc","nhw","nhc"]

    xs = np.pi / par['nkx0'] * np.arange(2*par['nkx0']) / par['kxmin']
    ys = 2*np.pi / par['nky0'] * np.arange(par['nky0']) / par['kymin']
    zs = 2*np.pi / par['nkz0'] * np.arange(par['nkz0']) / par['kzmin']
    
    modespecs = [lwk,lck,rwk,rck]

    fig,ax = plt.subplots(2,2)
    ax = ax.flatten()
    fig2,ax2 = plt.subplots(2,2)
    ax2 = ax2.flatten()

    fig3,ax3 = plt.subplots(2,2)
    ax3 = ax3.flatten()

    fig3,ax3 = plt.subplots(2,2)
    ax3 = ax3.flatten()

    fig4,ax4 = plt.subplots(2,2)
    ax4 = ax4.flatten()

    m = 10**4
    M = 0
    
    for I,ms in enumerate(modespecs):
        mm = convert_spec_to_real(lpath,ms/1j) # Divide by 1j because mode amplitudes are anti-Hermitian
            # assume axisymmetric for transverse struct fn - test this later
        # very small - lets rescale to check Parseval's theorem print(np.amax(mm**2)); sum(ms**2) = sum(mm**2)/N
        mm *= np.sqrt(np.sum(np.abs(ms[0,:,:])**2+2*np.abs(ms[1:,:,:])**2)*np.size(mm)/(np.sum(mm**2)))*np.sqrt(8*np.pi**3)/np.sqrt(par["kxmin"]*par["kymin"]*par["kzmin"])
        
        if par["splitx"]:
            nx = 2*par["nkx0"]
        else:
            nx = par["nkx0"]

        # perpendicular structure function
        str_perp = np.zeros(nx)
        for j in range(nx//2):
            if j != 0:
                str_perp[j] = np.average((np.roll(mm,j,0)-mm)**2)
            
        # parallel structure function
        str_par = np.zeros(par["nkz0"])
        for k in range(par["nkz0"]//2):
            if k != 0:
                str_par[k] = np.average((np.roll(mm,k,2)-mm)**2)

        pt = np.nonzero(str_perp)
        pz = np.nonzero(str_par)


        print("Maximum",np.amax(str_perp),np.amax(str_par))
        
        ax[I].plot(xs[pt],str_perp[pt],marker="s",markersize=1,color="tomato",label="Transverse",linestyle="")
        ax[I].plot(zs[pz],str_par[pz],marker="s",markersize=1,color="mediumturquoise",label="Parallel",linestyle="")
        if (I==2 or I == 3):
            ax[I].set_xlabel("r ($d_i$)")
        
        ax2[I].plot(str_perp[pt],marker="s",markersize=1,color="tomato",label="Transverse",linestyle="")
        ax2[I].plot(str_par[pz],marker="s",markersize=1,color="mediumturquoise",label="Parallel",linestyle="")

        ax3[I].plot(str_perp[pt],xs[pt],marker="s",markersize=1,color="tomato",label="Transverse",linestyle="")
        ax3[I].plot(str_par[pz],zs[pz],marker="s",markersize=1,color="mediumturquoise",label="Parallel",linestyle="")
        if (I==2 or I == 3):
            ax[I].set_xlabel("r ($d_i$)",size="large")
            ax2[I].set_xlabel("r (Grid Position)",size="large")
            ax3[I].set_xlabel("Structure Function",size="large")
            ax4[I].set_xlabel("x ($d_i$)",size="large")

        if (I == 0 or I == 2):
            ax[I].set_ylabel("$S_{"+pmsymbols[I]+"}^2(r)$")
            ax2[I].set_ylabel("$S_{"+pmsymbols[I]+"}^2(r)$")
            ax4[I].set_ylabel("$\chi=k_z/(k_\perp \sqrt{SF^2})$")
        
        strmin = max(np.amin(str_perp[pt]),np.amin(str_par[pz]))
        strmax = min(np.amax(str_perp[pt]),np.amax(str_par[pz]))

        # Set ax limits to be same
        mI = min(np.amin(str_perp[pt]),np.amin(str_par[pz]))
        MI = max(np.amax(str_perp[pt]),np.amax(str_par[pz]))

        if mI < m:
            m = mI
        if MI > M:
            M = MI

        # Find where GS parameter is well defined

        x1 = str_perp[pt]
        y1 = xs[pt]
        x2 = str_par[pz]
        y2 = zs[pz]

        a = np.argwhere((str_perp > strmin) * (str_par> strmin) * (str_perp<strmax) * (str_par<strmax))
        pgs = np.nonzero((str_perp > strmin) * (str_par> strmin) * (str_perp<strmax) * (str_par<strmax))
        if np.size(a) > 2:
            ax4[I].plot(xs[pgs],xs[pgs]/(zs[pgs]*np.sqrt(str_perp[pgs])),color="mediumturquoise",marker="s",markersize=1,linestyle="")
            
        ax[I].set_title(labels[I])
        ax[I].set_xscale("log")
        ax[I].set_yscale("log")
        ax2[I].set_title(labels[I])
        ax2[I].set_xscale("log")
        ax2[I].set_yscale("log")
        ax3[I].set_title(labels[I])
        ax3[I].set_xscale("log")
        ax3[I].set_yscale("log")
        ax4[I].set_title(labels[I])
        ax4[I].set_xscale("log")
        ax4[I].set_xlim(np.amin(xs),np.amax(xs))
        ax4[I].set_ylim(10**(-3),10)
        ax4[I].set_yscale("log")
        if (I < 2):
            ax[I].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False,labeltop=False)
            ax2[I].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False,labeltop=False)
            ax3[I].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False,labeltop=False)
            ax4[I].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False,labeltop=False)
        ax[I].legend(loc="lower right")
        ax2[I].legend(loc="lower right")

    for axi in ax:
        axi.set_ylim(m/2,M*2)
    for axi in ax2:
        axi.set_ylim(m/2,M*2)
        
    fig.suptitle("Structure Functions "+"t = %.2f $(\omega_c^{-1})$" % (time)) 
    fig2.suptitle("Structure Functions "+"t = %.2f $(\omega_c^{-1})$ "% (time))
    fig3.suptitle("Goldreich Sridhar Nonlinearity Parameter Calculation")
    fig4.suptitle("Goldreich Sridhar Nonlinearity Parameter")
    plt.tight_layout()
    plt.tight_layout()
    plt.tight_layout()
    plt.tight_layout()
    if lpath[-1] == '/':
        lpath = lpath[:-1]
    if not os.path.exists(lpath + '/eplots/'):
        os.mkdir(lpath + '/eplots/')    
    fig.savefig(lpath + "/eplots/stfns",bbox_inches="tight")
    fig2.savefig(lpath + "/eplots/stfns2",bbox_inches="tight")
    fig3.savefig(lpath + "/eplots/gs95pcalc",bbox_inches="tight")
    fig4.savefig(lpath+"/eplots/gs95p",bbox_inches="tight")
    plt.close()

    return(0)
        
def mode_nlparam(lpath,tt,dim,show=False,tmax=200000):

    read_parameters(lpath)
    kx,ky,kz = get_grids()

    t,lwk,lck,rwk,rck = modes_from_bv(lpath)

    kyy,kxx,kzz = np.meshgrid(ky[1:],kx[1:],kz[1:])
    kmags = np.sqrt(kxx**2 +kyy**2 + kzz**2)
    kperps = np.sqrt(kxx**2 + kyy**2)

    sq = np.sqrt(1+kmags**2 / 4)
    whist = kzz * (kmags/2 + sq)
    cyclo = kzz * (sq - kmags/2)

    mode = [lwk[tt,1:,1:,1:],lck[tt,1:,1:,1:],rwk[tt,1:,1:,1:],rck[tt,1:,1:,1:]]
    freq = [whist,cyclo,whist,cyclo]

    fmts = ['b--','b:','r--','r:']
    labels = ['+ Helicity Whistler','+ Helicity Cyclotron','- Helicity Whistler','- Helicity Cyclotron']
    plt.figure()
    for i in range(4):
        perp = kperps
        if dim == 1:
            perp,e = integrated_spectrum_1d(mode[i],lpath)
            whist = kzz[0,0,1] * (perp/2 + np.sqrt(1 + perp**2 / 4))
            cyclo = kzz[0,0,1] * (np.sqrt(1 + perp**2 / 4) -perp/2)
            freq = [whist,cyclo,whist,cyclo]
        else:
            perp = kperps[:,:,1]
            e = mode[i][:,:,1]
        rmse = np.sqrt(2*e)
        if dim == 3:
            xi = rmse*perp/freq[i][:,:,1]
        else:
            xi = rmse*perp/freq[i]

        if dim == 1:
            plt.plot(perp.flatten(),xi.flatten(),fmts[i],label=labels[i])
        else:
            plt.plot(kmags[:,:,1].flatten(),xi.flatten(),fmts[i],label=labels[i])
    
    plt.ylabel("Nonlinearity Parameter",size="large")
    plt.xlabel("|k| ($d_i^{-1}$)",size="large")
    plt.title(str(dim)+"D Nonlinearity Parameter Spectrum t = "+np.format_float_positional(t[tt],0)+" ($\omega_c^{-1}$)")
    plt.ylim(10**(-5),10**2)
    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    plt.savefig(lpath+'/eplots/nlpar'+str(dim)+'d'+str(int(t[tt])))
    
    return(0)

def patch_mhc(mhc):

    # Recovers magnetic helicity correction from a warm restart
    
    N = np.size(mhc)
    mhc1 = np.zeros(N)

    for i in range(N-1):
        # Only correct after a substantial drop
        if np.abs(mhc[i]) > np.abs(mhc[i+1]) and np.abs(mhc[i+1]) < 0.1 * np.abs(mhc[i]):
            mhc1[i+1:] += mhc[i]
        mhc1[i] += mhc[i]
    mhc1[-1] += mhc[-1]

    return(mhc1)

def threewaveenergy(lpath):

    read_parameters(lpath)

    timeen,enval = getenergy(lpath,tmax=20000000)
    
    if int(par["init_cond"]) >= 31:

        kxgrid,kygrid,kzgrid = get_grids()
        if lpath[-1] == "/":
            lpath = lpath[:-1]
            
        f = open(lpath+"/threewave_out.dat","rb")
        timeenergies = np.fromfile(f,dtype="float64",count=-1)

        kxs = []
        kys = []
        kzs = []
        ks = []

        for i in range(3):
            kx = kxgrid[par["wave"+str(i+1)+"x"]-1]
            kxs.append(kx)
            ky = kygrid[par["wave"+str(i+1)+"y"]-1]
            kys.append(ky)
            kz = kzgrid[par["wave"+str(i+1)+"z"]-1]
            kzs.append(kz)
            ks.append(np.sqrt(kx**2 + ky**2 + kz**2))

        wave1 = 1
        wave2 = 5
        wave3 = 9
        
        waves = [wave1,wave2,wave3]

        ii = np.argsort(ks)
        print(ii)

        waves_ii = [waves[ii[0]],waves[ii[1]],waves[ii[2]]]
        colors = ["mediumturquoise","tomato","olivedrab"]

        plt.figure()
        for i in range(3):
            plt.plot(timeenergies[::13],timeenergies[waves_ii[i]::13]/(4*np.pi**3.0),color=colors[i],label=(ff(kxs[ii[i]],3),ff(kys[ii[i]],3),ff(kzs[ii[i]],3)))
        plt.xlabel("Time ($\omega_c^{-1}$)")
        plt.ylabel("Wave Energy / Guide Field Energy")
        plt.ylim(10**(-7),10**1)
        if lpath == "/pscratch/sd/e/echansen/DNA2MHDruns/austinsherwoodperp2":
            plt.xlim(0,2000)
        plt.title("Resonance Condition Wave Interaction")
        plt.yscale("log")
        plt.legend(loc="upper right")
        plt.savefig(lpath+"/eplots/threewaves")

        for i in range(3):
            plt.figure()
            plt.plot(timeenergies[::13],timeenergies[4*i+1::13]/(4*np.pi**3.0),color="b",linestyle="--",label="Positive Whistler")
            plt.plot(timeenergies[::13],timeenergies[4*i+2::13]/(4*np.pi**3.0),color="b",linestyle=":",label="Positive Cyclotron")
            plt.plot(timeenergies[::13],timeenergies[4*i+3::13]/(4*np.pi**3.0),color="r",linestyle="--",label="Negative Whistler")
            plt.plot(timeenergies[::13],timeenergies[4*i+4::13]/(4*np.pi**3.0),color="r",linestyle=":",label="Negative Cyclotron")
            plt.legend()
            plt.xlabel("Time ($\omega_c^{-1}$)")
            plt.ylabel("Wave Energy / Guide Field Energy")
            plt.yscale("log")
            plt.ylim(10**(-7),10**1)
            plt.title("Wave Energies k "+ff(kxs[i],3)+" "+
                      ff(kys[i],3)+" "+ff(kzs[i],3))
            plt.legend(loc="lower right")
            plt.savefig(lpath+"/eplots/wavebreakdown"+str(i+1))
            plt.close()

        # Signal to Noise Ratio Plot

        plt.figure()
        energyratio = (timeenergies[wave1::13]+timeenergies[wave2::13]+timeenergies[wave3::13])/enval[:,1]
        plt.plot(timeenergies[::13],energyratio,"b8")
        plt.xlabel("Time ($\omega_c^{-1}$)")
        plt.ylabel("Triplet Energy / Wave Ensemble Energy")
        plt.yscale("log")
        plt.title("Triplet Signal to Noise Ratio")
        plt.savefig(lpath+"/eplots/tripletsignoise")
        plt.close()
            
        return(timeenergies)
    else:
        print("Not a three wave simulation")
        return(None)
