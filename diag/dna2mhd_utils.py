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
    if (1==1): #lpath[-4:-1].isnumeric() and int(lpath[-4:-1]) > 280: 
        kzgrid=np.arange(par['nkz0'])
        kzgrid=kzgrid*par['kzmin']
        kygrid=np.empty(par['nky0'])
        kxgrid=np.empty(par['nkx0'])
        for i in range(par['nky0']//2):
            kygrid[par['nky0']-1-i]=-float(i+1)*par['kymin']
            kygrid[i]=float(i)*par['kymin']
        kygrid[par['nky0']//2]=par['nky0']/2*par['kymin']
        for i in range(par['nkx0']//2):
            kxgrid[par['nkx0']-1-i]=-float(i+1)*par['kxmin']
            kxgrid[i]=float(i)*par['kxmin']
        kxgrid[par['nkx0']//2]=par['nkx0']//2*par['kxmin']
    else:
        kxgrid=np.arange((par['nkx0']))
        kxgrid=kxgrid*par['kxmin']
        kygrid=np.empty(par['nky0'])
        kzgrid=np.empty(par['nkz0'])
        herm_grid=np.arange(2)
        herm_grid=1.0*herm_grid
        for i in range(par['nky0']//2):
            kygrid[par['nky0']-1-i]=-float(i+1)*par['kymin']
            kygrid[i]=float(i)*par['kymin']
        kygrid[par['nky0']//2]=par['nky0']/2*par['kymin']
        for i in range(par['nkz0']//2):
            kzgrid[par['nkz0']-1-i]=-float(i+1)*par['kzmin']
            kzgrid[i]=float(i)*par['kzmin']
        kzgrid[par['nkz0']//2]=par['nkz0']//2*par['kzmin']
    return kxgrid,kygrid,kzgrid

def read_time_step_b(which_itime,swap_endian=False):
   """Reads a time step from b_out.dat.  Time step determined by \'which_itime\'"""
   file_name = par['diagdir'][1:-1]+'/b_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*3#par['nv0']
   mem_tot=ntot*16
   gt0=np.empty((3,par['nkz0'],par['nky0'],par['nkx0']))
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='complex128',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   #print sum(gt0)
   f.close()
   return gt0

def read_time_step_v(which_itime,swap_endian=False):
   """Reads a time step from v_out.dat.  Time step determined by \'which_itime\'"""
   file_name = par['diagdir'][1:-1]+'/v_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*3#par['nv0']
   mem_tot=ntot*16
   gt0=np.empty((3,par['nkz0'],par['nky0'],par['nkx0']))
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='complex128',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   #print sum(gt0)
   f.close()
   return gt0

def read_time_step_opt(which_itime,opt,swap_endian=False):
   """Reads a time step from opt_out.dat.  Time step determined by \'which_itime\'"""
   file_name = par['diagdir'][1:-1]+'/'+opt+'_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*3#par['nv0']    
   mem_tot=ntot*16
   gt0=np.empty((3,par['nkz0'],par['nky0'],par['nkx0']))
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='complex128',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   f.close()
   return gt0

def read_time_step_energy(which_itime,swap_endian=False):
   """Reads a time step from opt_out.dat.  Time step determined by \'which_itime\'"""
   file_name = par['diagdir'][1:-1]+'/energy_out.dat'
   f = open(file_name,'rb')
   gt0=np.empty((1))
   ntot = 5
   mem_tot = (ntot)*8
   gt0 = np.empty(5)
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='float64',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   #print sum(gt0)                                                                                      \
   f.close()
   return gt0

def read_time_step_energyspec(which_itime,swap_endian=False):
   """Reads a time step from opt_out.dat.  Time step determined by \'which_itime\'"""
   file_name = par['diagdir'][1:-1]+'/energyspec_out.dat'
   f = open(file_name,'rb')
   gt0=np.empty((1))
   ntot = par['nkx0']*par['nky0']*par['nkz0']*2
   mem_tot = (ntot)*8
   gt0 = np.empty((par['nkx0'],par['nky0'],par['nkz0'],2))
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='float64',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   #print sum(gt0)                                                                                      \                            
   f.close()
   return gt0

def get_time_from_bout(swap_endian=False):
   """Returns time array taken from b_out.dat"""
   file_name = par['diagdir'][1:-1]+ '/b_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*3#par['nv0']
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
   f.close()
   return time

def get_time_from_vout(swap_endian=False):
   """Returns time array taken from v_out.dat"""
   file_name = par['diagdir'][1:-1]+ '/v_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*3#par['nv0']
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
   f.close()
   # print(time)
   # work = input('Proceed? Y/N ')
   # if work == 'N':
   #    quit('Wrong itimes')
   return time

def get_time_from_optout(opt,swap_endian=False):
   """Returns time array taken from v_out.dat"""
   file_name = par['diagdir'][1:-1]+ '/'+opt+'_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*3
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

     if inp==0 or inp:
         time = np.append(time,inp)
     else:
         continue_read=0
   f.close()
   print(time)
   work = input('Proceed? Y/N ')
   if work == 'N':
       quit('Wrong times')
   return time

def get_time_from_energyout(swap_endian=False):
   """Returns time array taken from v_out.dat"""
   file_name = par['diagdir'][1:-1]+ '/energy_out.dat'
   f = open(file_name,'rb')
   ntot=5
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
   print(time)
   # work = input('Proceed? Y/N ')
   # if work == 'N':
   #    quit('Wrong times')
   f.close()
   return time

def get_time_from_energyspecout(swap_endian=False):
   """Returns time array taken from v_out.dat"""
   file_name = par['diagdir'][1:-1]+ '/energyspec_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*2
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

     if inp==0 or inp:
         time = np.append(time,inp)
     else:
         continue_read=0
   f.close()
   return time

def getb(lpath):
    """Saves b_out.dat (located in the directory specified by lpath) into a python-readable format b_xyz.dat
    which will also be located in the lpath directory.
    """
    #lpath='/scratch/04943/akshukla/hammet_dna_output/full/omt%g_nu%1.2f'%(omt,nu)
    #if lpath==None:
    #    lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    read_parameters(lpath)
    time = get_time_from_bout()
    #time=time[:1000]
    kx,ky,kz=get_grids()
    i_n=[0,1,2]
    savepath = lpath+'/b_xyz.dat'
    #g=np.zeros((len(time)-1,len(kx),len(ky,),len(kz),len(i_n)), dtype='complex64')
    #print('allocating array')
    g=np.memmap(savepath,dtype='complex64',mode='w+', shape=(len(time),len(kx),len(ky,),len(kz),len(i_n)) )
    np.save(lpath+'/bshape.npy',g.shape)
    np.save(lpath+'/timeb.npy',time)
    #g=np.zeros((len(time),len(kx),len(ky,),len(kz),len(i_n)), dtype='complex64')
    #print('starting loop')
    print(par)
    print('time length = ', len(time))
    for t in range(len(time)):
        if(t%1000==0):
            print(str(t))
        gt = read_time_step_b(t)
        gt = np.reshape(gt,(par['nkx0'],par['nky0'],par['nkz0'],3),order='F')
        g[t] = gt
    #np.save(lpath+'/g_allk_g04',g)
    #print('finished loop')
    return time, g

def getv(lpath):
    """Saves v_out.dat (located in the directory specified by lpath) into a python-readable format v_xyz.dat
    which will also be located in the lpath directory.
    """
    #lpath='/scratch/04943/akshukla/hammet_dna_output/full/omt%g_nu%1.2f'%(omt,nu)
    #if lpath==None:
    #    lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    read_parameters(lpath)
    time = get_time_from_vout()
    #time=time[:1000]
    kx,ky,kz=get_grids()
    i_n=[0,1,2]
    savepath = lpath+'/v_xyz.dat'
    #g=np.zeros((len(time)-1,len(kx),len(ky,),len(kz),len(i_n)), dtype='complex64')
    #print('allocating array')
    g=np.memmap(savepath,dtype='complex64',mode='w+', shape=(len(time),len(kx),len(ky,),len(kz),len(i_n)) )
    np.save(lpath+'/vshape.npy',g.shape)
    np.save(lpath+'/timev.npy',time)
    #g=np.zeros((len(time),len(kx),len(ky,),len(kz),len(i_n)), dtype='complex64')
    #print('starting loop')
    print(par)
    print('time length = ', len(time))
    for t in range(len(time)):
        if(t%1000==0):
            print(str(t))
        gt = read_time_step_v(t)
        gt = np.reshape(gt,(par['nkx0'],par['nky0'],par['nkz0'],3),order='F')
        g[t] = gt
    #np.save(lpath+'/g_allk_g04',g)
    #print('finished loop')
    return time, g

def getopt(lpath,opt):
    """Saves opt_out.dat (located in the directory specified by lpath) into a python-readable format opt_xyz.dat                                                                                              
    which will also be located in the lpath directory.                                                                                                                                                    
    """
    read_parameters(lpath)
    time = get_time_from_optout(opt)

    kx,ky,kz=get_grids()
    i_n=[0,1,2]
    savepath = lpath+'/'+opt+'_xyz.dat'

    g=np.memmap(savepath,dtype='complex64',mode='w+', shape=(len(time),len(kx),len(ky,),len(kz),len(i_n)) )
    np.save(lpath+'/'+opt+'shape.npy',g.shape)
    np.save(lpath+'/time'+opt+'.npy',time)

    print(par)
    print('time length = ', len(time))
    for t in range(len(time)):
        if(t%10==0):
            print(str(t))
        gt = read_time_step_opt(t,opt)
        gt = np.reshape(gt,(par['nkx0'],par['nky0'],par['nkz0'],3),order='F')
        g[t] = gt

    f = open(lpath+'/dum'+opt+'.txt','w')
    f.write('Finished loop')
    f.close()
     
    return time, g

def getenergy(lpath):
    time = get_time_from_energyout()
    #time=time[:1000] 

    kx,ky,kz=get_grids()
    i_n=[0,1,2]
    savepath = lpath+'/energy_xyz.dat'
    #g=np.zeros((len(time)-1,len(kx),len(ky,),len(kz),len(i_n)), dtype='complex64') 
   #print('allocating array') 
    g=np.memmap(savepath,dtype='float64',mode='w+', shape=(len(time),5))
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
        gt = np.reshape(gt,5,order='F')
        g[t] = gt
    #np.save(lpath+'/g_allk_g04',g)               

    #print('finished loop')

    f = open(lpath+'/dumen.txt','w')
    f.write('Finished loop')
    f.close()

    return time,g

def getenergyspec(lpath):
    """Saves energyspec_out.dat (located in the directory specified by lpath) into a python-readable format v_xyz.dat                         
    which will also be located in the lpath directory.                                                                               
    """
    read_parameters(lpath)
    time = get_time_from_energyspecout()
    kx,ky,kz=get_grids()
    i_n=[0,1,2]
    savepath = lpath+'/energyspec_xyz.dat'
    g=np.memmap(savepath,dtype='float64',mode='w+', shape=(len(time),len(kx),len(ky,),len(kz),2) )
    np.save(lpath+'/energyspecshape.npy',g.shape)
    np.save(lpath+'/timeenergyspec.npy',time)
    print(par)
    print('time length = ', len(time))
    for t in range(len(time)):
        if(t%1000==0):
            print(str(t))
        gt = read_time_step_energyspec(t)
        gt = np.reshape(gt,(par['nkx0'],par['nky0'],par['nkz0'],2),order='F')
        g[t] = gt
    f = open(lpath+'/dumenspec.txt','w')
    f.write('Finished loop')
    f.close()
    evk = np.reshape(g[:,:,:,:,0],(len(time),len(kx),len(ky),len(kz)))
    ebk = np.reshape(g[:,:,:,:,1],(len(time),len(kx),len(ky),len(kz)))
    return time, evk,ebk

def load_b(lpath):
    """
    This method can only be run after getb has been called at least once to save the b_xyz.dat file_name
    This quickly loads the b array which will have indices [time,kx,ky,kz, x/y/z]
    """
    read_parameters(lpath)
    if lpath==None:
        lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    time = np.load(lpath+'/timeb.npy')
    bload=np.memmap(lpath+'/b_xyz.dat',dtype='complex64',mode='r',shape=tuple(np.load(lpath+'/bshape.npy')))
    return time, bload

def load_v(lpath):
    """
    This method can only be run after getv has been called at least once to save the v_xyz.dat file_name
    This quickly loads the v array which will have indices [time,kx,ky,kz, x/y/z]
    """
    read_parameters(lpath)
    if lpath==None:
        lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    time = np.load(lpath+'/timev.npy')
    vload=np.memmap(lpath+'/v_xyz.dat',dtype='complex64',mode='r',shape=tuple(np.load(lpath+'/vshape.npy')))
    return time, vload

def load_opt(lpath,opt):
    """                                                                                                                                                                                                   
    This method can only be run after getv has been called at least once to save the v_xyz.dat file_name                                                                                                  
    This quickly loads the v array which will have indices [time,kx,ky,kz, x/y/z]                                                                                                                         
    """
    read_parameters(lpath)
    if lpath==None:
        lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    time = np.load(lpath+'/time'+opt+'.npy')
    optload=np.memmap(lpath+'/'+opt+'_xyz.dat',dtype='complex64',mode='r',shape=tuple(np.load(lpath+'/'+opt+'shape.npy')))
    return time, optload

def load_energy(lpath):
    """  This method can only be run after getv has been called at least once to save the v_xyz.dat file_name    
    This quickly loads the v array which will have indices [time,kx,ky,kz, x/y/z]"""
    read_parameters(lpath)
    if lpath==None:
        lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    time = np.load(lpath+'/timeenergy.npy')
    enload=np.memmap(lpath+'/energy_xyz.dat',dtype='float64',mode='r',shape=tuple(np.load(lpath+'/energyshape.npy')))
    return time, enload

def load_energyspec(lpath):
    """ This method can only be run after getv has been called at least once to save the v_xyz.dat file_name 
    This quickly loads the v array which will have indices [time,kx,ky,kz, x/y/z]  """
    read_parameters(lpath)
    if lpath==None:
        lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    time = np.load(lpath+'/timeenergyspec.npy')
    g = np.memmap(lpath+'/energyspec_xyz.dat',dtype='float64',mode='r',shape=tuple(np.load(lpath+'/energyspecshape.npy')))
    evk = np.reshape(g[:,:,:,:,0],(len(time),par['nkx0'],par['nky0'],par['nkz0']))
    ebk = np.reshape(g[:,:,:,:,1],(len(time),par['nkx0'],par['nky0'],par['nkz0']))
    return time, evk,ebk

def plot_bv(lpath,ix,iy,iz,ind,show=True,ask=True):
    """
    This is an example method that plots the timetraces of b and v at the specified wavevector (kx[ix],ky[iy],kz[iz]).
    ind specifies whether you want the x(0),y(1), or z(2) component.
    """
    ind_strings= ['x','y','z']
    ind_string=ind_strings[ind]
    read_parameters(lpath)
    timeb,b=load_b(lpath)
    timev,v=load_v(lpath)
    fig,ax=plt.subplots(2)
    ax[0].plot(timeb,b[:,ix,iy,iz,ind].real,label='Re')
    ax[0].plot(timeb,b[:,ix,iy,iz,ind].imag,label='Im')
    ax[0].set_ylabel('b_%s'%ind_string)
    ax[1].plot(timev,v[:,ix,iy,iz,ind].real,label='Re')
    ax[1].plot(timev,v[:,ix,iy,iz,ind].imag,label='Im')
    ax[1].set_ylabel('v_%s'%ind_string)
    ax[0].set_ylim(-3*np.median(np.abs(b[:,ix,iy,iz,ind])),3*np.median(np.abs(3*b[:,ix,iy,iz,ind])))
    ax[1].set_ylim(-3*np.median(np.abs(v[:,ix,iy,iz,ind])),3*np.median(np.abs(3*v[:,ix,iy,iz,ind])))
    ax[0].legend()
    ax[1].legend()
    kx,ky,kz=get_grids()
    fig.suptitle('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
    fig.supxlabel('time (1/wc)')
    if lpath[-1] != '/':
        lpath = lpath + '/'
    if not os.path.exists(lpath + 'bvs/'):
        os.mkdir(lpath + 'bvs/')
    plt.savefig(lpath+'bvs/bv_%s_%d_%d_%d'%(ind_string,ix,iy,iz))
    if show == True:
        plt.show()
    if ask:
        dum = input('Log Zoom? Y/N ')
    else:
        dum = 'N'
    if dum == 'Y':
        xl = int(input('Xmin? Integer '))
        xh = int(input('Xmax? Integer '))
        fig,ax=plt.subplots(2)
        ax[0].plot(timeb,b[:,ix,iy,iz,ind].real,label='Re')
        ax[0].plot(timeb,b[:,ix,iy,iz,ind].imag,label='Im')
        ax[0].set_ylabel('b_%s'%ind_string)
        ax[0].set_xlabel('time')
        ax[1].plot(timev,v[:,ix,iy,iz,ind].real,label='Re')
        ax[1].plot(timev,v[:,ix,iy,iz,ind].imag,label='Im')
        ax[1].set_ylabel('v_%s'%ind_string)
        ax[0].set_xlabel('time')
        ax[0].legend()
        ax[1].legend()       
        fig.suptitle('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
        ax[0].set_xlim(xl,xh)
        ax[1].set_xlim(xl,xh)
        ax[0].set_ylim(0.01,10**8)
        ax[0].set_yscale('log')
        ax[1].set_ylim(0.01,10**8)
        ax[1].set_yscale('log')
        plt.yscale('log')
        plt.savefig(lpath+'bvs/zoomed_bv_%s_%d_%d_%d'%(ind_string,ix,iy,iz))
        plt.show()
    plt.close()
    return timeb,b,timev,v

def plot_vreal_spectrum(lpath,ix,iy,iz,ind):
    """
    ix,iy,iz specifies the wavevector
    ind specifies x/y/z (0/1/2) component
    This is an example method that performs the fft on the real part of v and plots the result.
    It will return an array of the frequencies found.
    *** Right now it seems like freqs need to multiplied by 2*pi to get the right dispersion relation.
        I think this makes sense because w = 2*pi*f
    """
    ind_strings= ['x','y','z']
    ind_string=ind_strings[ind]
    read_parameters(lpath)
    time,v=load_v(lpath)
    v_k = v[:,ix,iy,iz,ind]
    #plt.plot(time,v_k)
    #plt.show()
    sp=fftshift(fft(v_k-np.mean(v_k.real)))
    freq = fftshift(fftfreq(time.shape[-1],d=.01))
    peaks,_ = find_peaks(np.abs(sp),threshold=10)
    print(freq[peaks])
    print(freq[peaks]*2*np.pi)
    plt.plot(np.abs(sp))
    plt.plot(peaks, sp[peaks], "x")
    plt.show()
    return 2*np.pi*freq[peaks]

def plot_vspectrum(lpath,ix,iy,iz,ind,show=True):
    """
    ix,iy,iz specifies the wavevector
    ind specifies x/y/z (0/1/2) component
    This is an example method that performs the fft on the real part of v and plots the result.
    It will return an array of the frequencies found.
    *** Right now it seems like freqs need to multiplied by 2*pi to get the right dispersion relation.
        I think this makes sense because w = 2*pi*f
    """
    ind_strings= ['x','y','z']
    ind_string=ind_strings[ind]
    read_parameters(lpath)
    kx,ky,kz=get_grids()
    time,v=load_v(lpath)
    v_k = v[:,ix,iy,iz,ind]
    #plt.plot(time,v_k)
    #plt.show()
    sp=fftshift(fft(v_k-np.mean(v_k)))
    freq = fftshift(fftfreq(time.shape[-1],d=.1))
    omega = 2*np.pi*freq
    peaks,_ = find_peaks(np.abs(sp),threshold=10)
    print(freq[peaks])
    print(freq[peaks]*2*np.pi)
    omega_plot = omega[(omega>-2*np.pi)&(omega<2*np.pi)]
    sp_plot= sp[(omega>-2*np.pi)&(omega<2*np.pi)]
    plt.plot(omega_plot, np.abs(sp_plot))
    rews = analytical_omega(lpath,ix,iy,iz)
    plt.plot([rews[0],rews[0]],[-1,1])
    plt.plot([rews[1],rews[1]],[-1,1])
    plt.plot([-rews[0],-rews[0]],[-1,1])
    plt.plot([-rews[1],-rews[1]],[-1,1])
    plt.ylim(-1.2*np.max(np.abs(sp_plot)),1.2*np.max(np.abs(sp_plot)))
    plt.xlim(max(-max(rews)*3,-6.3),min(max(rews)*3,6.3))
    #plt.plot(peaks, sp[peaks], "x")
    plt.ylabel('|FFT(v_%s)|'%ind_string )
    plt.xlabel('frequency')
    plt.title('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
    if lpath[-1] != '/':
        lpath =lpath +'/'
    if not os.path.exists(lpath + 'vspectra/'):
        os.mkdir(lpath + 'vspectra/')
    plt.savefig(lpath+'vspectra/vspectrum_%s_%d_%d_%d'%(ind_string,ix,iy,iz))
    if show == True:
        plt.show()
    plt.close()
    return freq[peaks]*2*np.pi

def plot_nls(lpath,ix,iy,iz,ind,show=True,ask=True):
    """This an example method that plots the timetraces of b and v at the specified wavevector (kx[ix],ky[iy],kz[iz]).
ind specifies whether you want the x(0),y(1), or z(2) component."""
    ind_strings= ['x','y','z']
    ind_string=ind_strings[ind]
    read_parameters(lpath)
    opts = ['vdb','bdcb','cbdb','vdv','bdb','db2','bdv']
    fmts = {'bdv':'om','vdb':'om','bdcb':'<g','cbdb':'sg','vdv':'Hr',
        'bdb':'xb','db2':'8b'}
 
    fig,ax = plt.subplots(2)
    i = 0

    optlist = []
    for opt in opts:
        if os.path.isfile(lpath+'/dum'+opt+'.txt'):
            optlist.append(load_opt(lpath,opt))
        else:
            optlist.append(getopt(lpath,opt))
        t = optlist[i][0]
        opty = np.array(optlist[i][1][:,ix,iy,iz,ind])
        fig,ax = plt.subplots(2)
        ax[0].plot(t,opty.real,fmts[opt],markersize=1)
        ax[1].plot(t,opty.imag,fmts[opt],markersize=1)
        i = i + 1
        ax[0].set_ylabel('real '+opt+' %s'%ind_string)
        ax[0].set_xlabel('time')
        ax[1].set_ylabel('imag '+opt+' %s'%ind_string)
        ax[0].set_xlabel('time')
        ax[0].set_ylabel('real '+opt+' %s'%ind_string)
        ax[0].set_xlabel('time')
        ax[1].set_ylabel('imag '+opt+' %s'%ind_string)
        ax[1].set_xlabel('time')
        kx,ky,kz=get_grids()
        fig.suptitle('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
        if lpath[-1] != '/':
            lpath = lpath + '/'
        if not os.path.exists(lpath + 'nls/'):
            os.mkdir(lpath + 'nls/')
        plt.savefig(lpath+'nls/'+opt+'_%s_%d_%d_%d'%(ind_string,ix,iy,iz))
        if show == True:
            plt.show()
        plt.close()

    if ask:
        duml = input('Log Zoom All? Y/N ')
    else:
        duml = 'N'
    
    if duml == 'Y':
        xl = int(input('Xmin? Integer '))
        xh = int(input('Xmax? Integer '))
        fig,ax = plt.subplots(2)
        i = 0
        for opt in opts:
            t = optlist[i][0]
            opty = np.array(optlist[i][1][(t>xl)*(t<xh),ix,iy,iz,ind])
            t = np.array(t[(t>xl)*(t<xh)])
            ax[0].plot(t,opty.real,fmts[opt],markersize=1,label=opt)
            ax[1].plot(t,opty.imag,fmts[opt],markersize=1,label=opt)
            i = i + 1
        ax[0].set_ylabel('real nonlinearity %s'%ind_string)
        ax[0].set_xlabel('time')
        ax[1].set_ylabel('imag nonlinearity %s'%ind_string)
        ax[0].set_xlabel('time')
        ax[0].legend()
        ax[1].legend()
        ax[0].set_ylabel('real nonlinearity %s'%ind_string)
        ax[0].set_xlabel('time')
        ax[1].set_ylabel('imag nonlinearity %s'%ind_string)
        ax[1].set_xlabel('time')
        ax[0].legend()
        ax[1].legend()
        ax[0].set_ylim(.01,10**8)
        ax[1].set_ylim(.01,10**8)
        ax[0].set_yscale('log')
        ax[1].set_yscale('log')
        fig.suptitle('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
        plt.savefig(lpath+'nls/nllogs_%s_%d_%d_%d'%(ind_string,ix,iy,iz))
        if show == True:
            plt.show()
        plt.close()

    if ask:
        dum = input('Bar Zoom All? Y/N ')
    else:
        dum = 'N'
    
    if dum == 'Y':
        Nb = 4
        bottom = np.zeros(Nb)
        upbottom = 0.0001*np.ones(Nb)
        lowbottom = -1*0.0001*np.ones(Nb)
        cs = {'bdv':'m','vdb':'m','bdcb':'g','cbdb':'g','vdv':'r',
            'bdb':'b','db2':'b'}
        hatchs = {'bdv':'++','vdb':'00','bdcb':'\\\\','cbdb':'////','vdv':'**','bdb':'xx','db2':'..'}

        # xl = int(input('Xmin? Integer '))
        # xh = int(input('Xmax? Integer '))
        xl = 0
        xh = 4*par['max_itime']
        x =  np.linspace(xl+(xh-xl)/(2*Nb),xh-(xh-xl)/(2*Nb),num=Nb)
        fig,ax=plt.subplots(2)
        i = 0 

        for opt in opts:
            y = np.zeros(Nb)
            t = optlist[i][0]
            opty = np.array(optlist[i][1][:,ix,iy,iz,ind])
            t = np.array(t) 
            dt = par['dt_max']
            N = np.size(t)
            for j in range(np.size(t)):
                if (j%4 == 0 or j % 4 == 3):
                    opty[j] = opty[j] * 1/6
                else:
                    opty[j] = opty[j] * 1/3
            w = 2/3*(xh-xl)/(Nb)
            for j in range(Nb):
                y[j] = np.sum(opty[j*(N//Nb):(j+1)*(N//Nb)],dtype='complex64')
                if y[j] >= 0:
                    bottom[j] = upbottom[j]
                if y[j] < 0:
                    bottom[j] = lowbottom[j]
            y0 = y.real.astype('int')
            y1 = y.imag.astype('int')
            ax[0].bar(x,y0,width=w,color=cs[opt],hatch=hatchs[opt],label=opt,bottom=bottom,align='center')        
            ax[1].bar(x,y1,width=w,color=cs[opt],hatch=hatchs[opt],label=opt,bottom=bottom,align='center')
            
            print(y0,y1)
            for j in range(Nb):
                if (j%10000) == 0:
                    print(j)
                if y[j] >= 0:
                    upbottom[j] += y[j]
                if y[j] < 0:
                    lowbottom[j] += y[j]
            i = i + 1

        ax[0].set_ylabel('real + nonlinearity')
        ax[0].set_xlabel('itime')
        ax[1].set_ylabel('real - nonlinearity')
        ax[0].set_xlabel('itime')
        ax[0].legend(loc='upper left')
        ax[1].legend(loc='upper left')
        fig.suptitle('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
        ax[0].set_xlim(xl,xh)
        ax[1].set_xlim(xl,xh)
        ax[0].grid(visible=True,which='major',axis='both')
        ax[1].grid(visible=True,which='major',axis='both')
        #ax[0].set_ylim(0.01,10**8)
        #ax[0].set_yscale('log')
        #ax[1].set_ylim(-10**8,-0.01)
        #ax[1].set_yscale('log')
        plt.savefig(lpath+'nls/bar_nls_%s_%d_%d_%d'%(ind_string,ix,iy,iz))
        plt.show()
        plt.close()

    return 0

def plot_energy(lpath,ntp,show=True,log=False):
    """ Plots Scalars Written in energy_out.dat """

    read_parameters(lpath)
    if os.path.isfile(lpath+'/dumen.txt'):
        timeen,enval = load_energy(lpath)
    else:
        timeen,enval = getenergy(lpath)

    shapes = {1:(1,1),2:(2,1),3:(2,2),4:(2,2),5:(2,3),6:(2,3),7:(3,3),8:(3,3),9:(3,3)}
    s = shapes[ntp+1]
    labels = {0:'Energy',1:'Magnetic Helicity',2:'Cross Helicity',3:'Entropy',4:'Next Param',5:'Next Param',6:'Next Param',7:'Next Param',8:'Next Param'}
    fnames = {0:'energy',1:'maghcty',2:'crosshcty',3:'entropy',4:'par4',5:'par5',6:'par6',7:'par7',8:'par8'}
    
    if not os.path.exists(lpath + '/eplots/'):
        os.mkdir(lpath + '/eplots/')
       
    for i in range(ntp+1):
        if not log:
            fig,ax = plt.subplots(1)
            ax.plot(timeen,enval[:,i])
            ax.set_xlabel('time (1/wc)')
            ax.set_ylabel(labels[i])
            fig.suptitle(labels[i])
            if np.amax(enval) > 10 ** 5:
                ax.set_yscale('log')
            plt.savefig(lpath+'/eplots/'+fnames[i])
            if show == True:
                plt.show()
            plt.close()
        else:
            fig,ax = plt.subplots(1)
            ax.plot(timeen,np.abs(enval[:,i]))
            ax.set_xlabel('time (1/wc)')
            ax.set_ylabel(labels[i])
            ax.set_yscale('log')
            fig.suptitle(labels[i])
            plt.savefig(lpath+'/eplots/'+fnames[i])
            if show == True:
                plt.show()
            plt.close()
    
    return timeen,enval

def plot_enspec(lpath,npt=1,zz=-1,show=True,log=False,linplot=False,newload=False):
    read_parameters(lpath)
    kx,ky,kz = get_grids()

    if not newload:
        t,bkt = load_b(lpath)
        t,vkt = load_v(lpath)
    else:
        if os.path.isfile(lpath+'/dumenspec.txt'):
            t,ekvf,ekbf = load_energyspec(lpath)
            ekvm = ekvf[:,:,:,1:]
            ekbm = ekbf[:,:,:,1:]
        else:
            t,ekvf,ekbf = getenergyspec(lpath)
            ekvm = ekvf[:,:,:,1:]
            ekbm = ekbf[:,:,:,1:]

    fmts = ["ks","mo","b^","g*","r8"]
    fig,ax = plt.subplots(2)
    kmag = np.sqrt(np.log(np.tensordot(np.tensordot(np.exp(kx**2),np.exp(ky**2),axes=0),np.exp(kz[1:]**2),axes=0)))
    if log:
        kmag = np.log(kmag)
    j = 0
    for i in range(0,np.size(t),np.size(t)//npt):
        if not newload:
            ekb = 0.5*(np.abs(bkt[i,:,:,1:,0])**2 + np.abs(bkt[i,:,:,1:,1])**2)
            ekv = 0.5*(np.abs(vkt[i,:,:,1:,0])**2 + np.abs(vkt[i,:,:,1:,1])**2)
        elif log:
            ekb = np.log(ekbm[i,:,:,:])
            ekv = np.log(ekvm[i,:,:,:])
        else:
            ekb = ekbm[i,:,:,:]
            ekv = ekvm[i,:,:,:]
        if (not newload) and log:
            ekb[:,:,:] = np.log(ekb[:,:,:])
            ekv[:,:,:] = np.log(ekv[:,:,:])
        if (zz == -1):
            x = np.reshape(kmag[:,:,:],np.size(kmag[:,:,:]))
            a = np.argsort(x)
            yb = np.reshape(ekb[:,:,:],np.size(kmag[:,:,:]))
            yv = np.reshape(ekv[:,:,:],np.size(kmag[:,:,:]))
        else:
            x = np.reshape(kmag[:,:,0],np.size(kmag[:,:,0]))
            a = np.argsort(x)
            yb = np.reshape(ekb[:,:,zz],np.size(kmag[:,:,0]))
            yv = np.reshape(ekv[:,:,zz],np.size(kmag[:,:,0]))

        if linplot:
            ax[0].plot(1/(x[a]**2),yb[a],fmts[np.mod(i,5)],label=np.format_float_positional(t[i],2),markersize=1/(j+1))
            ax[1].plot(1/(x[a]**2),yv[a],fmts[np.mod(i,5)],label=np.format_float_positional(t[i],2),markersize=1/(j+1))
        else:
            ax[0].plot(x[a],yb[a],fmts[np.mod(j,5)],label=np.format_float_positional(t[i],2),markersize=1/(j+1))
            ax[1].plot(x[a],yv[a],fmts[np.mod(j,5)],label=np.format_float_positional(t[i],2),markersize=1/(j+1))
        j = j + 1

    ax[0].legend(loc=1)
    ax[1].legend(loc=1)
    if log:
        label = " Log"
    else:
        label = ""
    if (zz == -1):
        fig.suptitle("Energy Spectrum")
        ax[0].set_ylabel("Magnetic Energy Spectrum"+label)
        ax[1].set_ylabel("Kinetic Energy Spectrum"+label)
        if linplot:
            fig.supxlabel("1/k^2")
        else:
            fig.supxlabel("k")
    else:
        fig.suptitle("Perpendicular Energy Spectrum kz = "+np.format_float_positional(kz[zz],2))
        ax[0].set_ylabel("Magnetic Energy Spectrum"+label)
        ax[1].set_ylabel("Kinetic Energy Spectrum"+label)
        fig.supxlabel("kperp")
    if not os.path.exists(lpath + '/eplots/'):
        os.mkdir(lpath + '/eplots/')
    plt.savefig(lpath+'/eplots/enspec'+str(zz)+'.png')
    if show == True:
        plt.show()
    plt.close()
    return(0)

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
    wp = kz[iz] * np.sqrt((1+0.5*k**2) + np.sqrt((1+0.5*k**2)**2 - 1))
    wm = kz[iz] * np.sqrt((1+0.5*k**2) - np.sqrt((1+0.5*k**2)**2 - 1))
    #wp = kz[iz]*(-np.sqrt(kx[ix]**2+ky[iy]**2+kz[iz]**2)/2 + np.sqrt(1+ (kx[ix]**2+ky[iy]**2+kz[iz]**2)/4))
    #wm = kz[iz]*(-np.sqrt(kx[ix]**2+ky[iy]**2+kz[iz]**2)/2 - np.sqrt(1+ (kx[ix]**2+ky[iy]**2+kz[iz]**2)/4))
    return wp,wm

def wherenezero(arr):
    zs = []
    sz = np.shape(arr)
    for i in range(sz[1]):
        for j in range(sz[2]):
            for k in range(sz[3]):
                for l in range(sz[4]):
                    if np.amax(np.abs(arr)) < 10**(-20):
                        zs.append([i,j,k,l])
    return(zs)

def plot_bspectrum(lpath,ix,iy,iz,ind,show=True):
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
    read_parameters(lpath)
    kx,ky,kz=get_grids()
    time,b=load_b(lpath)
    b_k = b[:,ix,iy,iz,ind]
    sp=fftshift(fft(b_k-np.mean(b_k)))    
    freq = fftshift(fftfreq(time.shape[-1],d=.1))
    omega = 2*np.pi*freq
    peaks,_ = find_peaks(np.abs(sp),threshold=10)
    print(freq[peaks])
    print(freq[peaks]*2*np.pi)
    omega_plot = omega[(omega>-2*np.pi)&(omega<2*np.pi)]
    sp_plot= sp[(omega>-2*np.pi)&(omega<2*np.pi)]
    plt.plot(omega_plot, np.abs(sp_plot))
    rews = analytical_omega(lpath,ix,iy,iz)
    plt.plot([rews[0],rews[0]],[-1,1])
    plt.plot([rews[1],rews[1]],[-1,1])
    plt.plot([-rews[0],-rews[0]],[-1,1])
    plt.plot([-rews[1],-rews[1]],[-1,1])
    plt.ylim(-1.2*np.max(np.abs(sp_plot)),1.2*np.max(np.abs(sp_plot)))
    plt.xlim(max(-max(rews)*3,-6.3),min(max(rews)*3,6.3))
    plt.ylabel('|FFT(b_%s)|'%ind_string )
    plt.xlabel('frequency')
    plt.title('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
    if lpath[-1] != '/':
        lpath =lpath +'/'
    if not os.path.exists(lpath + 'bspectra/'):
        os.mkdir(lpath + 'bspectra/')
    plt.savefig(lpath+'bspectra/bspectrum_%s_%d_%d_%d'%(ind_string,ix,iy,iz))
    if show == True:
        plt.show()
    plt.close()
    return freq[peaks]*2*np.pi

def errpct_energy(lpath):
    read_parameters(lpath)
    Nx = par['nkx0']
    Ny = par['nky0']
    Nz = par['nkz0']
    p = par['init_kolm']
    kx = par['kxmin']
    ky = par['kymin']
    kz = par['kzmin']
    E = 0
    Et = 0
    for i in range(0,1024):
        for j in range(0,1024):
            for k in range(1,1024):
                E += 1/((kx*i)**2 + (ky*j)**2 + (kz*k)**2)
                if ((i < Nx//2) and (j < Ny//2)) and (k < Nz):
                    Et += 1/((kx*i)**2 + (ky*j)**2 + (kz*k)**2)
        print(i)
    return (E-Et)/E

