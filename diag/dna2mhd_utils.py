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
from matplotlib.ticker import ScalarFormatter,MultipleLocator

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
    if (0==1): #lpath[-4:-1].isnumeric() and int(lpath[-4:-1]) > 280: 
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
   if opt == 'xi':
      mem_tot = 8 * (par['nkx0']-1)*(par['nky0']-1)*(par['nkz0']-1)*3
   gt0=np.empty((3,par['nkz0'],par['nky0'],par['nkx0']))
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='complex128',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
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

def read_time_step_energyspec(which_itime,swap_endian=False,mode=False):
   """Reads a time step from energyspec_out.dat.  Time step determined by \'which_itime\'"""
   file_name = par['diagdir'][1:-1]+'/energyspec_out.dat'
   if mode:
       file_name = par['diagdir'][1:-1]+'/mode_out.dat'
   f = open(file_name,'rb')
   gt0=np.empty((1))
   ntot = par['nkx0']*par['nky0']*par['nkz0']*2
   if mode:
       ntot *= 2
   mem_tot = (ntot)*8
   if not mode:
       gt0 = np.empty((par['nkx0'],par['nky0'],par['nkz0'],2))
   else:
       gt0 = np.empty((par['nkx0'],par['nky0'],par['nkz0'],4))
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='float64',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   #print sum(gt0)                                                                                      \                            
   f.close()
   return gt0

def read_time_step_xi(which_itime,swap_endian=False):
   """Reads a time step from xi_out.dat.  Time step determined by \'which_itime\'"""
   file_name = par['diagdir'][1:-1]+'/xi_out.dat'
   f = open(file_name,'rb')
   gt0=np.empty((1))
   ntot = (par['nkx0']-1)*(par['nky0']-1)*(par['nkz0']-1)
   mem_tot = (ntot)*8
   gt0 = np.empty((par['nkx0']-1,par['nky0']-1,par['nkz0']-1))
   f.seek(8+which_itime*(8+mem_tot))
   gt0=np.fromfile(f,dtype='float64',count=ntot)
   if swap_endian:
       gt0=gt0.newbyteorder()
   #print sum(gt0)                                                                                      \                                                                                                                                                                                                                      
   f.close()
   return gt0

def get_time_from_bout(swap_endian=False,tmax=2000000):
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
     if inp >= tmax:
         continue_read=0
     
   f.close()
   return time

def get_time_from_vout(swap_endian=False,tmax=2000000):
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
     if inp >= tmax:
         continue_read=0

   f.close()
   # print(time)
   # work = input('Proceed? Y/N ')
   # if work == 'N':
   #    quit('Wrong itimes')
   return time

def get_time_from_optout(opt,swap_endian=False,tmax=2000000):
   """Returns time array taken from v_out.dat"""
   file_name = par['diagdir'][1:-1]+ '/'+opt+'_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*3
   mem_tot=ntot*16
   if opt == 'xi':
       mem_tot =8 * (par['nkx0']-1)*(par['nky0']-1)*(par['nkz0']-1)*3
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
     if inp >= tmax:
         continue_read=0
   f.close()
   print(time)

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

def get_time_from_energyspecout(swap_endian=False,tmax=2000000,mode=False):
   """Returns time array taken from v_out.dat"""
   file_name = par['diagdir'][1:-1]+ '/energyspec_out.dat'
   if mode:
       file_name = par['diagdir'][1:-1]+'/mode_out.dat'
   f = open(file_name,'rb')
   ntot=par['nkx0']*par['nky0']*par['nkz0']*2
   if mode:
       ntot *= 2
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
     if inp >= tmax:
         continue_read=0

   f.close()
   return time

def get_time_from_xiout(swap_endian=False,tmax=2000000):
   """Returns time array taken from v_out.dat"""
   file_name = par['diagdir'][1:-1]+ '/xi_out.dat'
   f = open(file_name,'rb')
   ntot=(par['nkx0']-1)*(par['nky0']-1)*(par['nkz0']-1)
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
     if inp >= tmax:
         continue_read=0

   f.close()
   return time

def getb(lpath,tmax=2000000):
    """Saves b_out.dat (located in the directory specified by lpath) into a python-readable format b_xyz.dat
    which will also be located in the lpath directory.
    """
    #lpath='/scratch/04943/akshukla/hammet_dna_output/full/omt%g_nu%1.2f'%(omt,nu)
    #if lpath==None:
    #    lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    read_parameters(lpath)
    time = get_time_from_bout(tmax=tmax)
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

def getv(lpath,tmax=2000000):
    """Saves v_out.dat (located in the directory specified by lpath) into a python-readable format v_xyz.dat
    which will also be located in the lpath directory.
    """
    #lpath='/scratch/04943/akshukla/hammet_dna_output/full/omt%g_nu%1.2f'%(omt,nu)
    #if lpath==None:
    #    lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    read_parameters(lpath)
    time = get_time_from_vout(tmax=tmax)
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

def getopt(lpath,opt,tmax=2000000):
    """Saves opt_out.dat (located in the directory specified by lpath) into a python-readable format opt_xyz.dat                                                                                              
    which will also be located in the lpath directory.                                                                                                                                                    
    """
    read_parameters(lpath)
    time = get_time_from_optout(opt,tmax=tmax)

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

    f = open(lpath+'/dumen.txt','w')
    f.write('Finished loop')
    f.close()

    return time,g

def getenergyspec(lpath,tmax=2000000,mode=False):
    """Saves energyspec_out.dat (located in the directory specified by lpath) into a python-readable format v_xyz.dat                         
    which will also be located in the lpath directory.                                                                               
    """
    read_parameters(lpath)
    time = get_time_from_energyspecout(tmax=tmax,mode=mode)
    kx,ky,kz=get_grids()
    i_n=[0,1,2]

    if mode:
        savepath = lpath+'/modespec_xyz.dat'
        g=np.memmap(savepath,dtype='float64',mode='w+', shape=(len(time),len(kx),len(ky,),len(kz),4) )
        np.save(lpath+'/modespecshape.npy',g.shape)
        np.save(lpath+'/timemodespec.npy',time)
    else:
        savepath = lpath+'/energyspec_xyz.dat'
        g=np.memmap(savepath,dtype='float64',mode='w+', shape=(len(time),len(kx),len(ky,),len(kz),2) )
        np.save(lpath+'/energyspecshape.npy',g.shape)
        np.save(lpath+'/timeenergyspec.npy',time)        
    print(par)
    print('time length = ', len(time))
    for t in range(len(time)):
        if(t%100==0):
            print(str(t))
        gt = read_time_step_energyspec(t,mode=mode)
        if not mode:
            gt = np.reshape(gt,(par['nkx0'],par['nky0'],par['nkz0'],2),order='F')
            g[t] = gt
        else:
            gt = np.reshape(gt,(par['nkx0'],par['nky0'],par['nkz0'],4),order='F')
            g[t] = gt
    if not mode:
        f = open(lpath+'/dumenspec.txt','w')
        f.write('Finished loop')
        f.close()
        evk = np.reshape(g[:,:,:,:,0],(len(time),len(kx),len(ky),len(kz)))
        ebk = np.reshape(g[:,:,:,:,1],(len(time),len(kx),len(ky),len(kz)))
        return time, evk,ebk
    else:
        f = open(lpath+'/dummodespec.txt','w')
        f.write('Finished loop')
        f.close()
        lwk = np.reshape(g[:,:,:,:,0],(len(time),len(kx),len(ky),len(kz)))
        lck = np.reshape(g[:,:,:,:,1],(len(time),len(kx),len(ky),len(kz)))
        rwk = np.reshape(g[:,:,:,:,2],(len(time),len(kx),len(ky),len(kz)))
        rck = np.reshape(g[:,:,:,:,3],(len(time),len(kx),len(ky),len(kz)))
        return time, lwk,lck,rwk,rck

def getxi(lpath,tmax=2000000):
    """Saves xi_out.dat (located in the directory specified by lpath) into a python-readable format v_xyz.dat                                                                                                                                                                                                          
    which will also be located in the lpath directory.                                                                                                                                                                                                                                                                         
    """
    read_parameters(lpath)
    time = get_time_from_xiout(tmax=tmax)
    kx,ky,kz=get_grids()
    i_n=[0,1,2]
    savepath = lpath+'/xi_xyz.dat'
    g=np.memmap(savepath,dtype='float64',mode='w+', shape=(len(time),len(kx)-1,len(ky)-1,len(kz)-1) )
    np.save(lpath+'/xishape.npy',g.shape)
    np.save(lpath+'/timexi.npy',time)
    print(par)
    print('time length = ', len(time))
    for t in range(len(time)):
        if(t%1000==0):
            print(str(t))
        gt = read_time_step_xi(t)
        gt = np.reshape(gt,(par['nkx0']-1,par['nky0']-1,par['nkz0']-1),order='F')
        g[t] = gt
    f = open(lpath+'/dumxi.txt','w')
    f.write('Finished loop')
    f.close()
    return time, g

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

def load_energyspec(lpath,mode=False):
    """ This method can only be run after getv has been called at least once to save the v_xyz.dat file_name 
    This quickly loads the v array which will have indices [time,kx,ky,kz, x/y/z]  """
    read_parameters(lpath)
    if lpath==None:
        lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    if not mode:
        time = np.load(lpath+'/timeenergyspec.npy')
        g = np.memmap(lpath+'/energyspec_xyz.dat',dtype='float64',mode='r',shape=tuple(np.load(lpath+'/energyspecshape.npy')))
        evk = np.reshape(g[:,:,:,:,0],(len(time),par['nkx0'],par['nky0'],par['nkz0']))
        ebk = np.reshape(g[:,:,:,:,1],(len(time),par['nkx0'],par['nky0'],par['nkz0']))
        return time,evk,ebk
    else:
        time = np.load(lpath+'/timemodespec.npy')
        g = np.memmap(lpath+'/modespec_xyz.dat',dtype='float64',mode='r',shape=tuple(np.load(lpath+'/modespecshape.npy')))
        lwk = np.reshape(g[:,:,:,:,0],(len(time),par['nkx0'],par['nky0'],par['nkz0']))
        lck = np.reshape(g[:,:,:,:,1],(len(time),par['nkx0'],par['nky0'],par['nkz0']))
        rwk = np.reshape(g[:,:,:,:,2],(len(time),par['nkx0'],par['nky0'],par['nkz0']))
        rck = np.reshape(g[:,:,:,:,3],(len(time),par['nkx0'],par['nky0'],par['nkz0']))
        return time,lwk,lck,rwk,rck

def load_xi(lpath):
    """  This method can only be run after getv has been called at least once to save the v_xyz.dat file_name                                                                                                                                                                                                                  
    This quickly loads the v array which will have indices [time,kx,ky,kz, x/y/z]"""
    read_parameters(lpath)
    if lpath==None:
        lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    time = np.load(lpath+'/timexi.npy')
    xi=np.memmap(lpath+'/xi_xyz.dat',dtype='float64',mode='r',shape=tuple(np.load(lpath+'/xishape.npy')))
    return time, xi

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
    freq = fftshift(fftfreq(time.shape[-1],d=.01))
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
    plt.xlim(-6.3,6.3)
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

def plot_nls(lpath,ix,iy,iz,ind,show=True,ask=True,tmax=2000000):
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
            optlist.append(getopt(lpath,opt,tmax=tmax))
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

def plot_energy(lpath,ntp,show=True,log=False,rescale=True,xb=1,tmax=2000000):
    """ Plots Scalars Written in energy_out.dat """

    read_parameters(lpath)
    if os.path.isfile(lpath+'/dumen.txt'):
        timeen,enval = load_energy(lpath)
    else:
        timeen,enval = getenergy(lpath,tmax=tmax)

    shapes = {1:(1,1),2:(2,1),3:(2,2),4:(2,2),5:(2,3),6:(2,3),7:(3,3),8:(3,3),9:(3,3)}
    s = shapes[ntp+1]
    labels = {0:'Energy',1:'Magnetic Helicity',2:'MH Error',3:'MH Bound',4:'Canonical Helicity',5:'CH Error',6:'CH Bound',7:'Kinetic Energy',8:'Magnetic Energy',9:'MH Corr'}
    fnames = {0:'energy',1:'maghcty',4:'canhcty',7:'spliten',9:'nlp'}
    
    if not os.path.exists(lpath + '/eplots/'):
        os.mkdir(lpath + '/eplots/')
    xx = len(timeen)
    timeen = timeen[range(0,xx,xb)]
    enval = enval[range(0,xx,xb),:]
    plts = [0,1,4,7,9]
    for i in plts[:ntp+1]:
        if i == 7:
            fig,ax = plt.subplots(2)
            ax[0].plot(timeen,enval[:,i])
            ax[0].set_ylabel(labels[i])
            ax[1].plot(timeen,enval[:,i+1])
            ax[1].set_ylabel(labels[i+1])
            ax[0].set_xlabel('Time (1/$\omega_c$)')
            ax[1].set_xlabel('Time (1/$\omega_c$)')
            fig.suptitle('Kinetic and Magnetic Energies')
            plt.savefig(lpath+'/eplots/'+fnames[i])
            if show == True:
                plt.show()
            plt.close()
        elif not log:
            fig,ax = plt.subplots(1)

#            if i == 1:
#                rr = enval[:,2]
#            elif i == 4:
#                rr  = np.sqrt(enval[:,2]**2 + enval[:,5]**2)
            if i == 0 or i == 9:
                ev = enval[:,i]
            elif par['mhc']:
                ev = (enval[:,i]+enval[:,9])/enval[0,i+2]
                ax.plot(timeen,ev,"b")
            else:
                ev = enval[:,i]/enval[0,i+2]
            r = np.max(ev) - np.min(ev)
            # if r < 10.0**(-12.0):
            #    ax.plot(timeen,ev*10**(12.0),'b')
            #    ax.set_ylabel(labels[i]+" ($10^{-12})V_A^2$)")
            #    ax.set_ylim(bottom=10.0**(12.0)*(ev[0]-1.2*r),top=10.0**(12.0)*(ev[0]+1.2*r))
            # else:
            ax.plot(timeen,ev,'b')
            ax.set_ylabel(labels[i])
            ax.set_ylim(bottom=ev[0]-1.2*r,top=ev[0]+1.2*r)
            ax.yaxis.set_major_locator(MultipleLocator(base=r*0.4,offset=ev[0]))
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.set_xlabel('Time (1/$\omega_c$)')
            fig.suptitle(labels[i])
            fig.tight_layout()
            if np.amax(enval) > 10 ** 5:
                ax.set_yscale('log')
            plt.savefig(lpath+'/eplots/'+fnames[i])
            if show == True:
                plt.show()
            plt.close()
        else:
            fig,ax = plt.subplots(1)
            ax.plot(timeen,np.abs(enval[:,i]))
            ax.set_xlabel('Time (1/$\omega_c$)')
            ax.set_ylabel(labels[i])
            ax.set_yscale('log')
            fig.suptitle(labels[i])
            plt.savefig(lpath+'/eplots/'+fnames[i])
            if show == True:
                plt.show()
            plt.close()
    
    return timeen,enval

def plot_enspec(lpath,npt=1,zz=-1,show=True,log=False,linplot=False,newload=False,fullspec=False,old=True,tmaxfac=1,tmax=2000000):
    read_parameters(lpath)
    kx,ky,kz = get_grids()

    if not old:
        t,bkt = load_b(lpath)
        t,vkt = load_v(lpath)
    elif os.path.isfile(lpath+'/dumenspec.txt'):
        t,ekvf,ekbf = load_energyspec(lpath)
        ekvm = ekvf[:,1:,1:,1:]
        ekbm = ekbf[:,1:,1:,1:]
        ekm = ekvm+ekbm
    else:
        t,ekvf,ekbf = getenergyspec(lpath,tmax=tmax)
        ekvm = ekvf[:,1:,1:,1:]
        ekbm = ekbf[:,1:,1:,1:]
        ekm = ekvm+ekbm
    # Obtain needed phase information to correct initial spectrum
    if not old:
        bths = np.arctan2(np.abs(bkt[0,:,:,1:,1]),np.abs(bkt[0,:,:,1:,0]))
        vths = np.arctan2(np.abs(vkt[0,:,:,1:,1]),np.abs(vkt[0,:,:,1:,0]))
        bphixs = np.angle(bkt[0,:,:,1:,0])
        bphiys = np.angle(bkt[0,:,:,1:,1])
        vpsixs = np.angle(vkt[0,:,:,1:,0])
        vpsiys = np.angle(vkt[0,:,:,1:,1])
        badj = np.zeros(np.shape(bkt[0,:,:,1:,0]))
        vadj = np.zeros(np.shape(vkt[0,:,:,1:,0]))
        for i in range(0,par['nkx0']):
            for j in range(0,par['nky0']):
                for k in range(1,par['nkz0']):
                    badj = 1+(kx[i]**2 * np.cos(bths)**2 + ky[j]**2 * np.sin(bths)**2 + kx[i]*ky[j]*np.sin(2*bths)*np.cos(bphixs-bphiys))/(kz[k]**2)
                    vadj = 1+(kx[i]**2 * np.cos(vths)**2 + ky[j]**2 * np.sin(vths)**2 + kx[i]*ky[j]*np.sin(2*vths)*np.cos(vpsixs-vpsiys))/(kz[k]**2)
            print(i)
    fmts = ["r8","ks","ro","b^","m*"]
    fig1,ax1 = plt.subplots(1)
    fig2,ax2 = plt.subplots(1)
    fig3,ax3 = plt.subplots(1)
    kmag = np.sqrt(np.log(np.tensordot(np.tensordot(np.exp(kx[1:]**2),np.exp(ky[1:]**2),axes=0),np.exp(kz[1:]**2),axes=0)))
    j = 4

    ttp = np.linspace(0,(np.size(t)-1)/tmaxfac,num=npt)
    ttp = ttp.astype(int)

    for i in ttp[::-1]:
        if not newload:
            ekb = 0.5*(np.abs(bkt[i,1:,1:,1:,0])**2 + np.abs(bkt[i,1:,1:,1:,1])**2+np.abs(bkt[i,1:,1:,1:,2])**2)/badj
            ekv = 0.5*(np.abs(vkt[i,1:,1:,1:,0])**2 + np.abs(vkt[i,1:,1:,1:,1])**2+np.abs(vkt[i,1:,1:,1:,2])**2)/vadj
            if fullspec:
                ekb += 0.5*np.abs(bkt[i,:,:,1:,2])**2
                ekv += 0.5*np.abs(vkt[i,:,:,1:,2])**2
            ek = ekb+ekv
        # elif log:
        #    ekb = np.log10(ekbm[i,:,:,:])
        #    ekv = np.log10(ekvm[i,:,:,:])
        #    ek = np.log10(ekm[i,:,:,:])
        else:
            ekb = ekbm[i,:,:,:]
            ekv = ekvm[i,:,:,:]
            ek = ekm[i,:,:,:]
        # print(np.amax(np.abs(ekm[i,:,:,:]-ekm[0,:,:,:])),np.argmax(np.abs(ekm[i,:,:,:]-ekm[0,:,:,:])))

        #if (not newload) and log:
        #    ekb[:,:,:] = np.log10(ekb[:,:,:])
        #    ekv[:,:,:] = np.log10(ekv[:,:,:])
        #    ek = np.log10(ek[:,:,:])
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

        if linplot:
            ax1.plot(1/(x[a]**par["init_kolm"]),yb[a],fmts[np.mod(i,5)],label=np.format_float_positional(t[i],2),markersize=1)
            ax2.plot(1/(x[a]**par["init_kolm"]),yv[a],fmts[np.mod(i,5)],label=np.format_float_positional(t[i],2),markersize=1)
            ax3.plot(1/(x[a]**par["init_kolm"]),yt[a],fmts[np.mod(i,5)],label=np.format_float_positional(t[i],2),markersize=1)
        else:
            ax1.plot(x[a],yb[a],fmts[np.mod(j,5)],label=np.format_float_positional(t[i],2)+" $\omega_c^{-1}$",markersize=1)
            ax2.plot(x[a],yv[a],fmts[np.mod(j,5)],label=np.format_float_positional(t[i],2)+" $\omega_c^{-1}$",markersize=1)
            ax3.plot(x[a],yt[a],fmts[np.mod(j,5)],label=np.format_float_positional(t[i],2)+" $\omega_c^{-1}$",markersize=1)
        
        j = j - 1

    ax1.legend(loc=3)
    ax2.legend(loc=3)
    ax3.legend(loc=3)
    ax1.set_ylim(10**(-10),1)
    ax2.set_ylim(10**(-10),1)
    ax3.set_ylim(10**(-10),1)
    if log:
       ax1.set_yscale("log")
       ax1.set_xscale("log")
       ax2.set_yscale("log")
       ax2.set_xscale("log")
       ax3.set_yscale("log")
       ax3.set_xscale("log")

    label=" Energy Spectrum"
    if (zz == -1):
        fig1.suptitle("3D Energy Spectrum")
        fig2.suptitle("3D Energy Spectrum")
        fig3.suptitle("3D Energy Spectrum")
        ax1.set_ylabel("Magnetic"+label)
        ax2.set_ylabel("Kinetic"+label)
        ax3.set_ylabel("Total"+label)
        if linplot:
            fig1.supxlabel("1/k^{}".format(par["init_kolm"]))
            fig2.supxlabel("1/k^{}".format(par["init_kolm"]))
            fig3.supxlabel("1/k^{}".format(par["init_kolm"]))
        else:
            fig1.supxlabel("|k| ($d_i^{-1}$)")
            fig2.supxlabel("|k| ($d_i^{-1}$)")
            fig3.supxlabel("|k| ($d_i^{-1}$)")

    else:
        fig1.suptitle("Perpendicular Energy Spectrum kz = "+np.format_float_positional(kz[zz+1],2))
        fig2.suptitle("Perpendicular Energy Spectrum kz = "+np.format_float_positional(kz[zz+1],2))
        fig3.suptitle("Perpendicular Energy Spectrum kz = "+np.format_float_positional(kz[zz+1],2))
        ax1.set_ylabel("Magnetic Energy Spectrum"+label)
        ax2.set_ylabel("Kinetic Energy Spectrum"+label)
        ax3.set_ylabel("Total Energy Spectrum"+label)
        fig1.supxlabel(label+" k")
        fig2.supxlabel(label+" k")
        fig3.supxlabel(label+" k")
    if not os.path.exists(lpath + '/eplots/'):
        os.mkdir(lpath + '/eplots/')
    fig1.savefig(lpath+'/eplots/benspec'+str(zz+1)+'.png')
    fig2.savefig(lpath+'/eplots/venspec'+str(zz+1)+'.png')
    fig3.savefig(lpath+'/eplots/tenspec'+str(zz+1)+'.png')
    if show == True:
        fig1.show()
        fig2.show()
        fig3.show()
    else:
        plt.close("all")
    
    tt = np.arange(np.size(t))
    sts = ["ks","ro","b^","m*"]
    ttplot = np.floor(np.linspace(0,tt[-1],4)).astype("int")
    fig,ax = plt.subplots(1)
    for i in range(4):
        kperps,ekti = integrated_spectrum_1d(ekm[ttplot[i],:,:,:],lpath)
        ax.plot(kperps,ekti,sts[i],label=np.format_float_positional(t[ttplot[i]],2))
    fig.suptitle("Integrated Total Energy Spectrum")
    ax.set_ylabel("Total Energy Spectrum")
    ax.set_xlabel("$k_{\perp}$")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.legend()
    fig.savefig(lpath+'/eplots/t1denspec.png')
    if show == True:
        plt.show()
    else:
        plt.close()

    return(t,ekm)

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

def plot_bspectrum(lpath,ix,iy,iz,ind,show=True,opt='b'):
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
    print(par['dt_max'])
    if opt=='b':
        time,b=load_b(lpath)
        dt = time[1]-time[0]
        b_k = b[:,ix,iy,iz,ind]
        sp=fftshift(fft(b_k-np.mean(b_k)))    
        freq = fftshift(fftfreq(time.shape[-1],d=dt))
    elif opt=='v':
        time,b=load_v(lpath)
        dt = time[1]-time[0]
        b_k = b[:,ix,iy,iz,ind]
        sp=fftshift(fft(b_k-np.mean(b_k)))
        freq = fftshift(fftfreq(time.shape[-1],d=dt))
    else:
        time,b=load_opt(lpath,opt)
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
    
    ax.set_ylabel('|FFT('+opt+'_%s)|'%ind_string )
    ax.set_xlabel('frequency')
    fig.suptitle('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
    if lpath[-1] != '/':
        lpath =lpath +'/'
    if not os.path.exists(lpath + opt+'spectra/'):
        os.mkdir(lpath + opt+'spectra/')
    plt.savefig(lpath+opt+'spectra/'+opt+'spectrum_%s_%d_%d_%d'%(ind_string,ix,iy,iz))
    if show == True:
        plt.show()
    plt.close()
    return freq[peaks]*2*np.pi

def plot_profile_enspec(lpath,ix,iy,iz,tbv='t',show=True,navg=20):
    read_parameters(lpath)
    kx,ky,kz = get_grids()
    time,ebk,evk =load_energyspec(lpath)
    ebkt = ebk[:,ix,iy,iz]
    evkt = evk[:,ix,iy,iz]
    etkt = ebkt + evkt
    profs = {'t':etkt,'b':ebkt,'v':evkt}
    labels = {'t':'Total ','b':'Magnetic ','v':'Kinetic '}
    ekt = profs[tbv]
    sampled_profile = np.array([])
    tt = time[np.arange(0,np.size(time),navg)]
    mover = 0
    while mover < np.size(time):
        sampled_profile = np.append(sampled_profile,np.mean(ekt[mover:mover+navg]))
        mover += navg
    fig,ax = plt.subplots(1)
    ax.plot(time,ekt,'k:',label='Full Profile')
    ax.plot(tt,sampled_profile,'bs',label='Averaged')
    ax.set_ylabel(labels[tbv]+'Energy Spectrum')
    ax.set_ylim(0,1.2*np.amax(ekt))
    ts, cts = profile_chartime_numdiv(tt,sampled_profile)
    fig2,ax2 = plt.subplots(1)
    ax2.plot(ts,cts,label='Derivative Times')
    popt,pcov = profile_chartime_expfit(tt,sampled_profile)
    yy = func_exp(tt,*popt)
    ax.plot(tt,yy,'r',label="Exponential Fit Rate "+np.format_float_positional(np.abs(popt[1]),2))
    ax.legend()
    ax2.set_ylabel("Derivative Times")
    fig.suptitle('$ k_x, k_y, k_z$ = %1.2f,%1.2f,%1.2f ($d_i^{-1}$)'%(kx[ix],ky[iy],kz[iz]))
    fig2.suptitle('$ k_x, k_y, k_z$ = %1.2f,%1.2f,%1.2f ($d_i^{-1}$)'%(kx[ix],ky[iy],kz[iz]))
    fig.supxlabel('Time ($\omega_c^{-1}$)')
    fig2.supxlabel('Time ($\omega_c^{-1}$)')
    if lpath[-1] != '/':
        lpath = lpath + '/'
    if not os.path.exists(lpath + 'esps/'):
        os.mkdir(lpath + 'esps/')
    plt.savefig(lpath+'esps/esp_%s_%d_%d_%d'%(tbv,ix,iy,iz))
    if show == True:
        plt.show()
    plt.close()
    return tt,sampled_profile

def plot_profile_xi(lpath,ix,iy,iz,tbv='t',show=True,newload=False):
    read_parameters(lpath)
    kx,ky,kz = get_grids()
    if newload==False:
        wp,wm = analytical_omega(lpath,ix,iy,iz)
    time,xik =load_xi(lpath)
    xikt = xik[:,ix-1,iy-1,iz-1]
    fig,ax = plt.subplots(1)
    ax.plot(time,xikt,'ks',label='Alfven Whistler',markersize=1)
    ax.plot(time,xikt*wp,'rs',label='MHD',markersize=1)
    ax.plot(time,xikt*wp/wm,'bs',label='Alfven Cyclotron',markersize=1)
    ax.set_ylabel('Nonlinearity Parameter')
    ax.set_ylim(10**(-4),10)
    ax.set_yscale('log')
    ax.legend()
    fig.suptitle('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
    fig.supxlabel('time (1/wc)')
    if lpath[-1] != '/':
        lpath = lpath + '/'
    if not os.path.exists(lpath + 'xis/'):
        os.mkdir(lpath + 'xis/')
    plt.savefig(lpath+'xis/xi_%s_%d_%d_%d'%(tbv,ix,iy,iz))
    if show == True:
        plt.show()
    plt.close()
    return time,xikt

def profile_chartime_numdiv(time,ekt):
    de = ekt[2:]-ekt[:-2]
    dt = 2*(time[2:]-time[:-2])
    ee = -ekt[1:-1] + ekt[-1]*np.ones(np.size(ekt[1:-1]))
    ts = time[1:-1]
    # print(ee)
    # print(de/dt)
    cts = np.abs(1/ee * de/dt)
    return ts, cts

def func_exp(x,b,w):
    return( b * np.exp( -w*x))

def profile_chartime_expfit(time,ekt):
    tg = -np.log(ekt[-1]/ekt[0])/(time[-1]-time[0])
    # print(tg)
    popt,pcov = spo.curve_fit(func_exp,time,ekt,p0=[ekt[0],tg],method='dogbox')
    return(popt,pcov)

    
def plot_xispec(lpath,npt=1,zz=-1,show=True,log=True,tmaxfac=1,tmax=2000000,tt=1,newload=False):
    read_parameters(lpath)
    kx,ky,kz = get_grids()

    if os.path.isfile(lpath+'/dumxi.txt'):
        t,xik = load_xi(lpath)
    else:
        t,xik = getxi(lpath,tmax=tmax)
    
    fmts = ["k","m","b","g","r"]
    fig1,ax1 = plt.subplots(1)
    kmag = np.sqrt(np.log(np.tensordot(np.tensordot(np.exp(kx[1:]**2),np.exp(ky[1:]**2),axes=0),np.exp(kz[1:]**2),axes=0)))
    j = 0

    ttp = np.linspace(0,(np.size(t)-1)/tmaxfac,num=npt)
    ttp = ttp.astype(int)

    if newload==False:
        kxx,kyy,kzz = np.meshgrid(ky[1:],kx[1:],kz[1:])
        wp = kzz * np.sqrt((1+0.5*kmag**2) + np.sqrt((1+0.5*kmag**2)**2 - 1))
        wm = kzz * np.sqrt((1+0.5*kmag**2) - np.sqrt((1+0.5*kmag**2)**2 - 1))
    for i in ttp[tt-1:tt]:
        xit = xik[i,:,:,:]
        if (zz < 1):
            x = np.reshape(kmag,np.size(kmag))
            a = np.argsort(x)
            y = np.reshape(xit,np.size(kmag))
            if newload==False:
                wpp = np.reshape(wp,np.size(kmag))
                wmm = np.reshape(wm,np.size(kmag))
        else:
            x = np.reshape(kmag[:,:,zz],np.size(kmag[:,:,zz]))
            a = np.argsort(x)
            y = np.reshape(xit[:,:,zz],np.size(kmag[:,:,0]))
            if newload==False:
                wpp = np.reshape(wp[:,:,zz],np.size(kmag[:,:,0]))
                wmm = np.reshape(wm[:,:,zz],np.size(kmag[:,:,0]))

        ax1.plot(x[a],y[a],"ks",label="Alfven Whister "+np.format_float_positional(t[i],2),markersize=1)
        ax1.plot(x[a],y[a]*wpp/wmm,"bs",label="Alfven Cyclotron "+np.format_float_positional(t[i],2),markersize=1)
        ax1.plot(x[a],y[a]*wpp,"rs",label="MHD "+np.format_float_positional(t[i],2),markersize=1)
        print(np.argmax(y[a]*wpp/wmm))
        j = j - 1

    ax1.legend()
    ax1.set_ylim(10**(-4),10)

    if log:
       ax1.set_yscale("log")
       ax1.set_xscale("log")

    label=""
    if (zz < 1):
        fig1.suptitle("Nonlinearity Parameter Spectrum")
        ax1.set_ylabel("Nonlinearity Parameter")
        fig1.supxlabel("k")
    else:
        fig1.suptitle("Nonlinearity Parameter Spectrum kz = {}".format(np.format_float_positional(kz[zz],2)))
        ax1.set_ylabel("Nonlinearity Parameter")
        fig1.supxlabel("k")

    if not os.path.exists(lpath + '/xis/'):
        os.mkdir(lpath + '/xis/')
    fig1.savefig(lpath+'/xis/xispec'+str(max(zz,0))+'.png')
    if show == True:
        fig1.show()
    else:
        plt.close("all")

#    kperps,xit1 = integrated_spectrum_1d(xik[-1,:,:,:],lpath)
#
#    fig,ax = plt.subplots(1)
#    ax.plot(kperps,xit1)
#    fig.suptitle("Integrated Nonlinearity Parameter Spectrum")
#    ax.set_ylabel("Nonlinearity Spectrum")
#    ax.set_xlabel("$k_{\perp}$")
#    ax.set_yscale("log")
#    ax.set_xscale("log")
#    fig.savefig(lpath+'/xis/t1dxispec.png')
#    if show == True:
#        plt.show()
#    else:
#        plt.close()
#
    return(t,xik)

def integrated_spectrum_1d(spec,lpath,v1=False,v2=True):
    read_parameters(lpath)
    kx,ky,kz = get_grids()
    kxx = kx[1:]
    kyy = ky[1:]
    kzz = kz[1:]
    kygrid,kxgrid = np.meshgrid(kyy,kxx)
    kpgrid = np.sqrt(kxgrid**2 + kygrid**2)
    kps = kpgrid.flatten()
    kperps = np.unique(kps)

    if v1:
        spec2d = np.zeros([np.size(kperps),np.size(kzz)])
        for j in np.arange(np.size(kzz)):
            for i in np.arange(np.size(kperps)):
                specperp = spec[:,:,j]
                w = np.argwhere(kpgrid==kperps[i])
                spec2d[i,j] = kperps[i]*np.sum(specperp[w])
        spec1d = np.sum(spec2d,axis=1)

    elif v2:
        specperp = np.sum(spec,axis=2)
        a = np.argsort(kps)
        kps_sorted = kps[a]
        spiral = specperp.flatten()[a]

        cumulative_sum = []
        kperps = np.linspace(2*kx[1],np.amax(kperps),num=int(np.amax(kperps)//(2*kx[1])))

        for ring in kperps:
            test = np.nonzero(kps_sorted < ring)
            cumulative_sum.append(np.sum(spiral[test]))

        spec1d = np.array(cumulative_sum[1:]) - np.array(cumulative_sum[:-1])

    return(kperps[:-1],spec1d)

def nlparam(lpath,tt):

    read_parameters(lpath)
    kx,ky,kz = get_grids()    
    if os.path.isfile(lpath+'/dumenspec.txt'):
        t,ekvf,ekbf = load_energyspec(lpath)
        ekvm = ekvf[:,1:,1:,1:]
        ekbm = ekbf[:,1:,1:,1:]
        ekm = ekvm+ekbm
    else:
        t,ekvf,ekbf = getenergyspec(lpath,tmax=tmax)
        ekvm = ekvf[:,1:,1:,1:]
        ekbm = ekbf[:,1:,1:,1:]
        ekm = ekvm+ekbm

    kyy,kxx,kzz = np.meshgrid(ky[1:],kx[1:],kz[1:])
    kmags = np.sqrt(kxx**2 +kyy**2 + kzz**2)
    kperps = np.sqrt(kxx**2 + kyy**2)
    
    sq = np.sqrt(1+kmags**2 / 4)
    whist = kzz * (kmags/2 + sq)
    cyclo = kzz * (sq - kmags/2)
    mhd = kzz

    rmse = np.sqrt(2*ekm[tt,:,:,:])
    xi_mhd = rmse*kperps/mhd
    xi_whist = rmse*kperps/whist
    xi_cyclo = rmse*kperps/cyclo

    plt.plot(kmags.flatten(),xi_cyclo.flatten(),'bs',label="Cyclotron",markersize=1)
    plt.plot(kmags.flatten(),xi_mhd.flatten(),'ks',label="MHD",markersize=1)
    plt.plot(kmags.flatten(),xi_whist.flatten(),'rs',label="Whister",markersize=1)
    plt.ylabel("Nonlinearity Parameter")
    plt.xlabel("|k| ($d_i^{-1}$)")
    plt.title("3D Nonlinearity Parameter Spectrum t = "+np.format_float_positional(t[tt],1)+" ($\omega_c^{-1}$)")
    plt.ylim(10**(-5),10**1)
    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    plt.show()

    return(kperps,xi_cyclo,xi_mhd,xi_whist)

def mode_break(lpath,show=False,tmax=200000):

    read_parameters(lpath)
    kx,ky,kz = get_grids()

    if os.path.isfile(lpath+'/dummodespec.txt'):
        t,lwk,lck,rwk,rck = load_energyspec(lpath,mode=True)
        lwkm = lwk[:,1:,1:,1:]
        lckm = lck[:,1:,1:,1:]
        rwkm = rwk[:,1:,1:,1:]
        rckm = rck[:,1:,1:,1:]
    else:
        t,lwk,lck,rwk,rck = getenergyspec(lpath,tmax=tmax,mode=True)
        lwkm = lwk[:,1:,1:,1:]
        lckm = lck[:,1:,1:,1:]
        rwkm = rwk[:,1:,1:,1:]
        rckm = rck[:,1:,1:,1:]

    print(t)
    fmts = ['b--','b:','r--','r:']
    labels = ['+ Helicity Whistler','+ Helicity Cyclotron','- Helicity Whistler','- Helicity Cyclotron']

    mode_ks = np.stack((lwkm,lckm,rwkm,rckm))
    mode_profs = np.sum(mode_ks,axis=(2,3,4))
    plt.figure(1)
    for i in range(4):
        plt.plot(t,mode_profs[i,:],fmts[i],label=labels[i],markersize=1)
    plt.xlabel('Time ($\omega_c^{-1}$)')
    plt.ylabel('Mode Energy')
    plt.title('Energy Distribution in Hall MHD Modes vs Time')
    plt.yscale('log')
    plt.ylim(10**(-10),10**1)
    plt.legend()

    if not os.path.exists(lpath+'/eplots/'):
        os.mkdir(lpath+'/eplots/')
    plt.savefig(lpath+'/eplots/modeev.png')
    if show == True:
        plt.show()
    plt.close()

    # Other plots: 1D spectra

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
    plt.close()
    return (t,lwkm,lckm,rwkm,rckm)

def steadystate(lpath,graph=True):
    '''
    Return an Linf metric of the final relative change in the recorded spectrum
    '''

    if os.path.isfile(lpath+'/dumenspec.txt'):
        t,ekvf,ekbf = load_energyspec(lpath)
        ekvm = ekvf[:,1:,1:,1:]
        ekbm = ekbf[:,1:,1:,1:]
        ekm = ekvm+ekbm
    else:
        t,ekvf,ekbf = getenergyspec(lpath,tmax=tmax)
        ekvm = ekvf[:,1:,1:,1:]
        ekbm = ekbf[:,1:,1:,1:]
        ekm = ekvm+ekbm

    tspectrum = ekm
    prev = tspectrum[:-1,:,:,:]
    current = tspectrum[1:,:,:,:]

    ignore = np.nonzero(tspectrum[:-1,:,:,:])
    prev[ignore] = 1.0
    current[ignore] = 0.0

    met = np.abs(current/prev - 1.0)
    
    answer = np.amax(met,axis=(1,2,3))

    if graph:
        plt.plot(t[1:],answer,"b*")
        plt.xlabel("Time ($\omega_c^{-1}$)")
        plt.ylabel("Max Relative Change in Energy Spectrum")
        plt.title("Evolution of Relative Energy Spectrum Change")
        plt.savefig(lpath+"/eplots/steadystate.png")

    return(answer[-1])

def planeplotter(lpath,t,normvec=2,planenum=0,show=False):
    """
    Plots b and v fields in a specified plane normal to the x, y, or z directions
    Arguments:
    lpath - file of run to show
    t - time point to consider
    normvec - number (x 0, y 1, z 2) of the normal vector to the plane 
    planenum - grid point along normvec direction to show a plane for
    Example: normvec = 2, planenum = 0 will plot the bv fields in the plane z = 0
    """
    read_parameters(lpath)
    kx,ky,kz = get_grids()

    xx = np
    if os.path.isfile(lpath+'/v_xyz.dat'):
        time,b = load_b(lpath)
        time,v = load_v(lpath)
    else:
        time,b = getb(lpath)
        time,v = getv(lpath)

    # get real fields
    # note that irfftn assumes the last dimension is the one to be doubled, opposite of the current splitx and FFTW 90
    # so we need to transpose

    bt = b[t,:,:,:,:]
    vt = v[t,:,:,:,:]
    bint = np.transpose(bt,axes=(2,1,0,3))
    vint = np.transpose(vt,axes=(2,1,0,3))
    brealint = irfftn(bint,s=[par['nkz0'],par['nky0'],2*par['nkx0']],axes=(0,1,2))
    vrealint = irfftn(vint,s=[par['nkz0'],par['nky0'],2*par['nkx0']],axes=(0,1,2))
    br = np.transpose(brealint,axes=(2,1,0,3))
    vr = np.transpose(vrealint,axes=(2,1,0,3))

    bplane = np.take(br,indices=planenum,axis=normvec)
    vplane = np.take(vr,indices=planenum,axis=normvec)

    xs = np.pi / par['nkx0'] * np.arange(2*par['nkx0']) / par['kxmin']
    ys = 2*np.pi / par['nky0'] * np.arange(par['nky0']) / par['kymin']
    zs = 2*np.pi / par['nkz0'] * np.arange(par['nkz0']) / par['kzmin']

    print(np.size(xs))
    print(np.size(ys))
    
    inds = ['x','y','z']

    ii = np.mod(normvec+1,3)
    jj = np.mod(normvec+2,3)

    ggrid = {0:xs,1:ys,2:zs}
    print("xgrid ",np.size(ggrid[0]))
    print("ygrid ",np.size(ggrid[1]))

    if t == -1:
        ts = 'fin'
    else:
        ts = str(t)
    
    plt.figure()
    plt.streamplot(ggrid[ii],ggrid[jj],bplane[:,:,ii],bplane[:,:,jj],broken_streamlines=False)
    plt.xlabel(inds[ii]+" ($d_i$)")
    plt.ylabel(inds[jj]+" ($d_i$)")
    plt.xlim(0,ggrid[ii][0]+ggrid[ii][-1])
    plt.ylim(0,ggrid[jj][0]+ggrid[jj][-1])
    plt.title('Magnetic Field Lines Time '+time[t])
    if lpath[-1] != '/':
        lpath = lpath + '/'
    if not os.path.exists(lpath + 'bvs/'):
        os.mkdir(lpath + 'bvs/')        
    plt.savefig(lpath+'bvs/bfield'+ts+'.png')
    if show:
        plt.show()

    plt.figure()
    plt.streamplot(ggrid[ii],ggrid[jj],vplane[:,:,ii],vplane[:,:,jj],broken_streamlines=False)
    plt.xlabel(inds[ii]+" ($d_i$)")
    plt.ylabel(inds[jj]+" ($d_i$)")
    plt.title('Velocity Field Lines Time '+time[t])
    if not os.path.exists(lpath + 'bvs/'):
        os.mkdir(lpath + 'bvs/')
    plt.savefig(lpath+'bvs/vfield'+ts+'.png')
    if show:
        plt.show()

    return ggrid[ii],ggrid[jj],bplane,vplane

def enheldev(lpath):
    te,e = plot_energy(lpath,3,show=False)
    dt = te[-1]-te[0]
    de = e[:,0]-e[0,0]
    dmh = (e[:,1]+e[:,-3]) - (e[0,1]+e[0,-3])
    dch = (e[:,4]+e[:,-3]) - (e[0,4]+e[0,-3])
    
    print(dt)
    
    # print(de[0])
    # print(de[1]+de[-3])
    # print(de[4]+de[-3])
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

    if os.path.isfile(lpath+'/v_xyz.dat'):
        time,bk = load_b(lpath)
        time,vk = load_v(lpath)
    else:
        time,bk = getb(lpath)
        time,vk = getv(lpath)
    # decompose bk,vk into constituent modes

    Ky,Kx,Kz = np.meshgrid(ky,kx,kz)
    kmags = np.sqrt(Kx**2 + Ky**2 + Kz**2)
    alpha_lw = -kmags/2 - np.sqrt(1+kmags**2 /4)
    alpha_lc = -kmags/2 + np.sqrt(1+kmags**2 /4)

    Kvec = np.zeros([par['nkx0'],par['nky0'],par['nkz0'],3],dtype='complex64')
    Zvec = np.zeros([par['nkx0'],par['nky0'],par['nkz0'],3],dtype='complex64')
    Kvec[:,:,:,0] = Kx
    Kvec[:,:,:,1] = Ky
    Kvec[:,:,:,2] = Kz
    Zvec[:,:,:,2] = 1

    pceig = np.zeros([par['nkx0'],par['nky0'],par['nkz0'],3],dtype='complex64')
    for i in range(par['nkx0']):
        for j in range(par['nky0']):
            for k in range(par['nkz0']):
                pceig[i,j,k,:] = np.cross(Kvec[i,j,k,:],Zvec[i,j,k,:])+1/kmags[i,j,k] * 1.0j * np.cross(Kvec[i,j,k,:],np.cross(Kvec[i,j,k,:],Zvec[i,j,k,:]))
    pceig = pceig / (np.sqrt(2) * np.sqrt(Kx[:,:,:,None]**2 + Ky[:,:,:,None]**2))
    pceig[0,0,:,:] = 0
    bkt = bk[-1,:,:,:,:]
    vkt = vk[-1,:,:,:,:]

    lwk = np.zeros([par['nkx0'],par['nky0'],par['nkz0']],dtype='complex64')
    lck = np.zeros([par['nkx0'],par['nky0'],par['nkz0']],dtype='complex64')
    rwk = np.zeros([par['nkx0'],par['nky0'],par['nkz0']],dtype='complex64')
    rck = np.zeros([par['nkx0'],par['nky0'],par['nkz0']],dtype='complex64')

    for i in range(par['nkx0']):
        for j in range(par['nky0']):
            for k in range(par['nkz0']):
                lwk[i,j,k] = np.dot(np.conj(pceig[i,j,k,:]),alpha_lw[i,j,k]*bkt[i,j,k,:]+vkt[i,j,k,:])/np.sqrt(alpha_lw[i,j,k]**2+1)
                lck[i,j,k] = np.dot(np.conj(pceig[i,j,k,:]),alpha_lc[i,j,k]*bkt[i,j,k,:]+vkt[i,j,k,:])/np.sqrt(alpha_lc[i,j,k]**2+1)
                rwk[i,j,k] = np.dot((pceig[i,j,k,:]),-alpha_lw[i,j,k]*bkt[i,j,k,:]+vkt[i,j,k,:])/np.sqrt(alpha_lw[i,j,k]**2+1)
                rck[i,j,k] = np.dot((pceig[i,j,k,:]),-alpha_lc[i,j,k]*bkt[i,j,k,:]+vkt[i,j,k,:])/np.sqrt(alpha_lc[i,j,k]**2+1)
    
    labels = ['+ Helicity Whistler','+ Helicity Cyclotron','- Helicity Whistler','- Helicity Cyclotron']
    stname = ["phw","phc","nhw","nhc"]

    xs = np.pi / par['nkx0'] * np.arange(2*par['nkx0']) / par['kxmin']
    ys = 2*np.pi / par['nky0'] * np.arange(par['nky0']) / par['kymin']
    zs = 2*np.pi / par['nkz0'] * np.arange(par['nkz0']) / par['kzmin']
    
    modespecs = [lwk,lck,rwk,rck]

    fig,ax = plt.subplots(2,2)
    ax = ax.flatten()
    fig2,ax2 = plt.subplots(2,2)
    ax2 = ax2.flatten()
    
    for I,ms in enumerate(modespecs):
        mm = convert_spec_to_real(lpath,ms)
            # assume axisymmetric for transverse struct fn - test this later
        # very small - lets rescale to check Parseval's theorem print(np.amax(mm**2)); 8 pi^3 sum(ms**2) = sum(mm**2)/N
        mm *= np.sqrt(8*np.pi**3 * np.sum(np.abs(ms)**2)*np.size(mm)/(np.sum(mm**2)))
        
        if par["splitx"]:
            nx = 2*par["nkx0"]
        else:
            nx = par["nkx0"]
            
        str_perp = np.zeros(nx)
        for j in range(nx//2):
            if j != 0:
                str_perp[j] = np.average((mm[j:,:,:]-mm[:-j,:,:])**2)
            
            # parallel structure function
        str_par = np.zeros(par["nkz0"])
        for k in range(par["nkz0"]//2):
            if k != 0:
                str_par[k] = np.average((mm[:,:,k:]-mm[:,:,:-k])**2)

        pt = np.nonzero(str_perp)
        pz = np.nonzero(str_par)

        ax[I].plot(xs[pt],str_perp[pt],"--",label="Transverse")
        ax[I].plot(zs[pz],str_par[pz],":",label="Parallel")
        if (I==2 or I == 3):
            ax[I].set_xlabel("r (d_i)")

        
        ax2[I].plot(str_perp[pt],"--",label="Transverse")
        ax2[I].plot(str_par[pz],":",label="Parallel")
        if (I==2 or I == 3):
            ax2[I].set_xlabel("r (Grid Position)")
        ax[I].set_title(labels[I])
        ax[I].set_xscale("log")
        ax[I].set_yscale("log")
        ax2[I].set_title(labels[I])
        ax2[I].set_xscale("log")
        ax2[I].set_yscale("log")
        if (I < 2):
            ax[I].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False,labeltop=False)
            ax2[I].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False,labeltop=False)
        ax[I].legend(loc=4)
        ax2[I].legend(loc=4)
    fig.suptitle("Transverse and Longitudinal Structure Functions")
    fig2.suptitle("Transverse and Longitudinal Structure Functions")
    if lpath[-1] != '/':
        lpath = lpath + '/'
    if not os.path.exists(lpath + 'eplots/'):
        os.mkdir(lpath + 'eplots/')    
    fig.savefig(lpath + 'eplots/'+"stfns")
    fig2.savefig(lpath + 'eplots/'+"stfns2")
    plt.close()
        
def convert_spec_to_real(lpath,spectra):
    """
    Performs IFFTs needed to get real space values on one time b or mode
    """

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

def lastbv(lpath):
    if lpath[-1] != '/':
        lpath = lpath + '/'
    read_parameters(lpath)
    memstep = 8 + 16 * 3 * par['nkx0']*par['nky0']*par['nkz0']
    size = os.path.getsize(lpath+'b_out.dat')
    Nt = size//memstep
    print(Nt)
    
    f = open(lpath+'b_out.dat','rb')
    tlast = np.fromfile(f,dtype='float64',count=1,offset=(Nt-1)*memstep)
    blast_back = np.empty((3,par['nkz0'],par['nky0'],par['nkx0']))
    blast_back = np.fromfile(f,dtype='complex64',count=3 * par['nkx0']*par['nky0']*par['nkz0'],offset=(Nt-1)*memstep+8)
    f.close()

    print(tlast, blast_back)
    
    g = open(lpath+'v_out.dat','rb')
    vlast_back = np.empty((3,par['nkz0'],par['nky0'],par['nkx0']))
    vlast_back = np.fromfile(g,dtype='complex64',count=3 * par['nkx0']*par['nky0']*par['nkz0'],offset=(Nt-1)*memstep+8)
    g.close()

    print(vlast_back)

    blast = np.reshape(blast_back,(par['nkx0'],par['nky0'],par['nkz0'],3),order='F')
    vlast = np.reshape(vlast_back,(par['nkx0'],par['nky0'],par['nkz0'],3),order='F')

    np.save(lpath+'b_fin',blast)
    np.save(lpath+'v_fin',vlast)

    h = open(lpath+'dumlasts.txt','w')
    h.write("Wrote Lasts")
    h.close()
    
    return(tlast,blast,vlast)

def mode_nlparam(lpath,tt,dim,show=False,tmax=200000):

    read_parameters(lpath)
    kx,ky,kz = get_grids()

    if os.path.isfile(lpath+'/dummodespec.txt'):
        t,lwk,lck,rwk,rck = load_energyspec(lpath,mode=True)
        lwkm = lwk[:,1:,1:,1:]
        lckm = lck[:,1:,1:,1:]
        rwkm = rwk[:,1:,1:,1:]
        rckm = rck[:,1:,1:,1:]
    else:
        t,lwk,lck,rwk,rck = getenergyspec(lpath,tmax=tmax,mode=True)
        lwkm = lwk[:,1:,1:,1:]
        lckm = lck[:,1:,1:,1:]
        rwkm = rwk[:,1:,1:,1:]
        rckm = rck[:,1:,1:,1:]
    
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
    
    plt.ylabel("Nonlinearity Parameter")
    plt.xlabel("|k| ($d_i^{-1}$)")
    plt.title(str(dim)+"D Nonlinearity Parameter Spectrum t = "+np.format_float_positional(t[tt],1)+" ($\omega_c^{-1}$)")
    plt.ylim(10**(-5),10**1)
    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    plt.savefig(lpath+'/eplots/nlpar'+str(dim)+'d'+str(int(t[tt])))
    
    return(0)
