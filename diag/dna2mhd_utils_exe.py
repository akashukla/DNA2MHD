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
from scipy.fft import fft,fftfreq,fftshift


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
     input=np.fromfile(f,dtype='float64',count=1)
     if swap_endian:
         input=input.newbyteorder()
     #print input
     if input==0 or input:
         time = np.append(time,input)
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
     input=np.fromfile(f,dtype='float64',count=1)
     if swap_endian:
         input=input.newbyteorder()
     #print input
     if input==0 or input:
         time = np.append(time,input)
     else:
         continue_read=0

   f.close()
   return time


def getb(lpath=None):
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

def getv(lpath=None):
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

def load_b(lpath):
    """
    This method can only be run after getb has been called at least once to save the b_xyz.dat file_name
    This quickly loads the b array which will have indices [time,kx,ky,kz, x/y/z]
    """
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
    if lpath==None:
        lpath='/scratch/04943/akshukla/dna2mhd_output_0'
    time = np.load(lpath+'/timev.npy')
    vload=np.memmap(lpath+'/v_xyz.dat',dtype='complex64',mode='r',shape=tuple(np.load(lpath+'/vshape.npy')))
    return time, vload


def plot_bv(ix,iy,iz,ind,lpath=None):
    """
    This is an example method that plots the timetraces of b and v at the specified wavevector (kx[ix],ky[iy],kz[iz]).
    ind specifies whether you want the x(0),y(1), or z(2) component.
    """
    timeb,b=load_b(lpath)
    timev,v=load_v(lpath)
    fig,ax=plt.subplots(2)
    ax[0].plot(timeb,np.abs(b[:,ix,iy,iz,ind]))
    ax[0].set_title('b')
    ax[1].plot(timev,np.abs(v[::,ix,iy,iz,ind]))
    ax[1].set_title('v')
    kx,ky,kz=get_grids()
    fig.suptitle('kx,ky,kz = %1.2f,%1.2f,%1.2f'%(kx[ix],ky[iy],kz[iz]))
    plt.show()
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

def plot_vspectrum(lpath,ix,iy,iz,ind):
    """
    ix,iy,iz specifies the wavevector
    ind specifies x/y/z (0/1/2) component
    This is an example method that performs the fft on the real part of v and plots the result.
    It will return an array of the frequencies found.
    *** Right now it seems like freqs need to multiplied by 2*pi to get the right dispersion relation.
        I think this makes sense because w = 2*pi*f
    """
    time,v=load_v(lpath)
    v_k = v[:,ix,iy,iz,ind]
    #plt.plot(time,v_k)
    #plt.show()
    sp=fftshift(fft(v_k-np.mean(v_k)))
    freq = fftshift(fftfreq(time.shape[-1],d=.01))
    peaks,_ = find_peaks(np.abs(sp),threshold=10)
    print(freq[peaks])
    print(freq[peaks]*2*np.pi)
    plt.plot(np.abs(sp))
    plt.plot(peaks, sp[peaks], "x")
    plt.show()
    return 2*np.pi*freq[peaks]




if __name__ == '__main__':
    lpath = str(sys.argv[1])
    getb(lpath)
    getv(lpath)

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


