# -*- coding: utf-8 -*-
"""
** Cleaned up script based on lock-in-emulator-2022-10-27-withquasiCompressibility.py **

Created on Thu Apr 21 17:30:15 2022

Testing lock-in amplifier code for use with the torque mixing simulation script output. 

The lock-in still requires input of parameters from the simulation: dt and f_RF (as fref)
as well as the cut-off frequency desired for the low-pass filter

06may2022  adding energy and angle lock-in analyses, for independent calculation of cross-product torques

20may2022 fixing energy/angle analysis, and adding anisotropy fields and torques

How to use: 
    - set parameters at top of script (i.e. time step, reference frequency, path information)
    - Call plotting and saving functions at bottom of script


@author: AWGPC1
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.close('all')


"""  lock-in parameters:  """
dt = 25e-3  #ns
f_ref = 0.208#GHz

fo = 4 #filter order for the low-pass filter. 

""" Constant Frequency Settings:"""
#cut-off frequency for the lock-in. I have found it works best if the cutoff
#frequency is about 0.4*reference frequency
fc_frac = 0.4
f_cutoff = f_ref*fc_frac 

""" Swept Frequency Settings: """
#parameters for frequency sweeps
f1 = 0.15   #GHz (starting frequency)
f2 = 0      #GHz (ending frequency)

fc_freq = 1e-2 #cutoff frequeny for a swept frequency. This typically needs to 
#be played with until the lock-in outputs appropriately, typically should be quite small

"""   Loading data:    """
path = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Conferences and Papers/NRC paper/Simulation Summary/Macrospin High Field z Phase offset/'
filename = 'PyShroom_zDrive_208MHz_2p5us_dt25ps_ph-60'
dat = np.loadtxt(path+filename+'.csv',delimiter=',',skiprows=11)

#define variables from the datafile
t = dat[:,0]  # time in ns

Hx = dat[:,1]
Hy = dat[:,2]
Hz = dat[:,3]
DC = True #boolean to indicate whether (true) or not (false) to include the DC component of the signal


fslope = (f2-f1)/t[-1] 
freq = fslope*t + f1 #linear frequency sweep based on simulation time and initial/final frequencies as input above


""" Functions to be used in the lock-in emulator """
def lowPass(sig,fc,dt,FO):
    n = len(sig) #length of array (to be used in FFT)
    ft_sig = np.fft.fft(sig) #FFT of signal
    f = np.fft.fftfreq(n,d=dt) #frequencies associated with signal FFT
    H = (1/(1 + 1j*2*np.pi*f/(2*np.pi*fc)))**FO  #Transfer function to be applied to FFT of signal
    ft_filt = ft_sig*H #apply transfer function to signal FFT
    lp = np.fft.ifft(ft_filt) # perform IFFT to return filtered signal to time-domain. 
    return lp
    

def lockIn(t,sig,dt,fref,fc,FO,DC): 
    """
    Lock-In function for constant frequency
    ----------
    Parameters
    ----------
    t : time array associated with signal
    sig : signal array to be processed
    dt : time step size of signal
    fref : reference frequency for the lock-in modulation
    fc : cut-off frequency for low-pass filter used within the lock-in
    FO : filter order for the low-pass filter used within the lock-in
    """
    
    if not DC: #remove DC component of signal 
        sigfft = np.fft.fft(sig)
        sigfft[0]=0
        sig = np.fft.ifft(sigfft).real
    ref = (-1j)*np.sqrt(2)*np.exp(1j*2*np.pi*fref*t) #reference signal in complex plane
    mix = sig*ref #mixed signal with reference
    filt = lowPass(mix,fc,dt,FO) #apply low-pass filter to mixed signal to remove the f1 + f2 component

    X = filt.real 
    Y = filt.imag
    R = np.sqrt(X**2 + Y**2)
    Th = np.unwrap(np.arctan2(Y,X))
    return X,Y,R,Th

def freqSweepLockIn(t,sig,dt,fc,FO,DC):
    """
    Lock-In function for swept frequency
    ----------
    Parameters
    ----------
    t : time array associated with signal
    sig : signal array to be processed
    dt : time step size of signal
    fref : reference frequency for the lock-in modulation
    fc : cut-off frequency for low-pass filter used within the lock-in
    FO : filter order for the low-pass filter used within the lock-in
    """
    
    if not DC: #remove DC component of signal 
        sigfft = np.fft.fft(sig)
        sigfft[0]=0
        sig = np.fft.ifft(sigfft).real
    
    fint= (0.5*fslope)*t**2 + f1*t #integral over time of a *linear* frequency function (to make a sinusoid with frequency that varies like f(t))
    ref = (-1j)*np.sqrt(2)*np.exp(1j*2*np.pi*fint)
    mix = sig*ref #mixed signal with reference
    fcint = fint*fc_frac
    fcint[0] = fcint[0]+1e-6
    filt = lowPass(mix,fc,dt,FO) #apply low-pass filter to mixed signal to remove the f1 + f2 component

    X = filt.real 
    Y = filt.imag
    R = np.sqrt(X**2 + Y**2)
    Th = np.unwrap(np.arctan2(Y,X))
    return X,Y,R,Th


""" PLOT DEFINITIONS """

def EdHTorqueTime(fx,fy,fz,fc): #EdH torque plotted as a function of time
    xtorque = dat[:,28]
    ytorque = dat[:,29]
    ztorque = dat[:,30]
    #t,sig,dt,fref,fc,FO
    Xx,Yx,Rx,Tx = lockIn(t,xtorque,dt,fx,fc,fo,DC)
    Xy,Yy,Ry,Ty = lockIn(t,ytorque,dt,fy,fc,fo,DC)
    Xz,Yz,Rz,Tz = lockIn(t,ztorque,dt,fz,fc,fo,DC)
    
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(t,Rx)
    ax[1,0].plot(t,np.rad2deg(Tx))
    
    ax[0,1].plot(t,Ry)
    ax[1,1].plot(t,np.rad2deg(Ty))
    
    ax[0,2].plot(t,Rz)
    ax[1,2].plot(t,np.rad2deg(Tz))
    
    ax[0,0].set_ylabel('Magnitude [aNm]')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude [aNm]')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude [aNm]')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('Time [ns]')
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$\tau_{EdH}$ (x) '+'- demod {:.0f} MHz'.format(fx*1000))
    ax[0,1].set_title(r'$\tau_{EdH}$ (y) '+'- demod {:.0f} MHz'.format(fy*1000))
    ax[0,2].set_title(r'$\tau_{EdH}$ (z) '+'- demod {:.0f} MHz'.format(fz*1000))
    fig.suptitle('EdH torques - lock-in emulator - demod {:.0f} MHz'.format(f_ref*1000))

def EdHTorqueField(fx,fy,fz,fc,H): #EdH torque plotted as a function of field
    xtorque = dat[:,28]
    ytorque = dat[:,29]
    ztorque = dat[:,30]
    Xx,Yx,Rx,Tx = lockIn(t,xtorque,dt,fx,fc,fo,DC)
    Xy,Yy,Ry,Ty = lockIn(t,ytorque,dt,fy,fc,fo,DC)
    Xz,Yz,Rz,Tz = lockIn(t,ztorque,dt,fz,fc,fo,DC)
    
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(H*1e-3,Rx)
    ax[1,0].plot(H*1e-3,np.rad2deg(Tx))
    
    ax[0,1].plot(H*1e-3,Ry)
    ax[1,1].plot(H*1e-3,np.rad2deg(Ty))
    
    ax[0,2].plot(H*1e-3,Rz)
    ax[1,2].plot(H*1e-3,np.rad2deg(Tz))
    
    ax[0,0].set_ylabel('Magnitude [aNm]')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude [aNm]')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude [aNm]')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('$H_x^{DC}$ [kA/m]')
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$\tau_{EdH}$ (x) '+'- demod {:.0f} MHz'.format(fx*1000))
    ax[0,1].set_title(r'$\tau_{EdH}$ (y) '+'- demod {:.0f} MHz'.format(fy*1000))
    ax[0,2].set_title(r'$\tau_{EdH}$ (z) '+'- demod {:.0f} MHz'.format(fz*1000))
    fig.suptitle('EdH torques - lock-in emulator - demod {:.0f} MHz'.format(f_ref*1000))

def EdHTorqueFreq(fc): #EdH torque plotted as a function of frequency
    xtorque = dat[:,28]
    ytorque = dat[:,29]
    ztorque = dat[:,30]
    #t,sig,dt,fref,fc,FO
    Xx,Yx,Rx,Tx = freqSweepLockIn(t,xtorque,dt,fc,fo,DC)
    Xy,Yy,Ry,Ty = freqSweepLockIn(t,ytorque,dt,fc,fo,DC)
    Xz,Yz,Rz,Tz = freqSweepLockIn(t,ztorque,dt,fc,fo,DC)
    
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(freq*1000,Rx)
    ax[1,0].plot(freq*1000,np.rad2deg(Tx))
    
    ax[0,1].plot(freq*1000,Ry)
    ax[1,1].plot(freq*1000,np.rad2deg(Ty))
    
    ax[0,2].plot(freq*1000,Rz)
    ax[1,2].plot(freq*1000,np.rad2deg(Tz))
    
    ax[0,0].set_ylabel('Magnitude [aNm]')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude [aNm]')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude [aNm]')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('Frequency [MHz]')
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$\tau_{EdH}$ (x) '+'- demod {:.0f} MHz'.format(f_ref*1000))
    ax[0,1].set_title(r'$\tau_{EdH}$ (y) '+'- demod {:.0f} MHz'.format(f_ref*1000))
    ax[0,2].set_title(r'$\tau_{EdH}$ (z) '+'- demod {:.0f} MHz'.format(f_ref*1000))
    fig.suptitle('EdH torques - lock-in emulator - demod {:.0f} MHz'.format(f_ref*1000))
    
def mxBTorqueTime(fx,fy,fz,fc): #cross-product torque plotted as a function of time
    xtorque = dat[:,19] # applied field
    ytorque = dat[:,20]
    ztorque = dat[:,21]
    
    Xx,Yx,Rx,Tx = lockIn(t,xtorque,dt,fx,fc,fo,DC)
    Xy,Yy,Ry,Ty = lockIn(t,ytorque,dt,fy,fc,fo,DC)
    Xz,Yz,Rz,Tz = lockIn(t,ztorque,dt,fz,fc,fo,DC)
    
    i0 = 0
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(t,xtorque)
    ax[0,0].plot(t[i0:],Rx[i0:])
    ax[1,0].plot(t[i0:],np.rad2deg(Tx[i0:]))
    
    ax[0,1].plot(t,ytorque)
    ax[0,1].plot(t[i0:],Ry[i0:])
    ax[1,1].plot(t[i0:],np.rad2deg(Ty[i0:]))
    
    ax[0,2].plot(t,ztorque)
    ax[0,2].plot(t[i0:],Rz[i0:])
    ax[1,2].plot(t[i0:],np.rad2deg(Tz[i0:]))
    
    ax[0,0].set_ylabel('Magnitude [aNm]')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude [aNm]')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude [aNm]')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('Time [ns]')
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$\tau_{mxB}$ (x) '+'- demod {:.0f} MHz'.format(fx*1000))
    ax[0,1].set_title(r'$\tau_{mxB}$ (y) '+'- demod {:.0f} MHz'.format(fy*1000))
    ax[0,2].set_title(r'$\tau_{mxB}$ (z) '+'- demod {:.0f} MHz'.format(fz*1000))
    fig.suptitle('Precession (specify) torques - lock-in emulator - demod {:.0f} MHz'.format(f_ref*1000))

def mxBTorqueField(fx,fy,fz,fc,H): #cross-product torque plotted as a function of field
    xtorque = dat[:,19] # applied field
    ytorque = dat[:,20]
    ztorque = dat[:,21]
    
    Xx,Yx,Rx,Tx = lockIn(t,xtorque,dt,fx,fc,fo,DC)
    Xy,Yy,Ry,Ty = lockIn(t,ytorque,dt,fy,fc,fo,DC)
    Xz,Yz,Rz,Tz = lockIn(t,ztorque,dt,fz,fc,fo,DC)
    
    i0 = 2000
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(H[i0:]*1e-3,Rx[i0:])
    ax[1,0].plot(H[i0:]*1e-3,np.rad2deg(Tx[i0:]))
    
    ax[0,1].plot(H[i0:]*1e-3,Ry[i0:])
    ax[1,1].plot(H[i0:]*1e-3,np.rad2deg(Ty[i0:]))
    
    ax[0,2].plot(H[i0:]*1e-3,Rz[i0:])
    ax[1,2].plot(H[i0:]*1e-3,np.rad2deg(Tz[i0:]))
    
    ax[0,0].set_ylabel('Magnitude [aNm]')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude [aNm]')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude [aNm]')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('$H_x^{DC}$ [kA/m]')
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$\tau_{mxB}$ (x) '+'- demod {:.0f} MHz'.format(fx*1000))
    ax[0,1].set_title(r'$\tau_{mxB}$ (y) '+'- demod {:.0f} MHz'.format(fy*1000))
    ax[0,2].set_title(r'$\tau_{mxB}$ (z) '+'- demod {:.0f} MHz'.format(fz*1000))
    fig.suptitle('Precession (specify) torques - lock-in emulator - demod {:.0f} MHz'.format(f_ref*1000))

def mxBTorqueFreq(fc): #cross-product torque plotted as a function of frequency
    xtorque = dat[:,22] # anisotropy field
    ytorque = dat[:,23]
    ztorque = dat[:,24]
    
    Xx,Yx,Rx,Tx = freqSweepLockIn(t,xtorque,dt,fc,fo,DC)
    Xy,Yy,Ry,Ty = freqSweepLockIn(t,ytorque,dt,fc,fo,DC)
    Xz,Yz,Rz,Tz = freqSweepLockIn(t,ztorque,dt,fc,fo,DC)
    
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(freq*1000,Rx)
    ax[1,0].plot(freq*1000,np.rad2deg(Tx))
    
    ax[0,1].plot(freq*1000,Ry)
    ax[1,1].plot(freq*1000,np.rad2deg(Ty))
    
    ax[0,2].plot(freq*1000,Rz)
    ax[1,2].plot(freq*1000,np.rad2deg(Tz))
    
    ax[0,0].set_ylabel('Magnitude [aNm]')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude [aNm]')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude [aNm]')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('Frequency [MHz]')   
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$\tau_{mxB}$ (x) '+'- demod {:.0f} MHz'.format(f_ref*1000))
    ax[0,1].set_title(r'$\tau_{mxB}$ (y) '+'- demod {:.0f} MHz'.format(f_ref*1000))
    ax[0,2].set_title(r'$\tau_{mxB}$ (z) '+'- demod {:.0f} MHz'.format(f_ref*1000))
    fig.suptitle('Precession (specify) torques - lock-in emulator - demod {:.0f} MHz'.format(f_ref*1000))

def resTorqueTime(fx,fy,fz,fc): #resultant-product torque plotted as a function of time
    xtorque = dat[:,19]-dat[:,28] # applied field
    ytorque = dat[:,20]-dat[:,29]
    ztorque = dat[:,21]-dat[:,30]

    
    Xx,Yx,Rx,Tx = lockIn(t,xtorque,dt,fx,fc,fo,DC)
    Xy,Yy,Ry,Ty = lockIn(t,ytorque,dt,fy,fc,fo,DC)
    Xz,Yz,Rz,Tz = lockIn(t,ztorque,dt,fz,fc,fo,DC)
    
    i0 = 0
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    # ax[0,0].plot(t,xtorque)
    ax[0,0].plot(t[i0:],Rx[i0:])
    ax[1,0].plot(t[i0:],np.rad2deg(Tx[i0:]))
    
    # ax[0,1].plot(t,ytorque)
    ax[0,1].plot(t[i0:],Ry[i0:])
    ax[1,1].plot(t[i0:],np.rad2deg(Ty[i0:]))
    
    # ax[0,2].plot(t,ztorque)
    ax[0,2].plot(t[i0:],Rz[i0:])
    ax[1,2].plot(t[i0:],np.rad2deg(Tz[i0:]))
    
    ax[0,0].set_ylabel('Magnitude [aNm]')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude [aNm]')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude [aNm]')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('Time [ns]')
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$\tau_{res}$ (x) '+'- demod {:.0f} MHz'.format(fx*1000))
    ax[0,1].set_title(r'$\tau_{res}$ (y) '+'- demod {:.0f} MHz'.format(fy*1000))
    ax[0,2].set_title(r'$\tau_{res}$ (z) '+'- demod {:.0f} MHz'.format(fz*1000))
    fig.suptitle('Resultant torques - lock-in emulator - demod {:.0f} MHz'.format(f_ref*1000))
    
    plt.savefig(path+'lockin-output-ResTorques-'+filename+'.png',format='png',dpi=300,bbox_inches='tight')


def resTorqueField(fx,fy,fz,fc,H):#resultant-product torque plotted as a function of field

    xtorque = dat[:,19] - dat[:,28] # anisotropy field
    ytorque = dat[:,20] - dat[:,29]
    ztorque = dat[:,21] - dat[:,30]
    
    Xx,Yx,Rx,Tx = lockIn(t,xtorque,dt,fx,fc,fo,DC)
    Xy,Yy,Ry,Ty = lockIn(t,ytorque,dt,fy,fc,fo,DC)
    Xz,Yz,Rz,Tz = lockIn(t,ztorque,dt,fz,fc,fo,DC)
    
    i0 = 2000
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(H[i0:]*1e-3,Rx[i0:])
    ax[1,0].plot(H[i0:]*1e-3,np.rad2deg(Tx[i0:]))
    
    ax[0,1].plot(H[i0:]*1e-3,Ry[i0:])
    ax[1,1].plot(H[i0:]*1e-3,np.rad2deg(Ty[i0:]))
    
    ax[0,2].plot(H[i0:]*1e-3,Rz[i0:])
    ax[1,2].plot(H[i0:]*1e-3,np.rad2deg(Tz[i0:]))
    
    ax[0,0].set_ylabel('Magnitude [aNm]')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude [aNm]')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude [aNm]')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('$H_x^{DC}$ [kA/m]')
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$\tau_{res}$ (x) '+'- demod {:.0f} MHz'.format(fx*1000))
    ax[0,1].set_title(r'$\tau_{res}$ (y) '+'- demod {:.0f} MHz'.format(fy*1000))
    ax[0,2].set_title(r'$\tau_{res}$ (z) '+'- demod {:.0f} MHz'.format(fz*1000))
    fig.suptitle('Resultant (specify) torques - lock-in emulator - demod {:.0f} MHz'.format(f_ref*1000))

def resTorqueFreq(fc): #resultant-product torque plotted as a function of frequency
    xtorque = dat[:,19] - dat[:,28]
    ytorque = dat[:,20] - dat[:,29]
    ztorque = dat[:,21] - dat[:,30]
    
    Xx,Yx,Rx,Tx = freqSweepLockIn(t,xtorque,dt,fc,fo,DC)
    Xy,Yy,Ry,Ty = freqSweepLockIn(t,ytorque,dt,fc,fo,DC)
    Xz,Yz,Rz,Tz = freqSweepLockIn(t,ztorque,dt,fc,fo,DC)
    
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(freq*1000,Rx)
    ax[1,0].plot(freq*1000,np.rad2deg(Tx))
    
    ax[0,1].plot(freq*1000,Ry)
    ax[1,1].plot(freq*1000,np.rad2deg(Ty))
    
    ax[0,2].plot(freq*1000,Rz)
    ax[1,2].plot(freq*1000,np.rad2deg(Tz))
    
    ax[0,0].set_ylabel('Magnitude [aNm]')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude [aNm]')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude [aNm]')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('Frequency [MHz]')   
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$\tau_{res}$ (x) '+'- demod {:.0f} MHz'.format(f_ref*1000))
    ax[0,1].set_title(r'$\tau_{res}$ (y) '+'- demod {:.0f} MHz'.format(f_ref*1000))
    ax[0,2].set_title(r'$\tau_{res}$ (z) '+'- demod {:.0f} MHz'.format(f_ref*1000))
    fig.suptitle('Resultant torques - lock-in emulator - demod {:.0f}-{:.1f} MHz'.format(f1*1000,f2*1000))
    
    plt.savefig(path+'lockin-output-ResTorques-'+filename+'.png',format='png',dpi=300,bbox_inches='tight')

def EnergyTime(fx,fy,fz): #plot energy as a function of time
    Edemag = dat[:,31]
    EZeem = dat[:,32]
    Euniax = dat[:,33]
    #Ecubic = dat[:,34]
    
    Xx,Yx,Rx,Tx = lockIn(t,Edemag,dt,fx,0.45*fx,fo,DC)
    Xy,Yy,Ry,Ty = lockIn(t,EZeem,dt,fy,0.45*fy,fo,DC)
    Xz,Yz,Rz,Tz = lockIn(t,Euniax,dt,fz,0.45*fz,fo,DC)
    
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(t,Rx)
    ax[1,0].plot(t,np.rad2deg(Tx))
    
    ax[0,1].plot(t,Ry)
    ax[1,1].plot(t,np.rad2deg(Ty))
    
    ax[0,2].plot(t,Rz)
    ax[1,2].plot(t,np.rad2deg(Tz))
    
    ax[0,0].set_ylabel('Magnitude ')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude ')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude ')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('Time [ns]')
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$E_{\mathrm{demag}}$ '+'- demod {:.0f} MHz'.format(fx*1000))
    ax[0,1].set_title(r'$E_{\mathrm{Zeeman}}$ '+'- demod {:.0f} MHz'.format(fy*1000))
    ax[0,2].set_title(r'$E_{\mathrm{uniax}}$ '+'- demod {:.0f} MHz'.format(fz*1000))
    fig.suptitle('Energies - lock-in emulator - demod {:.0f} MHz'.format(f_ref*1000))

def AngleTime(fx,fy,fz): #plot magnetization as a function of time
    mx = dat[:,13]
    my = dat[:,14]
    mz = dat[:,15]
    
    Xx,Yx,Rx,Tx = lockIn(t,mx,dt,fx,0.45*fx,fo,DC)
    Xy,Yy,Ry,Ty = lockIn(t,my,dt,fy,0.45*fy,fo,DC)
    Xz,Yz,Rz,Tz = lockIn(t,mz,dt,fz,0.45*fz,fo,DC)
    
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(t,Rx)
    ax[1,0].plot(t,np.rad2deg(Tx))
    
    ax[0,1].plot(t,Ry)
    ax[1,1].plot(t,np.rad2deg(Ty))
    
    ax[0,2].plot(t,Rz)
    ax[1,2].plot(t,np.rad2deg(Tz))
    
    ax[0,0].set_ylabel('Magnitude ')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude ')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('Magnitude ')
    ax[1,2].set_ylabel('Phase (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('Time [ns]')
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$m_x$ '+'- demod {:.0f} MHz'.format(fx*1000))
    ax[0,1].set_title(r'$m_y$ '+'- demod {:.0f} MHz'.format(fy*1000))
    ax[0,2].set_title(r'$m_z$ '+'- demod {:.0f} MHz'.format(fz*1000))
    fig.suptitle('Angles - lock-in emulator - demod {:.0f} MHz'.format(f_ref*1000))

def tauz(fy):
    my = dat[:,14] # for angle vs. time, dither along y
    mz = dat[:,15] # for angle vs. time, dither along z
    #Edemag = dat[:,31] # for work done by torque vs. shape anisotropy
    Euniax = dat[:,33] # for work done by torque vs. uniax anisotropy
    
    X1,Y1,R1,T1 = lockIn(t,my,dt,fy,0.45*fy,fo,DC)
    X2,Y2,R2,T2 = lockIn(t,Euniax,dt,2*fy,0.45*fy,fo,DC)
        
    fig, ax = plt.subplots(2,3,figsize=(18,6))
    ax[0,0].plot(t,R1)
    ax[1,0].plot(t,np.rad2deg(T1))
    
    ax[0,1].plot(t,R2)
    ax[1,1].plot(t,np.rad2deg(T2))
    
    rmstorque = 2*np.sqrt(2)*np.divide(R2,R1)
    phasecheck = np.subtract(np.rad2deg(T2),np.rad2deg(T1))
    
    ax[0,2].plot(t,rmstorque)
    ax[1,2].plot(t,phasecheck)
    
    ax[0,0].set_ylabel('Magnitude ')
    ax[1,0].set_ylabel('Phase (deg)')
    ax[0,1].set_ylabel('Magnitude ')
    ax[1,1].set_ylabel('Phase (deg)')
    ax[0,2].set_ylabel('RMS Torque ')
    ax[1,2].set_ylabel('Phase diff (deg)')
    
    for i in range(3):
        ax[0,i].set_xticklabels([])
        ax[1,i].set_xlabel('Time [ns]')
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.25)
    fig.align_ylabels()
    
    ax[0,0].set_title(r'$m_y$ '+'- demod {:.0f} MHz'.format(fy*1000))
    ax[0,1].set_title(r'$E_{\mathrm{uniax}}$ '+'- demod {:.0f} MHz'.format(2*fy*1000))
    ax[0,2].set_title(r'$2\sqrt{2}E_{\mathrm{uniax}}/m_y$ ')
    fig.suptitle('Torque from energy and angles - lock-in emulator')


def saveLIoutput(f_ref,fc,comp):
    """
    Lock-in emulator applied to all torques, magnetization, energies that are
    calculated in the macroscale torque mixing simulation script, so the output 
    csv file has the same shape. 
    
    The comp variable indicates the spacing between subsequent saved points 
    (i.e. if comp = n, every nth data point will be saved to the csv file.)
    This allows much smaller output files to be saved, as lockin output is 
    typically a smooth function and doesn't need as many points as the "raw" 
    simulation output. 
    """
    xtorque_damp = dat[:,16]
    ytorque_damp = dat[:,17]
    ztorque_damp = dat[:,18]
    xtorque_mxB_app = dat[:,19] # applied field
    ytorque_mxB_app = dat[:,20]
    ztorque_mxB_app = dat[:,21]
    xtorque_mxB_anis = dat[:,22] # anisotropy field
    ytorque_mxB_anis = dat[:,23]
    ztorque_mxB_anis = dat[:,24]
    xtorque_mxB_eff = dat[:,25] # effective field
    ytorque_mxB_eff = dat[:,26]
    ztorque_mxB_eff = dat[:,27]
    xtorque_EdH = dat[:,28]
    ytorque_EdH = dat[:,29]
    ztorque_EdH = dat[:,30]
    xtorque_res = dat[:,19] - dat[:,28]
    ytorque_res = dat[:,20] - dat[:,29]
    ztorque_res = dat[:,21] - dat[:,30]
    mx = dat[:,13]
    my = dat[:,14]
    mz = dat[:,15]
    Edemag = dat[:,31]
    EZeem = dat[:,32]
    Euniax = dat[:,33]
    Ecub = dat[:,34]
    
    Xdampz,Ydampx,Rdampx,Tdampx = lockIn(t,xtorque_damp,dt,f_ref,fc,fo,DC)
    Xdampy,Ydampy,Rdampy,Tdampy = lockIn(t,ytorque_damp,dt,f_ref,fc,fo,DC)
    Xdampz,Ydampz,Rdampz,Tdampz = lockIn(t,ztorque_damp,dt,f_ref,fc,fo,DC)
    
    XmxBx_eff,YmxBx_eff,RmxBx_eff,TmxBx_eff = lockIn(t,xtorque_mxB_eff,dt,f_ref,fc,fo,DC)
    XmxBy_eff,YmxBy_eff,RmxBy_eff,TmxBy_eff = lockIn(t,ytorque_mxB_eff,dt,f_ref,fc,fo,DC)
    XmxBz_eff,YmxBz_eff,RmxBz_eff,TmxBz_eff = lockIn(t,ztorque_mxB_eff,dt,f_ref,fc,fo,DC)
    
    XmxBx_app,YmxBx_app,RmxBx_app,TmxBx_app = lockIn(t,xtorque_mxB_app,dt,f_ref,fc,fo,DC)
    XmxBy_app,YmxBy_app,RmxBy_app,TmxBy_app = lockIn(t,ytorque_mxB_app,dt,f_ref,fc,fo,DC)
    XmxBz_app,YmxBz_app,RmxBz_app,TmxBz_app = lockIn(t,ztorque_mxB_app,dt,f_ref,fc,fo,DC)
    
    XmxBx_anis,YmxBx_anis,RmxBx_anis,TmxBx_anis = lockIn(t,xtorque_mxB_anis,dt,f_ref,fc,fo,DC)
    XmxBy_anis,YmxBy_anis,RmxBy_anis,TmxBy_anis = lockIn(t,ytorque_mxB_anis,dt,f_ref,fc,fo,DC)
    XmxBz_anis,YmxBz_anis,RmxBz_anis,TmxBz_anis = lockIn(t,ztorque_mxB_anis,dt,f_ref,fc,fo,DC)
    
    XEdHx,YEdHx,REdHx,TEdHx = lockIn(t,xtorque_EdH,dt,f_ref,fc,fo,DC)
    XEdHy,YEdHy,REdHy,TEdHy = lockIn(t,ytorque_EdH,dt,f_ref,fc,fo,DC)
    XEdHz,YEdHz,REdHz,TEdHz = lockIn(t,ztorque_EdH,dt,f_ref,fc,fo,DC)
    
    Xresx,Yresx,Rresx,Tresx = lockIn(t,xtorque_res,dt,f_ref,fc,fo,DC)
    Xresy,Yresy,Rresy,Tresy = lockIn(t,ytorque_res,dt,f_ref,fc,fo,DC)
    Xresz,Yresz,Rresz,Tresz = lockIn(t,ztorque_res,dt,f_ref,fc,fo,DC)
    
    Xmx,Ymx,Rmx,Tmx = lockIn(t,mx,dt,f_ref,fc,fo,DC)
    Xmy,Ymy,Rmy,Tmy = lockIn(t,my,dt,f_ref,fc,fo,DC)
    Xmz,Ymz,Rmz,Tmz = lockIn(t,mz,dt,f_ref,fc,fo,DC)
    
    XEd,YEd,REd,TEd = lockIn(t,Edemag,dt,f_ref,fc,fo,DC)
    XEZ,YEZ,REZ,TEZ = lockIn(t,EZeem,dt,f_ref,fc,fo,DC)
    XEu,YEu,REu,TEu = lockIn(t,Euniax,dt,f_ref,fc,fo,DC)
    XEc,YEc,REc,TEc = lockIn(t,Ecub,dt,f_ref,fc,fo,DC)
    
    c = comp
    
    df = pd.DataFrame({
        'time (ns)':t[::c],
        'H applied (x) [A/m]':Hx[::c], #index 1
        'H applied (y) [A/m]':Hy[::c], #index 2
        'H applied (z) [A/m]':Hz[::c],
        'H anisotropy (x) [A/m]':dat[:,4][::c], #index 4
        'H anisotropy (y) [A/m]':dat[:,5][::c], #index 5
        'H anisotropy (z) [A/m]':dat[:,6][::c], #index 6
        'H effective (x) [A/m]':dat[:,7][::c], #index 7
        'H effective (y) [A/m]':dat[:,8][::c], #index 8
        'H effective (z) [A/m]':dat[:,9][::c], #index 9
        'RF Frequency (x) [GHz]':dat[:,10][::c], #index 10
        'RF Frequency (y) [GHz]':dat[:,11][::c], #index 11
        'RF Frequency (z) [GHz]':dat[:,12][::c], #index 12  
        'mx':Rmx[::c], #index 13
        'my':Rmy[::c], #index 14
        'mz':Rmz[::c], #index 15
        'tau damping (x) [aN m]':Rdampx[::c], #index 16
        'tau damping (y) [aN m]':Rdampy[::c], #index 17
        'tau damping (z) [aN m]':Rdampz[::c], #index 18
        'tau precession (H_app) (x) [aN m]':RmxBx_app[::c], #index 19
        'tau precession (H_app) (y) [aN m]':RmxBy_app[::c], #index 20
        'tau precession (H_app) (z) [aN m]':RmxBz_app[::c], #index 21
        'tau precession (H_anis) (x) [aN m]':RmxBx_anis[::c], #index 22
        'tau precession (H_anis) (y) [aN m]':RmxBy_anis[::c], #index 23
        'tau precession (H_anis) (z) [aN m]':RmxBz_anis[::c], #index 24
        'tau precession (H_eff) (x) [aN m]':RmxBx_eff[::c], #index 25
        'tau precession (H_eff) (y) [aN m]':RmxBy_eff[::c], #index 26
        'tau precession (H_eff) (z) [aN m]':RmxBz_eff[::c], #index 27
        'tau EdH (x) [aN m]':REdHx[::c], #index 28
        'tau EdH (y) [aN m]':REdHy[::c], #index 29
        'tau EdH (z) [aN m]':REdHz[::c], #index 30
        'tau res (x) [aN m]':Rresx[::c], #index 31
        'tau res (y) [aN m]':Rresy[::c], #index 32
        'tau res (z) [aN m]':Rresz[::c], #index 33
        'E_demag [aJ]':REd[::c], #index 34
        'E_Zeeman [aJ]':REZ[::c], #index 35
        'E_anisotropy_uniaxial [aJ]':dat[:,33][::c], #index 36
        'E_anisotropy_cubic [aJ]':REc[::c], #index 37
        'tau precession phase (H_app) (x) [rad]':TmxBx_app[::c],#index 38
        'tau precession phase (H_app) (y) [rad]':TmxBy_app[::c],#index 39
        'tau precession phase (H_app) (z) [rad]':TmxBz_app[::c],#index 40
        'tau precession phase (H_anis) (x) [rad]':TmxBx_anis[::c],#index 41
        'tau precession phase (H_anis) (y) [rad]':TmxBy_anis[::c],#index 42
        'tau precession phase (H_anis) (z) [rad]':TmxBz_eff[::c],#index 43
        'tau precession phase (H_eff) (x) [rad]':TmxBx_eff[::c],#index 44
        'tau precession phase (H_eff) (y) [rad]':TmxBy_eff[::c],#index 45
        'tau precession phase (H_eff) (z) [rad]':TmxBz_eff[::c],#index 46
        'tau EdH phase (x) [rad]':TEdHx[::c],#index 47
        'tau EdH phase (y) [rad]':TEdHy[::c],#index 48
        'tau EdH phase (z) [rad]':TEdHz[::c],#index 49
        'tau res phase (x) [rad':Tresx[::c], #index 50
        'tau res phase (y) [rad]':Tresy[::c], #index 51
        'tau res phase (z) [rad]':Tresz[::c], #index 52
        
        })
    # plt.figure()
    # plt.plot(Hx,RmxBy_anis)
    if comp != 1:
        file_path = path +'lockin-output-compressed{}pt-'.format(comp) + filename+'.csv'
    else:
        file_path = path+'lockin-output-'+filename+'.csv'
    
    if DC:
        dctxt = ' Yes'
    else:
        dctxt = ' No'
    t1 = 'reference frequency = {0} GHz; cut-off frequency = {1} GHz\n '.format(f_ref,f_cutoff)
    t2 = 'dt = {0} ns; filter order = {1}\n '.format(dt,fo)
    t3 = 'DC component included?' + dctxt + '\n'
    parametersString = ' """ Lock-in Emulator Parameters:\n' + t1 + t2 + t3
        
    #write file header with all the simulation parameter info first
    with open(file_path, 'w') as f: 
          f.write(parametersString)
    #Then write dataFrame to the same file. 
    df.to_csv(file_path, header=True, index=False,mode="a")

def saveFreqSweepLIoutput(fc,comp):
    """
    Lock-in emulator applied to all torques, magnetization, energies that are
    calculated in the macroscale torque mixing simulation script, so the output 
    csv file has the same shape. 
    
    The comp variable indicates the spacing between subsequent saved points 
    (i.e. if comp = n, every nth data point will be saved to the csv file.)
    This allows much smaller output files to be saved, as lockin output is 
    typically a smooth function and doesn't need as many points as the "raw" 
    simulation output. 
    """
    xtorque_damp = dat[:,16]
    ytorque_damp = dat[:,17]
    ztorque_damp = dat[:,18]
    xtorque_mxB_app = dat[:,19] # applied field
    ytorque_mxB_app = dat[:,20]
    ztorque_mxB_app = dat[:,21]
    xtorque_mxB_anis = dat[:,22] # anisotropy field
    ytorque_mxB_anis = dat[:,23]
    ztorque_mxB_anis = dat[:,24]
    xtorque_mxB_eff = dat[:,25] # effective field
    ytorque_mxB_eff = dat[:,26]
    ztorque_mxB_eff = dat[:,27]
    xtorque_EdH = dat[:,28]
    ytorque_EdH = dat[:,29]
    ztorque_EdH = dat[:,30]
    xtorque_res = dat[:,19] - dat[:,28]
    ytorque_res = dat[:,20] - dat[:,29]
    ztorque_res = dat[:,21] - dat[:,30]
    
    mx = dat[:,13]
    my = dat[:,14]
    mz = dat[:,15]
    Edemag = dat[:,31]
    EZeem = dat[:,32]
    Euniax = dat[:,33]
    Ecub = dat[:,34]

    Xdampz,Ydampx,Rdampx,Tdampx = freqSweepLockIn(t,xtorque_damp,dt,fc,fo,DC)
    Xdampy,Ydampy,Rdampy,Tdampy = freqSweepLockIn(t,ytorque_damp,dt,fc,fo,DC)
    Xdampz,Ydampz,Rdampz,Tdampz = freqSweepLockIn(t,ztorque_damp,dt,fc,fo,DC)
    
    XmxBz_eff,YmxBx_eff,RmxBx_eff,TmxBx_eff = freqSweepLockIn(t,xtorque_mxB_eff,dt,fc,fo,DC)
    XmxBy_eff,YmxBy_eff,RmxBy_eff,TmxBy_eff = freqSweepLockIn(t,ytorque_mxB_eff,dt,fc,fo,DC)
    XmxBz_eff,YmxBz_eff,RmxBz_eff,TmxBz_eff = freqSweepLockIn(t,ztorque_mxB_eff,dt,fc,fo,DC)
    
    XmxBz_app,YmxBx_app,RmxBx_app,TmxBx_app = freqSweepLockIn(t,xtorque_mxB_app,dt,fc,fo,DC)
    XmxBy_app,YmxBy_app,RmxBy_app,TmxBy_app = freqSweepLockIn(t,ytorque_mxB_app,dt,fc,fo,DC)
    XmxBz_app,YmxBz_app,RmxBz_app,TmxBz_app = freqSweepLockIn(t,ztorque_mxB_app,dt,fc,fo,DC)
    
    XmxBz_anis,YmxBx_anis,RmxBx_anis,TmxBx_anis = freqSweepLockIn(t,xtorque_mxB_anis,dt,fc,fo,DC)
    XmxBy_anis,YmxBy_anis,RmxBy_anis,TmxBy_anis = freqSweepLockIn(t,ytorque_mxB_anis,dt,fc,fo,DC)
    XmxBz_anis,YmxBz_anis,RmxBz_anis,TmxBz_anis = freqSweepLockIn(t,ztorque_mxB_anis,dt,fc,fo,DC)
    
    XEdHx,YEdHx,REdHx,TEdHx = freqSweepLockIn(t,xtorque_EdH,dt,fc,fo,DC)
    XEdHy,YEdHy,REdHy,TEdHy = freqSweepLockIn(t,ytorque_EdH,dt,fc,fo,DC)
    XEdHz,YEdHz,REdHz,TEdHz = freqSweepLockIn(t,ztorque_EdH,dt,fc,fo,DC)
    
    Xresx,Yresx,Rresx,Tresx = freqSweepLockIn(t,xtorque_res,dt,fc,fo,DC)
    Xresy,Yresy,Rresy,Tresy = freqSweepLockIn(t,ytorque_res,dt,fc,fo,DC)
    Xresz,Yresz,Rresz,Tresz = freqSweepLockIn(t,ztorque_res,dt,fc,fo,DC)
    
    Xmx,Ymx,Rmx,Tmx = freqSweepLockIn(t,mx,dt,fc,fo,DC)
    Xmy,Ymy,Rmy,Tmy = freqSweepLockIn(t,my,dt,fc,fo,DC)
    Xmz,Ymz,Rmz,Tmz = freqSweepLockIn(t,mz,dt,fc,fo,DC)
    
    XEd,YEd,REd,TEd = freqSweepLockIn(t,Edemag,dt,fc,fo,DC)
    XEZ,YEZ,REZ,TEZ = freqSweepLockIn(t,EZeem,dt,fc,fo,DC)
    XEu,YEu,REu,TEu = freqSweepLockIn(t,Euniax,dt,fc,fo,DC)
    XEc,YEc,REc,TEc = freqSweepLockIn(t,Ecub,dt,fc,fo,DC)
      
    c = comp
    
    df = pd.DataFrame({
        'time (ns)':t[::c],
        'H applied (x) [A/m]':Hx[::c], #index 1
        'H applied (y) [A/m]':Hy[::c], #index 2
        'H applied (z) [A/m]':Hz[::c],
        'H anisotropy (x) [A/m]':dat[:,4][::c], #index 4
        'H anisotropy (y) [A/m]':dat[:,5][::c], #index 5
        'H anisotropy (z) [A/m]':dat[:,6][::c], #index 6
        'H effective (x) [A/m]':dat[:,7][::c], #index 7
        'H effective (y) [A/m]':dat[:,8][::c], #index 8
        'H effective (z) [A/m]':dat[:,9][::c], #index 9
        'RF Frequency (x) [GHz]':dat[:,10][::c], #index 10
        'RF Frequency (y) [GHz]':dat[:,11][::c], #index 11
        'RF Frequency (z) [GHz]':dat[:,12][::c], #index 12  
        'mx':Rmx[::c], #index 13
        'my':Rmy[::c], #index 14
        'mz':Rmz[::c], #index 15
        'tau damping (x) [aN m]':Rdampx[::c], #index 16
        'tau damping (y) [aN m]':Rdampy[::c], #index 17
        'tau damping (z) [aN m]':Rdampz[::c], #index 18
        'tau precession (H_app) (x) [aN m]':RmxBx_app[::c], #index 19
        'tau precession (H_app) (y) [aN m]':RmxBy_app[::c], #index 20
        'tau precession (H_app) (z) [aN m]':RmxBz_app[::c], #index 21
        'tau precession (H_anis) (x) [aN m]':RmxBx_anis[::c], #index 22
        'tau precession (H_anis) (y) [aN m]':RmxBy_anis[::c], #index 23
        'tau precession (H_anis) (z) [aN m]':RmxBz_anis[::c], #index 24
        'tau precession (H_eff) (x) [aN m]':RmxBx_eff[::c], #index 25
        'tau precession (H_eff) (y) [aN m]':RmxBy_eff[::c], #index 26
        'tau precession (H_eff) (z) [aN m]':RmxBz_eff[::c], #index 27
        'tau EdH (x) [aN m]':REdHx[::c], #index 28
        'tau EdH (y) [aN m]':REdHy[::c], #index 29
        'tau EdH (z) [aN m]':REdHz[::c], #index 30
        'tau res (x) [aN m]':Rresx[::c], #index 31
        'tau res (y) [aN m]':Rresy[::c], #index 32
        'tau res (z) [aN m]':Rresz[::c], #index 33
        'E_demag [aJ]':REd[::c], #index 34
        'E_Zeeman [aJ]':REZ[::c], #index 35
        'E_anisotropy_uniaxial [aJ]':dat[:,33][::c], #index 36
        'E_anisotropy_cubic [aJ]':REc[::c], #index 37
        'tau precession phase (H_app) (x) [rad]':TmxBx_app[::c],#index 38
        'tau precession phase (H_app) (y) [rad]':TmxBy_app[::c],#index 39
        'tau precession phase (H_app) (z) [rad]':TmxBz_app[::c],#index 40
        'tau precession phase (H_anis) (x) [rad]':TmxBx_anis[::c],#index 41
        'tau precession phase (H_anis) (y) [rad]':TmxBy_anis[::c],#index 42
        'tau precession phase (H_anis) (z) [rad]':TmxBz_eff[::c],#index 43
        'tau precession phase (H_eff) (x) [rad]':TmxBx_eff[::c],#index 44
        'tau precession phase (H_eff) (y) [rad]':TmxBy_eff[::c],#index 45
        'tau precession phase (H_eff) (z) [rad]':TmxBz_eff[::c],#index 46
        'tau EdH phase (x) [rad]':TEdHx[::c],#index 47
        'tau EdH phase (y) [rad]':TEdHy[::c],#index 48
        'tau EdH phase (z) [rad]':TEdHz[::c],#index 49
        'tau res phase (x) [rad':Tresx[::c], #index 50
        'tau res phase (y) [rad]':Tresy[::c], #index 51
        'tau res phase (z) [rad]':Tresz[::c], #index 52
        })
    
    if comp != 1:
        file_path = path +'lockin-output-compressed{}pt-'.format(comp) + filename+'.csv'
    else:
        file_path = path+'lockin-output-'+filename+'.csv'
    
    if DC:
        dctxt = ' Yes'
    else:
        dctxt = ' No'
    t1 = 'reference frequency = {0} GHz; cut-off frequency = {1} GHz\n '.format(f_ref,f_cutoff)
    t2 = 'dt = {0} ns; filter order = {1}\n '.format(dt,fo)
    t3 = 'DC component included?' + dctxt + '\n'
    parametersString = ' """ Lock-in Emulator Parameters:\n' + t1 + t2 + t3
        
    #write file header with all the simulation parameter info first
    with open(file_path, 'w') as f: 
         f.write(parametersString)
    #Then write dataFrame to the same file. 
    df.to_csv(file_path, header=True, index=False,mode="a")



""" Call any plotting or saving functions that are desired """

# mxBTorqueField(f_ref,f_ref,f_ref,f_cutoff,Hx)
# EdHTorqueField(f_ref,f_ref,f_ref,f_cutoff,Hx)

# mxBTorqueFreq(fc_freq)
# EdHTorqueFreq(fc_freq)

resTorqueTime(f_ref,f_ref,f_ref,f_cutoff)
saveLIoutput(f_ref,f_cutoff,comp=50)

# resTorqueFreq(fc_freq)
# saveFreqSweepLIoutput(fc_freq,comp=500)

# resTorqueField(f_ref,f_ref,f_ref,f_cutoff,k)


