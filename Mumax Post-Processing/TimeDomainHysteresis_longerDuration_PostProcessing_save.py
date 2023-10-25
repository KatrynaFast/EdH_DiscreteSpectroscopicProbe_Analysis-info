# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 11:52:35 2023

@author: titanw_u
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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


Ms = 767e3
V = 5.42e-20
e = 1.602e-19
m_e = 9.11e-31
g = 2

AEdH = 1e18*2*m_e*Ms*V/(e*g)
AmxB = 1e18*Ms*V*1e-4 #1e-4 to account for conversion of G to T

path = '' #folder path for mumax simulation output 
tab = pd.read_csv(path+'table.txt',delimiter='\t')
savename = 'Processed-Torques_HighFields_218MHz_yDrive.csv'

ti = 9175 #number of times saved for each DC field value 
tiA=5000 #index at which to begin averaging (tiA:)

## Parameters for lock-in filtering
f0 = 218e6
fc = 0.3*f0
dt = 10e-12
FO = 4

m = tab[tab.columns[1:4]].to_numpy().T
B = tab[tab.columns[4:7]].to_numpy().T*1e4
t = tab[tab.columns[0]].to_numpy()

#calculate torques from sim output
mxb = AmxB*np.cross(m.T,B.T).T
edh = AEdH*np.gradient(m,dt,axis=1)
res = mxb + edh

""" Lock-in filtering for m, and all torques"""

Xmx,Ymx,Rmx,Tmx = lockIn(t,m[0],dt,f0,fc,FO,1)
Xmy,Ymy,Rmy,Tmy = lockIn(t,m[1],dt,f0,fc,FO,1)
Xmz,Ymz,Rmz,Tmz = lockIn(t,m[2],dt,f0,fc,FO,1)

XXx,YXx,RXx,TXx = lockIn(t,mxb[0],dt,f0,fc,FO,1)
XXy,YXy,RXy,TXy = lockIn(t,mxb[1],dt,f0,fc,FO,1)
XXz,YXz,RXz,TXz = lockIn(t,mxb[2],dt,f0,fc,FO,1)

XEx,YEx,REx,TEx = lockIn(t,edh[0],dt,f0,fc,FO,1)
XEy,YEy,REy,TEy = lockIn(t,edh[1],dt,f0,fc,FO,1)
XEz,YEz,REz,TEz = lockIn(t,edh[2],dt,f0,fc,FO,1)

XRx,YRx,RRx,TRx = lockIn(t,res[0],dt,f0,fc,FO,1)
XRy,YRy,RRy,TRy = lockIn(t,res[1],dt,f0,fc,FO,1)
XRz,YRz,RRz,TRz = lockIn(t,res[2],dt,f0,fc,FO,1)



Rm = np.array([Rmx,Rmy,Rmz])
RX = np.array([RXx,RXy,RXz])
RE = np.array([REx,REy,REz])
RR = np.array([RRx,RRy,RRz])

Tm = np.array([Tmx,Tmy,Tmz])
TX = np.array([TXx,TXy,TXz])
TE = np.array([TEx,TEy,TEz])
TR = np.array([TRx,TRy,TRz])


NB = np.unique(B[0]).size
BxDC = np.unique(B[0])

db = np.where(B[0]==BxDC[1])[0][0] - np.where(B[0]==BxDC[0])[0][0] #

#if the field sweep is decreasing in magnitude, BxDC needs to be reversed
if db<0:
    BxDC = BxDC[::-1]

Ravm = []
RavX = []
RavE = []
RavR = []
Tavm = []
TavX = []
TavE = []
TavR = []

for i in range(NB):
    #break up locked-in arrays into smaller arrays of each individual DC field
    RTm = Rm[:,i*ti:(i+1)*ti-1]
    RTX = RX[:,i*ti:(i+1)*ti-1]
    RTE = RE[:,i*ti:(i+1)*ti-1]
    RTR = RR[:,i*ti:(i+1)*ti-1]
    
    #take the "saturated" point (where the magnetization is settled)
    Rmsat = RTm[:,tiA:]
    RXsat = RTX[:,tiA:]
    REsat = RTE[:,tiA:]
    RRsat = RTR[:,tiA:]
    
    #grab the mean of the settled output
    Ravm.append(Rmsat.mean(axis=1))
    RavX.append(RXsat.mean(axis=1))
    RavE.append(REsat.mean(axis=1))
    RavR.append(RRsat.mean(axis=1))
    
    TTm = Tm[:,i*ti:(i+1)*ti-1]
    TTX = TX[:,i*ti:(i+1)*ti-1]
    TTE = TE[:,i*ti:(i+1)*ti-1]
    TTR = TR[:,i*ti:(i+1)*ti-1]
    
    #take the "saturated" point (where the magnetization is settled)
    Tmsat = TTm[:,tiA:]
    TXsat = TTX[:,tiA:]
    TEsat = TTE[:,tiA:]
    TRsat = TTR[:,tiA:]
    
    #grab the mean of the settled output
    Tavm.append(Tmsat.mean(axis=1))
    TavX.append(TXsat.mean(axis=1))
    TavE.append(TEsat.mean(axis=1))
    TavR.append(TRsat.mean(axis=1))
    
Ravm = np.asarray(Ravm).T
RavX = np.asarray(RavX).T
RavE = np.asarray(RavE).T
RavR = np.asarray(RavR).T

Tavm = np.asarray(Tavm).T
TavX = np.asarray(TavX).T
TavE = np.asarray(TavE).T
TavR = np.asarray(TavR).T


df = pd.DataFrame({'DC Field (G)':BxDC,
                  'R (mx)':Ravm[0],
                  'R (my)':Ravm[1],
                  'R (mz)':Ravm[2],
                  'R (mxB_x) (aNm)': RavX[0],
                  'R (mxB_y) (aNm)': RavX[1],
                  'R (mxB_z) (aNm)': RavX[2],
                  'R (edH_x) (aNm)': RavE[0],
                  'R (edH_y) (aNm)': RavE[1],
                  'R (edH_z) (aNm)': RavE[2],
                  'R (res_x) (aNm)':RavR[0],
                  'R (res_y) (aNm)':RavR[1],
                  'R (res_z) (aNm)':RavR[2],
                  'Phase (mx) (rad)':Tavm[0],
                  'Phase (my) (rad)':Tavm[1],
                  'Phase (mz) (rad)':Tavm[2],
                  'Phase (mxB_x) (rad)':TavX[0],
                  'Phase (mxB_y) (rad)':TavX[1],
                  'Phase (mxB_z) (rad)':TavX[2],
                  'Phase (edH_x) (rad)':TavE[0],
                  'Phase (edH_y) (rad)':TavE[1],
                  'Phase (edH_z) (rad)':TavE[2],
                  'Phase (res_x) (rad)':TavR[0],
                  'Phase (res_y) (rad)':TavR[1],
                  'Phase (res_z) (rad)':TavR[2]
                  })
df.to_csv(path+savename,index=False)


