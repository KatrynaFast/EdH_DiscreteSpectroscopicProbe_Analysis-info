#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 09:50:37 2023

@author: katrynafast
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.close('all')

plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
# plt.rc('font',size=20)
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
# plt.rc('legend',fontsize=SMALLEST_SIZE)
plt.rc('axes',labelsize=15,grid=False)
plt.rc('lines',linewidth=3)
# plt.rc('grid',c='0.5',ls='-',lw=0.5)
plt.rc('savefig',format='png',dpi=400,bbox='tight',transparent=True)

dt = 10e-12
f0 = 208e6
fc = 0.3*f0
FO = 4

ti = 14424
td0 = 14423
tiA = 10000

mu0 = 4*np.pi*1e-7 #N/A^2
g = 2.
m_e = 9.1083837e-31 #kg #electron mass
e = -1.60217663e-19 #C #electron charge
gamma =  1.76086e11 #rad/sT
gprime = 4-g

RF = 'y'

#material constants 
Ms = 767e3 #A/m 
V = 5.42e-20

AEdH = -1e18*2*m_e*Ms*V/(e*g)
AmxB = 1e18*Ms*V


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
    ref = np.sqrt(2)*(-1j)*np.exp(1j*2*np.pi*fref*t) #reference signal in complex plane (and to align with np.sin at zero phase)
    mix = sig*ref #mixed signal with reference
    filt = lowPass(mix,fc,dt,FO) #apply low-pass filter to mixed signal to remove the f1 + f2 component

    X = filt.real 
    Y = filt.imag
    R = np.sqrt(X**2 + Y**2)
    Th = np.arctan2(Y,X)
    return X,Y,R,Th

pathPk1 = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Research/mumax/PyShroom/DC Hysteresis/TimeDomain_Hysteresis_LtH_v5_FineRes_LongRunTime_20nmPinningAt50G_208MHzDrive_FirstPeak39-45G.out/table.txt'
pathPk2 = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Research/mumax/PyShroom/DC Hysteresis/TimeDomain_Hysteresis_LtH_v5_FineRes_LongRunTime_20nmPinningAt50G_208MHzDrive_SecondPeak57-65G.out/table.txt'
pathMainLoop = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Research/mumax/PyShroom/DC Hysteresis/TimeDomain_Hysteresis_LtH_v4_20nmPinningAt50G_1Gydrive_208MHz.out/table.txt'
trajPath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Research/mumax/PyShroom/DC Hysteresis/DCHysteresis_LtH_v5_FineRes_20nmPinning_ExtractedVortexCorePositions.csv'

savepath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Conferences and Papers/NRC paper/Data Analysis Writeup/Simulation Plots/'

tabPk1 = pd.read_csv(pathPk1,delimiter='\t')
tabPk2 = pd.read_csv(pathPk2,delimiter='\t')
tabML = pd.read_csv(pathMainLoop,delimiter='\t')

mPk1 = tabPk1[tabPk1.columns[1:4]].to_numpy().T
BPk1 = tabPk1[tabPk1.columns[4:7]].to_numpy().T*1e4
tPk1 = tabPk1[tabPk1.columns[0]].to_numpy()

mPk2 = tabPk2[tabPk2.columns[1:4]].to_numpy().T
BPk2 = tabPk2[tabPk2.columns[4:7]].to_numpy().T*1e4
tPk2 = tabPk2[tabPk2.columns[0]].to_numpy()

mML = tabML[tabML.columns[1:4]].to_numpy().T
BML = tabML[tabML.columns[4:7]].to_numpy().T*1e4
tML = tabML[tabML.columns[0]].to_numpy()

traj = pd.read_csv(trajPath)
trajB = traj[traj.columns[0]].to_numpy()
trajX = traj[traj.columns[1]].to_numpy()
trajY = traj[traj.columns[2]].to_numpy()


#delete the repeated time values in both arrays (same length of arrays;
#deleting the last table entry for each field step because the time is repeated for the first
#table entry of the next field step. )
# mPk1 = np.delete(mPk1,np.s_[td0::ti],axis=1)
# BPk1 = np.delete(BPk1,np.s_[td0::ti],axis=1)
# tPk1 = np.delete(tPk1,np.s_[td0::ti])

# mML = np.delete(mML,np.s_[td0::ti],axis=1)
# BML = np.delete(BML,np.s_[td0::ti],axis=1)
# tML = np.delete(tML,np.s_[td0::ti])


Pk1i = 39
Pk1f = 45.1
Pk2i = 57
Pk2f = 65.1

indML_a = np.where(BML[0]<Pk1i)[0]

indML_b = np.where((BML[0]<Pk2i)*(BML[0]>Pk1f))[0]

indML_c = np.where(BML[0]>Pk2f)[0]

# indML1 = np.concatenate([indML_a1,indML_b1])
# 

tMLa = tML[indML_a]
tPk1 = tPk1 - tPk1[0] + tMLa[-1] + dt
tMLb = tML[indML_b] - tML[indML_b][0] + tPk1[-1] + dt
tPk2 = tPk2 - tPk2[0] + tMLb[-1] + dt
tMLc = tML[indML_c] - tML[indML_c][0] + tPk2[-1] + dt


### combine the tables into unique matrices
# B = np.concatenate([BML[:,indML_a1],BPk1,BML[:,indML_c],BPk2,BML[:,indML_b2]],axis=1)
# m = np.concatenate([mML[:,indML_a1],mPk1,mML[:,indML_c],mPk2,mML[:,indML_b2]],axis=1)
# t = np.concatenate([tMLa1,tPk1,tMLb1,tPk2,tMLb2])
B = np.concatenate([BML[:,indML_a],BPk1,BML[:,indML_b],BPk2,BML[:,indML_c]],axis=1)
m = np.concatenate([mML[:,indML_a],mPk1,mML[:,indML_b],mPk2,mML[:,indML_c]],axis=1)
t = np.concatenate([tMLa,tPk1,tMLb,tPk2,tMLc])


mxb = AmxB*np.cross(m.T,(B*1e-4).T).T
edh = AEdH*np.gradient(m,dt,axis=1)

res = mxb+edh

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


# plt.plot(BML[0],mML[0])
# plt.plot(BPk1[0],mPk1[0])

BxDC = np.unique(B[0])
NB = BxDC.size


Ravm = []
RavX = []
RavE = []
RavR = []

Tavm = []
TavX = []
TavE = []
TavR = []

for i in range(NB):
    # if i==NB-1: #need to treat last one differently I think...
    #     RTm = Rm[:,i*ti:]
    #     RTX = RX[:,i*ti:]
    #     RTE = RE[:,i*ti:]
    # else:
        
    #break up locked-in arrays into smaller arrays of each individual DC field
    RTm = Rm[:,i*ti:(i+1)*ti-1]
    RTX = RX[:,i*ti:(i+1)*ti-1]
    RTE = RE[:,i*ti:(i+1)*ti-2]
    RTR = RR[:,i*ti:(i+1)*ti-2]
    
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
    
    #break up locked-in arrays into smaller arrays of each individual DC field
    TTm = Tm[:,i*ti:(i+1)*ti-1]
    TTX = TX[:,i*ti:(i+1)*ti-1]
    TTE = TE[:,i*ti:(i+1)*ti-2]
    TTR = TR[:,i*ti:(i+1)*ti-2]
    
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
csvSavePath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Research/mumax/PyShroom/DC Hysteresis/'
df.to_csv(csvSavePath+'LongDurationSims_combined_v4_v5.csv',index=False)
                  


b1 = 34.9
b2a = 47
b2 = 70
b1b = 54.9
b2b = 70


# fitPk1 = np.polyfit(trajB[trajB<39],trajY[trajB<39],1)
# fitPin = np.polyfit(trajB[(trajB>=43.5)*(trajB<57)],trajY[(trajB>=43.5)*(trajB<57)],1)
# fitPk2 = np.polyfit(trajB[trajB>=60],trajY[trajB>=60],1)

# # plt.figure()
# fig,ax = plt.subplots(3,1,figsize=(6,12))
# ax[0].plot(BxDC[(BxDC>=b1)*(BxDC<b2a)],Ravm[0][(BxDC>=b1)*(BxDC<b2a)],'.-',ms=5,lw=0.5)

# ax0t = ax[0].twinx()
# ax0t.plot(trajB[(trajB>=b1)*(trajB<b2a)],trajY[(trajB>=b1)*(trajB<b2a)],'k.-',lw=0.5,ms=2)
# ax0t.plot(trajB[(trajB>=b1)*(trajB<b2a)],trajB[(trajB>=b1)*(trajB<b2a)]*fitPk1[0] +fitPk1[1],'m-',lw=1)
# ax0t.plot(trajB[(trajB<b2a)*(trajB>=40)],trajB[(trajB<b2a)*(trajB>=40)]*fitPin[0] + fitPin[1],'g-',lw=1)

# ax[1].plot(BxDC[(BxDC>=b1)*(BxDC<b2)],Ravm[0][(BxDC>=b1)*(BxDC<b2)],'.-',ms=5,lw=0.5)
# ax1t = ax[1].twinx()
# # ax0t.plot(trajB[(trajB>=b1)*(trajB<b2a)],trajY[(trajB>=b1)*(trajB<b2a)],'k.-',lw=0.5,ms=2)
# ax1t.plot(trajB,trajY,'k.-',lw=0.5,ms=2)
# ax1t.plot(trajB,trajB*fitPin[0]+fitPin[1],'g-',lw=1)


# ax[2].plot(BxDC[(BxDC>=b1b)*(BxDC<b2b)],Ravm[0][(BxDC>=b1b)*(BxDC<b2b)],'.-',ms=5,lw=0.5)
# ax2t = ax[2].twinx()

# ax2t.plot(trajB[(trajB>=b1b)*(trajB<b2b)],trajY[(trajB>=b1b)*(trajB<b2b)],'k.-',lw=0.5,ms=2)
# ax2t.plot(trajB[(trajB>=b1b)*(trajB<b2b)],trajB[(trajB>=b1b)*(trajB<b2b)]*fitPk2[0]+fitPk2[1],'-',color='purple',lw=1)
# ax2t.plot(trajB[(trajB>=b1b)*(trajB<=63)],trajB[(trajB>=b1b)*(trajB<=63)]*fitPin[0] + fitPin[1],'g-',lw=1)

# plt.figure(figsize=(8,6))
# b1 = 34.9
# b2 = 67.5
# plt.plot(BxDC[(BxDC>=b1)*(BxDC<b2)]/(4*np.pi),RavE[1][(BxDC>=b1)*(BxDC<b2)],'.-',ms=5,lw=1)
# plt.xlabel('DC Field, $H_x^{DC}$ (kA/m)')
# plt.ylabel(r'$\tau^{EdH}_y$ (aNm)')
# plt.xlim(b1/(4*np.pi),b2/(4*np.pi))
# plt.twinx()
# plt.plot(trajB/(4*np.pi),trajY,'k.-',ms=3,lw=0.5)
# plt.ylabel('Vortex core position (nm)')



fig,ax =plt.subplots(1,1,figsize=(8,6))
# ax.fill_between(BxDC[(BxDC>=b1)*(BxDC<b2)]/(4*np.pi),0.25,3.2,alpha=1,facecolor='w',edgecolor='k',linestyle='--')
ax.plot(BxDC/(4*np.pi),RavE[1],'.-',ms=5,lw=1)
ax.set_xlabel('DC Field, $H_x^{DC}$ (kA/m)')
ax.set_ylabel(r'$\tau^{EdH}_y$ (aNm)')
# ax.set_xlim(BxDC[0]/(4*np.pi),BxDC[-1]/(4*np.pi))
ax.set_xlim(BxDC[0]/(4*np.pi),8)
ax.set_xlim(25/(4*np.pi),75/(4*np.pi))





# left, bottom, width, height = [0.64, 0.6, 0.25, 0.25]
# ax2 = fig.add_axes([left, bottom, width, height])
# # ax.set_xlim(b1/(4*np.pi),b2/(4*np.pi))


# ax2.plot(BxDC[(BxDC>=b1)*(BxDC<b2)]/(4*np.pi),RavE[1][(BxDC>=b1)*(BxDC<b2)],'.-',ms=5,lw=1)

# plt.savefig(savepath+'LongDurationSimulation_EdHTorque_y.png',format='png',dpi=400,bbox_inches='tight')

# plt.savefig(savepath+'FineFieldStep_EdHy_w_traj.png',format='png',dpi=400,bbox_inches='tight')
    
    