# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 09:33:55 2023

@author: AWGPC2
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf


plt.close('all')

def Lorentz(x,A,b,x0,gam):
    return b + (A*gam/2)/((x-x0)**2 + (gam/2)**2)
# 
# path = r'D:\Katryna\PyShroom\DC Hysteresis\Hysteresis_LtH_w20nmPinningAt50G_0p8Ms FromThreadripper\RingDowns_36-46G_50nsRuntime_v2.out\table.txt'
path = r'D:\Katryna\PyShroom\DC Hysteresis\Hysteresis_LtH_w20nmPinningAt50G_0p8Ms FromThreadripper\RingDowns_25-75G_100nsRuntime.out\table.txt'


pathNP = r'D:\Katryna\PyShroom\DC Hysteresis\Hysteresis_LtH_noPinning_take2.out\NoPinnning_RingDowns_25-75G_100nsRuntime.out\table.txt'
savepath =r'D:/Katryna/PyShroom/DC Hysteresis/'
tab = pd.read_csv(path,delimiter='\t')
tabNP = pd.read_csv(pathNP,delimiter='\t')

m = tab[tab.columns[1:4]].to_numpy().T
B = tab[tab.columns[4:7]].to_numpy().T
t = tab[tab.columns[0]].to_numpy()

mNP = tabNP[tabNP.columns[1:4]].to_numpy().T
BNP = tabNP[tabNP.columns[4:7]].to_numpy().T
tNP = tabNP[tabNP.columns[0]].to_numpy()


ti = np.where(t==0)[0]
# i = 4
Bi = B[0][ti]*1e4

f0 = []
fM = []
lw = []

f0NP = []
FMNP = []
lwNP = []
# plt.plot(t[i*ti[i]:(i+1)*ti[i]],m[0][i*ti[i]:(i+1)*ti[i+1]])
# plt.plot(t[ti[i]:ti[i+1]-101],m[0][ti[i]:ti[i+1]-101])

for i in range(ti.size):
    if i==ti.size-1:
        mm = m[0][ti[i]:]#[:3000]
        tt = t[ti[i]:]#[:3000]
        
        mmy = m[0][ti[i]:]
        
        mmNP = mNP[0][ti[i]:]#[:3000]
        ttNP = t[ti[i]:]#[:3000]
    else:
        mm = m[0][ti[i]:ti[i+1]][:-100]#[:3000]
        tt = t[ti[i]:ti[i+1]][:-100]#[:3000]
        
        mmy = m[1][ti[i]:ti[i+1]][:-100]
        
        mmNP = mNP[0][ti[i]:ti[i+1]][:-100]#[:3000]
        ttNP = tNP[ti[i]:ti[i+1]][:-100]#[:3000]
        
    FT = np.fft.fft(mm)
    FTNP = np.fft.fft(mmNP)
    FTf = np.fft.fftfreq(mm.size,d=10e-12)*1e-6
    n = int(mm.size/2)
    
    FTy = np.fft.fft(mmy)
    
    
    
    ifit = np.where((FTf[1:n]<800)*(FTf[1:n]>100))[0]
    # ifitNP = np.where
    # ffi = np.where(abs(FT[10:n])==abs(FT[10:n]).max())[0][0]
    ffi = np.where(abs(FT[1:n][ifit])==abs(FT[1:n][ifit]).max())[0][0]
    ffiNP = np.where(abs(FTNP[1:n][ifit])==abs(FTNP[1:n][ifit]).max())[0][0]

    
    fM.append(FTf[1:n][ffi])
    FMNP.append(FTf[1:n][ffiNP])
    
    # fit,cov = cf(Lorentz,FTf[10:n],abs(FT[10:n]),p0=[10,0,ffi,10])
    fit,cov = cf(Lorentz,FTf[1:n][ifit],abs(FT[1:n][ifit]),p0=[100,0,FTf[1:n][ifit][ffi],25])
    f0.append(fit[2])
    lw.append(fit[3])
    
    fitNP,covNP = cf(Lorentz,FTf[1:n][ifit],abs(FTNP[1:n][ifit]),p0=[100,0,FTf[1:n][ifit][ffiNP],25])
    f0NP.append(fitNP[2])
    lwNP.append(fitNP[3])
    """
    plt.figure()
    plt.plot(FTf[1:n],abs(FT[1:n]),'.-')
    # plt.plot(FTf[1:n],abs(FTy[1:n]),'.-')
    plt.xlim(0,800)
    plt.title(int(np.round(Bi[i],0)))
    plt.plot(FTf[1:n],Lorentz(FTf[1:n],*fit),'k-')
    plt.plot([208,208],[0,Lorentz(FTf[1:n],*fit).max()],'r--',lw=1)
    """
    
# plt.figure()
# plt.plot(Bi,fM,label='max location')
# plt.plot(Bi,f0,label='fit peak')
# plt.legend()
# plt.xlabel('$B_x^{DC} (G)')
# plt.ylabel('FFT peak frequency (MHz)')
Bin1 = 40
Bin2 = 62.1

Bpin1 = 43.3
Bpin2 = 58

# Bi = Bi/(4*np.pi)
# Bpin1 = Bpin1/(4*np.pi)
# Bpin2 = Bpin2/(4*np.pi)
# Bin1 = Bin1/(4*np.pi)
# Bin2 = Bin2/(4*np.pi)

lw = np.asarray(lw)
f0 = np.asarray(f0)
lwNP = np.asarray(lwNP)
f0NP = np.asarray(f0NP)

# plt.figure()
# plt.plot(Bi,lw,'.-',label='linewidth')
# plt.legend(loc='upper left')
# plt.fill_between([Bin1,Bin2],15,250,color='gray',alpha=0.3)
# plt.fill_between([Bpin1,Bpin2],15,250,color='gray',alpha=0.5)
# plt.ylim(15,249)
# plt.ylabel('FWHM (MHz)')
# plt.xlabel('DC Field, $B_x^{DC}$ (G)')
# plt.grid()
# plt.twinx()
# plt.ylabel('Frequency (MHz)')
# plt.plot(Bi,f0,'k.-',label='frequency')

# plt.fill_between(Bi,f0-lw/2,f0+lw/2)

plt.figure()
plt.xlabel('DC Field, $H_x^{DC}$ (kA/m)')
plt.ylabel('Gyro Frequency (MHz)')


plt.fill_between([Bin1,Bin2],100,700,color='rebeccapurple',alpha=0.3)
plt.fill_between([Bpin1,Bpin2],100,700,color='gray',alpha=0.85)

plt.plot(Bi,f0NP,'.-',color='maroon',label='Without pinning')
plt.fill_between(Bi,f0NP-lwNP/2,f0NP+lwNP/2,color='darkorange',alpha=0.7)

plt.plot(Bi,f0,'k.-',label='With pinning')
plt.fill_between(Bi,f0-lw/2,f0+lw/2,alpha=0.7)
plt.legend(loc='upper left',framealpha=0.2,fontsize=8)

plt.xlim(Bi[0],Bi[-1])
# plt.grid()
# plt.ylim(150,675)

plt.ylim(140,690)
plt.plot(Bi,208*np.ones_like(Bi),'k--',lw=1)




# plt.plot(Bi,np.ones_like(Bi)*208,'g-')

# Bf = np.concatenate(([0],Bi))
# plt.xlim(0,Bi[-1])
# fit = np.polyfit(Bi,f0NP,1)
# plt.plot(Bf,Bf*fit[0] + fit[1],'g-')


# df = pd.DataFrame({'DC Field (G)':Bi*4*np.pi,
#                    'Gyro Frequency (MHz)':f0NP,
#                    'FWHM (MHz)':lwNP})
# df.to_csv(savepath+'RingDown_FFTFreq_FWWHM_frequencies_noPinning.csv',index=False)


# plt.legend(loc='upper right')

# plt.plot(Bi,208*np.ones_like(Bi),'m-')

# plt.plot([40.15,40.15],[150,600],'k--')
# plt.plot([61.95,61.95],[150,600],'k--')

# plt.fill_between([40.05,43.2],175,300,color='gray',alpha=0.5)
# plt.fill_between([57.95,62.05],175,300,color='gray',alpha=0.5)

# plt.savefig(savepath+'RingDown_FFTFreq_FWHM_pinned_v3.png',format='png',dpi=400,bbox_inches='tight')