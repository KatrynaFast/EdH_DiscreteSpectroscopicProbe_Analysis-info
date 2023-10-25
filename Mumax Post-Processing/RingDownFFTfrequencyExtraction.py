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
path = r'' #file location for mumax simulation output table
pathNP = r'' #file location for mumax simulation (with no pinning) output table
savepath =r''

tab = pd.read_csv(path,delimiter='\t')
tabNP = pd.read_csv(pathNP,delimiter='\t')

m = tab[tab.columns[1:4]].to_numpy().T
B = tab[tab.columns[4:7]].to_numpy().T
t = tab[tab.columns[0]].to_numpy()

mNP = tabNP[tabNP.columns[1:4]].to_numpy().T
BNP = tabNP[tabNP.columns[4:7]].to_numpy().T
tNP = tabNP[tabNP.columns[0]].to_numpy()


ti = np.where(t==0)[0] #simulation resets time for every field step, can identify through locations where t=0
Bi = B[0][ti]*1e4

f0 = []
fM = []
lw = []

f0NP = []
FMNP = []
lwNP = []


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
        
    FT = np.fft.fft(mm) #calculate FFT of mx
    FTNP = np.fft.fft(mmNP)
    FTf = np.fft.fftfreq(mm.size,d=10e-12)*1e-6
    n = int(mm.size/2)
    
    FTy = np.fft.fft(mmy)
    
    
    
    ifit = np.where((FTf[1:n]<800)*(FTf[1:n]>100))[0]
    ffi = np.where(abs(FT[1:n][ifit])==abs(FT[1:n][ifit]).max())[0][0]
    ffiNP = np.where(abs(FTNP[1:n][ifit])==abs(FTNP[1:n][ifit]).max())[0][0]

    
    fM.append(FTf[1:n][ffi])
    FMNP.append(FTf[1:n][ffiNP])
    
    fit,cov = cf(Lorentz,FTf[1:n][ifit],abs(FT[1:n][ifit]),p0=[100,0,FTf[1:n][ifit][ffi],25])
    f0.append(fit[2])
    lw.append(fit[3])
    
    fitNP,covNP = cf(Lorentz,FTf[1:n][ifit],abs(FTNP[1:n][ifit]),p0=[100,0,FTf[1:n][ifit][ffiNP],25])
    f0NP.append(fitNP[2])
    lwNP.append(fitNP[3])
    
Bin1 = 40
Bin2 = 62.1

Bpin1 = 43.3
Bpin2 = 58


""" From Lorentz fits to FFT, compile list of resonance frequencies and linewidths for each applied DC field, 
then save to csv. """
lw = np.asarray(lw) 
f0 = np.asarray(f0)
lwNP = np.asarray(lwNP)
f0NP = np.asarray(f0NP)


df = pd.DataFrame({'DC Field (G)':Bi*4*np.pi,
                    'Gyro Frequency (MHz)':f0NP,
                    'FWHM (MHz)':lwNP})
df.to_csv(savepath+'RingDown_FFTFreq_FWWHM_frequencies_noPinning.csv',index=False)

