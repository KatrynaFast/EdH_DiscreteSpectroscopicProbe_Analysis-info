    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:56:20 2022

@author: katrynafast

This script is used to load and process multiple subsequent hysteresis loops, collected back-to-back. 
The files for each set contain net signal (5x) and background (5x), each of which have the same base 
file name, with an index counting up from 0. 

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf

plt.close('all')

plt.rcParams['lines.linewidth']=2
plt.rcParams['lines.markersize']=0
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['axes.labelsize']=22
plt.rcParams['xtick.labelsize']=15
plt.rcParams['ytick.labelsize']=15
plt.rcParams['axes.titlesize']=18
plt.rcParams['legend.fontsize']=15
plt.rc('savefig',format='png',dpi=300,bbox='tight')


cols = ['#2F327A','#3E848D','#C36C42','#F8D554','#556A2F','#CACAE2','#C2C2C2','#9E5F42']


savepath = '' #path to save folder
path = '' #path to data

bgName = '{0}MHz_180mVH{2}RF_30ms_bg1598_Chan 0_00{1}' #name of background data files
sigName = '{0}MHz_180mVH{2}RF_30ms_Chan 0_00{1}' #name of net signal data files
f = [208.15]
HRF = ['y']
Datind = [0,1,2,3,4]

RA = []
XA = []
YA = []
HA = []
PA = []

bRA = []
bXA = []
bYA = []
bHA = []
bPA = []

RB = []
XB = []
YB = []
HB = []
PB = []

bRB = []
bXB = []
bYB = []
bHB = []
bPB = []

dRA = []
dXA = []
dYA = []
dPA = []

dRB = []
dXB = []
dYB = []
dPB = []

dR = []
dX = []
dY = []
dP = []


def r2a(x,a,b,c,d):
        return a/(x+b)**2 + d/(x+b) - c

HfA = []
HfB = []

j = 3
print('Loading data ...')
for k in range(len(HRF)):
        for j in range(len(f)):
            print(f[j],' MHz (H{} drive)'.format(HRF[k]))
            """ Create empty lists in which to store loaded data """
            r = []
            x = []
            y = []
            h = []
            hr = []
            p = []
            hf = []
            
            br = []
            bx = []
            by = []
            bh = []
            bhr = []
            bp = []
            bhf = []

            L = 120000 
            
        """ Load all data files into the empty lists and convert units:
                Fields: G -> kA/m
                Signals: V -> uV
                Phase: degrees  
            Each file has variable length, so the sizes are kept track of to 
            truncate to smallest value to enable processing of the data together  
            """
        for i in range(len(Datind)):
            dat = np.loadtxt(path+sigName.format(f[j],Datind[i],HRF[k]),delimiter='\t') #load and separate each set of data
            hr.append(dat[:,0]/(4*np.pi))
            h.append(dat[:,1]/(4*np.pi))
            x.append(dat[:,2]*1e6)
            y.append(dat[:,3]*1e6)
            r.append(dat[:,4]*1e6)
            p.append(dat[:,5])
            
            bdat = np.loadtxt(path+bgName.format(f[j],Datind[i],HRF[k]),delimiter='\t') #load and separate each set of background data
            bhr.append(bdat[:,0]/(4*np.pi))
            bh.append(bdat[:,1]/(4*np.pi))
            bx.append(bdat[:,2]*1e6)
            by.append(bdat[:,3]*1e6)
            br.append(bdat[:,4]*1e6)
            bp.append(bdat[:,5])
            if dat.shape[0]<L:
                L = dat.shape[0]
            if bdat.shape[0]<L:
                L = bdat.shape[0]
            
                
        for i in range(len(Datind)):
            r[i] = r[i][:L] #cut off points beyond index L to force all arrays to be same length
            x[i] = x[i][:L]
            y[i] = y[i][:L]
            h[i] = h[i][:L]
            hr[i] = hr[i][:L]
            p[i] = p[i][:L]
            
            br[i] = br[i][:L]
            bx[i] = bx[i][:L]
            by[i] = by[i][:L]
            bh[i] = bh[i][:L]
            bhr[i] = bh[i][:L]
            bp[i] = bp[i][:L]
            
            l = int(L/2) #identify halfway point to break data into sweep down (:l) and sweep up (l:) components
            ci = 0
            ci2 = 3
            
            indsA = np.arange(0,h[i][:l].shape[0],1)
            indsB = np.arange(0,h[i][l:].shape[0],1)
           
            """ Fit field data in both directions to a 2nd order polynomial function 
            in order to clean up the data due to high noise levels in recorded fields"""
            HfitA,HcovA = cf(r2a,indsA,h[i][:l],p0=[55e8,10000,1.5,10000]) 
            HfitB,HcovB = cf(r2a,indsB,h[i][l:],p0=[55e8,-65000,1.5,10000])
           
            hfA = r2a(indsA,*HfitA) 
            hfB = r2a(indsB,*HfitB)
            hfB = hfB - hfB[0]
            hfA = hfA - hfA[-1]
            
            hf.append(np.concatenate((hfA,hfB)))
                        
            bindsA = np.arange(0,bh[i][:l].shape[0],1)
            bindsB = np.arange(0,bh[i][l:].shape[0],1)
           
            bHfitA,bHcovA = cf(r2a,bindsA,bh[i][:l],p0=[55e8,10000,1.5,10000])
            bHfitB,bHcovB = cf(r2a,bindsB,bh[i][l:],p0=[55e8,-65000,1.5,10000])
            
            bhfA = r2a(bindsA,*HfitA) 
            bhfB = r2a(bindsB,*HfitB)
            bhfB = bhfB - bhfB[0]
            bhfA = bhfA - bhfA[-1]
            
            bhf.append(np.concatenate((bhfA,bhfB)))
             
            
        
        """ Take average of each sweep, then fill empty lists with data separated
        into sweep down and sweep up """
        rm = np.asarray(r).mean(axis=0) 
        xm = np.asarray(x).mean(axis=0)
        ym = np.asarray(y).mean(axis=0)
        pm = np.asarray(p).mean(axis=0)
        hm = np.asarray(h).mean(axis=0)
        hrm = np.asarray(hr).mean(axis=0)
        
        brm = np.asarray(br).mean(axis=0)
        bxm = np.asarray(bx).mean(axis=0)
        bym = np.asarray(by).mean(axis=0)
        bpm = np.asarray(bp).mean(axis=0)
        bhm = np.asarray(bh).mean(axis=0)
        bhrm = np.asarray(bhr).mean(axis=0)
        
        rA = rm[:int(L/2)] #take first half of data for high-to-low field sweep
        xA = xm[:int(L/2)]
        yA = ym[:int(L/2)]
        pA = pm[:int(L/2)]
        hA = hm[:int(L/2)]
        
        rB = rm[int(L/2):] #take second half of data for low-to-high field sweep
        xB = xm[int(L/2):]
        yB = ym[int(L/2):]
        pB = pm[int(L/2):]
        hB = hm[int(L/2):]
        
        brA = brm[:int(L/2)] #split background
        bxA = bxm[:int(L/2)]
        byA = bym[:int(L/2)]
        bpA = bpm[:int(L/2)]
        bhA = bhm[:int(L/2)]
        
        brB = brm[int(L/2):]
        bxB = bxm[int(L/2):]
        byB = bym[int(L/2):]
        bpB = bpm[int(L/2):]
        bhB = bhm[int(L/2):]
        
        RA.append(rA)
        XA.append(xA)
        YA.append(yA)
        PA.append(pA)
        HA.append(hA)
        
        RB.append(rB)
        XB.append(xB)
        YB.append(yB)
        PB.append(pB)
        HB.append(hB)
        
        bRA.append(brA)
        bXA.append(bxA)
        bYA.append(byA)
        bPA.append(bpA)
        bHA.append(bhA)
        
        bRB.append(brB)
        bXB.append(bxB)
        bYB.append(byB)
        bPB.append(bpB)
        bHB.append(bhB)
        
        """ Perform phasor subtractions of background from net signal """
        dxa = xA - bxA #background subtractions
        dya = yA - byA
        dra = np.sqrt(dxa**2 + dya**2)
        dpa = np.arctan2(dya,dxa)
        
        dxb = xB - bxB
        dyb = yB - byB
        drb = np.sqrt(dxb**2 + dyb**2)
        dpb = np.arctan2(dyb,dxb)
        
        dRA.append(dra)
        dXA.append(dxa)
        dYA.append(dya)
        dPA.append(dpa)
    
        dRB.append(drb)
        dXB.append(dxb)
        dYB.append(dyb)
        dPB.append(dpb)
        
        dx = xm - bxm
        dy = ym - bym
        dr = np.sqrt(dx**2 + dy**2)
        dp = np.arctan2(dy,dx)
        
        dR.append(dr)
        dY.append(dy)
        dX.append(dx)
        dP.append(dp)
        
        """ Save processed data for individual hysteresis loop  """
        ARR = [hf,x,y,r,p,bhf,bx,by,br,bp]
        datSavePath = ''
        saveName = ''
        np.save(datSavePath+saveName,ARR) 
 

""" Save output for averaged, separated hysteresis loops"""

n = []
for i in range(len(HA)):
        n.append(len(HA[i]))
        n.append(len(HB[i]))

n = np.asarray(n)
m = n.min()

for i in range(len(HA)):
        hrf = int((i/4)//1)
        df = pd.DataFrame({
            'H [kA/m] (HtL)':HA[i][:m], #0
            'X [uV] (HtL)':dXA[i][:m], #1
            'Y [uV] (HtL)':dYA[i][:m], #2
            'R [uV] (HtL)':dRA[i][:m], #3
            'Phase [rad] (HtL)':dPA[i][:m], #4
            'H [kA/m] (LtH)':HB[i][:m], #5
            'X [uV] (LtH)':dXB[i][:m], #6
            'Y [uV] (LtH)':dYB[i][:m], #7
            'R [uV] (LtH)':dRB[i][:m], #8
            'Phase [rad] (LtH)':dPB[i][:m] #9
            })
        bgSub_saveName = ''
        df.to_csv(savepath+bgSub_saveName,header=True,index=False)

for i in range(len(HA)):
        hrf = int((i/4)//1)
        df = pd.DataFrame({
            'H [kA/m] (HtL)':bHA[i][:m], #0
            'X [uV] (HtL)':bXA[i][:m], #1
            'Y [uV] (HtL)':bYA[i][:m], #2
            'R [uV] (HtL)':bRA[i][:m], #3
            'Phase [rad] (HtL)':bPA[i][:m], #4
            'H [kA/m] (LtH)':bHB[i][:m], #5
            'X [uV] (LtH)':bXB[i][:m], #6
            'Y [uV] (LtH)':bYB[i][:m], #7
            'R [uV] (LtH)':bRB[i][:m], #8
            'Phase [rad] (LtH)':bPB[i][:m] #9
            })
        bg_saveName = ''
        df.to_csv(savepath+bg_saveName,header=True,index=False)
       
for i in range(len(HA)):
        hrf = int((i/4)//1)
        df = pd.DataFrame({
            'H [kA/m] (HtL)':HA[i][:m], #0
            'X [uV] (HtL)':XA[i][:m], #1
            'Y [uV] (HtL)':YA[i][:m], #2
            'R [uV] (HtL)':RA[i][:m], #3
            'Phase [rad] (HtL)':PA[i][:m], #4
            'H [kA/m] (LtH)':HB[i][:m], #5
            'X [uV] (LtH)':XB[i][:m], #6
            'Y [uV] (LtH)':YB[i][:m], #7
            'R [uV] (LtH)':RB[i][:m], #8
            'Phase [rad] (LtH)':PB[i][:m] #9
            })
        raw_saveName = ''
        df.to_csv(savepath+raw_saveName,header=True,index=False)
