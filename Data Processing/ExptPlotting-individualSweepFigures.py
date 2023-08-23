    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:56:20 2022

@author: katrynafast
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

# path = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/Shared drives/NRC Experiment/data/2021/0607 Hyst loops HyHz 125 MHz PD HyAmp laser power adj/'
savepath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Conferences and Papers/NRC paper/Data/208 MHz/Raw Hysteresis Data/'
# sigName = '{0}MHz_H{2}RF_5mW_Chan 0_00{1}'
# bgName = '{0}MHz_H{2}RF_BG1598nm_5.3mW_Chan 0_00{1}'
# # f = [3,20.82,64.4,114.1]
# f = [3]
# FI = int(np.round(f,0))

# Datind = [0,1,2,3,4]
# # Datind = [0,1,2,3,4]

cols = ['#2F327A','#3E848D','#C36C42','#F8D554','#556A2F','#CACAE2','#C2C2C2','#9E5F42']

# path = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/Shared drives/NRC Experiment/data/2021/0525 High field HzRF hysteresis /'
path = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/Shared drives/NRC Experiment/data/2021/0628 Hyst of 208MHz mode/'
# savepath = '//Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Research/Barclavian Device/NRC Apparatus/Plots/2023-02-24 208 MHz Hysteresis from 0628/'
# sigName = '{0}MHz_H{2}RF_Chan 0_00{1}'
# bgName = '{0}MHz_H{2}RF_BG1598nm_Chan 0_00{1}'
sigName = 'Hyst_HyRF_180mVdrive_{0}MHz_Chan 0_00{1}'
bgName = 'Hyst_HyRF_180mVdrive_{0}MHzBG1598nm_Chan 0_00{1}'
# HRF=['z']
f = [208]
FI = 208
# f = [2.987,20.73,64.31,114.66]
# f = [114.66]
# Datind = [[0,1,2,3,4],[0,1,2,4]]
Datind = np.array([0,1,2,3,4])



# path = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Research/Barclavian Device/NRC Apparatus/Plots/2022-12-08 Single Sweeps From 2021-0525 and 2021-0514/0514 HyRF Hysteresis/data/'
# savepath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Research/Barclavian Device/NRC Apparatus/Plots/2022-12-08 Single Sweeps From 2021-0525 and 2021-0514/0514 HyRF Hysteresis/'
# sigName = '{0}MHz_180mVHyRF_30msTC_Chan 0_00{1}'
# bgName = '{0}MHz_180mVHyRF_30msTC_BG1598nm_Chan 0_00{1}'
HRF=['y']
# # f = [2.987,20.76,64.34,114.3]
# f=[20.76]


savepath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Conferences and Papers/NRC paper/Data/208 MHz/0611 Hy Hz/'
# path = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/Shared drives/NRC Experiment/data/2021/0611 Hysteresis HyHz 208MHz/'
path = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/Shared drives/NRC Experiment/data/2021/0514 High field hysteresis with backgrounds/'
# bgName = '{0}MHz_H{2}RF_BG1598nm_Chan 0_00{1}'
# sigName = '{0}MHz_H{2}RF_Chan 0_00{1}'
bgName = '{0}MHz_180mVH{2}RF_30ms_bg1598_Chan 0_00{1}'
sigName = '{0}MHz_180mVH{2}RF_30ms_Chan 0_00{1}'
f = [208.15]
HRF = ['y']
Datind = [0,1,2,3,4]


# savepath = '/Volumes/GoogleDrive/My Drive/Research/Barclavian Device/NRC Apparatus/Plots/2022-12-08 Single Sweeps From 2021-0525/'
# savepath = '/Volumes/GoogleDrive/My Drive/Research/Barclavian Device/NRC Apparatus/Plots/2022-12-08 Single Sweeps From 2021-0525 and 2021-0514/0514 HyRF Hysteresis/'
# f = [64.4]


# HRF = ['y','z']
# HRF = ['z']




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

    
# bHfA = []
# bHfB = []
# for i in range(len(HA)):
#     indsA = np.arange(0,bHrA[i].shape[0],1)
#     indsB = np.arange(0,bHrB[i].shape[0],1)
#     # HfitA,HcovA = cf(r2a,indsA,HA[i],p0=[55e8,10000,1.5,10000])
#     # HfitB,HcovB = cf(r2a,indsB,HB[i],p0=[55e8,-65000,1.5,10000])
#     HfitA,HcovA = cf(r2a,indsA,bHrA[i],p0=[55e8,10000,1.5,10000])
#     HfitB,HcovB = cf(r2a,indsB,bHrB[i],p0=[55e8,-65000,1.5,10000])
#     hfA = r2a(indsA,*HfitA) 
#     hfB = r2a(indsB,*HfitB)
#     hfB = hfB - hfB[0]
#     hfA = hfA - hfA[-1]
#     bHfA.append(hfA)
#     bHfB.append(hfB)
j = 3
print('Loading data ...')
for k in range(len(HRF)):
        for j in range(len(f)):
            print(f[j],' MHz (H{} drive)'.format(HRF[k]))
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
            # plt.figure()
            L = 120000
            
            fig,ax = plt.subplots(4,len(Datind),figsize=(30,20))
            ax[0,0].set_ylabel('Signal Magnitude ($\mu$V)')
            ax[1,0].set_ylabel('Signal Phase ($^{\circ}$)')
            ax[2,0].set_ylabel('BG Magnitude ($\mu$V)')
            ax[3,0].set_ylabel('BG Phase ($^{\circ}$)')
            plt.setp(ax,xlim=(-0.25,10))
            # plt.suptitle('{0} MHz - H{1} drive'.format(f[j],HRF[k]),fontsize=28,y=0.925)
            plt.subplots_adjust(hspace=0.1,wspace=0.2)
            fig.align_ylabels()
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
            
            
            
            dx = dat[:int(L/2),2] - bdat[:int(L/2),2] #Phasor subtraction 
            dy = dat[:int(L/2),3] - bdat[:int(L/2),3]
            dh = dat[:int(L/2),0]/(4*np.pi)
            dr = np.sqrt(dx**2 + dy**2)
            # plt.plot(dh,dr*1e6)
            # plt.title('{0} MHz - H{1} drive'.format(f[j],HRF[k]))
            # plt.xlim(0,8)
            # print(L)
            # print(dat.shape)
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
            
            l = int(L/2)
            ci = 0
            ci2 = 3
            
            
            # for i in range(len(HA)):
            indsA = np.arange(0,h[i][:l].shape[0],1)
            indsB = np.arange(0,h[i][l:].shape[0],1)
            # # HfitA,HcovA = cf(r2a,indsA,HA[i],p0=[55e8,10000,1.5,10000])
            # # HfitB,HcovB = cf(r2a,indsB,HB[i],p0=[55e8,-65000,1.5,10000])
            HfitA,HcovA = cf(r2a,indsA,h[i][:l],p0=[55e8,10000,1.5,10000])
            HfitB,HcovB = cf(r2a,indsB,h[i][l:],p0=[55e8,-65000,1.5,10000])
            hfA = r2a(indsA,*HfitA) 
            hfB = r2a(indsB,*HfitB)
            hfB = hfB - hfB[0]
            hfA = hfA - hfA[-1]
            
            hf.append(np.concatenate((hfA,hfB)))
            
            
            bindsA = np.arange(0,bh[i][:l].shape[0],1)
            bindsB = np.arange(0,bh[i][l:].shape[0],1)
            # # HfitA,HcovA = cf(r2a,bindsA,HA[i],p0=[55e8,10000,1.5,10000])
            # # HfitB,HcovB = cf(r2a,bindsB,HB[i],p0=[55e8,-65000,1.5,10000])
            bHfitA,bHcovA = cf(r2a,bindsA,bh[i][:l],p0=[55e8,10000,1.5,10000])
            bHfitB,bHcovB = cf(r2a,bindsB,bh[i][l:],p0=[55e8,-65000,1.5,10000])
            bhfA = r2a(bindsA,*HfitA) 
            bhfB = r2a(bindsB,*HfitB)
            bhfB = bhfB - bhfB[0]
            bhfA = bhfA - bhfA[-1]
            
            bhf.append(np.concatenate((bhfA,bhfB)))
            
            # HfA.append(hfA)
            # HfB.append(hfB)
            # ax[0,i].plot(hfA,r[i][:l],color=cols[ci],label='High to Low Sweep')
            # ax[1,i].plot(hfA,(p[i][:l]),color=cols[ci])
            # # ax[1,i].plot(hfA,np.rad2deg(np.unwrap(np.deg2rad(p[i][:l]))))
            # ax[2,i].plot(bhfA,br[i][:l],color=cols[ci])
            # # ax[3,i].plot(bhfA,np.rad2deg(np.unwrap(np.deg2rad(bp[i][:l]))))
            # ax[3,i].plot(bhfA,bp[i][:l],color=cols[ci])
            
            ax[0,i].plot(hfB,r[i][l:],color=cols[ci],label='Low to High Sweep')
            ax[1,i].plot(hfB,(p[i][l:]),color=cols[ci])
            # ax[1,i].plot(hfB,np.rad2deg(np.unwrap(np.deg2rad(p[i][l:]))))
            ax[2,i].plot(bhfB,br[i][l:],color=cols[ci])
            # ax[3,i].plot(bhfB,np.rad2deg(np.unwrap(np.deg2rad(bp[i][l:]))))
            ax[3,i].plot(bhfB,bp[i][l:],color=cols[ci])
            
            # ax[0,i].legend(loc='upper left')
            
            # ax[0,i].plot(h[i][:l],r[i][:l],color=cols[ci])
            # ax[1,i].plot(h[i][:l],(p[i][:l]),color=cols[ci])
            # # ax[1,i].plot(h[i][:l],np.rad2deg(np.unwrap(np.deg2rad(p[i][:l]))))
            # ax[2,i].plot(bh[i][:l],br[i][:l],color=cols[ci])
            # # ax[3,i].plot(bh[i][:l],np.rad2deg(np.unwrap(np.deg2rad(bp[i][:l]))))
            # ax[3,i].plot(bh[i][:l],bp[i][:l],color=cols[ci])
            
            # ax[0,i].plot(h[i][l:],r[i][l:],color=cols[ci2])
            # ax[1,i].plot(h[i][l:],(p[i][l:]),color=cols[ci2])
            # # ax[1,i].plot(h[i][l:],np.rad2deg(np.unwrap(np.deg2rad(p[i][l:]))))
            # ax[2,i].plot(bh[i][l:],br[i][l:],color=cols[ci2])
            # # ax[3,i].plot(bh[i][l:],np.rad2deg(np.unwrap(np.deg2rad(bp[i][l:]))))
            # ax[3,i].plot(bh[i][l:],bp[i][l:],color=cols[ci2])
            
            #z-driven mag lims
            # ax[0,i].set_ylim(250,800) ##3 MHz
            # ax[0,i].set_ylim(150,450) ##21 MHz
            # ax[0,i].set_ylim(20,75) ##64 MHz
            # ax[0,i].set_ylim(30,56) ##114 MHz
            
            # ax[3,i].set_xlabel('DC Field, $H_x^{DC}$ (kA/m)')
            ### 114 MHz 0607 limits: (LtH)
            # ax[0,i].set_ylim(0,250)
            # ax[1,i].set_ylim(-100,100)
            # ax[2,i].set_ylim(30,50)
            # ax[3,i].set_ylim(-70,-30)
            
            ### 114 MHz 0607 limits: (LtH)
            # ax[0,i].set_ylim(0,45)
            # ax[1,i].set_ylim(-182,182)
            # ax[2,i].set_ylim(5.5,15.5)
            # ax[3,i].set_ylim(8,62)
            
            ### 64 MHz 0607 limits:
            # ax[0,i].set_ylim(0,62)
            # ax[1,i].set_ylim(15,65)
            # ax[2,i].set_ylim(0,5)
            # ax[3,i].set_ylim(-182,182)
            
            ### 21 MHz 0607 limits:
            # ax[0,i].set_ylim(0,100)
            # ax[1,i].set_ylim(-182,182)
            # ax[2,i].set_ylim(0,4)
            # ax[3,i].set_ylim(-182,182)
            
            ### 3 MHz 0607 limits:
            # ax[0,i].set_ylim(0,165)
            # ax[1,i].set_ylim(-182,182)
            # ax[2,i].set_ylim(0,5.5)
            # ax[3,i].set_ylim(-182,182)
            
            
            ### 114 MHz 0607 limits: (HtL)
            # ax[0,i].set_ylim(0,25)
            # ax[1,i].set_ylim(-182,182)
            # ax[2,i].set_ylim(5.5,15.5)
            # ax[3,i].set_ylim(8,62)
            
            ### 64 MHz 0607 limits:
            # ax[0,i].set_ylim(0,62)
            # ax[1,i].set_ylim(15,65)
            # ax[2,i].set_ylim(0,5)
            # ax[3,i].set_ylim(-182,182)
            
            ### 21 MHz 0607 limits:
            # ax[0,i].set_ylim(0,120)
            # ax[1,i].set_ylim(0,100)
            # ax[2,i].set_ylim(0,4)
            # ax[3,i].set_ylim(-182,182)
            
            ### 3 MHz 0607 limits:
            # ax[0,i].set_ylim(0,195)
            # ax[1,i].set_ylim(0,182)
            # ax[2,i].set_ylim(0,5.5)
            # ax[3,i].set_ylim(-182,182)
            
            # plt.savefig(savepath+'SingleSweeps_HtL_{0}MHz_H{1}RF.png'.format(int(np.round(f[j],0)),HRF[k]))
            # plt.draw()
            
            # R.append()
        rm = np.asarray(r).mean(axis=0) #take mean of the the sweeps
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
        
        ARR = [hf,x,y,r,p,bhf,bx,by,br,bp]
        datSavePath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Conferences and Papers/NRC paper/Data/208 MHz/0525 Hz 0514 Hy/'
        np.save(datSavePath+'IndividualSweeps_0514_{0}drive_{1}MHz'.format(HRF[k],FI),ARR)
 
# ARR = [hf,x,y,r,p,bhf,bx,by,br,bp]

# datSavePath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Conferences and Papers/NRC paper/Data/Hysteresis/Individual Sweeps/'
# datSavePath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Conferences and Papers/NRC paper/Data/208 MHz/0611 Hy Hz'

# np.save(datSavePath+'IndividualSweeps_0611_{0}drive_{1}MHz'.format(HRF[0],FI),ARR)
# plt.savefig(savepath+'IndividualSweeps_signal_background_{}MHz_zDrive_0607_LtH.png'.format(FI))

# print(FI,' MHz:')
hmin = []
hminstd = []

i1 = 56700
i2 = 57750
for i in range(len(hr)):
       print(np.round(1000*hr[i][i1:i2].mean(),3),np.round(1000*hr[i][i1:i2].std(),3))
       hmin.append(hr[i][i1:i2].mean())
       hminstd.append(hr[i][i1:i2].std())
hmin = np.asarray(hmin)     
print('**********')
print(np.round(1000*hrm[i1:i2].mean(),3),np.round(1000*hrm[i1:i2].std(),3))
 


def r2a(x,a,b,c,d):
        return a/(x+b)**2 + d/(x+b) - c

HfA = []
HfB = []
for i in range(len(HA)):
        indsA = np.arange(0,HA[i].shape[0],1)
        indsB = np.arange(0,HB[i].shape[0],1)
        HfitA,HcovA = cf(r2a,indsA,HA[i],p0=[55e8,10000,1.5,10000])
        HfitB,HcovB = cf(r2a,indsB,HB[i],p0=[55e8,-65000,1.5,10000])
        hfA = r2a(indsA,*HfitA) 
        hfB = r2a(indsB,*HfitB)
        hfB = hfB - hfB[0]
        hfA = hfA - hfA[-1]
        HfA.append(hfA)
        HfB.append(hfB)

""" SAVE OUTPUT """

n = []
for i in range(len(HA)):
        n.append(len(HA[i]))
        n.append(len(HB[i]))

n = np.asarray(n)
m = n.min()
# for i in range(len(HRF)):
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
    # df.to_csv(savepath+'0611-BGsubtractedData_20230221_H_{0}RF_{1}MHz.csv'.format(HRF[hrf],f[i%4]),header=True,index=False)
    # df.to_csv(savepath+'0628-BGsubtractedData_20230405_H_{0}RF_{1}MHz.csv'.format(HRF[hrf],f[i%1]),header=True,index=False)

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
    # df.to_csv(savepath+'0611-BGAvg_20230221_H_{0}RF_{1}MHz.csv'.format(HRF[hrf],f[i%4]),header=True,index=False)
    # df.to_csv(savepath+'0628-BGAvg_20230405_H_{0}RF_{1}MHz.csv'.format(HRF[hrf],f[i%1]),header=True,index=False)

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
    # df.to_csv(savepath+'0611-RawSigAvg_20230221_H_{0}RF_{1}MHz.csv'.format(HRF[hrf],f[i%4]),header=True,index=False)
    # df.to_csv(savepath+'0628-RawSigAvg_20230405_H_{0}RF_{1}MHz.csv'.format(HRF[hrf],f[i%1]),header=True,index=False)
   

# path2 = '/Volumes/GoogleDrive/Shared drives/NRC Experiment/data/2021/0611 Hysteresis HyHz 208MHz/'
# sigName2 = '208.15MHz_H{0}RF_Chan 0_00{1}'
# bgName2 = '208.15MHz_H{0}RF_BG1598nm_Chan 0_00{1}'


# R = np.asarray(R)
# X = np.asarray(X)
# Y = np.asarray(Y)
# P = np.asarray(P)
# H = np.asarray(H)


# plt.plot(HA[4],dRA[4]/s[0] - dRA[5][:-1]/s[1])
# plt.plot(HA[4],dRA[4]/s[0] - dRA[6]/s[2])
# plt.plot(HA[4],dRA[4]/s[0] - dRA[7]/s[3])


print('Plotting...')
# plt.figure()

# plt.figure(figsize=(8,6))
# plt.plot(HA[0],dRA[0],label='High-to-Low field sweep')
# plt.plot(HB[0],dRB[0],label='Low-to-High field sweep')
# plt.legend(loc='upper right')
# plt.xlabel('DC Field, $H_x^{DC}$ (kA/m)')
# plt.ylabel('Measured Torque (arb.)')
# plt.xlim(0,20)


# plt.plot(51/(4*np.pi)*np.array([1,1]),[0,350],'k--')
# plt.plot(36/(4*np.pi)*np.array([1,1]),[0,350],'r--')
# savepath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Conferences and Papers/NRC paper/Data Analysis Writeup/Simulation Plots/'
# plt.savefig(savepath+'208MHz_Hysteresis_lowFieldRegion.png',format='png',dpi=400,bbox_inches='tight')

"""
s = []
ds = []
sE = []
dsE = []
p0 = []

for i in range(len(f)):
    inds = np.where(HA[i+4]>45)[0]
    s.append(dRA[i+4][inds].mean())
    ds.append(dRA[i+4][inds].std())
    indsE = np.where(HA[i]>45)[0]
    sE.append(dRA[i][indsE].mean())
    dsE.append(dRA[i][indsE].std())
    
for i in range(2*len(f)):
    inds = np.where(HA[i]>45)[0]
    p0.append(dPA[i][inds].mean())
    
def RA(N,x):
    avewindow = np.ones((N,))/N
    avg = np.convolve(x,avewindow,mode='valid')
    return avg


def r2(x,a,b,c): #fit field to 1/r^2 shape to smooth x-data
    return a/(x+b)**2 - c



# fig,ax = plt.subplots(2,1)
# N = 800
# for i in range(len(f)):
#     ax[0].plot(RA(N,HA[i+4]),RA(N,dRA[i+4])/s[i],label=f[i])
#     ax[1].plot(RA(N,HA[i]),RA(N,dRA[i])/s[i],label=f[i])
# ax[0].legend()
# ax[0].set_xlim(15,47)
# ax[1].set_xlim(15,47)
# ax[0].set_ylim(0.9,1.025)
# ax[1].set_ylim(0.15,0.42)
"""
cols = ['steelblue','darkorange','forestgreen','firebrick']

""" Full-scale plots """

# savepath2 = '/Volumes/GoogleDrive/My Drive/Research/Barclavian Device/NRC Apparatus/Plots/2022-12-07 Raw BG XY plots/'

# for i in range(len(f)):
#     fig,ax = plt.subplots(2,2,figsize=(10,6))
#     ax[0,0].plot(HA[i+4],XA[i+4])
#     ax[1,0].plot(HA[i+4],YA[i+4])
#     ax[0,0].set_title('High-to-Low - Hz drive\n{} MHz - raw data'.format(f[i]))
#     ax[0,0].set_ylabel('X')
#     ax[1,0].set_ylabel('Y')
    
#     ax[0,1].plot(bHA[i+4],bXA[i+4])
#     ax[1,1].plot(bHA[i+4],bYA[i+4])
#     ax[0,1].set_title('High-to-Low - Hz drive\n{} MHz - BG data'.format(f[i]))
    
#     plt.savefig(savepath2+'HzRF_HtL_Raw-BG-XY-{}MHz.png'.format(f[i]))
#     plt.close()
    
#     fig,ax = plt.subplots(2,2,figsize=(10,6))
#     ax[0,0].plot(HB[i+4],XB[i+4])
#     ax[1,0].plot(HB[i+4],YB[i+4])
#     ax[0,0].set_ylabel('X')
#     ax[1,0].set_ylabel('Y')
#     ax[0,0].set_title('Low-to-High - Hz drive\n{} MHz - raw data'.format(f[i]))
    
#     ax[0,1].plot(bHB[i+4],bXB[i+4])
#     ax[1,1].plot(bHB[i+4],bYB[i+4])
#     ax[0,1].set_title('Low-to-High - Hz drive\n{} MHz - BG data'.format(f[i]))
    
#     plt.savefig(savepath2+'HzRF_LtH_Raw-BG-XY-{}MHz.png'.format(f[i]))
#     plt.close()
    
#     fig,ax = plt.subplots(2,2,figsize=(10,6))
#     ax[0,0].plot(HA[i],XA[i])
#     ax[1,0].plot(HA[i],YA[i])
#     ax[0,0].set_title('High-to-Low - Hy drive\n{} MHz - raw data'.format(f[i]))
#     ax[0,0].set_ylabel('X')
#     ax[1,0].set_ylabel('Y')
#     ax[0,1].plot(bHA[i],bXA[i])
#     ax[1,1].plot(bHA[i],bYA[i])
#     ax[0,1].set_title('High-to-Low - Hy drive\n{} MHz - BG data'.format(f[i]))
    
#     plt.savefig(savepath2+'HyRF_HtL_Raw-BG-XY-{}MHz.png'.format(f[i]))
#     plt.close()
    
#     fig,ax = plt.subplots(2,2,figsize=(10,6))
#     ax[0,0].plot(HB[i],XB[i])
#     ax[1,0].plot(HB[i],YB[i])
#     ax[0,0].set_title('Low-to-High - Hy drive\n{} MHz - raw data'.format(f[i]))
#     ax[0,0].set_ylabel('X')
#     ax[1,0].set_ylabel('Y')
#     ax[0,1].plot(bHB[i],bXB[i])
#     ax[1,1].plot(bHB[i],bYB[i])
#     ax[0,1].set_title('Low-to-High - Hy drive\n{} MHz - BG data'.format(f[i]))
#     plt.close()
    
    # plt.savefig(savepath2+'HyRF_LtH_Raw-BG-XY-{}MHz.png'.format(f[i]))

# for i in rang


# fig,ax = plt.subplots(2,1,figsize=(8,5),gridspec_kw={'height_ratios':[2,1]})
# for i in range(len(f)):
#     ax[0].plot(HfA[i],dRA[i]/s[i],color=cols[i])
#     ax[1].plot(HfA[3-i],dPA[3-i]-p0[3-i],color=cols[3-i])
# ax[1].set_xlabel('$H_x^{DC}$ [kA/m]')
# ax[0].set_xticklabels([])
# ax[0].set_ylabel('Magnitude (scaled)')
# ax[1].set_ylabel('Phase [rad]')
# ax[0].set_title('High-to-Low field sweep - $H_y^{RF}$')
# ax[0].grid()
# ax[1].grid()
# ax[1].set_ylim(-0.5,0.5)
# ax[0].set_xlim(0,12)
# ax[1].set_xlim(0,12)
# # plt.savefig(savepath+'LowFieldLimit_HyRF_HtL_scaled.png')


# fig2,ax2 = plt.subplots(2,1,figsize=(8,5),gridspec_kw={'height_ratios':[2,1]})
# for i in range(len(f)):
#     ax2[0].plot(HfA[i+4],dRA[i+4]/s[i],color=cols[i])
#     ax2[1].plot(HfA[7-i],dPA[7-i]-p0[7-i],color=cols[3-i])
# ax2[1].set_xlabel('$H_x^{DC}$ [kA/m]')
# ax2[0].set_xticklabels([])
# ax2[0].set_ylabel('Magnitude (scaled)')
# ax2[1].set_ylabel('Phase [rad]')
# ax2[0].set_title('High-to-Low field sweep - $H_z^{RF}$')
# ax2[0].grid()
# ax2[1].grid()
# ax2[1].set_ylim(-0.3,0.5)
# ax2[0].set_xlim(0,12)
# ax2[1].set_xlim(0,12)
# # plt.savefig(savepath+'LowFieldLimit_HzRF_HtL_scaled.png')

# fig3,ax3 = plt.subplots(2,1,figsize=(8,5),gridspec_kw={'height_ratios':[2,1]})
# for i in range(len(f)):
#     ax3[0].plot(HfB[i],dRB[i]/s[i],color=cols[i])
#     ax3[1].plot(HfB[3-i],dPB[3-i]-p0[3-i],color=cols[3-i])
# ax3[1].set_xlabel('$H_x^{DC}$ [kA/m]')
# ax3[0].set_xticklabels([])
# ax3[0].set_ylabel('Magnitude (scaled)')
# ax3[1].set_ylabel('Phase [rad]')
# ax3[0].set_title('Low-to-High field sweep - $H_y^{RF}$')
# ax3[0].grid()
# ax3[1].grid()
# ax3[1].set_ylim(-0.5,0.7)
# ax3[0].set_xlim(0,12)
# ax3[1].set_xlim(0,12)
# # plt.savefig(savepath+'LowFieldLimit_HyRF_LtH_scaled.png')

# fig2,ax4 = plt.subplots(2,1,figsize=(8,5),gridspec_kw={'height_ratios':[2,1]})
# for i in range(len(f)):
#     ax4[0].plot(HfB[i+4],dRB[i+4]/s[i],color=cols[i])
#     ax4[1].plot(HfB[7-i],dPB[7-i]-p0[3-i],color=cols[3-i])
# ax4[1].set_xlabel('$H_x^{DC}$ [kA/m]')
# ax4[0].set_xticklabels([])
# ax4[0].set_ylabel('Magnitude (scaled)')
# ax4[1].set_ylabel('Phase [rad]')
# ax4[0].set_title('Low-to-High field sweep - $H_z^{RF}$')
# ax4[0].grid()
# ax4[1].grid()
# ax4[1].set_ylim(-0.2,0.5)
# ax4[0].set_xlim(0,12)
# ax4[1].set_xlim(0,12)
# plt.savefig(savepath+'LowFieldLimit_HzRF_LtH_scaled.png')






"""
plt.figure()
for i in range(1,len(f)):
    plt.plot(HfA[i][:57227],dRA[i][:57227]/s[i] - dRA[0]/s[0],color=cols[i])
plt.xlim(0,17)
"""
"""

# fig,ax = plt.subplots(2,1,figsize=(8,5),gridspec_kw={'height_ratios':[2,1]})
# for i in range(len(f)):
#     ax[0].plot(HfA[i],dRA[i]/s[i],color=cols[i])
#     ax[1].plot(HfA[3-i],np.unwrap(dPA[3-i]-p0[3-i]),color=cols[3-i])
# ax[1].set_xlabel('$H_x^{DC}$ [kA/m]')
# ax[0].set_xticklabels([])
# ax[0].set_ylabel('Magnitude (scaled)')
# ax[1].set_ylabel('Phase [rad]')
# ax[0].set_title('High-to-Low field sweep - $H_y^{RF}$')
# ax[0].grid()
# ax[1].grid()
# # ax[1].set_ylim(-0.5,0.5)
# ax[0].set_xlim(0,10)
# ax[1].set_xlim(0,10)


# fig2,ax2 = plt.subplots(2,1,figsize=(8,5),gridspec_kw={'height_ratios':[2,1]})
# for i in range(len(f)):
#     ax2[0].plot(HfA[7-i],dRA[7-i]/s[3-i],color=cols[i])
#     ax2[1].plot(HfA[7-i],np.unwrap(dPA[7-i]-p0[7-i]),color=cols[3-i])
# ax2[1].set_xlabel('$H_x^{DC}$ [kA/m]')
# ax2[0].set_xticklabels([])
# ax2[0].set_ylabel('Magnitude (scaled)')
# ax2[1].set_ylabel('Phase [rad]')
# ax2[0].set_title('High-to-Low field sweep - $H_z^{RF}$')
# ax2[0].grid()
# ax2[1].grid()
# # ax2[1].set_ylim(-0.3,0.5)
# ax2[0].set_xlim(0,10)
# ax2[1].set_xlim(0,10)
# ax2[0].set_ylim(0,0.85)

# fig3,ax3 = plt.subplots(2,1,figsize=(8,5),gridspec_kw={'height_ratios':[2,1]})
# for i in range(len(f)):
#     ax3[0].plot(HfB[i],dRB[i]/s[i],color=cols[i])
#     ax3[1].plot(HfB[3-i],np.unwrap(dPB[3-i]-p0[3-i]),color=cols[3-i])
# ax3[1].set_xlabel('$H_x^{DC}$ [kA/m]')
# ax3[0].set_xticklabels([])
# ax3[0].set_ylabel('Magnitude (scaled)')
# ax3[1].set_ylabel('Phase [rad]')
# ax3[0].set_title('Low-to-High field sweep - $H_y^{RF}$')
# ax3[0].grid()
# ax3[1].grid()
# # ax3[1].set_ylim(-0.5,0.7)
# ax3[0].set_xlim(0,10)
# ax3[1].set_xlim(0,10)


# fig2,ax4 = plt.subplots(2,1,figsize=(8,5),gridspec_kw={'height_ratios':[2,1]})
# for i in range(len(f)):
#     ax4[0].plot(HfB[7-i],dRB[7-i]/s[3-i],color=cols[3-i])
#     ax4[1].plot(HfB[7-i],np.unwrap(dPB[7-i]-p0[3-i]),color=cols[3-i])
# ax4[1].set_xlabel('$H_x^{DC}$ [kA/m]')
# ax4[0].set_xticklabels([])
# ax4[0].set_ylabel('Magnitude (scaled)')
# ax4[1].set_ylabel('Phase [rad]')
# ax4[0].set_title('Low-to-High field sweep - $H_z^{RF}$')
# ax4[0].grid()
# ax4[1].grid()
# # ax4[1].set_ylim(-0.2,0.5)
# ax4[0].set_xlim(0,10)
# ax4[1].set_xlim(0,10)
# ax4[0].set_ylim(0,0.75)


"""
