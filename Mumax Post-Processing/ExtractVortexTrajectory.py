# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 17:01:51 2023

@author: Threadripper
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as patches


# nm = np.arange(0,50,1)
nm = np.arange(0,270,2)
# B = np.arange(-100,631,1)
# nm = np.arange(0,500,1)
XC = []
YC = []
XS = []
YS = []
X = []
Y = []

# for i in range(len(nm)):
    # mNamei = 
# i = np.where(B==150)[0][0]
grid = 512

# csx = 2
# csy = 2
csx = 1454/grid
csy = 1346/grid

vcx_nm = -173.23
vcy_nm = -162.99

vcx2_nm = -167.55
vcy2_nm = -99.9

vcx = (vcx_nm + grid*csx/2)/csx
vcy = (vcy_nm + grid*csy/2)/csy

vcx2 = (vcx2_nm + grid*csx/2)/csx
vcy2 = (vcy2_nm + grid*csy/2)/csy

# path = 'D:/Katryna/PyShroom/DC Hysteresis/Hysteresis_LtH_w20nmPinningAt50G_0p8Ms.out/'
# path = 'D:/Katryna/PyShroom/DC Hysteresis/Hysteresis_LtH_w20nmPinningAt50G_0p8Ms.out/TimeDomain_ydriven208MHz_saveDuringDrive.out/'
# tabpath = 'D:/Katryna/PyShroom/DC Hysteresis/Hysteresis_LtH_w20nmPinningAt50G_0p8Ms.out/TimeDomain_ydriven208MHz_saveDuringDrive.out/table.txt'
path = 'D:\Katryna\PyShroom\DC Hysteresis/Hysteresis_LtH_w20nmPinningAt50Gand28G_0p8Ms.out\TimeDomain_ydriven230MHz_saveDuringDrive_DoublePinning_v4_35-43G.out/'
path = 'D:/Katryna/PyShroom/DC Hysteresis/LowtoHighsim_states_noByBz_KRF_From0G.out/'
tabpath = path+'table.txt'
# tab 

# path = 'D:/Katryna/PyShroom/DC Hysteresis/Hysteresis_LtH_noPinning_take2.out/RingDown_SteppedFieldv3_100psSave_50nsRun_alpha0p01.out/'
tab = pd.read_csv(path+'table.txt',delimiter='\t')
B = tab[tab.columns[4:7]].to_numpy().T*1e4
m = tab[tab.columns[1:4]].to_numpy().T
t = tab[tab.columns[0]].to_numpy()
# Bx = Bx[109:] #for starting at 218
savepath = 'D:/Katryna/PyShroom/DC Hysteresis/Plots/'


# ti = np.where(t==0)[0]
# ti = np.concatenate((ti,[t.size]))

# n = []
# fft = []
# freqs = []

# for i in range(ti.size-1):
#     N = m[0][ti[i]:ti[i+1]].size
#     FT = abs(np.fft.fft(m[0][ti[i]:ti[i+1]]))
#     FTf = np.fft.fftfreq(N,d=10e-12)*1e-6
    
#     fft.append(FT[1:(N//2)])
#     freqs.append(FTf[1:(N//2)])
    # plt.figure()
    # plt.plot(FTf[1:(N//2)],FT[1:(N//2)],'.-',label='{:.0f} G'.format(B[0][ti[i]]),ms=5)

    # plt.xlim(0,1000)
    # plt.xlabel('Frequency (MHz)')
    # plt.ylabel('$m_x$ FFT')
    # # plt.legend()
    # plt.title('{:.0f} G'.format(B[0][ti[i]]))
    
    

# Bx = np.arange(32,48,2)
Bx = np.unique(B[0])
# Bx = np.arange(39,45,1)

# DX_B = []
# DY_B = []
# x0_B = []
# y0_B = []


# for i in range(Bx.size):
# DX = np.zeros_like(Bx)
# DY = np.zeros_like(Bx)
# x0 = np.zeros_like(Bx)
# y0 = np.zeros_like(Bx)
DX = []
DY = []
x0 = []
y0 = []


# n_p = 31
# nm = np.arange(n_p*i,n_p*(i+1),1)
nm = np.arange(0,1261,2)
for nmi in nm:
# nmi = 250
    # print(nmi)
    mName = '{}'.format(nmi).zfill(6)
#     # path = 'D:/Katryna/PyShroom/Joe old scripts/HightoLowsim_states_noByBz_2.out/Joe_PyShroom_yRFDrive_208MHz_0p75G_Hy_full_50nsRun.out/m{}.csv'.format(mName)
#     # path = 'D:/Katryna/PyShroom/Joe old scripts/HightoLowsim_states_noByBz_2.out/GyroHysteresis_SteppedDC_TimeDomainRF_208MHz_Hy_2Gsteps.out/m{}.csv'.format(mName)
#     # path = 'D:/Katryna/PyShroom/Joe old scripts/HightoLowsim_states_noByBz_2.out/GyroHysteresis_SteppedDC_TimeDomainRF_0p75G_208MHz_Hy_2Gsteps_with45nmPinning.out/Joe_PyShroom_spec_1ns_pulse_doubleVortex_0G_Hy_full_100nsRun.out/m{}.csv'.format(mName)
#     # path = 'D:/Katryna/PyShroom/DC Hysteresis/Hysteresis_LtH_w20nmPinningAt50G_0p8Ms.out/RingDown_500psSave_100nsRun.out/m{}.csv'.format(mName)
    m = np.loadtxt(path+'m{}.csv'.format(mName),delimiter=',')
    
    mx = m[:grid,:]
    my = m[grid:2*grid,:]
    mz = m[2*grid:,:]
    # plt.figure()
    # plt.imshow(mz)
    
    coreInds = np.where(abs(mz)>=0.6)
    
    if bool(coreInds[0].size>0):
        maxX = coreInds[1].max()
        minX = coreInds[1].min()
        maxY = coreInds[0].max()
        minY = coreInds[0].min()
        
        dX = maxX - minX
        dY = maxY - minY
        
        meanx = dX/2 + minX
        meany = dY/2 + minY
        
        DX.append(dX)
        DY.append(dY)
        x0.append(meanx)
        y0.append(meany)
    else:
        DX.append(np.nan)
        DY.append(np.nan)
        x0.append(np.nan)
        y0.append(np.nan)
    
#     DX_B.append(DX)
#     DY_B.append(DY)
#     x0_B.append(x0)
#     y0_B.append(y0)
   
#     plt.figure()
#     plt.plot(x0,y0)
    
DX = np.asarray(DX)
DY = np.asarray(DY)
x0 = np.asarray(x0)
y0 = np.asarray(y0)
    
# DX_B = np.asarray(DX_B)
# DY_B = np.asarray(DY_B)
# x0_B = np.asarray(x0_B)
# y0_B = np.asarray(y0_B)


x0_nm = x0*csx - grid*csx/2
y0_nm = y0*csy - grid*csy/2


df = pd.DataFrame({'DC Field (x) (G)':Bx,
                    'x-positions (nm)':x0_nm,
                    'y-positions (nm)':y0_nm,
                    'x-positions (px)':x0,
                    'y-positions (px)':y0})

df.to_csv(path+'ExtractedVortexCorePositions.csv',index=False)

# P1 = patches.Circle((vcx_nm,vcy_nm),10,color='gray',alpha=0.7)
# P2 = patches.Circle((vcx2_nm,vcy2_nm),10,color='gray',alpha=0.7)
# fig,ax = plt.subplots(1,1,figsize=(6,6))


# dx = []
# dy = []
# ar = []

# for i in range(Bx.size):
#     ax.plot(x0_nm[i],y0_nm[i],'.-',ms=2,lw=1,label='{:.0f} G'.format(Bx[i]))
#     dx.append(x0_nm[i].max()-x0_nm[i].min())
#     dy.append(y0_nm[i].max()-y0_nm[i].min())
#     # ar.append((y0_nm[i].max()-y0_nm[i].min())/(x0_nm[i].max()-x0_nm[i].min()))
# ax.add_patch(P1)
# ax.add_patch(P2)
# plt.axis('equal')
# ax.legend()

# plt.ylabel('y-position (nm)')
# plt.xlabel('x-position (nm)')

# # plt.savefig(savepath+'driven-vortexCore-trajectories-vs-field-w-pinning.png',format='png',dpi=400,bbox_inches='tight')

# B = np.reshape(np.repeat(Bx,n_p),x0_nm.shape)

# xyB = np.array([B,x0_nm,y0_nm])

# np.save(path+'PinningRegion-VortexCoreTrajectories_driven',xyB)

#     # plt.figure()
#     # plt.plot(DX,DY)
    
#     # max_xC = np.where(abs(mz[:,:500])==abs(mz[:,:500]).max())[1][0]
#     # max_yC = np.where(abs(mz[:,:500])==abs(mz[:,:500]).max())[0][0]
#     # max_x = np.where(abs(mz)==abs(mz).max())[1][0]
#     # max_y = np.where(abs(mz)==abs(mz).max())[0][0]
    
#     max_x = np.where(abs(mz[:,:315])==abs(mz[:,:315]).max())[1][0] #taking the x-index up to 315 only captures the cap of the mushroom; not the stalk (where a secondary vortex core likes to show up )
#     max_y = np.where(abs(mz[:,:315])==abs(mz[:,:315]).max())[0][0]
    
#     # VL = np.where(abs(mz[:,:int(2*grid/3)])>0.4)
#     # max_y,max_x = VL[0].mean(),VL[1].mean()   
    
#     # max_xS = np.where(abs(mz[:,500:])==abs(mz[:,500:]).max())[1][0]
#     # max_yS = np.where(abs(mz[:,500:])==abs(mz[:,500:]).max())[0][0]
    
#     X.append(max_x)
#     Y.append(max_y)
#     # XC.append(max_xC)
#     # YC.append(max_yC)
#     # XS.append(max_xS)
#     # YS.append(max_yS)

# X = np.asarray(X)
# Y = np.asarray(Y)

# x = csx*X - grid*csx/2
# y = csy*Y - grid*csy/2

# # XC = np.asarray(XC)
# # YC = np.asarray(YC)
# # YS = np.asarray(YS)
# # XS = np.asarray(XS)
# plt.figure()
# plt.plot(x,y,'.-')




# plt.title('Cap')
# plt.figure()
# plt.plot(XS,YS)
# plt.title('Stalk')

# plt.imshow(mx[::-1])

# x = X*2 - 800
# y = Y*2 - 800



# X.append(max_x*2 - 800)
# Y.append(max_y*2 - 800)
# X = max_x*csx - grid*csx/2
# Y = (max_y)*csy -  grid*csy/2 
# plt.plot(X,Y)

# print('vortex core ({0} G) @ ({1:.2f},{2:.2f}) nm (mumax coord.)'.format(B[i],X,Y))
# print('vortex core ({0} G) @ ({1:.2f},{2:.2f}) px'.format(B[i],max_x,max_y))
# 
# plt.figure()
# plt.plot(X)

# plt.figure()
# plt.plot(Y)
# X = np.asarray(X)
# Y = np.asarray(Y)
# print('dX = {}'.format(X.max()-X.min()))
# print('dY = {}'.format(Y.max()-Y.min()))

# x1 = 300
# x2 = 335
# y1 = 540
# y2 = 570

# plt.imshow(abs(mz[y1:y2,x1:x2]))
# plt.colorbar()
# plt.plot(12,14.5,'rx')
# plt.plot(21,14.5,'rx')

# mi =np.array([mx[max_y,max_x],my[max_y,max_x],mz[max_y,max_x]])
# N = np.array([0.409,0.409,0.182])

# Ms = 767e3

# def Ed(Ms):
#     return 0.5*4*np.pi*1e-7*Ms**2*np.multiply(np.multiply(N,mi),mi)
# # print(Ed(767e3))
# # print(Ed(767e3*0.9))
# Ms0 = 767e3
# MsP = 0.8*Ms0

# print((Ed(Ms0)-Ed(MsP))/Ed(Ms0))