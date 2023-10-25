# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 17:01:51 2023

@author: Threadripper
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as patches

nm = np.arange(0,1261,2) #array of indices of magnetization arrays saved from mumax

XC = []
YC = []
XS = []
YS = []
X = []
Y = []

grid = 512

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

path = ''
tabpath = path+'table.txt'

tab = pd.read_csv(path+'table.txt',delimiter='\t') #load table from mumax

B = tab[tab.columns[4:7]].to_numpy().T*1e4 #convert field to G
m = tab[tab.columns[1:4]].to_numpy().T
t = tab[tab.columns[0]].to_numpy()

savepath = ''

Bx = np.unique(B[0])

DX = []
DY = []
x0 = []
y0 = []


for nmi in nm:

    mName = '{}'.format(nmi).zfill(6)
    m = np.loadtxt(path+'m{}.csv'.format(mName),delimiter=',') #load single saved m (converted from .ovf to .csv)
    
    mx = m[:grid,:] #identify mx, my, mz components
    my = m[grid:2*grid,:]
    mz = m[2*grid:,:]
    
    
    coreInds = np.where(abs(mz)>=0.6) #identify where the magnitude of mz is non-negligible 
    
    if bool(coreInds[0].size>0):
        maxX = coreInds[1].max()
        minX = coreInds[1].min()
        maxY = coreInds[0].max()
        minY = coreInds[0].min()
        
        dX = maxX - minX #find position of center of vortex core
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
    
    
DX = np.asarray(DX)
DY = np.asarray(DY)
x0 = np.asarray(x0)
y0 = np.asarray(y0)
    
x0_nm = x0*csx - grid*csx/2 #convert pixels to nm using simulation parameters
y0_nm = y0*csy - grid*csy/2


df = pd.DataFrame({'DC Field (x) (G)':Bx,
                    'x-positions (nm)':x0_nm,
                    'y-positions (nm)':y0_nm,
                    'x-positions (px)':x0,
                    'y-positions (px)':y0})

df.to_csv(path+'ExtractedVortexCorePositions.csv',index=False)
