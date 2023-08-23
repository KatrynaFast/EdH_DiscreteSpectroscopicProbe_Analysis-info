# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 14:55:00 2022 (Katryna)

27feb2022 mods:
    adding EdH torque calculation based on dm/dt (equilibrium model)
        - put -ve sign on electron charge

03mar2022 mods:
    annotation fixed (Katryna)
    add inspection of ratio of EdH to precession torque 
    explore dependence on damping and on prolate cylinder aspect ratios (off resonance, still)
    
04mar2022:
    explore dependence of torque ratios on oblate cylinder aspect ratios
    
11mar2022:
    adding YIG anisotropy term from Michael

24mar2022: Katryna adds precession torque calculation with full effective field

21apr2022: Katryna adds field sweep

22apr2022: jig for Py
    disk thickness, diameter: 40 nm, 1.5 um
    DC field along x, RF drives along y, z (separately, and together)
    
20may2022: adding six columns to the saved data frame, to include the anisotropy field
and the anisotropy cross-product torque (a.k.a the mxB 'reaction' torque, mxB_anis = mx(B_eff - B_applied))

27may2022: adding Ku sweep option.

How to use: 
    - Lines 455 - 540 in the version of the script (accurate as of January 5, 2023) 
    are where the simulation settings can be adjusted. Above these lines is 
    definitions of the functions used throughout the script. 
    - Output saves once computation is complete (make sure path and filename are
    appropriately set before running script every time)
    - Plots that are drawn and saved can be adjusted at the bottom of the script 
    (starting at line  891)
    

@author: AWGPC1
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.fft import fft, fftfreq
import scipy.special as sp
import pandas as pd
import time
import os

startTime = time.time()
# plt.close('all')
plt.rcParams['axes.formatter.useoffset']=False
plt.rcParams['axes.labelsize']=12
plt.rcParams['axes.grid']=True
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

print('Did you change the file path? (if not, abort and do that!)') #reminder to adjust savename and savepath to avoid overwriting or putting saved data in the wrong spot. 

""" BEGIN DEFINITIONS """ ####################################################

def mag(the,phi): #magnetization vector in Cartesian coordinates
    mx = np.sin(the)*np.cos(phi)
    my = np.sin(the)*np.sin(phi)
    mz = np.cos(the)
    return [mx,my,mz]  

def H_demag(the,phi): #demag field vector
    return -Ms*np.dot(N,mag(the,phi))

def H_anis_uniax(the,phi,t): 
    return Ku(t)*np.dot(K,mag(the,phi))*2/(mu0*Ms)

def H_anis_cubic(the,phi):
    m = mag(the,phi)
    Hacx = 2*m[0]*(Kc1*(m[1]**2 + m[2]**2) + Kc2*m[1]*m[2])
    Hacy = 2*m[1]*(Kc1*(m[0]**2 + m[2]**2) + Kc2*m[0]*m[2])
    Hacz = 2*m[2]*(Kc1*(m[0]**2 + m[1]**2) + Kc2*m[0]*m[1])
    Hac = np.array([Hacx,Hacy,Hacz])
    return Hac

def H_sweep(t,H1,H2,hystMode):
    """
    Parameters
    ----------
    t : time -> must be an array for the hysteresis fields to calculate correctly
    H1 : initial field
    H2 : final field
    hystMode : type of hysteresis loop to perform:
        0: start at H1,end at H2, linear time-dependent field switching that persists through the entire simulation duration
        1: start at H1, linear switch to H2, reverse direction for linear switch to H1
        2: start at 0, adjust field to H1. Reverse direction, adjust field to H2
        3: start at 0, adjust field to H1. Reverse direction, adjust field to H2. Reverse direction, adjust field to H1. 

    Returns
    -------
    H : TYPE
        DESCRIPTION.

    """
    if H1==H2:
        H = H1
    elif hystMode==0:
        slope = (H2 - H1)/(t2-t1)
        H = H1 + slope*t
    elif hystMode==1:
        slope = 2*(H2 - H1)/(t2-t1)
        th = np.heaviside(t-t2/2,1)
        H = (H1 + slope*t)*(1-th) + (H2 - slope*(t-t2/2))*th
    elif hystMode==2:
        slope = (2*H1 - H2)/(t2-t1)
        t_H1 = t[np.where(abs(t*slope)>=abs(H1))[0][0]] #the time at which the field equals H1x
        th1 = np.heaviside(t - t_H1,1)
        H = slope*t*(1-th1) + (H1 - slope*(t-t_H1))*th1
    elif hystMode==3:
        slope = (3*H1 - 2*H2)/(t2-t1)
        t_H1 = t[np.where(abs(t*slope)>=abs(H1))[0][0]]
        if H1>0:
            t_H2 = t[np.where(H1 - slope*(t - t_H1)<=H2)[0][0]]
        elif H1<0:
            t_H2 = t[np.where(H1 - slope*(t - t_H1)>=H2)[0][0]]
        th1 = np.heaviside(t - t_H1,1)
        th2 = np.heaviside(t - t_H2,1)
        H = slope*t*(1-th1)+(H1-slope*(t-t_H1))*th1*(1-th2) + (H2 + slope*(t-t_H2))*th2
    return H

def freqSweep_v3_freq(t,f1,f2): #the frequency for a frequency sweep 
    """ 
    To define a frequency sweep, one needs to use the integral of the function 
    f(t) as the argument for the sine term (as done in function freqSweep_v4) 
    rather than f(t). 
    
    This function outputs the frequency values of a linear sweep given
    t: time array
    f1: starting frequency
    f2: ending frequency
    """
    if f1==f2:
        f = np.ones_like(t)*f1
    else:
        s = (f2 - f1)/(t2-t1)
        f = (s*t + f1)
    return f

def freqSweep_v4(t,f1,f2,p):
    if f1==f2:
        h = np.sin(2*np.pi*f1*t + p)
    else:
        s = (f2-f1)/(t2-t1)
        arg = 2*np.pi*(0.5*s*t**2 + f1*t) + p
        h = np.sin(arg)
    return h


def H_app(t): #applied field vector
    """ 
    Calculate the applied field as the combination of the DC field component
    and the RF field component. The RF field component can be turned off
    at an arbitrary time within the simulation through the variables tox,toy,toz
    """
    Hax = H_sweep(t,HxDCi,HxDCf,HxDCM) + (1 - np.heaviside(t - tox,1))*H_sweep(t,HxRFi,HxRFf,HxRFM)*freqSweep_v4(t,fx1RF,fx2RF,pxRF)
    Hay = H_sweep(t,HyDCi,HyDCf,HyDCM) + (1 - np.heaviside(t - toy,1))*H_sweep(t,HyRFi,HyRFf,HyRFM)*freqSweep_v4(t,fy1RF,fy2RF,pyRF)
    Haz = H_sweep(t,HzDCi,HzDCf,HzDCM) + (1 - np.heaviside(t - toz,1))*H_sweep(t,HzRFi,HzRFf,HzRFM)*freqSweep_v4(t,fz1RF,fz2RF,pzRF)
    return (Hax,Hay,Haz)

def H_eff(t,the,phi): #define the effective field
    #based on Chapter 1 Section 2.2 in 2002 book: Spin Dynamics in Confined Magnetic Structures 1
    HD = H_demag(the,phi)
    HK = H_anis_uniax(the,phi,t) + H_anis_cubic(the,phi)
    HA = H_app(t)
    
    Heff = HA + HD + HK
    return Heff

def H_anis(t,the,phi): #define the anisotropy field
    HD = H_demag(the,phi)
    HK = H_anis_uniax(the,phi,t) + H_anis_cubic(the,phi)
        
    Hanis = HD + HK 
    return Hanis

def alpha(t): 
    """
    function to sweep Gilbert damping constant alpha. 
    If alpha is *not* defined to be swept (alpha2==0), the variable remains a 
    float rather than becoming an array to reduce computation time
    """
    if alpha2==0:
        alpha = alpha1
    else:
        s = (alpha2-alpha1)/(t2-t1)
        alpha = s*t + alpha1
    return alpha

def Ku(t):
    """
    function to sweep uniaxial anisotropy constant. 
    If alpha is *not* defined to be swept (Ku_i == Ku_f), the variable Ku remains a 
    float rather than becoming an array to reduce computation time
    """
    if Ku_f==Ku_i:
        Ku = Ku_i
    else:
        s = (Ku_f - Ku_i)/(t2-t1)
        Ku = s*t + Ku_i
    return Ku

def wD(the,phi): #demag energy density #units J/m^3
    AD = 0.5*mu0*Ms**2
    wDx = Nx*(np.sin(the)*np.cos(phi))**2
    wDy = Ny*(np.sin(the)*np.sin(phi))**2
    wDz = Nz*(np.cos(the))**2
    wd = AD*(wDx + wDy + wDz)
    return wd

def dthewD(the,phi): #demag energy density theta derivative
    AtD = mu0*Ms**2*np.sin(the)*np.cos(the)
    dtwDx = Nx*(np.cos(phi))**2
    dtwDy = Ny*(np.sin(phi))**2
    dtwDz = -Nz
    dtwD = AtD*(dtwDx + dtwDy + dtwDz)
    return dtwD

def dphiwD(the,phi): #demag energy density phi derivative
    ApD = -mu0*np.cos(phi)*np.sin(phi)*(Ms*np.sin(the))**2
    dpwDx = Nx
    dpwDy = -Ny
    dpwDz = 0
    dpwD = ApD*(dpwDx + dpwDy + dpwDz)
    return dpwD

def wZ(the,phi,t): #Zeeman energy density #units: J/m^3
    AZ = -mu0*Ms
    wZx = H_app(t)[0]*np.sin(the)*np.cos(phi)
    wZy = H_app(t)[1]*np.sin(the)*np.sin(phi)
    wZz = H_app(t)[2]*np.cos(the)
    wz = AZ*(wZx + wZy + wZz)
    return wz

def dthewZ(the,phi,t): #zeeman energy density theta derivative
    AtZ = -mu0*Ms
    dtwZx = H_app(t)[0]*np.cos(the)*np.cos(phi)
    dtwZy = H_app(t)[1]*np.cos(the)*np.sin(phi)
    dtwZz = -H_app(t)[2]*np.sin(the)
    dtwZ = AtZ*(dtwZx + dtwZy + dtwZz)
    return dtwZ 

def dphiwZ(the,phi,t): #zeeman energy density phi derivative
    ApZ = mu0*Ms*np.sin(the)
    dpwZx = H_app(t)[0]*np.sin(phi)
    dpwZy = -H_app(t)[1]*np.cos(phi)
    dpwZz = 0
    dpwZ = ApZ*(dpwZx + dpwZy + dpwZz)
    return dpwZ

def wA_uniax(the,phi,t): #uniaxial anisotropy energy density #units J/m^3
    wAx = Kx*(np.sin(the)*np.cos(phi))**2
    wAy = Ky*(np.sin(the)*np.sin(phi))**2
    wAz = Kz*(np.cos(the))**2
    wa = -Ku(t)*(wAx + wAy + wAz)
    return wa

def dthewA_uniax(the,phi,t): #uniaxial anisotropy energy density theta derivative
    AtA = -2*np.sin(the)*np.cos(the)
    dtwAx = Kx*Ku(t)*(np.cos(phi))**2
    dtwAy = Ky*Ku(t)*(np.sin(phi))**2
    dtwAz = -Kz*Ku(t)
    dtwA = AtA*(dtwAx + dtwAy + dtwAz)
    return dtwA

def dphiwA_uniax(the,phi,t): #uniaxial anisotropy energy density phi derivative
    ApA = 2*np.cos(phi)*np.sin(phi)*(np.sin(the))**2
    dpwAx = Kx*Ku(t)
    dpwAy = -Ky*Ku(t)
    dpwAz = 0
    dpwA = ApA*(dpwAx + dpwAy + dpwAz)
    return dpwA
    
def wA_cub(the,phi): #cubic anisotropy energy density #units J/m^3
    wAc1 = (Kc1 + Kc2*np.cos(the)**2)*(np.sin(the)**2*np.cos(phi)*np.sin(phi))**2
    wAc2 = Kc1*(np.cos(the)*np.sin(the))**2
    wAc = wAc1 + wAc2
    return wAc

def dthewA_cub(the,phi): #cubic anisotropy energy density theta derivative
    A = 2*np.sin(the)*np.cos(the)
    dtwA1 = Kc1*(2*(np.sin(the)*np.cos(phi)*np.sin(phi))**2 + np.cos(the)**2 - np.sin(the)**2)
    dtwA2 = Kc2*(np.sin(the)*np.sin(phi)*np.cos(phi))**2*(2*np.cos(the)**2 - np.sin(the)**2)
    dtwAc = A*(dtwA1 + dtwA2)
    return dtwAc

def dphiwA_cub(the,phi): #cubic anisotropy energy density phi derivative
    dpwA1 = Kc1 + Kc2*np.cos(the)**2
    dpwA2 = np.sin(the)**4
    dpwA3 = 2*np.cos(phi)*np.sin(phi)*(np.cos(phi)**2 - np.sin(phi)**2)
    dpwAc = dpwA1*dpwA2*dpwA3
    return dpwAc

def wA(the,phi,t): #total anisotropy energy density (sum of uniaxial and cubic) # units J/m^3
    wAu = wA_uniax(the,phi,t)
    wAc = wA_cub(the,phi)
    wA_tot = wAu + wAc
    return wA_tot
    
def dthewA(the,phi,t): #total anisotropy energy density theta derivative
    dtwAu = dthewA_uniax(the,phi,t)
    dtwAc = dthewA_cub(the,phi)
    dtwA = dtwAu + dtwAc
    return dtwA

def dphiwA(the,phi,t): #total anisotropy energy density phi derivative
    dpwAu = dphiwA_uniax(the,phi,t)
    dpwAc = dphiwA_cub(the,phi)
    dpwA = dpwAu + dpwAc
    return dpwA

def w(the,phi,t): #total energy density
    ww = wD(the,phi) + wZ(the,phi,t) + wA(the,phi,t)
    return ww

def dthew(the,phi,t): #total energy density theta derivative
    dtw = dthewD(the,phi) + dthewZ(the,phi,t) + dthewA(the,phi,t)
    return dtw

def dphiw(the,phi,t): #total energy density phi derivative
    dpw = dphiwD(the,phi) + dphiwZ(the,phi,t) + dphiwA(the,phi,t)
    return dpw

def ode(t,y): #ode to be solved
    the = y[0]
    phi = y[1]
    
    A_ode = -gamma/(Ms*(1+alpha(t)**2))
    tdot = A_ode*(alpha(t)*dthew(the,phi,t) + (1/np.sin(the))*dphiw(the,phi,t))
    pdot = (A_ode/np.sin(the))*(-dthew(the,phi,t) + (alpha(t)/np.sin(the))*dphiw(the,phi,t))
    return (tdot,pdot)

### Define functions for the Jacobian matrix
def dthedthew(the,phi,t): #d/dthe (dw/dthe)
    A1 = (np.cos(the)**2 - np.sin(the)**2)
    a = -2*A1*(Kx*np.cos(phi)**2 + Ky*np.sin(phi)**2 - Kz)*Ku(t)
    b = A1*mu0*Ms**2*(Nx*np.cos(phi)**2 + Ny*np.sin(phi)**2 - Nz)
    c = mu0*Ms*(H_app(t)[0]*np.sin(the)*np.cos(phi) + H_app(t)[1]*np.sin(the)*np.sin(phi) + H_app(t)[2]*np.cos(the))
    
    #cubic anisotropy terms (added later, also very long):
    wcA = (np.sin(the)*np.sin(phi)*np.cos(phi))**2
    wc1 = (np.cos(the)**2 - np.sin(the)**2)*(4*Kc1 + 12*Kc2*np.cos(the)**2 - Kc2*np.sin(the)**2)
    wc2 = 8*Kc1*np.cos(the)**2
    wc3 = 2*Kc1*(np.cos(the)**2 - np.sin(the)**2)**2
    
    wc = wcA*(wc1 + wc2) + wc3
    #total derivative
    tot = a + b + c + wc
    return tot

def dphidthew(the,phi,t):#d/dphi (dw/dthe)
    A1 = 2*np.sin(the)*np.cos(the)*np.sin(phi)*np.cos(phi)
    a = A1*2*(Kx - Ky)*Ku(t)
    b = -A1*mu0*Ms**2*(Nx - Ny)
    c = mu0*Ms*np.cos(the)*(H_app(t)[0]*np.sin(phi) - H_app(t)[1]*np.cos(phi))
    
    #cubic anisotropy terms
    wc1 = 4*np.sin(the)**3*np.cos(the)*np.sin(phi)*np.cos(phi)*(np.cos(phi)**2 - np.sin(phi)**2)
    wc2 = 2*Kc1 + Kc2*(2*np.cos(the)**2 - np.sin(the)**2)
    wc = wc1*wc2
    
    #total derivative
    tot = a + b + c + wc
    return tot

def dphidphiw(the,phi,t): #d/dphi (dw/dphi)
    A1 = np.sin(the)**2*(np.cos(phi)**2 - np.sin(phi)**2)
    a = 2*A1*(Kx - Ky)*Ku(t)
    b = -A1*mu0*Ms**2*(Nx - Ny) 
    c = mu0*Ms*np.sin(the)*(H_app(t)[0]*np.cos(phi) + H_app(t)[1]*np.sin(phi))
    
    #cubic anisotropy terms
    wc1 = 2*np.sin(the)**4*(Kc1 + Kc2*np.cos(the)**2)
    wc2 = np.cos(phi)**4 - 6*(np.sin(phi)*np.cos(phi))**2 + np.sin(the)**4
    wc = wc1*wc2
    
    #total derivative
    tot = A1*(a + b) + c + wc
    return tot

def dthedphiw(the,phi,t): #d/dthe (dw/dphi)
    A1 = 2*np.sin(the)*np.cos(the)*np.sin(phi)*np.cos(phi)
    a = 2*A1*(Kx - Ky)*Ku(t)
    b = -A1*mu0*Ms**2*(Nx - Ny)
    c = mu0*Ms*np.cos(the)*(H_app(t)[0]*np.sin(phi) - H_app(t)[1]*np.cos(phi))
    
    #cubic aniostropy terms
    wc1 = 4*np.sin(the)**3*np.cos(the)*np.sin(phi)*np.cos(phi)*(np.cos(phi)**2 - np.sin(phi)**2)
    wc2 = Kc2*(2*np.cos(the)**2 - np.sin(the)**2) + 2*Kc1
    wc = wc1*wc2
    
    #total derivative
    tot = A1*(a + b) + c + wc
    return tot

def jacobian(t,y): #define jacobian matrix for solve_ivp
    the = y[0]
    phi = y[1]
    
    A = -gamma/(Ms*(1 + alpha(t)**2))
    
    a = alpha(t)*dthedthew(the,phi,t) - (np.cos(the)/(np.sin(the)**2))*dphiw(the,phi,t) + (1/np.sin(the))*dthedphiw(the,phi,t)
    b = alpha(t)*dphidthew(the,phi,t) + (1/np.sin(the))*dphidphiw(the,phi,t)
    c = (np.cos(the)/(np.sin(the)**2))*dthew(the,phi,t) - (1/np.sin(the))*dthedthew(the,phi,t) - (2*alpha(t)*np.cos(the)/(np.sin(the)**3))*dphiw(the,phi,t) + (alpha(t)/np.sin(the)**2)*dthedphiw(the,phi,t)
    d = (-1/np.sin(the))*dphidthew(the,phi,t) + (alpha(t)/np.sin(the)**2)*dphidphiw(the,phi,t)
    
    
    j = np.zeros(shape = (2,2))
    j[0,0] = A*a
    j[0,1] = A*b
    j[1,0] = A*c
    j[1,1] = A*d
    
    return j

"""Mark's code for calculating demag factor"""
#For cylinder:
def kappa(a):
    #here, a is the height/diameter aspect ratio
     return 1/np.sqrt(1 + a**2)

def Nzocyl(a): #Nz for oblate cylinder
    return 1 - (4/(3*np.pi*a))*(np.sqrt(1 + a**2)*(a**2*sp.ellipk(kappa(a)**2) + (1 - a**2)*sp.ellipe(kappa(a)**2)) - 1)

def Nzpcyl(a): #Nz for prolate cylinder
    return 1 - (4*a/(3*np.pi))*(np.sqrt(1 + (1/a)**2)*((1/a)**2*sp.ellipk(kappa(1/a)**2) + (1 - (1/a)**2)*sp.ellipe(kappa(1/a)**2)) - 1)

#for ellipsoid:
def eps_sq(x):
    #here, x is the aspect ratio (between 0 and 1, for both oblate and prolate)
     return (1-x**2)
 
def Nzoell(a): #Nz for oblate ellipsoid of rotation
    return sp.hyp2f1(1, 1, 5/2, eps_sq(a))/3.0

def Nzpell(a): #Nz for prolate ellipsoid of rotation
    return (1 - eps_sq(a))*sp.hyp2f1(1, 3/2, 5/2, eps_sq(a))/3.0

#Radial demag factor
def Nr_ce(Nz): #Nr (radial demag factor) for a given Nz as calculated above
    return (1 - Nz)/2

""" END DEFINITIONS """ #######################################################

###############################################################################
""""""""""""""""""""" RELEVANT WORKING CODE STARTS HERE! """""""""""""""""""""
#constants:
mu0 = 4*np.pi*1e-7 #N/A^2
g = 2.
m_e = 9.1083837e-31 #kg #electron mass
e = -1.60217663e-19 #C #electron charge
# gamma = abs(e)*g/(2*m_e) #rad/sT -- gamma for g given above, (Py value, based on Wallis 2006 g' value)
gamma =  1.76086e11 #rad/sT
# gprime = 1.83
gprime = 4-g

#material constants 
Ms = 767e3 #A/m - nominal Py magnetization (thin film value, used in mumax; 763 in many sibling scripts)
alpha1= 0.01 #damping constant
alpha2 = 0 #if no sweep of alpha is desired, set alpha2=0 and the value of alpha1 will be constant. 
h = 4e-8 #height of disk  #oblate geometry!
# D = 1.3e-6 #diameter of disk 
D = 0.47e-6 #equivalent spherical diameter for Py island, as calculated by Mark
# V = h*(np.pi*(D/2)**2) #volume of the disk, m^3
V = (4*np.pi/3)*(D/2)**3

# Na = Nzocyl(h/D) #demag factors for the height/diameter ratio given (oblate)
# Na = Nzpcyl(D/h) #demag factors for the height/diameter ratio given (prolate; send parameter D/h for function as-written)
# Na = Nzoell(h/D)
# Na = Nzpell(D/h)
# Nr = Nr_ce(Na)

# prefactor to convert dm/dt (from normalized LLG m) to EdH torque, per equilibrium model
A_EdH = -2*m_e*Ms*V/(e*g) #use g in EdH calculation so I can post-process it later to account for g/g' scaling factor between EdH and EdH from LLG

# Initial conditions
dp = 0.000001
the0 = np.pi/2  # pi/2 for x-y plane; 2pi for along z, but precisely = problematic;
phi0 = np.pi*0  # 0 along x, pi/2 along y
#phi0 = 0.048 #initial 50 mrad tilt along y, to look at ringdown

# simulation parameters
Nx,Ny,Nz = [0.0423,0.046,0.9117] #est. for PyShroom from Joe's mumax with Ms = 767 kA/m
# Nx,Ny,Nz = [Nr,Nr,Na] #demag factors as calculated with Mark's functions above
# Nx,Ny,Nz = [1/3,1/3,1/3] #demag factors for a sphere
# 
HxDCi,HyDCi,HzDCi = np.asarray([65e3,0,0])  #A/m -> H1 for DC field calculation
HxDCf,HyDCf,HzDCf = np.asarray([5e3,0,0])   #A/m -> H2 for DC field calculation
HxDCM,HyDCM,HzDCM = np.asarray([0,0,0])     #mode of hysteresis to use for the  DC field (see H_sweep function for what each mode represents)

HxRFi,HyRFi,HzRFi = np.asarray([0,0,10])    #A/m - initial RF field amplitudes
HxRFf,HyRFf,HzRFf = np.asarray([0,0,10])    #A/m - final RF field amplitudes
HxRFM,HyRFM,HzRFM = np.asarray([0,0,0])     #"mode" of RF field sweep (see H_sweep function for what each mode represents)

fRF = 3e6
fx1RF,fy1RF,fz1RF = np.asarray([0,fRF,fRF])  #starting frequency in frequency sweep
fx2RF,fy2RF,fz2RF = np.asarray([0,fRF,fRF]) #ending frequency in frequency sweep
pxRF,pyRF,pzRF = np.asarray([0,0,-np.pi/4])

Ku_i = 0   #uniaxial anisotropy constant
Ku_f = 0    # for  Ku sweep; if Ku_f==Ku_i, no sweep is performed, the value is simply Ku_i
Kc1 = 0     #cubic anisotropy constant (-610 for YIG)
Kc2 = 0     #cubic anisotropy constant (-26 for YIG)
Kx,Ky,Kz = [1,0,0] #define easy axis of anisotropy (for uniaxial anisotropy function) (negative value indicates hard axis)

t1 = 0 #s # starting time of simulation 
t2 = 150e-6 #s # ending time of simulation
dt = 1500e-12 #s #time step 
n = int((t2-t1)/dt)
t_vals = np.arange(t1,t2,dt) #time array to be used in solution

# Set times at which to turn off the RF fields 
tox = t2 + dt#1000e-9
toy = t2 + dt
toz = t2 + dt

T1 = t1*1e9 #to make plot limits easier, re-express t1,t2 in ns
T2 = t2*1e9 #ns


N = [[Nx,0,0],[0,Ny,0],[0,0,Nz]] #demag factor tensor
K = [[Kx,0,0],[0,Ky,0],[0,0,Kz]] #uniaxial energy direction definition
U = [[1,0,0],[0,1,0],[0,0,1]] #unit tensor for adding M to the effective field


""" Save settings"""

# savepath = '/Volumes/GoogleDrive/My Drive/Research/Torque Mixing - Kittel Mode/2022-12-05 Script Cleanup Testing/'
# savepath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Research/Torque Mixing - Kittel Mode/2023-06-27 High Field Sweeps/'
savepath = '/Users/katrynafast/Library/CloudStorage/GoogleDrive-kfast@ualberta.ca/My Drive/Conferences and Papers/NRC paper/Simulation Summary/Macrospin High Field z Phase offset/'
nametosave = 'PyShroom_zDrive_3MHz_150us_dt1500ps_ph-45'
savename = nametosave+'.csv'

print('demag factors: Nx = {:.3f}, Ny = {:.3f}, Nz = {:.3f}'.format(Nx,Ny,Nz))

"""""" 
sol = solve_ivp(ode,[t1,t2],[the0,phi0],method='Radau',t_eval=t_vals,dense_output=False,rtol=1e-9,atol=1e-12,jac=jacobian)
t = sol.t*1e9 #time in ns
m = np.array(mag(sol.y[0],sol.y[1])) #magnetization vector from solution
B = mu0*np.array(H_app(sol.t)) #B field from solution
H = np.array(H_app(sol.t))  #applied H field
Heff = H_eff(sol.t,sol.y[0],sol.y[1]) #effective field from solution
Hanis = H_anis(sol.t,sol.y[0],sol.y[1]) #effective field terms responsible for the 'reaction' cross-product torque
f_app = (freqSweep_v3_freq(sol.t,fx1RF,fx2RF),freqSweep_v3_freq(sol.t,fy1RF,fy2RF),freqSweep_v3_freq(sol.t,fz1RF,fz2RF)) #frequency throughout simulation

""""""

""" Calculate torques """
dtm = np.gradient(m,dt,axis=1) #time-derivative of the magnetization to be used to calculate torques.
tau_prec = Ms*V*np.cross(m.transpose(),B.transpose()) #cross-product torque (Nm)
tau_prec_LLG = Ms*V*mu0*np.cross(m.transpose(),Heff.transpose()) #first term in the LLG equation, using the effective field 
tau_prec_anis = Ms*V*mu0*np.cross(m.transpose(),Hanis.transpose()) #'reaction' cross-product torque, from anisotropies
# ('reaction' cross-product torque = internal mag. response to applied field torque, which can then transfer to mechanical torque)

if alpha2==0: #calculate damping torque based on solution (calculation depends on whether or not alpha was swept)
    tau_damp = Ms*V*(alpha(sol.t)/gamma)*np.cross(m.transpose(),dtm.transpose()) #damping torque (Nm)
else:
    a2 = np.asarray([alpha(sol.t),alpha(sol.t),alpha(sol.t)])
    cross_mdtm = np.cross(m.transpose(),dtm.transpose())
    tau_damp = (Ms*V/gamma)*a2.transpose()*cross_mdtm

tau_p = 1e18*tau_prec.transpose()               # aN m
tau_d = 1e18*tau_damp.transpose()               # aN m
tau_p_LLG = 1e18*tau_prec_LLG.transpose()       # aN m
tau_p_anis = 1e18*tau_prec_anis.transpose()     # aN m


""" Troubleshooting checks / print some variables """
calcTime = time.time() - startTime
cTh = int(calcTime/3600)
cTm = int((calcTime-cTh*3600)/60)
cTs = int((calcTime - cTh*3600 - cTm*60))

print('Too late!    (Calculation Time = {0} h {1} m {2:.1f} s)'.format(cTh,cTm,cTs))

m0 = mag(the0,phi0) #mostly for troubleshooting
thef,phif = sol.y[0][-1],sol.y[1][-1] #final angles from solution
mf = mag(thef,phif) #final magnetization (mostly for troubleshooting)


""" equilibrium EdH Torque calculation """  
tau_EdH = 1e18*A_EdH*dtm  # the EdH vector torque, in aN m

tau_res = tau_p - tau_EdH # resultant torque # (- sign in front of EdH torque to account for sign error in A_EdH)


"""and ratio of precession : EdH torques  """
ratio = np.divide(np.absolute(tau_p[0]),np.absolute(tau_EdH[0]+1e-7)) #tiny offset to avoid divide by zero


""" Set of strings to be used for the header in csv file: """
### csv file includes details of the simulation
### (Msat,alpha, Ku,Kc,starting angles, time details, geometry, demag factors, frequencies)
ps0 = '### Simulation Parameters:\n'
ps1 = ' Msat = {0}\n alpha = {1} alpha_2 = {5}\n Ku1 = {2}\n Kc1 = {3}; Kc2 = {4}\n'.format(Ms,alpha1,Ku_i,Kc1,Kc2,alpha2)
ps2 = ' theta_0 = {0}; phi_0 = {1}\n'.format(the0,phi0)
ps3 = ' starting time = {0} ns; ending time = {1} ns; time step = {2} ns\n'.format(t1*1e9,t2*1e9,dt*1e9)
ps4 = ' cylinder geometry: height = {0} um; diameter = {1} um\n'.format(h*1e6,D*1e6)
#ps5 = ' demag factors: Na = {0}; Nr = {1}\n'.format(Na,Nr)
ps5 = ' demag factors: Nx = {0}; Ny = {1}; Nz = {2}\n'.format(Nx,Ny,Nz)
ps6 = ' starting frequencies: fx_1 = {0} MHz; fy_1 = {1} MHz; fz_1 = {2} MHz\n'.format(fx1RF*1e6,fy1RF*1e6,fz1RF*1e6)
ps7 = ' starting frequencies: fx_2 = {0} MHz; fy_2 = {1} MHz; fz_2 = {2} MHz\n'.format(fx2RF*1e6,fy2RF*1e6,fz2RF*1e6)

parametersString = ps0 + ps1 + ps2 + ps3 +ps4 + ps5 + '######\n'

""" dataFrame with all the parameters to save in the csv file: """
df = pd.DataFrame({
    'time (ns)':t,
    'H applied (x) [A/m]':H[0], #index 1
    'H applied (y) [A/m]':H[1], #index 2
    'H applied (z) [A/m]':H[2], #index 3
    'H anisotropy (x) [A/m]':Hanis[0], #index 4
    'H anisotropy (y) [A/m]':Hanis[1], #index 5
    'H anisotropy (z) [A/m]':Hanis[2], #index 6
    'H effective (x) [A/m]':Heff[0], #index 7
    'H effective (y) [A/m]':Heff[1], #index 8
    'H effective (z) [A/m]':Heff[2], #index 9
    'RF Frequency (x) [GHz]':f_app[0]/1e9, #index 10
    'RF Frequency (y) [GHz]':f_app[1]/1e9, #index 11
    'RF Frequency (z) [GHz]':f_app[2]/1e9, #index 12   
    'mx':m[0], #index 13
    'my':m[1], #index 14
    'mz':m[2], #index 15
    'tau damping (x) [aN m]':tau_d[0], #index 16
    'tau damping (y) [aN m]':tau_d[1], #index 17
    'tau damping (z) [aN m]':tau_d[2], #index 18
    'tau precession (H_app) (x) [aN m]':tau_p[0], #index 19
    'tau precession (H_app) (y) [aN m]':tau_p[1], #index 20
    'tau precession (H_app) (z) [aN m]':tau_p[2], #index 21
    'tau precession (H_anis) (x) [aN m]':tau_p_anis[0], #index 22
    'tau precession (H_anis) (y) [aN m]':tau_p_anis[1], #index 23
    'tau precession (H_anis) (z) [aN m]':tau_p_anis[2], #index 24
    'tau precession (H_eff) (x) [aN m]':tau_p_LLG[0], #index 25
    'tau precession (H_eff) (y) [aN m]':tau_p_LLG[1], #index 26
    'tau precession (H_eff) (z) [aN m]':tau_p_LLG[2], #index 27
    'tau EdH (x) [aN m]':tau_EdH[0], #index 28
    'tau EdH (y) [aN m]':tau_EdH[1], #index 29
    'tau EdH (z) [aN m]':tau_EdH[2], #index 30
    'tau result (x) [aN m]':tau_res[0], #index 31
    'tau result (y) [aN m]':tau_res[1], #index 32
    'tau result (z) [aN m]':tau_res[2], #index 33
    'E_demag [aJ]':1e18*V*wD(sol.y[0],sol.y[1]), #index 34
    'E_Zeeman [aJ]':1e18*V*wZ(sol.y[0],sol.y[1],sol.t), #index 35
    'E_anisotropy_uniaxial [aJ]':1e18*V*wA_uniax(sol.y[0],sol.y[1],sol.t), #index 36
    'E_anisotropy_cubic [aJ]':1e18*V*wA_cub(sol.y[0],sol.y[1]) #index 37
    })

file_path = savepath+savename #savepath, savename defined above 

#write file header with all the simulation parameter info first
with open(file_path, 'w') as f: 
     f.write(parametersString)

#Then write dataFrame to the same file. 
df.to_csv(file_path, header=True, index=False,mode="a")

# df.to_pickle(savepath+nametosave+'.gz') #pickling data saves space when running large simulations.

"""

PLOT DEFINITIONS!

"""

"""
#CAPTION TEXT label includes:
- Python script name
- Save data name
- field sweep (range + axis)
- RF field amp (+ axis)
- freq sweep range 
- Ms
- dimensions + geometry
- damping constant, time steps

caption text included to easily identify which parameters were used in a given simulation.
"""
scriptName = os.path.basename(__file__)
fileName = file_path.replace('.csv','')
if Ku_i == Ku_f:
    Kutxt = 'Ku {} J/m^3; '.format(Ku_i)
else:
    Kutxt = 'Sweep Ku {0}-{1} J/m^3; '.format(Ku_i,Ku_f)

if HxDCi==HxDCf:
    if HxDCi==0:
        HxDCtxt = ''
    else:
        HxDCtxt = 'HxDC {} kA/m; '.format(np.round(HxDCi/1000,1))
else:
    HxDCtxt = 'Sweep HxDC {0}-{1} kA/m; '.format(np.round(HxDCi/1000,1),np.round(HxDCf/1000,1))
if HyDCi==HyDCf:
    if HyDCi==0:
        HyDCtxt = ''
    else:
        HyDCtxt = 'HyDC {} kA/m; '.format(np.round(HyDCi/1000,1))
else:
    HyDCtxt = 'Sweep HyDC {0}-{1} kA/m; '.format(np.round(HyDCi/1000,1),np.round(HyDCf/1000,1))
if HzDCi==HzDCf:
    if HzDCi==0:
        HzDCtxt = ''
    else:
        HzDCtxt = 'HzDC {} kA/m; '.format(np.round(HzDCi/1000,1))
else:
    HzDCtxt = 'Sweep HzDC {0}-{1} kA/m; '.format(np.round(HzDCi/1000,1),np.round(HzDCf/1000,1))

if HxRFi==HxRFf:
    if HxRFi==0:
        HxRFtxt = ''
    else:
        if fx1RF==fx2RF:
            HxRFtxt = 'RF x {0} A/m @ {1} MHz; '.format(HxRFi,fx1RF)
        else:
            HxRFtxt = 'RF x {0} A/m @ {1}-{2} MHz Sweep; '.format(HxRFi,fx1RF,fx2RF)
else:
    if fx1RF==fx2RF:
        HxRFtxt = 'RF x Sweep {0}-{1} A/m @ {2} MHz; '.format(HxRFi,HxRFf,fx1RF)
    else:
        HxRFtxt = 'RF x {0}-{1} A/m @ {2}-{3} MHz Sweep; '.format(HxRFi,HxRFf,fx1RF,fx2RF)

if HyRFi==HyRFf:
    if HyRFi==0:
        HyRFtxt = ''
    else:
        if fy1RF==fy2RF:
            HyRFtxt = 'RF y {0} A/m @ {1:.1f} MHz; '.format(HyRFi,fy1RF*1e-6)
        else:
            HyRFtxt = 'RF y {0} A/m @ {1:.1f}-{2:.1f} MHz Sweep; '.format(HyRFi,fy1RF*1e-6,fy2RF*1e-6)
else:
    if fy1RF==fy2RF:
        HyRFtxt = 'RF y Sweep {0}-{1} A/m @ {2:.1f} MHz; '.format(HyRFi,HyRFf,fy1RF*1e-6)
    else:
        HyRFtxt = 'RF y {0}-{1} A/m @ {2:.1f}-{3:.1f} MHz Sweep; '.format(HyRFi,HyRFf,fy1RF*1e-6,fy2RF*1e-6)

if HzRFi==HzRFf:
    if HzRFi==0:
        HzRFtxt = ''
    else:
        if fz1RF==fz2RF:
            HzRFtxt = 'RF z {0} A/m @ {1:.1f} MHz; '.format(HzRFi,fz1RF*1e-6)
        else:
            HzRFtxt = 'RF z {0} A/m @ {1:.1f}-{2:.1f} MHz Sweep; '.format(HzRFi,fz1RF*1e-6,fz2RF*1e-6)
else:
    if fz1RF==fz2RF:
        HzRFtxt = 'RF z Sweep {0}-{1} A/m @ {2:.1f} MHz; '.format(HzRFi,HzRFf,fz1RF*1e-6)
    else:
        HzRFtxt = 'RF z {0}-{1} A/m @ {2:.1f}-{3:.1f} MHz Sweep; '.format(HzRFi,HzRFf,fz1RF*1e-6,fz2RF*1e-6)

Mstxt = 'M_s = {} kA/m, '.format(np.round(Ms/1000,1))
alphatxt = 'damping constant {}, '.format(alpha1)
demagTxt = 'Nx = {0:.3f} Ny = {1:.3f} Nz = {2:.3f}'.format(Nx,Ny,Nz)
timeTxt = '{0:.1f} ps time steps, {1:.3f} us duration'.format(dt*1e12,(t2-t1)*1e6)

caption_text = scriptName+'\n'\
                +fileName+'\n'\
                +Kutxt+HxDCtxt+HyDCtxt+HzDCtxt+HxRFtxt+HyRFtxt+HzRFtxt+'\n'\
                +Mstxt + alphatxt + demagTxt+'\n'\
                +timeTxt
                
""" Plot magnetization """
def plotMag(i0):
    fig,ax = plt.subplots(3,1,figsize=(4,6))
    ax[0].plot(t[i0:],m[0][i0:])
    ax[1].plot(t[i0:],m[1][i0:])
    ax[2].plot(t[i0:],m[2][i0:])
    
    ax[2].set_xlabel('Time (ns)')
    ax[0].set_ylabel('$m_x$')
    ax[1].set_ylabel('$m_y$')
    ax[2].set_ylabel('$m_z$')
    
    ax[0].set_title('Magnetization')
    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    fig.subplots_adjust(hspace=0.05)
    fig.align_ylabels()
    plt.setp(ax,xlim=(t[i0],t[-1]))
    # plt.savefig('02-25-troubleshooting-Radau-working-poorly')

def plotTorque(i0):
    """ 
    Plot precession, anisotropy, and demag torques on a set of subplots
    """
    fig,ax = plt.subplots(3,2,figsize=(12,6))
    ax[0,0].plot(t[i0:],tau_p[0][i0:])
    ax[1,0].plot(t[i0:],tau_p[1][i0:])
    ax[2,0].plot(t[i0:],tau_p[2][i0:])
    
    ax[0,0].plot(t[i0:],tau_p_anis[0][i0:])
    ax[1,0].plot(t[i0:],tau_p_anis[1][i0:])
    ax[2,0].plot(t[i0:],tau_p_anis[2][i0:])
    
    ax[0,1].plot(t[i0:],tau_d[0][i0:])
    ax[1,1].plot(t[i0:],tau_d[1][i0:])
    ax[2,1].plot(t[i0:],tau_d[2][i0:])
    plt.setp(ax,xlim=(t[i0],t[-1]))
    
    ax[0,0].set_xticklabels([])
    ax[1,0].set_xticklabels([])
    ax[0,1].set_xticklabels([])
    ax[1,1].set_xticklabels([])
    
    ax[0,0].set_ylabel(r"$\tau_x^{\mathrm{m}\times\mathrm{B}}$ (aNm)")
    ax[1,0].set_ylabel(r"$\tau_y^{\mathrm{m}\times\mathrm{B}}$ (aNm)")
    ax[2,0].set_ylabel(r"$\tau_z^{\mathrm{m}\times\mathrm{B}}$ (aNm)")
    
    ax[0,1].set_ylabel(r"$\tau_x^{\mathrm{damp}}$ (aNm)")
    ax[1,1].set_ylabel(r"$\tau_y^{\mathrm{damp}}$ (aNm)")
    ax[2,1].set_ylabel(r"$\tau_z^{\mathrm{damp}}$ (aNm)")
    
    ax[2,0].set_xlabel('Time (ns)')
    ax[2,1].set_xlabel('Time (ns)')
    
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.35)
    fig.align_ylabels()
    ax[0,0].set_title('Precession Torques (applied field and anisotropy)')
    ax[0,1].set_title('Damping Torque')
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.125, top=0.875, hspace=0.06)
    
    plt.text(0,-0.075, caption_text,fontsize=10,transform=fig.transFigure)

def plotResTorque(i0): #plot resultant torque 
    
    fig,ax = plt.subplots(3,1,figsize=(4,6))
    ax[0].plot(t[i0:],(-tau_EdH[0]+tau_p[0])[i0:])
    ax[1].plot(t[i0:],(-tau_EdH[1]+tau_p[1])[i0:])
    ax[2].plot(t[i0:],(-tau_EdH[2]+tau_p[2])[i0:])
    
    ax[2].set_xlabel('Time (ns)')
    ax[0].set_ylabel(r"$\tau_x^{\mathrm{EdH}}$ (aNm)")
    ax[1].set_ylabel(r"$\tau_y^{\mathrm{EdH}}$ (aNm)")
    ax[2].set_ylabel(r"$\tau_z^{\mathrm{EdH}}$ (aNm)")
    
    ax[0].set_title(r'Resultant Torque ($\tau^{mxB}$ + $\tau^{EdH}$)')
    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    fig.subplots_adjust(hspace=0.05)
    fig.align_ylabels()
    plt.setp(ax,xlim=(t[i0],t[-1]))
    
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.125, top=0.875, hspace=0.06)

    plt.text(0,-0.075, caption_text,fontsize=10,transform=fig.transFigure)


def plotEdH(i0): #plot EdH torque only
    fig,ax = plt.subplots(3,1,figsize=(4,6))
    ax[0].plot(t[i0:],tau_EdH[0][i0:])
    ax[1].plot(t[i0:],tau_EdH[1][i0:])
    ax[2].plot(t[i0:],tau_EdH[2][i0:])
    
    ax[2].set_xlabel('Time (ns)')
    ax[0].set_ylabel(r"$\tau_x^{\mathrm{EdH}}$ (aNm)")
    ax[1].set_ylabel(r"$\tau_y^{\mathrm{EdH}}$ (aNm)")
    ax[2].set_ylabel(r"$\tau_z^{\mathrm{EdH}}$ (aNm)")
    
    ax[0].set_title('EdH Torque')
    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    fig.subplots_adjust(hspace=0.05)
    fig.align_ylabels()
    plt.setp(ax,xlim=(t[i0],t[-1]))
    
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.125, top=0.875, hspace=0.06)

    plt.text(0,-0.075, caption_text,fontsize=10,transform=fig.transFigure)

def plotRatio(i0):
    plt.figure(figsize=(6,4))
    plt.plot(t[i0:],ratio[i0:],linestyle='none',marker='.',markersize=0.5)
        
    plt.ylim(0.95,1.05)
    #plt.xlim(2,3)
    
    plt.xlabel('Time (ns)')
    plt.ylabel(r"$\tau_x^{\mathrm{prec}}/\tau_x^{\mathrm{EdH}}$")
    
    plt.title('precession / EdH torque ratio')
    
    plt.text(0,0.91, caption_text,fontsize=10)


""" CALL PLOTS """

# plotMag(0) #plot magnetization (mx,my,mz components)
# plt.draw()

plotTorque(0)  # plot torque components (precession, anisotropy, demag)
plt.draw()

plotEdH(0)  # plot EdH torques
plt.draw()

plotResTorque(0) # plot resultant torque 
plt.draw()
# plt.savefig(savepath+nametosave+'.png', format='png', dpi=600, transparent=False, bbox_inches='tight') #save resultant torque plot



