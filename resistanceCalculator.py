# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 14:14:43 2018

@author: u0112721
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import datetime
import pygfunction as gt
from scipy.constants import pi
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.special import j0, j1, y0, y1, exp1
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd


def shortTermCorrection(time, gFunc, r_b, aSoi):
    
    def _CHS(u, Fo, p):
        CHS_integrand = 1./(u**2*pi**2)*(np.exp(-u**2*Fo) - 1.0) / (j1(u)**2 + y1(u)**2) * (j0(p*u)*y1(u) - j1(u)*y0(p*2))
        return CHS_integrand
    
    def _ILS(t, aSoi, dis):
        ILS = exp1(dis**2/(4*aSoi*t))
        return ILS
    
    for i in range(len(time)):
        ILS = _ILS(time[i], aSoi, r_b)
        CHS, err = quad(
            _CHS, 1e-12, 100., args=(aSoi*time[i]/r_b**2, 1.))
        gFunc[i] = gFunc[i] + 2*pi*CHS - 0.5*ILS

    return gFunc

def resistanceCalculator(time, gFunc, kSoi):
    f = interp1d(time, gFunc)
    Rh = f(6*3600)/(2*pi*kSoi)
    Rm = (f(30*24*3600+6*3600) - f(6*3600))/(2*pi*kSoi)
    Ra = (f(20*365*24*3600+30*24*3600+6*3600) - f(30*24*3600+6*3600))/(2*pi*kSoi)

    return Rh,Rm,Ra

def computeResistances(qh,qm,qa,kSoi, aSoi='sentinel' ,Rb=0.144,Tg=10,Tf=18):
    # -------------------------------------------------------------------------
    # Simulation parameters
    # -------------------------------------------------------------------------

    # Load parameters
    
    #qh =           # Peak load (W)
    #qm =           # Monthly load (W)
    #qa =           # Annual unbalance (W)

    # Borehole parameters
    D = 0.0             # Borehole buried depth (m)
    H = 100             # Borehole length (m)
    r_b = 0.075         # Borehole radius (m)
    B = 6.0             # Borehole spacing (m)
    #Rb = 0.144          # Borehole effective thermal resistance ((mK)/W)

    # Soil thermal properties
    if aSoi is 'sentinel':
        print("WARNING: No diffusivity value was provided, assuming a volumetric thermal capacity CSoi = ", 1200*1800,"J/(m3K)")
        aSoi = kSoi/1800/1200
    #kSoi = 1.8 				# Ground thermal conductivity (W/(mK))
    #cSoi = 1200				# Specific heat capacity of the soil (J/(kgK))
    #dSoi = 1800				# Density of the soil (kg/m3)
    CSoi = kSoi/aSoi            # Volumetric thermal capacity of the soil (J/(m3K))
    #aSoi = kSoi/cSoi/dSoi   # Ground thermal diffusivity (m2/s)
    #Tg = 10                 # Undisturbed ground temperature (degC)

    # Temperature constraint
    #Tf = 18                 # Check whether it is for cooling or heating! (degC)

    # Number of segments per borehole
    nSegments = 1

    # Geometrically expanding time vector. 

    nt = 75						   # Number of time steps
    ts = H**2/(9.*aSoi)            # Bore field characteristic time
    ttsMax = np.exp(5)
    dt = 300.		                # Time step (s)
    tmax = ttsMax*ts                # Maximum time

    time = gt.utilities.time_geometric(dt, tmax, nt)

    print("IMPORTANT: Heat injected to the ground is (+)")
    print("IMPORTANT: Heat extracted from the ground is (-)")
    print("\n-------- Introduced arguments ---------")
    print("Thermal conductivity of the ground kSoi = ", kSoi,"W/(mK)")
    print("Thermal diffusivity of the ground aSoi = ", aSoi, "m2/s")
    print("Undisturbed ground temperature Tg = ", Tg, "degC")
    print("Fluid temperature constraint Tf = ", Tf, "degC")
    print("Borehole effective resistance Rb =", Rb,"(mK)/W")
    print("Volumetric thermal capacity of the ground CSoi =", CSoi, "J/(m3.K)")


    # -------------------------------------------------------------------------
    # Borehole fields
    # -------------------------------------------------------------------------
    print("\nCalculations for one borehole in progress...")
    #Calculation of the thermal response of a 100m borehole
    boreHole = gt.boreholes.Borehole(H, D, r_b, 0, 0)
    gFunc = gt.gfunction.uniform_temperature([boreHole], time, aSoi, nSegments=nSegments, disp=True)
    gFunc = shortTermCorrection(time, gFunc, r_b, aSoi)
    #Adding zero as the first element
    time = np.insert(time, 0, 0)
    gFunc = np.insert(gFunc, 0, 0) 

    #gFunc = gFunc / (2*np.pi*kSoi*H*nBor)

    #Initial guess using a g-function of H=100 borehole
    Rh, Rm, Ra = resistanceCalculator(time,gFunc,kSoi)

    L = (qh*(Rb+Rh)+qm*Rm+qa*Ra)/(Tf-Tg)

    print("Required length for a borehole with no interferences: ", L,"m")

    N_b = int(np.ceil(L/100))         # Number of boreholes of a borefield with 50m each borehole, rounded up  
    N_1 = int(np.ceil(np.sqrt(N_b))) # Number of boreholes in one direction
    N_2 = int(np.ceil(N_b/N_1))      # Number of boreholes in opposite direction   

    L = N_1*N_2*100                      # Redefinition of the length 

    L_old = 0
    i = 0 
    print("Starting iterations with a field length of ", L, "m and ", N_1*N_2, " boreholes in a ", N_1,"x",N_2," field")

    while abs((L-L_old)/L*100.) > 1:
        #Iteration number
        i+=1
        print("\n +++ Iteration number: ", i, "+++")
        # Redefinition of time vector
        time = gt.utilities.time_geometric(dt, tmax, nt)
        # Redefinition of the borefield g-function
        #boreField = gt.boreholes.circle_field(N_b, B, H, D, r_b)
        boreField = gt.boreholes.rectangle_field(N_1, N_2, B, B, H, D, r_b)
        gFunc = gt.gfunction.uniform_temperature(boreField, time, aSoi, nSegments=nSegments, disp=True)
        gFunc = shortTermCorrection(time, gFunc, r_b, aSoi)
        #Adding zero as the first element
        time = np.insert(time, 0, 0)
        gFunc = np.insert(gFunc, 0, 0) 

        #Calculation of the new thermal resistances
        Rh, Rm, Ra = resistanceCalculator(time,gFunc,kSoi)
        print("Computed 6h peak resistance: ", Rh, "(mK)/W")
        print("Computed monthly resistance: ", Rm, "(mK)/W")
        print("Computed yearly resistance: ", Ra, "(mK)/W")

        L_old = L 
        L = (qh*(Rb+Rh)+qm*Rm+qa*Ra)/(Tf-Tg)
        H = L/N_b

        #print("New field length: ", L, "m")
        print("Convergence criterion: ", abs((L-L_old)/L*100.), "%")


    print("\n########## ITERATION PROCESS FINISHED ##########")
    #print("The final length of the field is", L, "m")
    #print("The final length per borehole is", L/N_1/N_2, "m")
    print("Final 6h peak resistance: ", Rh, "(mK)/W")
    print("Final monthly resistance: ", Rm, "(mK)/W")
    print("Final yearly resistance: ", Ra, "(mK)/W")

    # -------------------------------------------------------------------------
    # Figure
    # -------------------------------------------------------------------------

    gt.boreholes.visualize_field(boreField)

    plt.rc('figure')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # Axis labels
    ax1.set_xlabel(r'$ln(t/t_s)$')
    ax1.set_ylabel(r'$g(t/t_s)$')
    # Axis limits
    #ax1.set_xlim([-10.0, 5.0])
    #ax1.set_ylim([0., 20.])
    # Show minor ticks
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    # Adjust to plot window
    plt.tight_layout()
    # Draw g-function
    ax1.plot(np.log(time[1:]/ts), gFunc[1:], 'k-', lw=1.5, label='g-function')
    ax1.legend(loc='upper left')

    plt.show()

    return