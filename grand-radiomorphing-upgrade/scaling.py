#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 15:28:02 2021

@author: chiche
"""
import matplotlib.pyplot as plt
from ModuleScale import EnergyScale, GeomagneticScale, DensityScale, CerenkovStretch
   
#def myscale(sim_file, primary, energy, zenith, azimuth):
def myscale(RefShower, TargetShower):
     
    Nant = RefShower.nant
    
    # Translation in the shower plane
    TargetShower.pos, TargetShower.traces = RefShower.GetinShowerPlane()
           
    # Energy scaling
    TargetShower.traces[:,Nant:], kE = EnergyScale(RefShower, TargetShower)
    
    # Geomagnetic scaling
    TargetShower.traces[:,2*Nant:3*Nant], kgeo = GeomagneticScale(RefShower, TargetShower)
    
    # Density scaling
    TargetShower.traces[:,2*Nant:3*Nant], TargetShower.xmaxpos, krho = DensityScale(RefShower, TargetShower)
    
    # Layout and traces stretching
    TargetShower.pos, TargetShower.traces[:,Nant:], kstretch = CerenkovStretch(RefShower, TargetShower)
    
    # Back in the geographic plane
    TargetShower.pos, TargetShower.traces = TargetShower.GetinGeographicFrame()
    
    plt.scatter(TargetShower.pos[:,1], TargetShower.pos[:,2])
    # TODO: include magnetic field scaling 
    
    print(RefShower.distplane, TargetShower.distplane)
    
    return TargetShower
        