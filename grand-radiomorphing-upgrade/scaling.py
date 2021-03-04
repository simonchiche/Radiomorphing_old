#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 15:28:02 2021

@author: chiche
"""


import numpy as np
from ModuleScale import GetinShowerPlane, EnergyScale, GeomagneticScale, DensityScale, CerenkovStretch
import matplotlib.pyplot as plt
   
#def myscale(sim_file, primary, energy, zenith, azimuth):
def myscale(sim_file, RefShower, TargetShower):
     
    Nant = RefShower.nant
    
    # Translation in the shower plane
    TargetShower.pos, TargetShower.traces = GetinShowerPlane(RefShower)
          
    # Energy scaling
    TargetShower.traces[:,Nant:], kE = EnergyScale(RefShower, TargetShower)
    
    # Geomagnetic scaling
    TargetShower.traces[:,2*Nant:3*Nant], kgeo = GeomagneticScale(RefShower, TargetShower)
    
    # Density scaling
    TargetShower.traces[:,2*Nant:3*Nant], TargetShower.xmaxpos, krho = DensityScale(RefShower, TargetShower)
    
    # Layout and traces stretching
    TargetShower.pos, TargetShower.traces[:,Nant:], kstretch = CerenkovStretch(RefShower, TargetShower)
    
    # Back in the geographic plane
    #TargetShower.pos = GetinGeographicFrame(TargetShower)
    
    # TODO: include magnetic field scaling 
    
    return TargetShower
        