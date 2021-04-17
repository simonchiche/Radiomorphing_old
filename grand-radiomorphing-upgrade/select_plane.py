#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 19:08:09 2021

@author: chiche
"""


import numpy as np
import glob

def select_zenith(target_zenith):
    
    target_zenith = 180 - target_zenith # cosmic ray convention
    
    #print(target_zenith)
    
    zenith_sim  = np.array([67.8, 74.8, 81.3, 83.9, 86.5])
   # print(zenith_sim)
    
    min_index = np.argmin(abs(zenith_sim - target_zenith))
    
    selected_zenith =  zenith_sim[min_index]
    
    return selected_zenith

def select_azimuth(target_azimuth):
    
    if(target_azimuth>180.0):
        
        diff = target_azimuth - 180
        target_azimuth = target_azimuth - 2*diff
        
    #print(target_azimuth)
    
    azimuth_sim  = np.array([0, 90, 180])
    
    min_index = np.argmin(abs(azimuth_sim - target_azimuth))
    
    selected_azimuth =  azimuth_sim[min_index]
    
    return selected_azimuth

def select_path(path, dplane):
    
    sim = glob.glob(path)
    n = len(sim)
    dsim = np.zeros(n)
    
    for i in range(n):
        
        dsim[i] =  float(sim[i].split("_")[-1][:-5])
        
    min_index = np.argmin(abs(dsim - dplane))
        
    return sim[min_index], dsim[min_index]
        
        
def select_plane(primary, energy, zenith, azimuth, injection, altitude, fluctuations, dplane=5000):
    
    
    target_zenith = select_zenith(zenith)
    target_azimuth = select_azimuth(azimuth)
    
    path = "./Simulations/SelectedPlane/theta_%.1f/*.hdf5" \
                      %(target_zenith)


    selected_plane, dsim = select_path(path, dplane)    
    
    return selected_plane

def print_plane(RefShower, TargetShower):
    
    print("-----------------------")
    print("Target shower: Energy = %.2f, Azimuth = %.2f, Zenith = %.2f, \
          Dxmax = %.2d" %(TargetShower.energy, 180 -TargetShower.zenith,\
          TargetShower.azimuth, TargetShower.distplane))
    print("")
    print("Ref shower: Energy = %.2f, Azimuth = %.2f, Zenith = %.2f, \
          Dxmax = %.2d" %(RefShower.energy, 180 -RefShower.zenith,\
          RefShower.azimuth, RefShower.distplane))
    print("-----------------------")
    

    
    
    