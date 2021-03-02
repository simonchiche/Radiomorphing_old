#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 02:45:28 2021

@author: chiche
"""
from scaling import myscale
    
def process(sim_dir, shower,  out_dir):
    """Rescale and interpolate the radio traces for all antennas 
        - start the Radio Morphing procedure

    Parameters:
    ----------
        sim_dir: str
            path to the simulated traces
        shower: dict
            properties of the requested shower
        antennas: str
            path the requested antenna positions
        out_dir: str
            path where the output traces should be dumped
    """
    # Rescale the simulated showers to the requested one
    # print "ATTENTION scaling commented"
    
    for i in range(len(sim_dir)):
        myscale(sim_dir[i],  **shower)

        # interpolate the traces.
        #interpolate(antennas, sim_dir[i], out_dir,
                #shower["zenith"], shower["azimuth"], shower["injection_height"])

