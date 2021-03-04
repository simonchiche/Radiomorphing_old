#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 02:45:28 2021

@author: chiche
"""
import h5py
import numpy as np
import hdf5fileinout as hdf5io
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
        
        # TODO: add a function to select the right reference shower (and right planes? or do we use all the planes?)
        RefShower = extractData(sim_dir[i])
        
        TargetShower = Shower(shower['primary'], shower['energy'], shower['zenith'], \
        shower['azimuth'], shower['injection'], RefShower.nant, RefShower.inclination,\
        shower['altitude'], RefShower.pos, RefShower.traces , \
        RefShower.dplane, RefShower.xmaxdist, RefShower.xmaxpos)
        
        myscale(sim_dir[i],  RefShower, TargetShower)

        # interpolate the traces.
        # interpolate(antennas, sim_dir[i], out_dir,
                #shower["zenith"], shower["azimuth"], shower["injection_height"])
   
class Shower:
    
    def __init__(self, primary, energy, zenith, azimuth, injection_height, \
                 nantennas, inclination, GroundAltitude, Positions, \
                 Traces, Distplane, XmaxDistance, XmaxPosition):
    
        self.primary = primary
        self.energy = energy
        self.zenith = zenith
        self.azimuth = azimuth
        self.injection = injection_height
        self.nant = nantennas
        self.inclination = inclination
        self.pos = Positions
        self.traces = Traces
        self.dplane = Distplane
        self.glevel = GroundAltitude
        self.xmaxdist = XmaxDistance
        self.xmaxpos = XmaxPosition
       

def extractData(sim_file):
    
    simu_path = './' + sim_file
    InputFilename = simu_path
    filehandle = h5py.File(InputFilename, 'r')

    #Shower event access
    RunInfo = hdf5io.GetRunInfo(InputFilename)         
    NumberOfEvents = hdf5io.GetNumberOfEvents(RunInfo)
    EventNumber = NumberOfEvents-1 
    EventName = hdf5io.GetEventName(RunInfo,EventNumber)       
    #All simulations info

    #Shower parameters
    Zenith = hdf5io.GetEventZenith(RunInfo,EventNumber)
    Azimuth = hdf5io.GetEventAzimuth(RunInfo,EventNumber)
    Primary = "electron"
    Energy = hdf5io.GetEventEnergy(RunInfo,EventNumber)
    XmaxDistance = hdf5io.GetEventXmaxDistance(RunInfo,EventNumber)
    #SlantXmax = hdf5io.GetEventSlantXmax(RunInfo,EventNumber) # Slant ?
    Energy = hdf5io.GetEventEnergy(RunInfo,EventNumber)
    #HadronicModel = hdf5io.GetEventHadronicModel(RunInfo,EventNumber)    
    
    #Shower info
    EventInfo = hdf5io.GetEventInfo(InputFilename,EventName)
    XmaxPosition = hdf5io.GetEventXmaxPosition(EventInfo)
    BFieldIncl = hdf5io.GetEventBFieldIncl(EventInfo) #BfieldIncl/Bfieldecl
    #BFieldDecl = hdf5io.GetEventBFieldDecl(EventInfo)
    GroundAltitude = hdf5io.GetGroundAltitude(EventInfo)    
    #Antannas info
    AntennaInfo = hdf5io.GetAntennaInfo(InputFilename,EventName)


    NumberOfAntennas = hdf5io.GetNumberOfAntennas(AntennaInfo)
    IDs_bis = hdf5io.GetAntIDFromAntennaInfo(AntennaInfo)
    X = hdf5io.GetXFromAntennaInfo(AntennaInfo)
    Y = hdf5io.GetYFromAntennaInfo(AntennaInfo)
    Z = hdf5io.GetZFromAntennaInfo(AntennaInfo)
    Positions = np.transpose([X,Y,Z])

    #Traces
    AntennaID = IDs_bis[0]
    Efield_trace = hdf5io.GetAntennaEfield(InputFilename,EventName,AntennaID,OutputFormat="numpy")
    time_sample = len(Efield_trace[:,0])
    
    trace_x = np.zeros([time_sample,NumberOfAntennas])
    trace_y = np.zeros([time_sample,NumberOfAntennas])
    trace_z = np.zeros([time_sample,NumberOfAntennas])
    Time = np.zeros([time_sample,NumberOfAntennas])
    
    for i in range(NumberOfAntennas):
    
        AntennaID = IDs_bis[i]
        path = str(EventName)+"/AntennaTraces/"+str(AntennaID)+"/efield"
        Efield_trace = filehandle[path]
        Time[:,i] = Efield_trace['Time']
        trace_x[:,i] = Efield_trace['Ex']
        trace_y[:,i] = Efield_trace['Ey']
        trace_z[:,i] = Efield_trace['Ez']
               
    Traces = np.transpose(np.concatenate((np.transpose(Time), np.transpose(trace_x), \
                                          np.transpose(trace_y), np.transpose(trace_z))))
              
    Distplane = 50000
    Nant = len(X)
    Injection = 1e5 # TODO: get it from the hdf5
    
    
    RefShower = Shower(Primary, Energy, Zenith, Azimuth, Injection, Nant, BFieldIncl, GroundAltitude,
                        Positions, Traces, Distplane, XmaxDistance, XmaxPosition[0])

    return RefShower  
    