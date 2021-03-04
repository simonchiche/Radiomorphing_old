#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 14:51:26 2021

@author: chiche
"""


import os
import h5py
import numpy as np
import hdf5fileinout as hdf5io
import sys
import glob


# =============================================================================
#                            Get in shower plane 
# =============================================================================

def get_center(shower):

    azimuth = shower.azimuth
    zenith = shower.zenith
    distground = shower.xmaxdist
    xmaxpos = shower.xmaxpos
    x_Xmax, y_Xmax, z_Xmax = xmaxpos[0], xmaxpos[1], xmaxpos[2]
    GroundLevel = shower.glevel
    
    pi = np.pi
    zenith = zenith*pi/180.0
    azimuth = azimuth*pi/180.0
    uv = np.array([np.sin(zenith)*np.cos(azimuth), np.sin(zenith)*np.sin(azimuth) , np.cos(zenith)]) # direction of the shower

    distground = np.sqrt(x_Xmax**2 + y_Xmax**2 + (z_Xmax-GroundLevel)**2)
    distplane = np.mean(get_distplane(shower))

    dist_plane_ground = distground - distplane
    core = uv*(dist_plane_ground)
    core[2] = core[2] - GroundLevel
    
    return core

def get_distplane(shower):
 
    #function that returns "w" at each antenna, i.e. the angle between the 
    #direction that goes from Xmax to the core and the direction that 
    #goes from Xmax to a given antenna
    
    azimuth = shower.azimuth
    zenith = shower.zenith
    pos =  shower.pos
    x, y, z = pos[:,0], pos[:,1], pos[:,2]
    xmaxpos = shower.xmaxpos
    x_Xmax, y_Xmax, z_Xmax = xmaxpos[0], xmaxpos[1], xmaxpos[2]
    
    pi = np.pi
    zenith = zenith*pi/180.0
    azimuth = azimuth*pi/180.0
    
    x_antenna = x - x_Xmax # distance along the x-axis between the antennas postions and Xmax
    y_antenna = y - y_Xmax
    z_antenna = z - z_Xmax
    
    uv = np.array([np.sin(zenith)*np.cos(azimuth), np.sin(zenith)*np.sin(azimuth) , np.cos(zenith)]) # direction of the shower
    u_antenna = np.array([x_antenna, y_antenna, z_antenna]) # direction of the unit vectors that goes from Xmax to the position of the antennas
    distplane = np.dot(np.transpose(u_antenna), uv)
    
    
    return distplane

def GetinShowerPlane(shower):
    
    # function that returns the trcaes in the shower plane (v, vxb, vxvxb) 
    #from the traces in the geographic plane (x, y, z)

    zenith = shower.zenith*np.pi/180.0
    azimuth = shower.azimuth*np.pi/180.0
    inclination = shower.inclination*np.pi/180.0
    
    pos =  shower.pos
    x, y, z = pos[:,0], pos[:,1], pos[:,2]
    n = len(x) # number of antennas
    
    # We move the core position in (0,0,0) before changing the reference frame
    core = get_center(shower)
    
    x = x + core[0]
    y = y + core[1] 
    z = z + core[2]
    
    Traces = shower.traces
    time_sample = len(Traces[:,0])
    
    # antennas positions in the  shower reference frame (v, vxB, vxvxB)
    v = np.zeros(n)   
    vxb = np.zeros(n)
    vxvxb = np.zeros(n)
    
    #Traces in the shower reference frame
    Traces_Ev = np.zeros([time_sample,n])
    Traces_Evxb = np.zeros([time_sample,n])
    Traces_Evxvxb = np.zeros([time_sample,n])
    Time = np.zeros([time_sample, n])
    
    # unit vectors 
    uv = np.array([np.sin(zenith)*np.cos(azimuth), np.sin(zenith)*np.sin(azimuth) , np.cos(zenith)]) # direction of the shower
    uB = np.array([np.cos(inclination), 0, -np.sin(inclination)]) # direction of the magnetic field 
    
    uv_x_uB = np.cross(uv, uB) # unit vector along the vxb direction
    uv_x_uB /= np.linalg.norm(uv_x_uB) # normalisation
    
    uv_x_uvxB  = np.cross(uv, uv_x_uB) # unit vector along the vxvxb direction
    uv_x_uvxB /= np.linalg.norm(uv_x_uB) # normalisation
    
    P = np.transpose(np.array([uv, uv_x_uB, uv_x_uvxB])) # matrix to go from the shower reference frame to the geographic reference frame
    
    P_inv = np.linalg.inv(P) # matrix to go from the geographic reference frame to the shower reference frame
    
    # We calculate the positions in the shower plane
    Position_geo = np.array([x,y,z]) # position in the geographic reference frame
    Position_shower = np.dot(P_inv, Position_geo) # position in the shower reference frame
    
    # We deduce the different components
    v = Position_shower[0, :] 
    vxb = Position_shower[1, :]
    vxvxb =  Position_shower[2, :]
    
    # We calulate the traces in the shower plane
    Traces_geo = np.zeros([time_sample,3])
    Traces_shower_temp = np.zeros([3, time_sample])
    
    for i in range(n):
        
        Traces_geo = np.array([Traces[:,i + n], Traces[:, i + 2*n], Traces[:, i + 3*n]])
        
        Traces_shower_temp = np.dot(P_inv, Traces_geo)
        
        Traces_Ev[:,i] = np.transpose(Traces_shower_temp[0,:]) # Ev component of the traces
        Traces_Evxb[:,i] = np.transpose(Traces_shower_temp[1,:]) # Evxb component of the traces
        Traces_Evxvxb[:,i] = np.transpose(Traces_shower_temp[2,:]) # Evxvxb component of the traces
        
        Time[:,i] = Traces[:,i]
    
    # We derive the traces in the shower plane
    
    Traces_sp = np.transpose(np.concatenate((np.transpose(Time), np.transpose(Traces_Ev), np.transpose(Traces_Evxb), np.transpose(Traces_Evxvxb))))

    Positions_sp = np.zeros([n,3]) # To chck
    Positions_sp[:,0], Positions_sp[:,1], Positions_sp[:,2] = v, vxb, vxvxb

    return Positions_sp, Traces_sp

# =============================================================================
#                           Energy scale
# =============================================================================
    
def EnergyScale(RefShower, TargetShower):
    
    Nant = RefShower.nant
    ref_energy = RefShower.energy
    target_energy = TargetShower.energy
    
    kE = target_energy/ref_energy
    
    scaled_traces = TargetShower.traces[:,Nant:]*kE
    
    return scaled_traces, kE

# =============================================================================
#                       Geomagnetic angle scale
# =============================================================================

def get_alpha(shower):
    
    # function that returns the angle between the direction of the shower and the direction of the magnetic field 
    
    zenith = shower.zenith
    azimuth = shower.azimuth
    inclination = shower.inclination
    
    pi = np.pi
    zenith = zenith*pi/180.0
    inclination = inclination*pi/180.0
    azimuth = azimuth*pi/180.0
    
    # unit vectors    
    uv = np.array([np.sin(zenith)*np.cos(azimuth), np.sin(zenith)*np.sin(azimuth) , np.cos(zenith)]) # direction of the shower
    uB = np.array([np.cos(inclination), 0, -np.sin(inclination)]) # direction of the magnetic field
    cos_alpha = np.dot(uv,uB)
    alpha = np.arccos(cos_alpha) # angle between the direction of the shower and the direction of the magnetic field
    
    return alpha


def get_w(shower):
 
    #function that returns "w" at each antenna, i.e. the angle between the direction that goes from Xmax to the core and the direction that goes from Xmax to a given antenna
    
    zenith = shower.zenith
    azimuth = shower.azimuth
    inclination = shower.inclination
    xmaxpos = shower.xmaxpos
    pos = shower.pos
    
    x_Xmax, y_Xmax, z_Xmax =  xmaxpos[0], xmaxpos[1], xmaxpos[2]
    x, y, z = pos[:,0], pos[:,1], pos[:,2] 
    
    
    pi = np.pi
    zenith = zenith*pi/180.0
    inclination = inclination*pi/180.0
    azimuth = azimuth*pi/180.0
    
    x_antenna = x - x_Xmax # distance along the x-axis between the antennas postions and Xmax
    y_antenna = y - y_Xmax
    z_antenna = z - z_Xmax
    
    uv = np.array([np.sin(zenith)*np.cos(azimuth), np.sin(zenith)*np.sin(azimuth) , np.cos(zenith)]) # direction of the shower
    u_antenna = np.array([x_antenna, y_antenna, z_antenna]) # direction of the unit vectors that goes from Xmax to the position of the antennas
    u_antenna /= np.linalg.norm(u_antenna, axis =0)
    w = np.arccos(np.dot(np.transpose(u_antenna), uv))
    w = w*180.0/pi # we calculte w in degrees
    
    return w  
    

def GeomagneticScale(RefShower, TargetShower):
    
    Nant = RefShower.nant
    
    ref_alpha = get_alpha(RefShower)
    target_alpha = get_alpha(TargetShower)
    
    kgeo  = np.sin(target_alpha)/np.sin(ref_alpha)
    
    if(kgeo<0.9): # TODO: refine this condition, rather a condition over theta
    
        vxb, vxvxb = TargetShower.pos[:,1], TargetShower.pos[:,2]
        cos_eta = vxb/np.sqrt(vxb**2 + vxvxb**2)
    
        Evxb, Evxvxb = TargetShower.traces[:,2*Nant:3*Nant], TargetShower.traces[:,3*Nant:4*Nant] 
        
        w = get_w(RefShower)
        
        w_kxkxb = []
        Ece = []
        Egeo = []
        
        for j in range(len(w)):
            if(abs(vxb[j])<1):
                w_kxkxb.append(w[j])
                Ece.append(Evxvxb[:,j])
                Egeo.append(Evxb[:,j])
          
        
        Evxb_scaled = np.zeros(np.shape(Evxb))
        
        for i in range(len(w)):
            
            diff = abs(w_kxkxb - w[i])
            minimum = np.argmin(diff)
            Evxb_scaled[:,i] = Egeo[minimum]*kgeo + Ece[minimum]*cos_eta[i]
            
    else: 
        Evxb_scaled = TargetShower.traces[:,2*Nant:3*Nant]*kgeo
    
    
    return Evxb_scaled, kgeo
    
 
# =============================================================================
#                          Density scaling
# =============================================================================

def Xmax_param(primary, energy):

    #input energy in EeV
    
    if(primary == 'electron'):
        a =65 
        c =270
        return a*np.log10(energy*1e6) +c
    else:
        print("missing primary")  
        sys.exit()
        

def _get_CRzenith(zen,injh,GdAlt):
    ''' Corrects the zenith angle for CR respecting Earth curvature, zenith seen by observer
        ---fix for CR (zenith computed @ shower core position
    
    Arguments:
    ----------
    zen: float
        GRAND zenith in deg
    injh: float
        injection height wrt to sealevel in m
    GdAlt: float
        ground altitude of array/observer in m (should be substituted)
    
    Returns:
    --------
    zen_inj: float
        GRAND zenith computed at shower core position in deg
        
    Note: To be included in other functions   
    '''

    #Note: To be included in other functions
    Re= 6370949 # m, Earth radius

    a = np.sqrt((Re + injh)**2. - (Re+GdAlt)**2 *np.sin(np.pi-np.deg2rad(zen))**2) - (Re+GdAlt)*np.cos(np.pi-np.deg2rad(zen))
    zen_inj = np.rad2deg(np.pi-np.arccos((a**2 +(Re+injh)**2 -Re**2)/(2*a*(Re+injh))))
    
    return zen_inj
    
    
    
def _dist_decay_Xmax(zen, injh2, Xmax_primary): 
    ''' Calculate the height of Xmax and the distance injection point to Xmax along the shower axis
    
    Arguments:
    ----------
    zen: float
        GRAND zenith in deg, for CR shower use _get_CRzenith()
    injh2: float
        injectionheight above sealevel in m
    Xmax_primary: float
        Xmax in g/cm2 
        
    Returns:
    --------
    h: float
        vertical Xmax_height in m
    ai: float
        Xmax_distance injection to Xmax along shower axis in m
    '''
    
    #% Using isothermal Model
    rho_0 = 1.225*0.001#; % kg/m3 to 0.001g/cm3: 1g/cm3=1000kg/m3, since X given in g/cm2
    M = 0.028966#;  %kg/mol - 1000g/mol
    g = 9.81#; %ms-2
    T = 288.#; % K
    R = 8.32#; J/K/mol , J=kg m2/s2

    zen2 = np.deg2rad(zen)
    
    hD=injh2
    step=10 #m
    if hD>10000:
        step=100 #m
    Xmax_primary= Xmax_primary#* 10. # g/cm2 to kg/m2: 1g/cm2 = 10kg/m2
    gamma=np.pi-zen2 # counterpart of where it goes to
    Re= 6370949 # m, Earth radius
    X=0.
    i=0.
    h=hD
    ai=0
    while X< Xmax_primary:
        i=i+1
        ai=i*step #100. #m
        hi= -Re+np.sqrt(Re**2. + ai**2. + hD**2. + 2.*Re*hD - 2*ai*np.cos(gamma) *(Re+hD))## cos(gamma)= + to - at 90dg
        deltah= abs(h-hi) #(h_i-1 - hi)= delta h
        h=hi # new height
        rho = _getAirDensity(hi, "linsley")
        X=X+ rho * step*100. #(deltah*100) *abs(1./np.cos(np.pi-zen2)) # Xmax in g/cm2, slanted = Xmax, vertical/ cos(theta); density in g/cm3, h: m->100cm, np.pi-zen2 since it is defined as where the showers comes from, abs(cosine) so correct for minus values
        
    return h, ai # Xmax_height in m, Xmax_distance in m    
    
    
def _getAirDensity(_height,model):
    
    '''Returns the air density at a specific height, using either an isothermal model or the Linsley atmoshperic model as in ZHAireS

    Parameters:
    ---------
        h: float
            height in meters

    Returns:
    -------
        rho: float
            air density in g/cm3
    '''

    if model == "isothermal":
            #Using isothermal Model
            rho_0 = 1.225*0.001    #kg/m^3
            M = 0.028966    #kg/mol
            g = 9.81        #m.s^-2
            T = 288.        #
            R = 8.32        #J/K/mol , J=kg m2/s2
            rho = rho_0*np.exp(-g*M*_height/(R*T))  # kg/m3

    elif model == "linsley":
        #Using Linsey's Model
        bl = np.array([1222., 1144., 1305.5948, 540.1778,1])*10  # g/cm2  ==> kg/cm3
        cl = np.array([9941.8638, 8781.5355, 6361.4304, 7721.7016, 1e7])  #m
        hl = np.array([4,10,40,100,113])*1e3  #m
        if (_height>=hl[-1]):  # no more air
            rho = 0
        else:
            hlinf = np.array([0, 4,10,40,100])*1e3  #m
            ind = np.logical_and([_height>=hlinf],[_height<hl])[0]
            rho = bl[ind]/cl[ind]*np.exp(-_height/cl[ind])
            rho = rho[0]*0.001
    else:
        print("#### Error in GetDensity: model can only be isothermal or linsley.")
        return 0

    return rho

def getXmaxHeight(XmaxPosition):
    
    XmaxPosition = XmaxPosition
    Rearth = 6370949 
    XmaxHeight = np.sqrt((Rearth + XmaxPosition[2])**2 + XmaxPosition[0]**2 + XmaxPosition[1]**2) - Rearth

    return XmaxHeight
    
def getXmaxPosition(TargetShower):
    
    Xmax = Xmax_param(TargetShower.primary, TargetShower.energy)
    
    print("Xmax" , Xmax)
    InjectionHeight = TargetShower.injection
    zenith_inj = _get_CRzenith(TargetShower.zenith, InjectionHeight, 0) # TODO: insert injection height
    print(zenith_inj)
    XmaxHeight, DistDecayXmax = _dist_decay_Xmax(zenith_inj, InjectionHeight, Xmax)
    
    zenith = zenith_inj*np.pi/180
    azimuth = TargetShower.azimuth*np.pi/180
    uv = np.array([np.sin(zenith)*np.cos(azimuth), np.sin(zenith)*np.sin(azimuth) , np.cos(zenith)])
    
    Injection = [-InjectionHeight*uv[0]/uv[2], -InjectionHeight*uv[1]/uv[2], InjectionHeight] # TODO: modify !!!!!
    XmaxPosition = Injection - uv*DistDecayXmax
    
    print(Injection[0]/1e3)
    
    return XmaxPosition, XmaxHeight

def DensityScale(RefShower, TargetShower):
    
    Nant = RefShower.nant
    
    xmax_target, XmaxHeight_target = getXmaxPosition(TargetShower)
    
    XmaxHeight_ref = getXmaxHeight(RefShower.xmaxpos)
        
    rho_ref = _getAirDensity(XmaxHeight_ref, "linsley")
    rho_target = _getAirDensity(XmaxHeight_target, "linsley")
    
    
    krho =  np.sqrt(rho_ref/rho_target)
    
    scaled_traces = TargetShower.traces[:,2*Nant:3*Nant]*krho
    
    return scaled_traces, xmax_target, krho
    
# =============================================================================
#                       Cerenkov Stretch   
# =============================================================================

def get_cerenkov_angle(hXmax):
    
    # TODO: use an averaged refractive index --> oui
    
    n_refraction =  1 + 0.000325*np.exp(-0.1218e-3*hXmax)
    cer_ang = np.rad2deg(np.arccos(1/n_refraction))
    
    return cer_ang


def CerenkovStretch(RefShower, TargetShower):
    
    cerangle_ref = get_cerenkov_angle(RefShower.xmaxpos[2])
    cerangle_target = get_cerenkov_angle(TargetShower.xmaxpos[2])
    
    kstretch = cerangle_ref/cerangle_target
    
    scaled_traces = TargetShower.traces[:,176:]*kstretch
    scaled_pos = TargetShower.pos/kstretch
    
    # TODO: use the w angle 
    
    return scaled_pos, scaled_traces, kstretch
    

    
   




