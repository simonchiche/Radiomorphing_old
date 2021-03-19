import numpy as np
import glob
#import grand_radiomorphing
from coreRadiomorphing import process


def run():
    # Settings of the radiomorphing
    
    # reference path
    #data_dir = join(root_dir, "examples", "data")
    # folder containing your refernence shower simulations
    sim_dir = glob.glob("./Simulations/*")
    # folder which will contain radio morphed traces afterwards
    out_dir = glob.glob("./OutputDirectory")
    # list of antenna positions you would like to simulate, stored in out_dir in the best case
    #antennas = glob.glob("./DesiredPositions/AntennasCoordinates.txt") 


    # definition of target shower parameters
    shower = {
        "primary" : "Iron",        # primary (electron, pion)
        "energy" : 3.9811,               # EeV
        "zenith" : 117.02,               # deg (GRAND frame)
        "azimuth" : 270.0,                # deg (GRAND frame)
        "injection" : 1e5,    # m (injection height in the local coordinate system)
        "altitude" : 1000. }   # m (alitude oj injection with respect to sealevel, 
                               #not necessarily eqivalent to injection height)
    # Perform the radiomorphing
    #grand_radiomorphing.process(sim_dir, shower, antennas, out_dir)
    process(sim_dir, shower, out_dir)



if __name__ == "__main__":
    run()
