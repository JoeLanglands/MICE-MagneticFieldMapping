import sys
import os

import numpy as np

import typeChecking as check
import cartesianMeasurement as xyz
import polarMeasurement as rphiz



def readFile(fileName):
    """
    Reads a tab-separated data file such as those created by readDataFromMapper.py or calculateField.py

    Depends on:
        - measurement.py
        - fieldFormat.py

    Sets:
        - Global array of measurements (co-ordinate system corrections applied or not)

    Inputs:
        - fileName: The name and location of the data file to be read
        - getMeasurements: (Optional) Determines whether or not the function immediately returns the read values (as well as storing them in the global variables)
            - True: Returns values and sets global variables
            - False: Only sets global variables
            - Default value is False

    Returns:
        - Returns measurements in cartesian or polar data format.
            - If data file was in Cartesian co-ordinates, returns variables in Cartesian co-ordinates
            - If data file was in cylindrical co-ordinates, returns variables in cylindrical co-ordinates.
            - NB: The returned values require the class cartesianMeasurement.py (Cartesian) or polarMeasurement.py (cylindrical) to interpret
    """
    check.readFile(fileName)

    measurements = []

    f = open(fileName, 'r')
    header = f.readline()
    if 'x' in header:
        # Data in Cartesian co-ordinates
        z, x, y, Bx, By, Bz, B, sensor, date, time = np.loadtxt(fileName, skiprows=1, unpack=True)

        for entry in range(0, z.size):
            measurements.append(xyz.Measurement(x[entry], y[entry], z[entry],
                                Bx[entry], By[entry], Bz[entry], sensor[entry], date[entry], time[entry]))



    elif 'phi' in header:
        # Data in cylindrical polar co-ordinates
        z, r, phi, Br, Bphi, Bz, B, sensor, date, time = np.loadtxt(fileName, skiprows=1, unpack=True)

        for entry in range(0, z.size):
            # no phi components in these files, so set phi = 0 and Bphi = 0
            measurements.append(rphiz.Measurement(r[entry], phi[entry], z[entry], Br[entry],
                                Bphi[entry], Bz[entry], sensor[entry], date[entry], time[entry]))
    else:
        print 'Unknown file headers. Acceptable headers are:'
        print '\t(1) #z (m)\tx (m)\ty (m) \tBx (T)\tBy (T)\tBz (T)\tB (T)\tprobeID\tDate (DDMMYYY)\tTime (24-HH:MM:SS)\n'
        print '\t(2) #z (m)\tr (m)\tphi (deg) \tBr (T)\tBphi (T)\tBz (T)\tB (T)\tprobeID\tDate (DDMMYYY)\tTime (24-HH:MM:SS)\n'
        print 'Header input was:'
        print '\t', header

    f.close()

    # return the measurements after sorting them in ascending z
    measurements.sort()
    return measurements
