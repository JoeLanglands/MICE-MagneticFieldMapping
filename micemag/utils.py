import os.path
import time
import sys

"""
Stores useful paths and constants.

These can then be easily accessed by importing this module. 
"""

#Define the data paths for ease of access
data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
fb_pickle_path = os.path.join(data_path, 'FBpickles')
resid_field_path = os.path.join(data_path, 'residFields')
SSD_data_path = os.path.join(data_path, 'SSD')
SSU_data_path = os.path.join(data_path, 'SSU')
geofit_field_path = os.path.join(data_path, 'geoFitFields')
maus_field_path = os.path.join(data_path, 'MAUS')

ssu_centre_dict = {'survey': {'M1': 0, 'M2': 0, 'E1': 0, 'CC': 0, 'E2': 0}, \
                   'mapper': {'M1':3.448, 'M2': 3.01, 'E1': 2.615, 'CC': 1.868, 'E2': 1.121, 'ECE': 1.868}}

#SSD M1 is broken
ssd_centre_dict = {'survey': {'M1': 0, 'M2': 0, 'E1': 0, 'CC': 0, 'E2': 0}, \
                   'mapper': {'M1': 0, 'M2': 1.841, 'E1': 2.2413, 'CC': 2.9811, 'E2': 3.734, 'ECE': 2.9811}}

#This is for the raw unflipped data.  SSD data should be passed through the flip_SSD_data function
#to make its coordinate system  parallel with the MICE coord system
ssd_centre_dict_RAW = {'survey': {'M1': 0, 'M2': 0, 'E1': 0, 'CC': 0, 'E2': 0}, \
                       'mapper': {'M1': 0, 'M2': 3.159, 'E1': 2.7587, 'CC': 2.0189, 'E2': 1.266}}

centres_dict = {'SSU': ssu_centre_dict, 'SSD': ssd_centre_dict}

#Default grid for g4bl field format
default_grid = {'x': {'start': 0, 'end': 0.18, 'step': 0.01}, \
                'y': {'start': 0, 'end': 0.18, 'step': 0.01}, \
                'z': {'start': 0, 'end': 4.99, 'step': 0.01}}


#These functions have no true home and are purely utility

def timeit(func): #Decorator so that functions can be timed
    def timed(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()

        _mins, secs = divmod(end_time - start_time, 60)
        hours, mins = divmod(_mins, 60)
        
        print '%r took %d hrs, %d mins, %d secs'%(func.__name__, int(hours), \
                                                  int(mins), int(secs))
        
        return result
    return timed


def progressBar(n, total, start_time=None, end_time=None):
    pcDone = int((float(n)/total)*100)/5
    progress = ''
 
    remaining = total - n
    time_for_one = (end_time - start_time)/n
    time_left = time_for_one*remaining
    etamin, etasec = divmod(time_left, 60)
    if pcDone > 0:
        for i in range(0, pcDone):
            progress += '#'
    if start_time == None:
        progString = ' [%s] %d / %d\r'%(progress.ljust(20), n, total)
    else:
        elapsed = end_time - start_time
        mins, secs = divmod(elapsed, 60)
        progString = ' [%s] %d / %d -- Elapsed: %02d:%02d  ETA: %02d:%02d\r'%(progress.ljust(20), n, \
                                                                          total, int(mins),int(secs), \
                                                                          etamin, etasec)
    if n < total:
        sys.stdout.write(progString)
        sys.stdout.flush()
    elif n == total:
        print progString


"""
The nested dictionaries here contain all the basic 'as-built' physical information about the coils:
layers/turns of conductor, length, rInner, rOuter etc...
These parameters are used as starting points for the geometrical/coil-model fit and are used to 
calculate the bracketing fields.

They do not contain information about the centre position since that can change between mice and
mapper coordinates. See dictionaries above.
"""


coil_datacards = {
    'SSU': {
        'M1': {
            'length': 0.2006,
            'rInner': 0.25704,
            'rOuter': 0.30159,
            'nLayers': 42,
            'nTurns': 115,
            '30A_data': os.path.join(SSU_data_path, 'run07_polarCoordinates.dat')
        },
        'M2': {
            'length': 0.1989,
            'rInner': 0.25698,
            'rOuter': 0.28667,
            'nLayers': 28,
            'nTurns': 114,
            '30A_data': os.path.join(SSU_data_path, 'run05_polarCoordinates.dat')
        },
        'E1': {
            'length': 0.1102,
            'rInner': 0.25709,
            'rOuter': 0.3165,
            'nLayers': 56,
            'nTurns': 64,
            '30A_data': None #New data has only ECE powered together
        },
        'CC': {
            'length': 1.3101,
            'rInner': 0.25712,
            'rOuter': 0.27835,
            'nLayers': 20,
            'nTurns': 768,
            '30A_data': os.path.join(SSU_data_path, 'run03_polarCoordinates.dat') #ECE run
        },
        'E2': {
            'length': 0.1102,
            'rInner': 0.25712,
            'rOuter': 0.3229,
            'nLayers': 62,
            'nTurns': 64,
            '30A_data': None #New data has only ECE powered together
        },
        'ECE': {
            '30A_data': os.path.join(SSU_data_path, 'run03_polarCoordinates.dat')
        }
    }, #End SSU
    'SSD': {
         'M1': {
            'length': 0.201268,
            'rInner': 0.258,
            'rOuter': 0.304483,
            'nLayers': 42,
            'nTurns': 115,
            '30A_data': None #Broken of course -> No data
        },
        'M2': {
            'length': 0.199492,
            'rInner': 0.258,
            'rOuter': 0.288608,
            'nLayers': 28,
            'nTurns': 114,
            '30A_data':  os.path.join(SSD_data_path, 'run05_polarCoordinates.dat') 
        },
        'E1': {
            'length': 0.110642,
            'rInner': 0.258,
            'rOuter': 0.319638,
            'nLayers': 56,
            'nTurns': 64,
            '30A_data':  None #same as above for SSU
        },
        'CC': {
            'length': 1.3143,
            'rInner': 0.258,
            'rOuter': 0.280416,
            'nLayers': 20,
            'nTurns': 768,
            '30A_data': os.path.join(SSD_data_path, 'run03_polarCoordinates.dat') 
        },
        'E2': {
            'length': 0.110642,
            'rInner': 0.258,
            'rOuter': 0.326220,
            'nLayers': 62,
            'nTurns': 64,
            '30A_data':  None #same as above for SSU
        },
        'ECE': {
            '30A_data': os.path.join(SSD_data_path, 'run03_polarCoordinates.dat')
        }
    } #End SSD
}

