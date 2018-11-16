import os.path
import time
import sys

"""
Stores useful paths and constants.

These can then be easily accessed by importing this module. 
e.g. import micemag.utils
or   from micemag.utils import *

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
                   'mapper': {'M1':3.448, 'M2': 3.01, 'E1': 2.615, 'CC': 1.868, 'E2': 1.121}}
#SSD M1 is broken
ssd_centre_dict = {'survey': {'M1': 0, 'M2': 0, 'E1': 0, 'CC': 0, 'E2': 0}, \
                   'mapper': {'M1': 0, 'M2': 3.1586, 'E1': 2.759, 'CC': 2.019, 'E2': 1.266}}

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


def progressBar(n, total):
    pcDone = int((float(n)/total)*100)/5
    progress = ''
    if pcDone > 0:
        for i in range(0, pcDone):
            progress += '#'
    progString = ' [%s] %d / %d\r'%(progress.ljust(20), n, total)
    sys.stdout.write(progString)
    if n != total:
        sys.stdout.flush()
