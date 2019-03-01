import multiprocessing as mp
import itertools
import gc
import sys

import scipy.interpolate
import numpy as np

from micemag.fieldmanip import polarMeasurement as rphiz
from micemag.fieldmanip import cartesianMeasurement as xyz
import micemag.fieldmanip.fieldManipulation as fm
import micemag.utils as utils


class interpolatePoint:
    """Class to be used to obtain the components of the field at any given point by interpolation.


    """
    def __init__(self, field):
        if field[0].identifier() == 'Polar Data':
            field = fm.convertToCartesian(field)
        self._Bx, self._By, self._Bz = [], [], []
        self._points = []

        for f in field:
            self._Bx.append(f.Bx)
            self._By.append(f.By)
            self._Bz.append(f.Bz)
            self._points.append([f.x, f.y, f.z])

        self._Bx = np.array(self._Bx)
        self._By = np.array(self._By)
        self._Bz = np.array(self._Bz)
        self._points = np.array(self._points)

        
    def getBxyz(self, x, y, z):
        #use better intepolation algorithm?
        print 'Doing Bx'
        Bx = scipy.interpolate.griddata(self._points, self._Bx, np.array([x, y, z]), \
                                        method='linear', fill_value=0, rescale=True)
        print 'Doing By'
        By = scipy.interpolate.griddata(self._points, self._By, np.array([x, y, z]), \
                                        method='linear', fill_value=0, rescale=True)
        print 'Doing Bz'
        Bz = scipy.interpolate.griddata(self._points, self._Bz, np.array([x, y, z]), \
                                        method='linear', fill_value=0, rescale=True)

        return Bx, By, Bz

    def getBxyz2(self, x, y, z):
        jobs = []
        q = mp.Queue()
        for comp, j in zip(['Bx', 'By', 'Bz'], [self._Bx, self._By, self._Bz]):
            _work = _interpWorker(j, self._points, [x, y, z], q, comp)
            jobs.append(_work)
            _work.start()

        results = [q.get() for job in jobs]

        for i in jobs:
            i.join()

        
        _bx, _by, _bz = None, None, None
        for i in results:
            if 'Bx' in i.keys():
                _bx = i['Bx']
            if 'By' in i.keys():
                _by = i['By']
            if 'Bz' in i.keys():
                _bz = i['Bz']

        return _bx[0], _by[0], _bz[0]

        

        
@utils.timeit
def interpolate_field_2(field_values, field_points):
    interp_class = interpolatePoint(field_values)

    if field_points[0].identifier() == 'Polar Data':
        field_points = fm.convertToCartesian(field_points)
    xi = []
    for i in field_points:
        xi.append([i.x, i.y, i.z])

    result = []
    print 'Doing interpolation'
    
    for point in xi:
        print point[0], point[1], point[2]
        _bx, _by, _bz = interp_class.getBxyz2(point[0], point[1], point[2])
        result.append(xyz.Measurement(point[0], point[1], point[2], _bx, _by, _bz, -1))
    
    #_bx, _by, _bz = interp_class.getBxyz(xi[25][0], xi[25][1], xi[25][2])
    #result.append(xyz.Measurement(xi[25][0], xi[25][1], xi[25][2], _bx, _by, _bz))
    
    result = fm.convertToPolar(result)
    return result

@utils.timeit
def interpolate_field(field_values, field_points):
    """
    This function takes two fields in the form of Cartesian lists of measurements and interpolates
    field_values's magnetic field onto field_points's grid so the chi_square can be taken.  It then converts
    the result back to polar coordinates so that the chi-square difference can be taken.
        Input:
            -> field_points: The field map of the data of which contains the coordinates
                          where you would like to know the field of the field_values
            -> field_values: The rotated field map to interpolate onto field_points's grid
                           so you can compare
        Output:
            ->outputField: A list of POLAR measurements that has field_values's Bx,By,Bz,B
                           interpolated onto field_points's coordinates (x,y,z). (Then converted to (r,phi,z)
    It doesn't touch field_points's field components, only the coordinates.
    """
    if field_points[0].identifier() == 'Polar Data':
        print 'Converting field_points to cartesians'
        field_points = fm.convertToCartesian(field_points)
    if field_values[0].identifier() == 'Polar Data':
        print 'Converting field_values to cartesians'
        field_values = fm.convertToCartesian(field_values)

    
    Bx, By, Bz = [], [], []

    #points and field values for field_values
    points = []
    for f in field_values:
        Bx.append(f.Bx)
        By.append(f.By)
        Bz.append(f.Bz)
        points.append([f.x, f.y, f.z])

    #Points we want field_values at
    xi = []
    for i in field_points:
        xi.append([i.x, i.y, i.z])


    xi = np.array(xi)
    points = np.array(points)
    Bx = np.array(Bx)
    By = np.array(By)
    Bz = np.array(Bz)


    print 'Interpolating field with %d points onto field with %d points'%(len(field_values), len(field_points))
    print 'N points -> xi: %d, data: %d'%(len(xi), len(field_points))

    del field_points
    del field_values
    
    print 'Doing interpolation...'
    grid_Bx = scipy.interpolate.griddata(points, Bx, xi, method='linear', fill_value=99, rescale=True)
    print '-> Finished Bx'
    grid_By = scipy.interpolate.griddata(points, By, xi, method='linear', fill_value=99, rescale=True)
    print '-> Finished By'
    grid_Bz = scipy.interpolate.griddata(points, Bz, xi, method='linear', fill_value=99, rescale=True)
    print '-> Finished Bz'

    if len(xi) != len(grid_Bx):
        print 'Something went wrong!'
        return

    outputField = []
    for XYZ, Bx, By, Bz in itertools.izip(xi, grid_Bx, grid_By, grid_Bz):
        outputField.append(xyz.Measurement(XYZ[0], XYZ[1], XYZ[2], Bx, By, Bz, -1))
    
    outputField = fm.convertToPolar(outputField)
    outputField.sort()
    return outputField

@utils.timeit
def interpolate_field_parallel(field_values, field_points):
    """Interpolates field_values onto the coordinate system of field_points.

    Uses a subclass of multiprocessing.Process to parallelise the interpolation to speed things up.
    Note: DO NOT attempt to use this on Windows because it relies on copy-on-write to not use a 
          metric crap tonne of RAM (it already uses quite a lot).

    Args:
        field_values: A list of polarMeasurements. This is the field that you would like to interpolate
                      onto the coordinates given by field_points.

        field_points: A list of polarMeasurements. The field containing the coordinates that you would
                      like field_values interpolated onto. (These points should be coarser)
    
    Returns:
        The field described by field_values but on the data points of field_points. Returned as
        an array of polarMeasurements.

    """
    if field_points[0].identifier() == 'Polar Data':
        field_points = fm.convertToCartesian(field_points)
    if field_values[0].identifier() == 'Polar Data':
        field_values = fm.convertToCartesian(field_values)

    Bx = []
    By = []
    Bz = []
    fv_coords = []
    for f in field_values:
        Bx.append(f.Bx)
        By.append(f.By)
        Bz.append(f.Bz)
        fv_coords.append([f.x, f.y, f.z])

    points = []
    for i in field_points:
        points.append([i.x, i.y, i.z])

    Bx = np.array(Bx)
    By = np.array(By)
    Bz = np.array(Bz)
    fv_coords = np.array(fv_coords)

    del field_points
    del field_values
    gc.collect()
    
    jobs = []
    q = mp.Queue()
    for comp, j in zip(['Bx', 'By', 'Bz'], [Bx, By, Bz]):
        _work = _interpWorker(j, fv_coords, points, q, comp)
        jobs.append(_work)
        _work.start()

    try:
        results = [q.get() for job in jobs]

        for i in jobs:
            i.join()

        _bx, _by, _bz = None, None, None
        for i in results:
            if 'Bx' in i.keys():
                _bx = i['Bx']
            if 'By' in i.keys():
                _by = i['By']
            if 'Bz' in i.keys():
                _bz = i['Bz']

        del results

    except KeyboardInterrupt:
        print 'Terminating processes and exiting...'
        for i in jobs:
            print 'Terminating process: ', i.pid
            i.terminate()
        gc.collect()
        sys.exit()
                
    

    _lengths = [len(points), len(_bx), len(_by), len(_bz)]
    if len(np.unique(_lengths)) != 1:
        raise Exception('Points array is not the same length as (Bx, By, Bz) arrays!')
    
    outputField = []
    for XYZ, _Bx, _By, _Bz in itertools.izip(points, _bx, _by, _bz):
        outputField.append(xyz.Measurement(XYZ[0], XYZ[1], XYZ[2], _Bx, _By, _Bz, -1))
    
    outputField = fm.convertToPolar(outputField)
    outputField.sort()

    return outputField
    
        
class _interpWorker(mp.Process):
    def __init__(self, B, fv_coords, points, queue, comp):
        super(_interpWorker, self).__init__()
        self.queue = queue
        self.B = B
        self.fv_coords = fv_coords
        self.points = points
        self.comp = comp

    def run(self):
        print 'Interpolating ', self.comp, ' -- PID: ', self.pid
        newB = scipy.interpolate.griddata(self.fv_coords, self.B, self.points, method='linear', \
                                          fill_value=0, rescale=True)
        print 'Finished ', self.comp
        self.queue.put({self.comp: newB})
    

def form_residual_field():
    #This is the function used to interpolate the data field to make the residual field for the
    #fourier bessel fit.
    pass
