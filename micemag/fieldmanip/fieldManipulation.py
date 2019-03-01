
import numpy as np
import scipy as sp

import cartesianMeasurement as xyz
import polarMeasurement as rphiz


def scaleField(field, lengthScale, fieldScale):
    """Scale the field strength/coordinates by two separate factors.

    Scales the coordinates and components of a  field, consisting to an array of the Measurement
    class (can be cartesian or polar).
    Scales each coordinate by the argument lengthScale and each component by fieldScale.
    Example x -> x*lengthScale, Bx -> Bx*fieldScale.
    If polar coordinates are used, phi is left untouched.
    Preserves and probe number information.
    
    Args:
        field (list of Measurement class or string): Array of measurement points that are to be scaled.
        lengthScale (float): Scale factor for the coordinates.
        fieldScale (float): Scale factor for the field components.

    Returns:
        A scaled field as a list of measurements. Returns whatever type of field was passed to it.

    """
    
    scaledField = []
    if field.identifier=='Polar Data':
        for f in field:
            z = f.z * lengthScale
            r = f.r * lengthScale
            phi = f.phi
            Br = f.Br * fieldScale
            Bphi = f.Bphi * fieldScale
            Bz = f.Bz * fieldScale

            scaledField.append(rphiz.Measurement(r, phi, z, Br, Bphi, Bz, f.sensorNumber))

    elif field.identifier=='Cartesian Data':
        for f in field:
            x = f.x * lengthScale
            y = f.y * lengthScale
            z = f.z * lengthScale
            Bx = f.Bx * fieldScale
            By = f.By * fieldScale
            Bz = f.Bz * fieldScale

            scaledField.append(xyz.Measurement(x, y, z, Bx, By, Bz, f.sensorNumber))
        
    return scaledField

def shiftField(field, dz):
    """Shifts the z-coordinate of the field by dz"""
    for f in field:
        if f.ID == 'Polar Data':
            f.set_RPhiZ(f.r, f.phi, f.z + dz)
        elif f.ID == 'Cartesian Data':
            f.set_XYZ(f.x, f.y, f.z + dz)

    return field #Don't need this really because it changes the actual values

def convertToCartesian(fieldMeasurements):
    """
    Convert a list of polarMeasurement's to cartesianMeasurements
    """

    newCoordinates = []
    for item in fieldMeasurements:
        r = item.r
        phi = np.radians(item.phi)
        z = item.z

        Br = item.Br
        Bphi = item.Bphi


        x = r*np.cos(phi)
        y = r*np.sin(phi)

        rotate = np.array([[np.cos(phi), -1.0*np.sin(phi)], [np.sin(phi), np.cos(phi)]])
        R = np.array([Br, Bphi])
        X = rotate.dot(R)

        Bx = X[0]
        By = X[1]
        Bz = item.Bz

        newCoordinates.append(xyz.Measurement(x, y, z, Bx, By, Bz, item.sensorNumber, item.Date, item.Time))

    newCoordinates.sort()
    return newCoordinates

def convertToPolar(fitted_fieldMeasurements):
    """
    Convert a list of cartesianMeasurements to polarMeasurements
    """
    polarData = []
    for item in fitted_fieldMeasurements:
        x = item.x
        y = item.y
        z = item.z

        Bx = item.Bx
        By = item.By
        Bz = item.Bz

        phi = np.arctan2(y, x)
        r = np.sqrt(x**2 + y**2)

        rotate = np.array([[np.cos(phi), -1.0*np.sin(phi)], [np.sin(phi), np.cos(phi)]])
        rotate = np.transpose(rotate)
        X = np.array([Bx, By])
        R = rotate.dot(X)

        Br = R[0]
        Bphi = R[1]

        phi_in_degrees = np.degrees(phi)
        if phi_in_degrees < 0.0:
            phi_in_degrees = 360.0 + phi_in_degrees
        elif phi_in_degrees > 360.0:
            phi_in_degrees = 360.0 - phi_in_degrees

        polarData.append(rphiz.Measurement(r, phi_in_degrees, z, Br, Bphi, Bz, item.sensorNumber))
    polarData.sort()
    return polarData


def flip_SSD_data(fieldMeasurements):
    """Function that flips the SSD data so its coordinate system is orientated like the MICE/SSU data.

    Since the field mapping machine went downstream to upstream for the SSD field map data, the y-axis
    and z-axis are reversed.  This function rectifies this so that when the field is fitted and made
    the coordinate system is the same as the MICE coordinate system.
    """
    if fieldMeasurements[0].ID == 'Polar Data':
        fieldMeasurements = convertToCartesian(fieldMeasurements)

    resultField = []
        
    for f in fieldMeasurements:
        #with a (z+5.0) to make it on the positive Z axis
        resultField.append(xyz.Measurement(f.x, -1.0*f.y, -1.0*f.z + 5.0, f.Bx, \
                                           -1.0*f.By, -1.0*f.Bz, f.sensorNumber))
    
    resultField = convertToPolar(resultField)

    return resultField

def print_field_from_list(listOfMeasurements, saveName, description=None):
    with open(saveName, 'w') as f:
        if listOfMeasurements[0].identifier() == 'Polar Data':
            f.write('#z (m)\tr (m)\tphi (deg) \tBr (T)\tBphi (T)\tBz (T)\tB (T)\tprobeID\tDate (DDMMYYY)\tTime (24-HH:MM:SS)\n')
        else:
            f.write('#z (m)\tx (m)\ty (m) \tBx (T)\tBy (T)\tBz (T)\tB (T)\tprobeID\tDate (DDMMYYY)\tTime (24-HH:MM:SS)\n')

        listOfMeasurements.sort()

        for m in listOfMeasurements:
            f.write(m.asFileLine())

        if description != None:
            # add the description to the bottom of the file
            f.write('\n')
            f.write(description)
