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
