import numpy as np

class Measurement:
    def __init__(self, r, phi, z, Br, Bphi, Bz, sensorNumber=None, date=None, time=None):
        """Class that stores cylindrical polar field data.

        This class stores the three cylindrical polar coordinates (r,phi,z) along with the
        corresponding field components (Br,Bphi,Bz) at a point in space.  It calculates |B| from these
        field components and stores this value.
        If real field mapper data is used, this class can also hold the date and time the data were
        recorded along with the sensor number the data is recorded by.
        If the phi coordinate is negative, it will be added to 360.0 so that phi always runs from 0
        to 360 degrees.
        Coordinates passed to this class should conventionally be in metres!
        The phi coordinate must be in degrees!
        Conventionally, the field components are in units of Tesla.

        Args:
            r    (float): The r coordinate of the point.
            phi  (float): The phi coordinate of the point (in degrees).
            z    (float): The z coordinate of the point.
            Br   (float): The Br component of the field.
            Bphi (float): The Bphi component of the field.
            Bz   (float): The Bz component of the field.
            sensorNumber (int, optional): The number of the Hall sensor that measured the data point.
                Defaults to None but is set to -1 if this is the case. Should be between 0-6.
            date (int, optional): The date that this data point was recorded. This is superfluous
                for normal work but is included here for continuity. Format is an integer - DDMMYYYY.
            time (int, optional): The time that this data point was recorded. Similar arguments to
                date. The format is HHMMSS.

        """
        if phi < 0.0:
            self.phi = np.float64(360.0 + phi)
        else:
            self.phi = np.float64(phi)
        self.r = np.float64(r)
        self.z = np.float64(z)

        self.Br = np.float64(Br)
        self.Bphi = np.float64(Bphi)
        self.Bz = np.float64(Bz)
        
        self.B = np.sqrt(self.Br**2 + self.Bphi**2 + self.Bz**2)
        self.ID = "Polar Data"

        if sensorNumber is None:
            self.sensorNumber = -1
        else:
            self.sensorNumber = sensorNumber

        if date is None:
            self.Date = 22052013
        else:
            self.Date = date

        if time is None:
            self.Time = 120000
        else:
            self.Time = time


    def __lt__(self, other):
        return (np.around(self.z, 3), np.around(self.r, 2), np.around(self.phi, 0)) < (np.around(other.z, 3), np.around(other.r, 2), np.around(other.phi, 0))
        
    
    def __add__(self, other):
        return Measurement(self.r, self.phi, self.z, self.Br + other.Br, self.Bphi + other.Bphi, self.Bz + other.Bz)

    
    def __str__(self):
        return '(z = %5.10f, r = %5.10f, phi = %5.10f, Br = %5.10f, Bphi = %5.10f, Bz = %5.10f, |B| = %5.10f, sensor = %d)' % (self.z, self.r, self.phi, self.Br, self.Bphi, self.Bz, self.B, self.sensorNumber)

    def __repr__(self):
        return '<' + self.__str__() + '>\n'

    def asFileLine(self):
        """Express this data point as a string that is suitable to print to a file."""
        return '%5.10f\t%5.10f\t%5.10f\t%5.10f\t%5.10f\t%5.10f\t%5.10f\t%d\t%d\t%s\n' % (self.z, self.r, self.phi, self.Br, self.Bphi, self.Bz, self.B, self.sensorNumber, self.Date, self.Time)

    def identifier(self):
        """Returns the identifier for this class.

        This class will return the string -- Polar Data.  This function is used to distinguish
        between Cartesian Measurements and polar Measurements.

        """
        return '%s' % (self.ID)

    def setData(self, date):
        """Sets the date of this instance."""
        self.Date = date

    def setTime(self, time):
        """Sets the time of this instance."""
        self.Time = time

    def date(self):
        """Returns the date of this instance."""
        return '%d' % (self.Date)

    def time(self):
        """Returns the time of this instance."""
        return '%s' % (self.Time)

    def set_RPhiZ(self, r, phi, z):
         """Set the r,phi,z coordinates of this instance.
        
         Args:
            r   (float): The desired r coordinate
            phi (float): The desired phi coordinate
            z   (float): The desired z coordinate

         """
         self.r = r
         self.phi = phi
         self.z = z
         
