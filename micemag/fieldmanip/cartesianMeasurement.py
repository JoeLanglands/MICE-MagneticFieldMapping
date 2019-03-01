import numpy as np

class Measurement:
    __slots__ = ['x', 'y', 'z', 'Bx', 'By', 'Bz', 'B', 'sensorNumber', 'Date', 'Time', 'ID']
    
    def __init__(self, x, y, z, Bx, By, Bz, sensorNumber=None, date=None, time=None):
        """Class that stores cartesian field data.

        This class stores the three coordinates (x,y,z) along with the corresponding field components 
        (Bx, By, Bz) at a point in space.  It calculates |B| from these field components and stores
        this value.
        If real field mapper data is used, this class can also hold the date and time the data were
        recorded along with the sensor number the data is recorded by.
        Coordinates passed to this class should conventionally be in metres!
        Conventionally, the field components are in units of Tesla.

        Args:
            x  (float): The x coordinate of the point.
            y  (float): The y coordinate of the point.
            z  (float): The z coordinate of the point.
            Bx (float): The Bx component of the magnetic field.
            By (float): The By component of the magnetic field.
            Bz (float): The Bz component of the magnetic field.
            sensorNumber (int, optional): The number of the Hall sensor that measured the data point.
                Defaults to None but is set to -1 if this is the case. Should be between 0-6.
            date (int, optional): The date that this data point was recorded. This is superfluous
                for normal work but is included here for continuity. Format is an integer - DDMMYYYY.
            time (int, optional): The time that this data point was recorded. Similar arguments to
                date. The format is HHMMSS.        

        """

        self.x = np.float64(x)
        self.y = np.float64(y)
        self.z = np.float64(z)

        self.Bx = np.float64(Bx)
        self.By = np.float64(By)
        self.Bz = np.float64(Bz)

        self.B = np.sqrt(self.Bx**2 + self.By**2 + self.Bz**2)

        if sensorNumber is None:
            self.sensorNumber = -1
        else:
            self.sensorNumber = sensorNumber

        self.ID = "Cartesian Data"

        if date is None:
            self.Date = 00000000
        else:
            self.Date = date

        if time is None:
            self.Time = 120000
        else:
            self.Time = time

    def __lt__(self, other):
        return (self.z, self.x, self.y) < (other.z, other.x, other.y)

    def __str__(self):
        return '(z = %5.10f, x = %5.10f, y = %5.10f, Bx = %5.10f, By = %5.10f, Bz = %5.10f, |B| = %5.10f, sensor = %d)' % (self.z, self.x, self.y, self.Bx, self.By, self.Bz, self.B, self.sensorNumber)

    def __repr__(self):
        return '<' + self.__str__() + '>\n'

    def asFileLine(self):
        """Express this data point as a string that is suitable to print to a file."""
        return '%5.10f\t%5.10f\t%5.10f\t%5.10f\t%5.10f\t%5.10f\t%5.10f\t%d\t%d\t%s\n' % (self.z, self.x, self.y, self.Bx, self.By, self.Bz, self.B, self.sensorNumber, self.Date, self.Time)

    def identifier(self):
        return self.ID

    def set_XYZ(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
