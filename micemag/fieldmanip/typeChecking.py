import numpy as np



def fieldFromCurrentLoop(current, radius, R, Z):
    """
    Checks inputs to fieldFromCurrentLoop() for TypeErrors etc.
    """
    if type(current) != type(0.):
        raise TypeError("Current should be a float, "+str(type(current))+" detected.")

    if type(radius) != type(0.):
        raise TypeError("Radius should be a float, "+str(type(radius))+" detected.")

    if R.ndim != 2:
        raise IndexError("R should be a 2D gridded array, "+str(R.ndim)+" dimensions detected.")

    if Z.ndim != 2:
        raise IndexError("Z should be a 2D gridded array, "+str(Z.ndim)+" dimensions detected.")



def makeCurrentLayer(numLoops, separation, startLoopPosition, current, radius, R, Z):
    """
    Checks inputs to makeCurrentLayer() for TypeErrors etc.
    """
    if type(numLoops) != type(0):
        raise TypeError("numLoops should be an int, "+str(type(numLoops))+" detected.")

    if type(separation) != type (0.):
        raise TypeError("separation should be a float, "+str(type(separation))+" detected.")

    if type(startLoopPosition) != type (0.):
        raise TypeError("startLoopPosition should be a float, "+str(type(startLoopPosition))+" detected.")

    if type(current) != type(0.):
        raise TypeError("current should be a double, "+str(type(current))+" detected.")

    if type(radius) != type(0.):
        raise TypeError("radius should be a double, "+str(type(radius))+" detected.")

    if R.ndim != 2:
        raise IndexError("R should be a 2D gridded array, "+str(R.ndim)+" dimensions detected.")

    if Z.ndim != 2:
        raise IndexError("Z should be a 2D gridded array, "+str(Z.ndim)+" dimensions detected.")



def makeCoil(numLayers, numLoopsPerLayer, layerSeparation, loopSeparation, startPosition, current, minRadius, R, Z):
    """
    Checks inputs to makeCoil() for TypeErrors etc.
    """
    if type(numLayers) != type(0):
        raise TypeError("numLayers should be an int, "+str(type(numLayers))+" detected.")

    if type(numLoopsPerLayer) != type(0):
        raise TypeError("numLoopsPerLayer should be an int, "+str(type(numLoops))+" detected.")

    if type(layerSeparation) != type (0.):
        raise TypeError("layerSeparation should be a float, "+str(type(separation))+" detected.")

    if type(loopSeparation) != type (0.):
        raise TypeError("loopSeparation should be a float, "+str(type(separation))+" detected.")

    if type(startPosition) != type (0.):
        raise TypeError("startPosition should be a float, "+str(type(startPosition))+" detected.")

    if type(current) != type(0.):
        raise TypeError("current should be a double, "+str(type(current))+" detected.")

    if type(minRadius) != type(0.):
        raise TypeError("minRadius should be a double, "+str(type(radius))+" detected.")

    if R.ndim != 2:
        raise IndexError("R should be a 2D gridded array, "+str(R.ndim)+" dimensions detected.")

    if Z.ndim != 2:
        raise IndexError("Z should be a 2D gridded array, "+str(Z.ndim)+" dimensions detected.")






def makeMagnet(numCoils, numLayers, numLoopsPerLayer, layerSeparation, loopSeparation, startPosition, current, minRadius, R, Z):
    """
    Checks inputs to makeMagnet() for TypeErrors etc.
    """
    if type(numCoils) != type(0):
        raise TypeError("numLayers should be an int, "+str(type(numCoils))+" detected.")

    if type(numLayers[0]) != type(0):
        raise TypeError("numLayers should be a list of ints, "+str(type(numLayers[0]))+" detected.")

    if len(numLayers) != numCoils:
        raise IndexError("numLayers should be a list of length numCoils ("+str(numCoils)+"), but numLayers has length "+str(len(numLayers))+".")

    if type(numLoopsPerLayer[0]) != type(0):
        raise TypeError("numLoopsPerLayer should be a list of ints, "+str(type(numLoopsPerLayer[0]))+" detected.")

    if len(numLoopsPerLayer) != numCoils:
        raise IndexError("numLoopsPerLayer should be a list of length numCoils ("+str(numCoils)+"), but numLoopsPerLayer has length "+str(len(numLoopsPerLayer))+".")

    if type(layerSeparation[0]) != type (0.):
        raise TypeError("layerSeparation should be a list of floats, "+str(type(layerSeparation[0]))+" detected.")

    if len(layerSeparation) != numCoils:
        raise IndexError("layerSeparation should be a list of length numCoils ("+str(numCoils)+"), but layerSeparation has length "+str(len(layerSeparation))+".")

    if type(loopSeparation[0]) != type (0.):
        raise TypeError("loopSeparation should be a list of floats, "+str(type(loopSeparation[0]))+" detected.")

    if len(loopSeparation) != numCoils:
        raise IndexError("loopSeparation should be a list of length numCoils ("+str(numCoils)+"), but loopSeparation has length "+str(len(loopSeparation))+".")

    if type(startPosition[0]) != type (0.):
        raise TypeError("startPosition should be a list of floats, "+str(type(startPosition[0]))+" detected.")

    if len(startPosition) != numCoils:
        raise IndexError("startPosition should be a list of length numCoils ("+str(numCoils)+"), but startPosition has length "+str(len(startPosition))+".")

    if type(current[0]) != type(0.):
        raise TypeError("current should be a list of floats, "+str(type(current[0]))+" detected.")

    if len(current) != numCoils:
        raise IndexError("current should be a list of length numCoils ("+str(numCoils)+"), but current has length "+str(len(current))+".")

    if type(minRadius[0]) != type(0.):
        raise TypeError("minRadius should be a list of floats, "+str(type(radius[0]))+" detected.")

    if len(minRadius) != numCoils:
        raise IndexError("minRadius should be a list of length numCoils ("+str(numCoils)+"), but minRadius has length "+str(len(minRadius))+".")

    if R.ndim != 2:
        raise IndexError("R should be a 2D gridded array, "+str(R.ndim)+" dimensions detected.")

    if Z.ndim != 2:
        raise IndexError("Z should be a 2D gridded array, "+str(Z.ndim)+" dimensions detected.")





def calcFieldOnAxis(numLayers, numLoopsPerLayer, layerSeparation, loopSeparation, startPosition, current, minRadius, z):
    """
    Checks inputs to calcFieldOnAxis() for TypeErrors etc.
    """
    if type(numLayers) != type(0):
        raise TypeError("numLayers should be an int, "+str(type(numLayers))+" detected.")

    if type(numLoopsPerLayer) != type(0):
        raise TypeError("numLoopsPerLayer should be an int, "+str(type(numLoops))+" detected.")

    if type(layerSeparation) != type (0.):
        raise TypeError("layerSeparation should be a float, "+str(type(separation))+" detected.")

    if type(loopSeparation) != type (0.):
        raise TypeError("loopSeparation should be a float, "+str(type(separation))+" detected.")

    if type(startPosition) != type (0.):
        raise TypeError("startPosition should be a float, "+str(type(startPosition))+" detected.")

    if type(current) != type(0.):
        raise TypeError("current should be a double, "+str(type(current))+" detected.")

    if type(minRadius) != type(0.):
        raise TypeError("minRadius should be a double, "+str(type(radius))+" detected.")

    if z.ndim != 1:
        raise IndexError("z should be a 1D array (i.e. not gridded with r), "+str(z.ndim)+" dimensions detected.")




def calcMagnetFieldOnAxis(numCoils, numLayers, numLoopsPerLayer, layerSeparation, loopSeparation, startPosition, current, minRadius, z):
    """
    Checks inputs to calcMagnetFieldOnAxis() for TypeErrors etc.
    """
    if type(numCoils) != type(0):
        raise TypeError("numLayers should be an int, "+str(type(numCoils))+" detected.")

    if type(numLayers[0]) != type(0):
        raise TypeError("numLayers should be an int, "+str(type(numLayers[0]))+" detected.")

    if len(numLayers) != numCoils:
        raise IndexError("numLayers should be a list of length numCoils ("+str(numCoils)+"), but numLayers has length "+str(len(numLayers))+".")

    if type(numLoopsPerLayer[0]) != type(0):
        raise TypeError("numLoopsPerLayer should be an int, "+str(type(numLoops[0]))+" detected.")

    if len(numLoopsPerLayer) != numCoils:
        raise IndexError("numLoopsPerLayer should be a list of length numCoils ("+str(numCoils)+"), but numLoopsPerLayer has length "+str(len(numLoopsPerLayer))+".")

    if type(layerSeparation[0]) != type (0.):
        raise TypeError("layerSeparation should be a float, "+str(type(layerSeparation[0]))+" detected.")

    if len(layerSeparation) != numCoils:
        raise IndexError("layerSeparation should be a list of length numCoils ("+str(numCoils)+"), but layerSeparation has length "+str(len(layerSeparation))+".")

    if type(loopSeparation[0]) != type (0.):
        raise TypeError("loopSeparation should be a float, "+str(type(loopSeparation[0]))+" detected.")

    if len(loopSeparation) != numCoils:
        raise IndexError("loopSeparation should be a list of length numCoils ("+str(numCoils)+"), but loopSeparation has length "+str(len(loopSeparation))+".")

    if type(startPosition[0]) != type (0.):
        raise TypeError("startPosition should be a float, "+str(type(startPosition[0]))+" detected.")

    if len(startPosition) != numCoils:
        raise IndexError("startPosition should be a list of length numCoils ("+str(numCoils)+"), but startPosition has length "+str(len(startPosition))+".")

    if type(current[0]) != type(0.):
        raise TypeError("current should be a float, "+str(type(current[0]))+" detected.")

    if len(current) != numCoils:
        raise IndexError("current should be a list of length numCoils ("+str(numCoils)+"), but current has length "+str(len(current))+".")

    if type(minRadius[0]) != type(0.):
        raise TypeError("minRadius should be a float, "+str(type(minRadius[0]))+" detected.")

    if len(minRadius) != numCoils:
        raise IndexError("minRadius should be a list of length numCoils ("+str(numCoils)+"), but minRadius has length "+str(len(minRadius))+".")

    if z.ndim != 1:
        raise IndexError("z should be a 1D array (i.e. not gridded with r), "+str(z.ndim)+" dimensions detected.")




def printField(R, Z, Br, Bz, saveName, description):
    """
    Checks inputs to printField() for TypeErrors etc.
    """
    if R.ndim != 2:
        raise IndexError("R should be a 2D gridded array, "+str(R.ndim)+" dimensions detected.")

    if Z.ndim != 2:
        raise IndexError("Z should be a 2D gridded array, "+str(Z.ndim)+" dimensions detected.")

    if Br.ndim != 2:
        raise IndexError("Br should be a 2D gridded array, "+str(Br.ndim)+" dimensions detected.")

    if Bz.ndim != 2:
        raise IndexError("Bz should be a 2D gridded array, "+str(Bz.ndim)+" dimensions detected.")

    if type(saveName) != type("I am a string"):
        raise TypeError("saveName should be a string, "+str(type(saveName))+" detected.")

    if description != None:
        if type(description) != type("I am a string"):
            raise TypeError("description should be a string, "+str(type(description))+"detected.")




def readOriginalFiles(fileList, sensorList=None, surveyedOffsets=None, surveyedAngles=None):
    if type(fileList[0]) != type("I am a string"):
        raise TypeError("fileList should be a python-like list of strings, "+str(type(fileList[0]))+" detected.")


def readFile(fileName):
    if type(fileName) != type("I am a string"):
        raise TypeError("fileName should be a string, "+str(type(fileName))+" detected.")


def setSensorPosition(sensorNumber, xPosition, yPosition, phiRotation):
    if type(sensorNumber) != type(0):
        raise TypeException("sensorNumber must be an integer between 0 and 6, type "+str(type(sensorNumber))+" detected with value "+str(sensorNumber)+".")
    if type(xPosition) != type(0.):
        raise TypeException("xPosition should be a float, "+str(type(xPosition))+" detected.")
    if type(yPosition) != type(0.):
        raise TypeException("yPosition should be a float, "+str(type(yPosition))+" detected.")
    if type(phiRotation) != type(0.):
        raise TypeException("phiRotation should be a float, "+str(type(phiRotation))+" detected.")


def getSensorPosition(sensorNumber):
    if type(sensorNumber) != type(0):
        raise TypeError("sensorNumber should be an integer between 0 and 6, type "+string(type(sensorNumber))+" detected.")


def rotateMapperCoordinates(rotationAngle, sensorNumber, x_local, B_local):
    if type(rotationAngle) != type(0.0):
        raise TypeError("rotationAngle should be a float, "+str(type(rotationAngle))+" detected.")
    if type(sensorNumber) != type(0):
        raise TypeError("x_local should be an int between 0 and 6, "+str(type(x_local))+" detected.")
    if type(x_local) != type(0.0):
        raise TypeError("x_local should be a float, "+str(type(x_local))+" detected.")
    if type(B_local[0]) != type(0.0):
        raise TypeError("B_local should be a list of floats, "+str(type(B_local[0]))+" detected.")
    if type(B_local[1]) != type(0.0):
        raise TypeError("B_local should be a list of floats, "+str(type(B_local[1]))+" detected.")
    if type(B_local[2]) != type(0.0):
        raise TypeError("B_local should be a list of floats, "+str(type(B_local[2]))+" detected.")
    if len(B_local) != 3:
        raise TypeError("B_local should be a list of length 3, length "+str(len(B_local))+" detected.")


def rotateToSurveyCoordinates(x_mapper, B_mapper, offsets, angles):
    test = np.array([0.0, 0.0, 0.0])
    if type(x_mapper.dtype) != type(test.dtype):
        raise TypeError("x_mapper should have type numpy.float64, "+str(type(x_mapper.dtype))+" detected")
    if x_mapper.size != 3:
        raise TypeError("x_mapper should have three components, "+str(x_mapper.size)+" components detected")

    if type(B_mapper.dtype) != type(test.dtype):
        raise TypeError("B_mapper should have type numpy.float64, "+str(type(B_mapper.dtype))+" detected")
    if B_mapper.size != 3:
        raise TypeError("B_mapper should have three components, "+str(B_mapper.size)+" components detected")

    if type(offsets.dtype) != type(test.dtype):
        raise TypeError("offsets should have type numpy.float64, "+str(type(offsets.dtype))+" detected")
    if offsets.size != 3:
        raise TypeError("offsets should have three components, "+str(offsets.size)+" components detected")

    if type(angles.dtype) != type(test.dtype):
        raise TypeError("angles should have type numpy.float64, "+str(type(angles.dtype))+" detected")
    if angles.size != 3:
        raise TypeError("angles should have three components, "+str(angles.size)+" components detected")



def plotVariables(data, xAxisVariable, yAxisVariable, zAxisVariable, cutVariable, HallProbeList, xRange, yRange, zRange, cutRange):
    polarVariables = ['r', 'phi', 'z', 'Br', 'Bphi', 'Bz', 'B', 'probe', 't', 'date']
    cartesianVariables = ['x', 'y', 'z', 'Bx', 'By', 'Bz', 'B', 'probe', 't', 'date']

    # 1. Make sure axis variables are strings:
    if type(xAxisVariable) != type('string'):
        raise TypeError('x-axis variable should be a string, e.g. \'x\' or \'r\'. Type '+str(type(xAxisVariable))+' detected.')
    if type(yAxisVariable) != type('string') and yAxisVariable != None:
        raise TypeError('y-axis variable should be a string, e.g. \'x\' or \'r\'. Type '+str(type(yAxisVariable))+' detected.')
    if type(zAxisVariable) != type('string') and zAxisVariable != None:
        raise TypeError('z-axis variable should be a string, e.g. \'x\' or \'r\'. Type '+str(type(zAxisVariable))+' detected.')
    if type(cutVariable) != type('string') and cutVariable != None:
        raise TypeError('cut-variable should be a string, e.g. \'x\' or \'r\'. Type '+str(type(cutVariable))+' detected.')

    # 2. Make sure they're the *right strings for the data type*:
    if data[0].identifier() == 'Polar Data':
        if xAxisVariable not in polarVariables:
            raise TypeError("x-axis variable must be a valid Polar co-ordinate or field component: "+xAxisVariable+" was requested.")
        if yAxisVariable not in polarVariables and yAxisVariable != None:
            raise TypeError("y-axis variable must be a valid Polar co-ordinate or field component: "+yAxisVariable+" was requested.")
        if zAxisVariable not in polarVariables and zAxisVariable != None:
            raise TypeError("z-axis variable must be a valid Polar co-ordinate or field component: "+zAxisVariable+" was requested.")
        if cutVariable not in polarVariables and cutVariable != None:
            raise TypeError("cut-variable must be a valid Polar co-ordinate or field component: "+cutVariable+" was requested.")

    if data[0].identifier() == 'Cartesian Data':
        if xAxisVariable not in cartesianVariables:
            raise TypeError("x-axis variable must be a valid Cartesian co-ordinate or field component: "+xAxisVariable+" was requested.")
        if yAxisVariable not in cartesianVariables and yAxisVariable != None:
            raise TypeError("y-axis variable must be a valid Cartesian co-ordinate or field component: "+yAxisVariable+" was requested.")
        if zAxisVariable not in cartesianVariables and zAxisVariable != None:
            raise TypeError("z-axis variable must be a valid Cartesian co-ordinate or field component: "+zAxisVariable+" was requested.")
        if cutVariable not in cartesianVariables and cutVariable != None:
            raise TypeError("cut-variable must be a valid Cartesian co-ordinate or field component: "+cutVariable+" was requested.")

    # 3. Make sure we don't have a spurious number of hall probes being requested, that they're in the correct range and are all ints:
    if HallProbeList != None:
        if len(HallProbeList) > 7:
            raise TypeError("Too many entries in list of Hall probes. Maximum considered = 7, "+str(len(HallProbeList))+" requested.")
        for probe in range(0, len(HallProbeList)):
            if type(HallProbeList[probe]) != type(0):
                raise TypeError("Hall probe identifiers should be ints, probe "+str(probe)+" in list is of type "+str(type(HallProbeList[probe]))+".")
            if probe > 6 or probe < 0:
                raise TypeError("Hall probe "+str(probe)+" in list is out of range. Valid probe ID's are 0..6, but "+str(HallProbeList[probe])+" was given.")


    # 4. Make sure xRange, yRange, zRange are all sensible:
    if xRange != None:
        if len(xRange) != 2:
            raise TypeError("xRange is specified as [min, max], list of length "+str(len(xRange))+" detected.")
        for x in range(0, len(xRange)):
            if xAxisVariable != 'probe' and xAxisVariable != 't' and xAxisVariable != 'date':
                # should be using floats:
                if type(xRange[x]) != type(0.0):
                    raise TypeError("xRange should be a list of floats, xRange["+str(x)+"] is of type "+str(type(xRange[x]))+".")
            else:
                # should be using ints:
                if type(xRange[x]) != type(0):
                    raise TypeError("xRange should be a list of ints, xRange["+str(x)+"] is of type "+str(type(xRange[x]))+".")

        # Finally, check that min < max:
        if xRange[0] > xRange[1]:
            raise TypeError("xRange should be specified as [min, max], but xRange[0] > xRange[1]: ("+str(xRange[0])+", "+str(xRange[1])+") given.")


    if yRange != None:
        if len(yRange) != 2:
            raise TypeError("yRange is specified as [min, max], list of length "+str(len(yRange))+" detected.")
        for y in range(0, len(yRange)):
            if yAxisVariable != 'probe' and xAxisVariable != 't' and xAxisVariable != 'date':
                # should be using floats:
                if type(xRange[y]) != type(0.0):
                    raise TypeError("yRange should be a list of floats, yRange["+str(y)+"] is of type "+str(type(yRange[y]))+".")
            else:
                # should be using ints:
                if type(yRange[y]) != type(0):
                    raise TypeError("yRange should be a list of ints, yRange["+str(y)+"] is of type "+str(type(yRange[y]))+".")

        # Finally, check that min < max:
        if yRange[0] > yRange[1]:
            raise TypeError("yRange should be specified as [min, max], but yRange[0] > yRange[1]: ("+str(yRange[0])+", "+str(yRange[1])+") given.")


    if zRange != None:
        if len(zRange) != 2:
            raise TypeError("zRange is specified as [min, max], list of length "+str(len(zRange))+" detected.")
        for z in range(0, len(zRange)):
            if zAxisVariable != 'probe' and zAxisVariable != 't' and zAxisVariable != 'date':
                # should be using floats:
                if type(zRange[z]) != type(0.0):
                    raise TypeError("zRange should be a list of floats, zRange["+str(z)+"] is of type "+str(type(zRange[z]))+".")
            else:
                # should be using ints:
                if type(zRange[z]) != type(0):
                    raise TypeError("zRange should be a list of ints, zRange["+str(z)+"] is of type "+str(type(zRange[z]))+".")

        # Finally, check that min < max:
        if zRange[0] > zRange[1]:
            raise TypeError("zRange should be specified as [min, max], but zRange[0] > zRange[1]: ("+str(zRange[0])+", "+str(zRange[1])+") given.")


    if cutRange != None:
        if len(cutRange) != 2:
            raise TypeError("cutRange is specified as [min, max], list of length "+str(len(cutRange))+" detected.")
        for z in range(0, len(cutRange)):
            if cutVariable != 'probe' and cutVariable != 't' and cutVariable != 'date':
                # should be using floats:
                if type(cutRange[z]) != type(0.0):
                    raise TypeError("cutRange should be a list of floats, cutRange["+str(z)+"] is of type "+str(type(cutRange[z]))+".")
            else:
                # should be using ints:
                if type(cutRange[z]) != type(0):
                    raise TypeError("cutRange should be a list of ints, cutRange["+str(z)+"] is of type "+str(type(cutRange[z]))+".")

        # Finally, check that min < max:
        if cutRange[0] > cutRange[1]:
            raise TypeError("cutRange should be specified as [min, max], but cutRange[0] > cutRange[1]: ("+str(cutRange[0])+", "+str(cutRange[1])+") given.")










