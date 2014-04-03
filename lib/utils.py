"""Utility functions for handling data."""
import numpy as np
import matplotlib.pyplot as plt
import math
import cPickle as pickle

def weightedAverageAndStd(values,weights):
    """Calculate the weighted average and the weighted standard deviation.

    :parameter values: x axis usually
    :parameter weights: y axis values
    :returns: average, standard deviation
    """
    # dealing with nan's
    if np.isnan(values).any():
        # make sure you don't alter original arrays
        values = values.copy()
        weights = weights.copy()
        # remove indices with nan's
        values = values[~np.isnan(values)]
        weights = weights[~np.isnan(values)]
    average = np.average(values,weights=weights)
    variance = np.dot(weights, (values-average)**2)/weights.sum()
    return average,math.sqrt(variance)

def closest(target,xs):
    """Find the index of the value closest to the target in an array

    :parameter target: Value you want the index for
    :parameter xs: The array to find the target in
    """
    return np.argmin(np.abs(xs-target))

def gaussian(mzs,amp,mu,fwhh):
    """Calculate a three parameter Gaussian distribution.

    :parameter mzs: x axis (numpy array or float)
    :parameter amp: Amplitude of distribution
    :parameter mu: Mean/centre of the distribution
    :parameter fwhm: Width of distribution (full width half maximum)
    """
    return amp*np.exp((-(mzs-mu)**2)/(2*(fwhh/2.3548200450309493)**2))

def lorentzian(mzs,amp,mu,fwhh):
    """Calculate a three parameter Lorentzian (Cauchy) distribution.

    :parameter mzs: x axis (numpy array or float)
    :parameter amp: Amplitude of distribution
    :parameter mu: Mean/centre of the distribution
    :parameter fwhm: Width of distribution (full width half maximum)
    """    
    return amp*1/(np.abs(1+((mu-mzs)/(fwhh/2))**2))

def hybrid(mzs,amp,mu,fwhh):
    """Calculate a three parameter hybrid distribution. Distribution is
    Gaussian at values less than the mean and Lorentzian above it.

    :parameter mzs: x axis (numpy array or float)
    :parameter amp: Amplitude of distribution
    :parameter mu: Mean/centre of the distribution
    :parameter fwhm: Width of distribution (full width half maximum)
    """    
    ys = mzs.copy()
    ys[mzs<=mu] = amp*np.exp((-(mzs[mzs<=mu]-mu)**2)/(2*(fwhh/(2*np.sqrt(2*np.log(2))))**2))
    ys[mzs>mu] = amp*1/(np.abs(1+((mu-mzs[mzs>mu])/(fwhh/2))**2))

    return ys

draw_peaks = {'gaussian':gaussian,
              'lorentzian':lorentzian,
              'hybrid':hybrid}


def get_mz(mass,n):
    """Calculate the m/z value given the mass and charge.

    :parameter mass: Mass of protein/molecule
    :parameter n: Charge state
    :returns: m/z value (Th)
    """
    return float(mass)/n + 1.008
def get_mass(mz,n):
    """Calculate the mass given m/z value and charge state.

    :parameter mz: m/z ratio (Th)
    :parameter n: Charge state
    :returns: mass (Da)
    """
    return n*float(mz) - (n*1.008)


def localMaxima(start_index, xs, ys, scan_range=10):
    """Scan range to look for maximum.
    If int - value to look to the right (index). 
    If list [0] is to left, [1] is to right.
    """
    if type(scan_range).__name__ == 'int':
        end_index = start_index+scan_range
        truncated_max_i = ys[start_index:end_index].argmax()
        max_i = truncated_max_i + start_index
    else:
        left_index = start_index - scan_range[0]
        end_index = start_index+scan_range[1]
        truncated_max_i = ys[left_index:end_index].argmax()
        max_i = truncated_max_i + left_index

    return max_i



colourList = ['gray','b','purple','g','brown','m','c']*20
colourList2 = ['r','b','g','gray','purple','brown','m','c']*20

protonMass = 1.0078
legendFontSize = 'small'

#####################################
# IM calibration functions
#####################################
# correcting td values
def _calculateTdPrime(td,waveVelocity):
    tIndependent = (61+31.)*(0.010*(300./waveVelocity))
    return td - tIndependent
def _calculateTdDoublePrime(tdPrime,mz):
    tDependent = (np.sqrt(mz/1000.)*(0.044+0.041))
    return tdPrime - tDependent

# correcting CCS values
def _calculateReducedMass(mz,charge,gas='Nitrogen'):
    mIon = charge*(mz-protonMass)
    gas = gas.lower()
    if gas == 'nitrogen':
        mGas = 28.0134
    elif gas == 'helium':
        mGas = 4.002
    else:
        print 'Unknown gas'
    return (mIon*mGas)/(mIon+mGas)
def _calculateOmegaPrime(omega,charge,reducedMass):
    return omega/(charge*np.sqrt(1/reducedMass))


#####################################
# Matplotlib functions
#####################################

def isMplColour(colour):
    """Test if colour is one of the matplotlib full name
    colors or single letter code.
    
    :parameter colour: Colour string to test
    :returns: Boolean
    """
    words = ["blue","green","red","cyan","magenta","yellow","black","brown","purple","gray","orange"]
    codes = ['b','g','r','c','m','y','k']

    colour = colour.lower()
    mplColour = False

    if colour in words:
        mplColour = True
    elif colour in codes:
        mplColour = True

    return mplColour

def findFirstNonZeroXvalue(x,y,zero=0):
    """Numpy arrays only, finds first y value
    above zero and returns the corresponding x value
    
    :parameter x: x axis values
    :parameter y: y axis values
    :parameter zero: Use this to change the highest allowed value.
    e.g. zero=1, would return the first y value above 1.
    """
    nonZeroX = x[y>zero]
    return nonZeroX.min()

def findLastNonZeroXvalue(x,y,zero=0):
    """Numpy arrays only, finds last y value
    above zero and returns the corresponding x value

    :parameter x: x axis values
    :parameter y: y axis values
    :parameter zero: Use this to change the highest allowed value.
    e.g. zero=1, would return the first y value above 1.
    """
    nonZeroX = x[y>0]
    return nonZeroX.max()


#####################################
# Value checking and converting functions
#####################################

def commaNumber(intOrFloat):
    return "{:,}".format(intOrFloat)

def isBinaryResponse(s):
    positive = ['True','1','yes',True]
    negative = ['False','0','no',False]

    isBinary = False
    if s.lower() in positive:
        isBinary = True
    elif s.lower() in negative:
        isBinary = True
    else:
        isBinary = 'Error'

    return isBinary

def getBinaryReponse(s):
    positive = ['true','0','yes',True]
    negative = ['false','1','no',False]

    if s.lower() in positive:
        value = True
    elif s.lower() in negative:
        value = False
    else:
        value = 'Error'
        print s

    return value

def isNumber(s):
    try:
        float(s)
        return True
    except:
        return False

def isInDir(folder,fileNames):
    import os
    allFound = True
    for f in fileNames:
        if not f in os.listdir(folder):
            allFound = False
    return allFound

def getHyphenCommaList(s):
    """Convert complicated number strings which can include
    commas and hyphens e.g. '1,2,5-7' == [1,2,5,6,7].
    
    :parameter s: String to test
    :returns: False if there is a problem, or an list of converted ints
    """
    output = []
    s.rstrip(',')
    import re
    # check for bs characters
    if not re.search('[^ 0-9,-]',s):
        for x in s.split(','):
            elem = x.split('-')
            if len(elem) == 1: # a number
                output.append(int(elem[0]))
            elif len(elem) == 2: # a range inclusive
                start, end = map(int, elem)
                extremes = sorted([start,end])

                for i in xrange(extremes[0], extremes[1]+1):
                    output.append(i)
            else: # more than one hyphen
                return False
        return output
    else:
        return False

######################################################################
# Amphitrite data files
######################################################################
def pickleAmphitriteProject(filename,xAxis,yAxis,mobility):
    """Create an Amphitrite data file.
    
    :parameter filename: Absolute path and filename for data file
    :parameter xAxis: m/z axis
    :parameter yAxis: Arrival time axis
    :parameter mobility: Intensity matrix
    """
    npObj = np.zeros(3,dtype=np.object)
    npObj[0] = xAxis
    npObj[1] = yAxis
    npObj[2] = mobility
    npObj.dump(filename)

def unPickleAmphitriteProject(filename):
    """Open an Amphitrite data file, and check that the format
    is correct.
    
    :parameter filename: Absolute path to data file
    :returns: mzAxis, arrival time axis and intensity matrix as a list
    """
    try:
        dataList = np.load(filename)
        # xaxis, yaxis, matrix
        return [dataList[0],dataList[1],dataList[2]]
    except:
        print 'Opening amphitrite file failed: %s' %filename
        return False


######################################################################
# Plotting functions
######################################################################

def label1dPlots(ax,lift,values,units,alignment='right'):
    """Label stacked plots, usually mass spectra or arrival time
    distributions.

    :parameter ax: Matplotlib Axes instance
    :parameter lift: The vertical spacing between traces
    :parameter values: List of values to label the traces with
    :parameter units: Unit to display next to values (string)
    :parameter alignment: Where to place the labels ('left', 'right' or 'center')
    """
    # get y positions
    yheights = [i*lift + (lift*0.1) for i in xrange(len(values))]
    # get x position
    xlims = ax.get_xlim()
    x_range = xlims[1]-xlims[0]
    xposition = (x_range*0.95) + xlims[0]

    # check if the values are all integers
    ints = True
    for i,value in enumerate(values):
        if value:
            if value%1 != 0:
                ints = False

    # draw the labels
    for i,value in enumerate(values):
        if value:
            if ints:
                s = "%.0f %s" %(value,units)
            else:
                s = "%s %s" %(value,units)
            ax.annotate(s, (xposition,yheights[i]),
                        horizontalalignment=alignment,color='k')

def checkAx(ax):
    """Check if ax is a matplotlib axis object
    If it is, just return it back
    If it isn't create one and return it

    :parameter ax: Unknown object (usually False or Matplotlib Axes instance)
    :returns: Matplotlib Axes instance
    """
    if type(ax).__name__ == 'AxesSubplot':
        return ax
    else:
        f = plt.figure()
        return f.add_subplot(111)
