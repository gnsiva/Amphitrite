import numpy as np
import matplotlib.pyplot as plt
from lib import utils
import lib.SG as SG
from collections import OrderedDict

class TwoDdata():

    def __init__(self):
        self.xvals = np.array([])
        self.yvals = np.array([])

        self.rawyvals = []
        self.gradient = []
        self.gPeaks = OrderedDict()

        self.normalisationType = 'none'

    #===========================================================================
    # Data manipulation
    #===========================================================================
    def readFile(self,filename,x_range=0,grain=1):
        '''Reads in x y coordinate pairs from text file
        ' ' separator as in copy spectrum list in MassLynx

        x_range - allows you to select lower and upper bounds
        in the format of [lower,upper]

        grain - allows the missing of data to speed up processing
        a grain of 2 means that every second value will be used'''
        raw_data = open(filename,'r').readlines()
        count = 0
        self.xvals = []
        self.yvals = []
        for x in raw_data:
            count += 1
            if count == grain:
                temp = x.rstrip('\r\n')
                vals = map(float, temp.split('\t'))
                if not x_range:
                    self.xvals.append(vals[0])
                    self.yvals.append(vals[1])
                else:
                    if vals[0] > x_range[0] and vals[0] < x_range[1]:
                        self.xvals.append(vals[0])
                        self.yvals.append(vals[1])
                count = 0
        self.xvals = np.array(self.xvals)
        self.yvals = np.array(self.yvals)

        # so that it isn't overwritten by smoothing
        self.rawyvals = self.yvals.copy()
        self._normalisePreset()

    def normalisationBpi(self):
        '''Normalise to base peak intensity (0-100)'''
        self.yvals = self.yvals/self.yvals.max()*100
        self.setNormalisationType(type='bpi')
    def normalisationArea(self):
        self.yvals = self.yvals/np.sum(self.yvals)
        self.setNormalisationType(type='area')

    def setNormalisationType(self,type):
        self.normalisationType = type

    def smoothingSG(self,window_len=3,smoothes=2,poly_order=1):
        '''Should only really be used on equally spaced data
        Actual window length used is 2*window_len+1 to avoid breakage'''
        window_len = 2*window_len + 1
        self.restoreRawYvals()
        for i in xrange(smoothes):
            self.yvals = SG.sg(self.yvals,window_size=window_len,order=poly_order)
        self._normalisePreset()

    def _normalisePreset(self):
        if self.normalisationType == 'bpi':
            self.normalisationBpi()
        elif self.normalisationType == 'area':
            self.normalisationArea()
        elif self.normalisationType == 'none':
            pass

    def restoreRawYvals(self):
        self.yvals = self.rawyvals.copy()
        self._normalisePreset()

    def getAxesWithoutNans(self):
        """The CCS calibration can cause some xvals
        to become NaNs, this function returns x and yvals
        truncated to remove the NaNs
        """
        xvals = self.xvals[np.invert(np.isnan(self.xvals))]
        yvals = self.yvals[np.invert(np.isnan(self.xvals))]
        return xvals,yvals

    #===========================================================================
    # Calculations
    #===========================================================================
    def calculateWeightedMeanStandardDeviation(self):
        yvals = self.yvals/self.yvals.max()*100
        average,stdev = utils.weightedAverageAndStd(self.xvals,yvals)
        return average,stdev
    def calculateAreaUnderCurve(self):
        """Integrates data using trapezium method"""
        xvals,yvals = self.getAxesWithoutNans()
        yvals = yvals/yvals.max()*100
        return np.trapz(yvals,xvals)

    #===========================================================================
    # Peak finding
    #===========================================================================

    def _calculateGradient(self):
        """when reconstructing the data, make data[0] the start value,
        skip the first gradient value then append on
        gradient[i] * (ys[i+1] - ys[i]) (actually gradient[i+1] ...)"""
        self.gradient = [0]
        for i,x in enumerate(self.xvals):
            if i+2 <= len(self.xvals):
                try:
                    gr = (float(self.yvals[i+1])-float(self.yvals[i]))/(float(self.xvals[i+1]) - float(x))
                except:
                    print 'Gradient calculation: divide by 0 replaced by 0.000001'
                    gr = 0.000001
                self.gradient.append(gr)
        self.gradient = np.array(self.gradient)

    def findPeaks(self,limit=0):
        """limit allows you to ignore slow peaks (remove noise)
        percentage of BPI e.g. 5 % cutoff should be 5"""

        # get gradient for peak picking
        self._calculateGradient()

        gPeaks = OrderedDict()
        count = 0
        for i,v in enumerate(self.gradient):
            if i+1 < len(self.gradient):
                if v > 0:
                    if self.gradient[i+1] <= 0:
                        gPeaks[count] = []
                        gPeaks[count].append(self.xvals[i])
                        gPeaks[count].append(self.yvals[i])
                        count += 1
        if limit:
            gPeaks_out = OrderedDict()
            lim = max([gPeaks[x][1] for x in gPeaks.keys()]) * float(limit)/100
            count = 0
            for i,(k,v) in enumerate(gPeaks.items()):
                if v[1] > lim:
                    gPeaks_out[count] = []
                    gPeaks_out[count].append(gPeaks[k][0])
                    gPeaks_out[count].append(gPeaks[k][1])
                    count += 1
            self.gPeaks = gPeaks_out

        else:
            self.gPeaks = gPeaks

    def addPeak(self,mz):
        """This function allows you to add additional peaks not found using the
        findPeaks method"""
        try:
            keys = sorted(self.gPeaks.keys())
            id = 1+keys[-1]
            self.gPeaks[id] = {}
        except:
            id = 0
            self.gPeaks = {}
            self.gPeaks[id] = {}
        self.gPeaks[id] = [[],[]]
        self.gPeaks[id][0] = mz
        index = utils.closest(mz,self.xvals)
        self.gPeaks[id][1] = self.yvals[index]

    #===========================================================================
    # Plotting
    #===========================================================================
    def plot(self,ax,**kwargs):
        """Plot 2D data (e.g. MS and ATDs)
        Can take matplotlib axes object, as well as any standard
        inputs for matplotlib.pyplot.plot().
        """
        ax = utils.checkAx(ax)
        if not 'color' in kwargs:
            kwargs['color'] = 'black'
        if not 'lw' in kwargs:
            kwargs['lw'] = 0.8
        ax.plot(self.xvals,self.yvals,**kwargs)
        ax.set_ylabel('Intensity')
        ax.set_xlabel('$m/z$')

    def plotgPeaks(self,ax,labels=0,**kwargs):
        """if labels=0 then only peak id's are displayed
        otherwise x-values are shown as well
        """
        if not 'color' in kwargs:
            kwargs['color'] = 'gray'
        if not 'alpha' in kwargs:
            kwargs['alpha'] = 0.5

        for i,k in enumerate(self.gPeaks):
            #print k
            #print self.gPeaks.keys()
            ax.axvline(self.gPeaks[k][0], **kwargs)
            if labels:
                label = str(i)+':'+str(self.gPeaks[k][0])
            else:
                label = str(i)
            ax.annotate(label,[self.gPeaks[k][0],self.gPeaks[k][1]])
