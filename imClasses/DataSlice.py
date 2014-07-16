"""Class for handling a portion of an IM dataset. Usually is used for holding
the data within a signle mass spectrum charge state peak.
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import ImData
import numpy as np
import matplotlib.pyplot as plt
from imClasses import Atd, CcsD
from lib import utils


class DataSlice(ImData.ImData):
    
    def __init__(self):
        ImData.ImData.__init__(self)
        self.atd = Atd()
        self.ccsAxisGrid = []
        self.ccsAxis = []
        self.charge = None

        self.ccsd = None
        """
        ---------------------
        Inherited:
        ---------------------
        self.xaxis = []
        self.yaxis = []
        self.matrix = []
        self.matrixUnaltered = []
        self.xlims = []
        self.ylims = []

        self.ccsMatrix # dno where from
        
        self.atd = Atd.Atd()
        self.massSpectrum = MassSpectrum.MassSpectrum()
        """
        
    def getAtdApex(self):
        """Get the arrival time value for the point of highest intensity.
        Get the highest point of the arrival time distribution.

        :returns: Arrival time of maximum intensity (float)
        """
        index = self.atd.yvals.argmax()
        return self.atd.xvals[index]

    def getMsApexI(self):
        """Get index of the m/z axis for the point of highest intensity.

        :returns: index
        """
        self.generateMassSpectrum()
        maxI = self.massSpectrum.yvals.argmax()
        return maxI
        
    def getMsApex(self):
        """Get the m/z value for the point of highest intensity.

        :returns: m/z value (float)
        """
        maxI = self.getMsApexI()
        maxMz = self.massSpectrum.xvals[maxI]
        return maxMz

    def getMsApexProportion(self):
        """Relative to width of the slice e.g. 0.5 for in the middle

        :returns: proportion of axis (float)
        """
        maxI = self.getMsApexI()
        entries = len(self.massSpectrum.xvals)
        return float(maxI)/entries

    def getAtdPeaks(self,smoothes,window_len,limit=0):
        """Find peaks using gradient method.

        :returns: atdPeaks (list of arrival time values)
        """
        self.atd.smoothingSG(window_len,smoothes)
        self.atd.findPeaks(limit=limit)
        atdPeaks = [v[0] for v in self.atd.gPeaks.values()]
        return atdPeaks

    def getCcsPeaks(self,smoothes,window_len,mz,calibrationOb,limit=0):
        """Smooth data then find peaks using the gradient method in the arrival time
        distribution. Convert the arrival time values to CCS.

        :parameter smoothes: Number of rounds of Savitzky-Golay smoothing 
        :parameter window_len: Savitzky-Golay smoothing window size (size used is window_len*2 + 1)
        :parameter mz: m/z value to use for CCS conversion
        :parameter calibrationObj: imClasses.Calibration() object
        :parameter limit: Percentage of base peak under which to ignore peaks found
        :returns: ccsPeaks (list of CCS values)
        """
        atdPeaks = self.getAtdPeaks(smoothes,window_len,limit)
        ccsPeaks = calibrationOb.apply1dCalibration(mz,atdPeaks,self.charge)
        return ccsPeaks
    
    
    def getCcsDistribution(self,calibrationObj):
        """Convert the matrix into CCS, sum it and return the intensity
        values for the CCS distribution.

        :parameter calibrationObj: imClasses.Calibration() object
        :returns: ccsLine (CCS distribution y axis values)
        """
        # TODO(gns) - test first but replace the code in this function with
        # generateCcsDistribution()
        '''PHASE THIS OUT, REPLACED BY WHATS SEEN IN generateCcsDistribution
        USE IN CCSD OBJECT FORM'''
        ccsMatrix, mzs, ccss = self.getCcsMatrix(calibrationObj)
        ccsLine = np.sum(ccsMatrix, axis=1)
        return ccsLine

    def generateCcsDistribution(self,calibrationObj):
        """Create a CcsD() out of this object's data.

        :parameter calibrationObj: imClasses.Calibration() object
        :returns: imClasses.CcsD() object
        """
        self.ccsd = CcsD(calibrationObj,self)

    
    def getData(self):
        """Get the intensity matrix, m/z axis and arrival time axis.

        :returns: matrix, x axis, y axis (all numpy arrays)
        """
        return self.matrix,self.xaxis,self.yaxis
    
    def getCcsMatrix(self,calibrationObj):
        """Use calibration to create a grid of CCS values associated with
        each entry in the intensity matrix.

        :parameter calibrationObj: imClasses.Calibration() object (IM calibration)
        :returns: CCS matrix, x axis and y axis (all numpy arrays)
        """
        self.generateCcsAxisAndGrid(calibrationObj)
        ccsMatrix = np.zeros([len(self.ccsAxis),len(self.xaxis)])
        for i,mz in enumerate(self.xaxis):
            ccsMatrix[:,i] = np.interp(x=self.ccsAxis, xp=self.ccsAxisGrid[:,i], fp=self.matrix[:,i])
        return ccsMatrix, self.xaxis, self.ccsAxis


    def generateCcsAxisAndGrid(self,calibrationObj,ccsInterval=1):
        """Create CCS axis which corresponds to the intensity information of an arrival
        time distribution, also create CCS grid which provides the CCS value for each
        datapoint in the intensity matrix.

        :parameter calibrationObj: imClasses.Calibration() object (IM calibration)
        :parameter ccsInterval: Spacing between each CCS value in Angstrom**2
        """
        '''calculates the equivalent CCS values for each
        td mz combination (used as an axis)'''
        mzs = self.xaxis
        tds = self.yaxis
        self.ccsAxisGrid = calibrationObj.getCcsAxisGrid(mzs,tds,self.charge)
        self.ccsAxis = np.arange(int(self.ccsAxisGrid.min()),
                                 int(self.ccsAxisGrid.max())+1,ccsInterval,dtype='float64')
        

    def plotCcsDistribution(self,ax,calibrationObj,**kwargs):
        """Plot the CCS distribution for this DataSlice.

        :parameter ax: matplotlib axes instance or False
        :parameter calibrationObj: imClasses.Calibration() object (IM calibration)
        :parameter \*\*kwargs: matplotlib.pyplot.plot() arguments
        """
        ccsLine = self.getCcsDistribution(calibrationObj)
        if not 'label' in kwargs:
            kwargs['label'] = '+%d' %self.charge
        ax.plot(self.ccsAxis,ccsLine,**kwargs)
        ax.xlabel('CCS ($\AA$)')
        ax.ylabel('Intensity')

