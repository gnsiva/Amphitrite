"""Class for collision cross section distributions (CCSDs) - plotting, peak identification and calculations."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import Im,Atd
from msClasses import TwoDdata
import numpy as np
from lib import utils

class CcsD(TwoDdata):
    """
    :parameter calibrationOb: Calibration() object for converting ATD to CCSD
    :parameter dataSlice: DataSlice() object containing the data for the CCSD
    """
    def __init__(self,calibrationOb,dataSlice):
        # TODO(gns) - Shouldn't you call the TwoDdata init here?
        self.xvals = None
        self.yvals = None

        self.rawyvals = None
        self.normalisationType = 'bpi'
        
        self.createCcsD(calibrationOb,dataSlice)

    def _createCcsD(self,calibrationOb,dataSlice):
        """Generate a CCS distribution from a DataSlice(), updates self.xvals and self.yvals.

        :parameter calibrationOb: Calibration() object for converting ATD to CCSD
        :parameter dataSlice: DataSlice() object containing the data for the CCSD
        """
        # TODO(gns) - this function used to be called self.createCcsD(), so check for problems.
        ds = dataSlice
        ccsMatrixAxis,mzs,ccsAxis = ds.getCcsMatrix(calibrationOb)
        ccsSignal = np.sum(ccsMatrixAxis,axis=1)

        # the range of values is too large so interpolate back to
        # original size

        # calibrate the td axis directly
        charge = ds.charge
        mz = ds.getMsApex()
        tds = ds.atd.xvals

        tdsConvertedToCcss = calibrationOb.apply1dCalibration(mz,tds,charge)

        # interpolate
        yInterp = np.interp(tdsConvertedToCcss,ccsAxis,ccsSignal)

        self.setXvals(tdsConvertedToCcss)
        self.setYvals(yInterp/yInterp.max()*100)

        
    def setYvals(self,yvals):
        """Set the y axis values for the object.
        
        :parameter yvals: Intensity values for the Atd() object
        :type numpy array:
        :returns: None
        """
        # TODO(gns) - This should be in TwoDdata, same problem in Atd()
        self.rawyvals = yvals
        self.yvals = self.rawyvals.copy()

    def setXvals(self,xvals):
        """Set the x axis values for the object.
        
        :parameter yvals: Arrival time values for the Atd() object
        :type numpy array:
        :returns: None
        """
        self.xvals = xvals
        
    
    def plot(self,ax,**kwargs):
        """Overwrites TwoDdata()'s function to get axis labels right
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
        ax.set_xlabel('CCS ($\AA^2$)')
