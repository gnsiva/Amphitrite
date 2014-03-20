import os, CcsCalculatorListCtrl, PlottingOptionsListCtrl, AtdPanel
from imClasses import Calibration
import cPickle as pickle
import numpy as np

class ApplyCalibrationGuiSettings():
    def __init__(self):
        self.calibrationPath = ''
        self.coordinatesPath = ''
        
        self.calibrationLoaded = False
        self.coordinatesLoaded = False
        
        self.mz = ''
        self.charge = ''
        self.ccs = ''
        self.tds = ''
        
        self.calibration = Calibration.Calibration()
        self.panelCcsCalculator = ''# CcsCalculatorListCtrl.CcsCalculatorListCtrl()
        self.panelPlottingOptions = ''#PlottingOptionsListCtrl.PlottingOptionsListCtrl()
        self.atdPanel = '' #AtdPanel.AtdPanel()
        
        self.displayedView = ''

    def setCalibrationPath(self,path):
        """Set the path to the pickled imClasses.Calibration() object
        :parameter path: Absolute path
        """
        self.calibrationPath = path
    def setCoordinatesPath(self,path):
        """Set the path to the arrival time distribution coordinate file.
        Data should be in the format of ArrivalTime\tIntensity\n.
        :parameter path: Absolute path
        """
        self.coordinatesPath = path
    
    def setCalibration(self):
        """Open pickled imClasses.Calibration() object.
        """
        self.calibration = pickle.load(open(self.calibrationPath,'rb'))
        self.calibrationLoaded = True

    def readyForCcs(self,mz,charge):
        """Check if the required information has been provided to calculate
        a CCS distribution.
        :parameter mz: m/z value
        :parameter charge: Charge state
        """
        ready = True
        if not os.path.isfile(self.coordinatesPath):
            ready = False
            print 'Co-ordinates file not found'
        elif not os.path.isfile(self.calibrationPath):
            ready = False
            print 'Calibration file not found'
        else:
            if mz == '' or charge == '':
                ready = False
            else:
                try: self.mz = float(mz)
                except: ready = False
                
                try: self.charge = int(charge)
                except: ready = False
        return ready

    def getCcsAxis(self):
        """Calculate CCS axis by converting each value in the arrival
        time axis (x) of the distribution. Then store it in self.ccs and return it.
        :returns: ccs axis 
        """
        self.tds = self.atdPanel.atd.xvals.copy()
        self.ccs = self.calibration.apply1dCalibration(self.mz, self.tds, self.charge)
        self.ccs[np.isnan(self.ccs)] = 0
        return self.ccs
