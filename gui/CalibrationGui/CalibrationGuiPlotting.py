"""Plotting panel for CalibrationGui."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt
import numpy as np
import wx
matplotlib.rc('xtick', labelsize=8) 
matplotlib.rc('ytick', labelsize=8) 

import imClasses.Calibration as Calibration

class CalGuiFigure():
    """
    :parameter panel: The wx panel in which this graph should be displayed
    :parameter calibrants: imClasses.Calibrant() objects to use in calibration
    :parameter settings: CalibrationGuiSettings() object to hold GUI settings
    :parameter yourself: 'self' from the Gui object
    """
    def __init__(self,panel,calibrants,settings,yourself):
        self.dpi = 80
        self.fig = Figure((3.0, 2.0), dpi=self.dpi)
        self.canvas = FigCanvas(panel, -1, self.fig)
        self.axMs = self.fig.add_subplot(311)
        self.axTds = self.fig.add_subplot(312)
        self.axCal = self.fig.add_subplot(313)
        
        self.fig.canvas.mpl_connect('pick_event',self.on_pick)
        
        self.toolbar = NavigationToolbar(self.canvas)  
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)
        panel.SetSizer(self.vbox)
        self.canvas.draw() 
        
        self.axMs.plot([1,2,3],[1,2,3])
        self.axTds.plot([3,2,1],[2,2,2])
        self.axCal.plot([1,2,3],[3,2,1])
        
        self.pickedValue = None
        self.pickingActive = False
        
        self.calibrants = calibrants
        self.settings = settings
        self.gui = yourself
        
    def on_pick(self, event):
        """When peak picking is active, this registers clicks, stores an
        arrival time value for that and redraws the mass spectrum with a
        vertical line showing the arrival time value.
        """
        if self.pickingActive:
            try:
                for i, line in enumerate(self.axTds.lines):
                    if i:
                        self.axTds.lines.pop(i)
                        line.remove()
            except:
                pass
            self.pickedValue = event.mouseevent.xdata
            self.lines = self.axTds.axvline(self.pickedValue,color='k')
            proName = self.settings.list_ctrl.list_ctrl.getProName()
            self.calibrants[proName].setTdValue(self.pickedValue,self.settings.list_ctrl.list_ctrl.getCharge())
            self.calibrants[proName]._updateSpecies(updateApex=0)
            self.createCalibration()
            self.canvas.draw()

    def createCalibration(self):
        """Create a travelling wave ion mobility calibration and plot the
        calibration curve.
        """
        gas = self.settings.getGas()
        waveVelocity = self.settings.getWaveVelocity()

        self.calibration = Calibration.Calibration()
        # This checks to make sure there are at least 2 peaks which can be used
        # in the least squared fitting (2 charges or multiple proteins)
        n = 0
        for name,calibrant in self.calibrants.items():
            calibrant.generateCorrectedTdsAndCcss(waveVelocity,gas)
            self.calibration.addCalibrant(calibrant)

        self.axCal.cla()
        if self.gui.panelCalibrantListCtrl.getNumberOfCheckedCharges() > 1:
            self.calibration.createCalibration(waveVelocity=waveVelocity, gas=gas)

            # Plotting the calibration curve
            self.calibration.plotCalibrationCurve(self.axCal)


            
