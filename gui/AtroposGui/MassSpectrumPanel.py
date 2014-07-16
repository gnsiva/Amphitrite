"""Right hand side panel for displaying matplotlib graphs."""

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
from msClasses import  MassSpectrum,Species
from imClasses import Im
from lib import utils
matplotlib.rc('xtick', labelsize=8) 
matplotlib.rc('ytick', labelsize=8) 

class MassSpectrumPanel():
    
    def __init__(self,panel,settings,yourself):
        self.dpi = 80
        self.fig = Figure((3.0, 2.0), dpi=self.dpi)
        self.canvas = FigCanvas(panel, -1, self.fig)
        self.axMs = self.fig.add_subplot(111)
        
        self.fig.canvas.mpl_connect('pick_event',self.on_pick)
        
        self.toolbar = NavigationToolbar(self.canvas)  
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)
        panel.SetSizer(self.vbox)
        self.canvas.draw() 
        
        self.axMs.plot([1,2,3],[1,2,3])
        
        self.pickedValue = None
        self.pickingActive = False
        self.pickedAlready = False
        self.tempAnnotation = None

        self.pickedPeakIdLimit = 5001
        
        self.settings = settings
        
        self.ms = MassSpectrum()
        self.ms.setNormalisationType('bpi')

        self.gui = yourself
        
    def plotMassSpectrum(self):
        """Draw mass spectrum in panel."""
        self.refresh_plot()
    
    def plotMsWithGpeaks(self,limit):
        """Draw mass spectrum with automatically found and manually added
        peak positions as vertical lines, labelled with peak identifiers.

        :parameter limit: percentage of BPI under which peaks found are ignored
        """
        self.refresh_plot()        
        try: limit = float(limit)
        except: limit = 0.0
        self.ms.findPeaks(limit=limit)

        
        if len(self.ms.gPeaks) > 2000:
            return len(self.ms.gPeaks)

        for mz in self.settings.addedPeakMzs:
            self.ms.addPeak(mz)
        self.ms.plotgPeaks(self.axMs)
        self.draw()
    
    def plotMsWithTheoZs(self,speciesName):
        """Draw mass spectrum with theoretical charge states of a species
        overlaid and labelled.

        :parameter speciesName: Name of the molecular species to use
        """
        self.refresh_plot()
        self.ms.plotTheoreticalMzs(self.axMs, speciesName)
        self.axMs.legend()
        self.draw()           
    
    def plotLeastSquaresSimulation(self,speciesNames):
        """Draw the result of the deconvolution in the panel.

        :parameter speciesNames: List of species which were deconvoluted
        """
        self.axMs.clear()
        self.ms.normalisationBpi()
        self.ms.plotAllComponents(self.axMs, self.ms.xvals, self.ms.yvals, 
                                  liftPercentage=20, speciesNames=speciesNames)

        # needs full integration with refresh_plot()
        self.draw()
    
    def setSettings(self,settings):
        """Set the AtroposGuiSettings object, so its functions and
        attributes can be accessed more readily.

        :parameter settings: AtroposGuiSettings() object
        """
        self.settings = settings
        
    def smoothing(self,wlen,smoothes,poly):
        """Smooth the mass spectrum (Savitzky-Golay) and redraw it.

        :parameter wlen: Window length for smoothing (actual size = wlen*2 + 1 as the number has to be odd)
        :parameter smoothes: How many rounds of smoothing to do
        :parameter poly: Polynomial to use for smoothing
        
        """
        self.ms.smoothingSG(wlen, smoothes, poly)
        self.refresh_plot()
    
    def loadTextFile(self,filename,grain=0):
        """Load mass spectrum from text file and plot it.

        :parameter filename: Text file filename
        :parameter grain: Can be used to reduce the data size. Grain of 2 means only use 1 in every 2 datapoints
        """
        if grain:
            self.ms.readFile(filename, grain=grain)
        else:
            self.ms.readFile(filename)
        self.refresh_plot()

    def loadAmphiFile(self,filename):
        """Load mass spectrum from .a file and display it.

        :parameter filename: Amphitrite data file (.a)
        """
        dataList = utils.unPickleAmphitriteProject(filename)
        imOb = Im()
        imOb.setDataFromAmphiExtract(dataList)

        import copy
        species = copy.deepcopy(self.ms.species)
        simulatedSpecies = copy.deepcopy(self.ms.simulatedSpecies)
        
        self.ms = imOb.getMassSpectrum()
        self.ms.setNormalisationType('bpi')

        self.ms.species = species
        self.ms.simulatedSpecies = simulatedSpecies

        
        self.refresh_plot()

    def refresh_plot(self):
        """Update plot panel."""
        # TODO(gns) - Much simpler than in other GUI elements
        # Just handles the final plotting information
        # Should later be used to control the other
        # plotting functions
        self.ms.normalisationBpi()
        self.axMs.clear()
        self.ms.plot(self.axMs)
        self.axMs.set_ylim(0,103)
        self.draw()
        
    
    def on_pick(self, event):
        """Check if peak picking is active. If it is draw vertical line at the
        peak's m/z value and label with peak ID.
        """
        if self.pickingActive:
            if self.pickedAlready:
                line = self.axMs.lines.pop(len(self.axMs.lines)-1)
                self.tempAnnotation.set_visible(False)
            self.pickedValue = event.mouseevent.xdata
            yval = self.ms.yvals[utils.closest(self.pickedValue,self.ms.xvals)]

            peakId = self.pickedPeakIdLimit #value controlled by gui
            
            self.axMs.axvline(self.pickedValue,color='k')
            self.tempAnnotation = self.axMs.annotate(str(peakId),[self.pickedValue,yval])
            self.draw()
            
            self.pickedAlready = True

    def temporaryTheoMzs(self,mass):
        """For use when only mass is given (no peakIds)."""
        self.refresh_plot()
        self.gui.radioBoxDisplay.SetSelection(2)
        tempSp = Species('Temporary Species')
        tempSp.plotTheoreticalMzs(self.axMs,self.ms.xvals,self.ms.yvals,mass=mass)
        self.gui.msPanel.draw()        
      
    def draw(self):
        """Convenience function for calling self.canvas.draw()
        (updates the matplotlib canvas).
        """
        self.canvas.draw()
