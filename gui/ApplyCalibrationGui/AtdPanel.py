"""Plotting panel for ApplyCalibrationGui()."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import PlotPanel
from imClasses import Atd
import ApplyCalibrationGuiSettings

class AtdPanel(PlotPanel.PlotPanel):
    """
    :parameter panel: wx panel to contain the plotting area
    """
    def __init__(self,panel):
        PlotPanel.PlotPanel.__init__(self, panel)
        self.settings = ApplyCalibrationGuiSettings.ApplyCalibrationGuiSettings()
        
        self.atd = Atd()
    
    def setSettings(self,settings):
        """Set the Gui Settings object.
        :parameter settings: ApplyCalibrationGuiSettings() object
        """
        self.settings = settings
    
    def openAtd(self):
        """Open coordinates file and create an imClasses.Atd() object
        with it.
        """
        path = self.settings.coordinatesPath
        print path
        self.atd.readFile(path)        
    
    def plotAtd(self,**kwargs):
        """Clear the plot, normalise arrival time distribution to
        base peak intensity and plot.
        """
        self.ax.clear()
        self.atd.normalisationBpi()
        self.atd.plot(self.ax,**kwargs)
        self.ax.set_xlabel('t$_d$')
        self.draw()
    
    def autoAxisX(self,ccs=False):
        """Turn on automatic x axis limits.
        :parameter ccs: Boolean for whether CCS is being used (as opposed to
        arrival time)
        """
        if ccs:
            ccs = self.settings.ccs        
            self.atd.autoAxisX(self.ax, limit=0.1,ccsAxisValues=ccs)
        else:
            self.atd.autoAxisX(self.ax, limit=0.1)
        self.draw()
    def autoAxisXoff(self):
        """Turn off automatic x axis limits.
        """
        ccs = 0
        if self.settings.displayedView == 'CCS':
            ccs = self.settings.ccs
        self.atd.autoAxisXoff(self.ax,ccs)
        self.draw()

    def plotCcs(self,**kwargs):
        """Clear the plot, normalise intensity values and
        plot the CCS distribution.
        """
        self.ax.clear()
        self.atd.normalisationBpi()
        self.atd.plotCcs(self.ax, self.settings.ccs,**kwargs)
        self.draw()
    
    def clear(self):
        """Clear the plot area.
        """
        self.ax.clear()

