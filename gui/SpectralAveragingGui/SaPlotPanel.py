"""Plotting area for SpectralAveragingGui()."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt
import numpy as np
import wx,os
import classes.MassSpectrum as MassSpectrum
from lib import utils
matplotlib.rc('xtick', labelsize=8) 
matplotlib.rc('ytick', labelsize=8) 


from collections import OrderedDict
import matplotlib.gridspec as gridspec
from AmphitriteEnums import *
import SaSettings


import matplotlib.font_manager
prop = matplotlib.font_manager.FontProperties(size='small')

class SaPlotPanel():
    
    def __init__(self,panel):
        self.dpi = 80
        self.fig = Figure((3.0, 2.0), dpi=self.dpi)
        self.canvas = FigCanvas(panel, -1, self.fig)
        self.gs = gridspec.GridSpec(1,1)
        self.axes = [None]
                
        #self.fig.canvas.mpl_connect('pick_event',self.on_pick)
        
        self.toolbar = NavigationToolbar(self.canvas)  
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)
        panel.SetSizer(self.vbox)
        self.canvas.draw() 

        self.pickedValue = None
        self.pickingActive = False
        self.pickedAlready = False
        self.tempAnnotation = None

        self.settings = None
        self.gui = None
        
        self.lift = 30

        self.picker = None
        self.releaser = None
        self.pickedValue = None
        self.releasedValue = None
        self.selectedAxis = None

        self.atroposLeft = None
        self.atroposRight = None
        
    #====================
    # Adding other panel components        
    def setSettings(self,settings):
        """Set the Gui settings object.
        :parameter settings: SaSettings() object
        """
        self.settings = settings
    def setGui(self,gui):
        """Set the plotting area object.
        :parameter plotPanel: SaPlotPanel() object
        """
        self.gui = gui
    #====================
        
    def _preparePlottingSections(self,rows,columns):
        """Clear currently plotted information and recreate the plot area.
        :parameter rows: Rows of Matplotlib Axes instances
        :parameter columns: Columns of Matplotlib Axes instances
        """
        self.fig.clf(keep_observers=True)
        self.gs = gridspec.GridSpec(rows,columns)

    def getSingleAxis(self):
        """Recreate plotting area with a single set of axes.
        :returns: Matplotlib Axes instance
        """
        self._preparePlottingSections(1,1)
        self.axes = [None]
        self.axes[0] = self.fig.add_subplot(self.gs[0,0])
        self.axes[0].plot([],[])
        return self.axes[0]
        
    def getDoubleColumnAxes(self):
        """Recreate plotting area with a two sets of axes (in a column).
        :returns: List of two Matplotlib Axes instances
        """
        self._preparePlottingSections(1,2)
        self.axes = [None,None]
        self.axes[0] = self.fig.add_subplot(self.gs[0])
        self.axes[1] = self.fig.add_subplot(self.gs[1])

        return self.axes[0],self.axes[1]

    def getDoubleRowAxes(self):
        """Recreate plotting area with a two sets of axes (in a row).
        :returns: List of two Matplotlib Axes instances
        """
        self._preparePlottingSections(2,1)
        self.axes = [None,None]
        self.axes[0] = self.fig.add_subplot(self.gs[0])
        self.axes[1] = self.fig.add_subplot(self.gs[1])

        return self.axes[0],self.axes[1]

    def getQuadAxes(self):
        """Recreate plotting area with 4 sets of axes.
        :returns: List of 4 Matplotlib Axes instances
        """
        self._preparePlottingSections(2,3)
        self.axes = [None,None,None,None]
        self.axes[0] = self.fig.add_subplot(self.gs[0,:2]) # MS
        self.axes[1] = self.fig.add_subplot(self.gs[1,:2]) # ATD
        self.axes[2] = self.fig.add_subplot(self.gs[0,2])   # MS Line
        self.axes[3] = self.fig.add_subplot(self.gs[1,2])   # ATD Line
        
        return self.axes[0],self.axes[1],self.axes[2],self.axes[3]
        
    
    def plotMassSpectra(self,ax):
        """Plot mass spectra and the currently selected variability
        representation.
        :parameter ax: Matplotlib Axes instance
        """
        xaxis,yaxes,fns = self.settings.getMsXaxisYaxes()
        self.plotAverageLine(ax,xaxis,yaxes,color='k',lw=0.5)

        ax.set_xlabel('$m/z$')
        ax.set_ylabel('Intensity')
        ax.set_yticks([])

        # boundaries
        if self.settings.boundaryAlgorithm == SaBndAlg.MIN_MAX:
            self.boundaryAlgorithmMinMax(ax,xaxis,yaxes)
        elif self.settings.boundaryAlgorithm == SaBndAlg.STD:
            self.boundaryAlgorithmStd(ax,xaxis,yaxes)
        elif self.settings.boundaryAlgorithm == SaBndAlg.NONE:
            pass

        # Click and drag selection
        if self.gui.checkboxDisplaySelection.IsChecked():
            if self.pickedValue and self.releasedValue:
                self.plotClickAndDragSelection(ax)
            elif self.atroposLeft and self.atroposRight:
                self.plotAtroposSelection(ax)


    def plotMsLines(self,ax):
        """Plot lines for each set of mass spectrum intensity axes values.
        :parameter ax: Matplotlib Axes instance
        """
        xaxis,yaxes,fns = self.settings.getMsXaxisYaxes()
        for i,yaxis in enumerate(yaxes):
            fn = os.path.basename(fns[i])
            ax.plot(xaxis,yaxis,color=utils.colourList[i],label=fn,lw=0.3)
        if self.settings.showLegend:
            ax.legend(prop=prop,loc='best')

        # Click and drag selection
        if self.gui.checkboxDisplaySelection.IsChecked():
            if self.pickedValue and self.releasedValue:
                self.plotClickAndDragSelection(ax)
            elif self.atroposLeft and self.atroposRight:
                self.plotAtroposSelection(ax)                
                
        ax.set_xlabel('$m/z$')
        ax.set_ylabel('Intensity')
        ax.set_yticks([])

    
    
    def plotClickAndDragSelection(self,ax):
        """Draw shaded region for click and release regions (mass
        spectrum use).
        :parameter ax: Matplotlib Axes instance
        """
        xpoints = [self.pickedValue,self.releasedValue]
        xpoints = sorted(xpoints)
        ylims = ax.get_ylim()
        ax.fill_betweenx(ylims,xpoints[0],xpoints[1],
                         color='red',alpha=0.2)

    def plotAverageLine(self,ax,xaxis,yaxes,**kwargs):
        """Plot the mean of the instensity values.
        :parameter ax: Matplotlib Axes instance
        :parameter xaxis: X values
        :parameter yaxes: Numpy array of Y axes values
        :parameter \*\*kwargs: Matplotlib.pyplot.plot compatible arguments
        """
        yaxis = np.average(yaxes,axis=0)
        ax.plot(xaxis,yaxis,**kwargs)
        
    def boundaryAlgorithmMinMax(self,ax,xaxis,yaxes):
        """Calculate the minimum and maximum value for the set
        of intensity values per m/z value.
        :parameter ax: Matplotlib Axes instance
        :parameter xaxis: Mass spectrum m/z axis
        :parameter xaxis: Mass spectrum intensity axes
        """
        ymin = np.min(yaxes,axis=0)
        ymax = np.max(yaxes,axis=0)
        ax.fill_between(xaxis,ymin,ymax,color='blue',alpha=0.2)

    def boundaryAlgorithmStd(self,ax,xaxis,yaxes):
        """Calculate the standard deviation for the set of intensity
        values per m/z value.
        :parameter ax: Matplotlib Axes instance
        :parameter xaxis: Mass spectrum m/z axis
        :parameter xaxis: Mass spectrum intensity axes
        """
        yaverage = np.average(yaxes,axis=0)
        ystd = np.std(yaxes,axis=0)
        ymin = yaverage-ystd
        ymax = yaverage+ystd
        ax.fill_between(xaxis,ymin,ymax,color='green',alpha=0.2)
        
    def getAtdsFromMzLims(self,mzLow,mzHigh):
        """Extract arrival times from mass spectrum limits.
        :parameter mzLow: Lower limit for extraction
        :parameter mzHigh: Upper limit for extraction
        :returns: Dictionary of imClasses.Atd() objects (d[filename] = Atd())
        """
        atds = OrderedDict()        
        for fn,imOb in self.settings.loadedFiles.items():
            mzLowIndex = utils.closest(
                mzLow,self.settings.massSpectra[fn].xvals)
            mzHighIndex = utils.closest(
                mzHigh,self.settings.massSpectra[fn].xvals)
            ds = imOb._generateSlice([mzLowIndex,mzHighIndex])
            ds.generateAtd()
            atds[fn] = ds.atd
        return atds
        
    def plotAtds(self,ax):
        """Plot averaged arrival time distribution.
        :parameter ax: Matplotlib Axes instance
        """
        if self.pickedValue and self.releasedValue:
            print 'using picked values'
            mzLow = sorted([self.pickedValue,self.releasedValue])[0]
            mzHigh = sorted([self.pickedValue,self.releasedValue])[1]
            atds = self.getAtdsFromMzLims(mzLow,mzHigh)
        elif self.atroposLeft and self.atroposRight:
            print 'using atropos values'
            atds = self.getAtdsFromMzLims(self.atroposLeft,
                                          self.atroposRight)
        else:
            print 'using full set of values'
            atds = OrderedDict()
            for fn,imOb in self.settings.loadedFiles.items():
                mzLeft = 0
                mzRight = len(self.settings.massSpectra[fn].xvals)-1
                ds = imOb._generateSlice([mzLeft,mzRight])
                ds.generateAtd()
                atds[fn] = ds.atd

        
        xaxis,yaxes = self.settings.getAtdXaxisYaxes(
            [atd.xvals for atd in atds.values()],
            [atd.yvals for atd in atds.values()])

        self.plotAverageLine(ax,xaxis,yaxes,color='k',lw=0.8)
        
        ax.set_xlabel('t$_d$')
        ax.set_ylabel('Intensity')
        ax.set_yticks([])
        
        # boundaries
        if self.settings.boundaryAlgorithm == SaBndAlg.MIN_MAX:
            self.boundaryAlgorithmMinMax(ax,xaxis,yaxes)
        elif self.settings.boundaryAlgorithm == SaBndAlg.STD:
            self.boundaryAlgorithmStd(ax,xaxis,yaxes)
        elif self.settings.boundaryAlgorithm == SaBndAlg.NONE:
            pass

        return atds

    
    def plotAtdLines(self,ax,atds):
        """Plot the arrival time distributions for each input
        data file.
        :parameter ax: Matplotlib Axes instance
        :parameter atds: Dictionary of imClasses.Atd() objects
        (d[filename] = Atd())
        """
        if len(atds):
            for i,(fn,atd) in enumerate(atds.items()):
                ax.plot(atd.xvals,atd.yvals,
                        color=utils.colourList[i],label=fn)
            if self.gui.checkboxLegend.IsChecked():
                ax.legend(loc=0,prop=prop)
            
        ax.set_xlabel('t$_d$')
        ax.set_ylabel('Intensity')
        ax.set_yticks([])

        
    def plotAtroposSelection(self,ax):
        """Draw a shaded region on the mass spectrum for use
        in extracting arrival time data.
        """
        species,z = self.settings.getSpeciesAndCharge()
        limits = self.settings._getMzLimits(species,z)
        ylims = ax.get_ylim()
        ax.fill_betweenx(ylims,limits[0],
                                   limits[1],color='red',alpha=0.2)
        
    def refresh_plot(self):
        """Update the plot using the current settings of the Gui.
        """
        ps = self.settings.plotState
        s = SaPlotState
        # MS main plot (axes[0])
        if ps == s.MS or ps == s.MS_ATD or ps == s.MS_LINES:
            if ps == s.MS_ATD:
                ax,ax2 = self.getDoubleRowAxes()
                self.plotAtds(ax2)
            elif ps == s.MS_LINES:
                ax,ax2 = self.getDoubleColumnAxes()
                self.plotMsLines(ax2)
            else:
                ax = self.getSingleAxis()
            self.plotMassSpectra(ax)
        # ATD main plot (axes[0])
        elif ps == s.ATD or ps == s.ATD_LINES:
            if ps == s.ATD_LINES:
                ax,ax2 = self.getDoubleColumnAxes()
                atds = self.plotAtds(ax)
                self.plotAtdLines(ax2,atds)
            else:
                ax = self.getSingleAxis()
                self.plotAtds(ax)
        # All plots
        elif ps == s.MS_ATD_LINES:
            ax,ax1,ax2,ax3 = self.getQuadAxes()
            self.plotMassSpectra(ax)
            self.plotMsLines(ax2)
            atds = self.plotAtds(ax1)
            self.plotAtdLines(ax3,atds)
        self.draw()
            
    def draw(self):
        """Update the plot panel.
        """
        self.canvas.draw()

             
    #===================================================================
    # Peak picking
    #===================================================================
 
    def toggleClickAndDrag(self,onoff,yourself):
        """Turn on or off click and drag to select a region of the
        mass spectrum (horizontally only).
        :parameter onoff: Boolean for toggle position
        :parameter yourself: SpectralAveragingGui() object
        """
        self.yourself = yourself
        if onoff:
            # reset picking data
            self.pickedValue = None
            self.releasedValue = None
            self.selectedAxis = None
            self.atroposLeft = None
            self.atroposRight = None

            self.refresh_plot()
            self.picker = self.fig.canvas.mpl_connect(
                'button_press_event',self._click)
        else:
            self.fig.canvas.mpl_disconnect(self.picker)
            try: self.fig.canvas.mpl_disconnect(self.releaser)
            except: pass

    def _click(self,event):
        """Check that the click is within the mass spectrum panel, and
        register the mouse click release event (to get position).
        """
        xlabel = event.inaxes.get_xlabel()
        import re
        if re.search('m/z', xlabel):
            self.pickedValue = event.xdata
            self.selectedAxis = event.inaxes
            self.releaser = self.fig.canvas.mpl_connect(
                'button_release_event',self._release)
        else:
            # TODO(gns) - consider removing or changing to warning dialog
            print 'Wrong axis! \nOnly works on mass spectra'

    def _release(self,event):
        if event.inaxes == self.selectedAxis:
            self.releasedValue = event.xdata
            print self.pickedValue,self.releasedValue
        else:
            print 'Release click on same axis!\nBug averted'
            
        # Reset everything
        self.fig.canvas.mpl_disconnect(self.picker)
        self.fig.canvas.mpl_disconnect(self.releaser)
        self.gui.buttonClickAndDrag.SetValue(False)

        self.refresh_plot()


