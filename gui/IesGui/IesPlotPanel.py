"""Plotting area for IesGui()."""

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
import msClasses.MassSpectrum as MassSpectrum
from lib import utils
matplotlib.rc('xtick', labelsize=8) 
matplotlib.rc('ytick', labelsize=8) 

#specific to iesplotpanel
from collections import OrderedDict
import matplotlib.gridspec as gridspec
import IesSettings
from AmphitriteEnums import *

class IesPlotPanel():
    
    def __init__(self,panel,yourself):
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

        self.settings = IesSettings.IesSettings('')
        self.currentPlot = IesPlotState.MASS_SPECTRA

        self.lift = 30

        self.picker = None
        self.pickedValue = None

        self.gui = yourself
        
    def setSettings(self,settings):
        """Set the Gui settings class.
        :parameter settings: IesSettings() object
        """
        self.settings = settings

    def getSingleAxis(self):
        """Recreate plotting area with a single set of axes.
        :returns: Matplotlib Axes instance
        """
        self._preparePlottingSections(1,1)
        self.axes = [None]
        self.axes[0] = self.fig.add_subplot(self.gs[0,0])
        self.axes[0].plot([],[])
        return self.axes[0]
    
        
    def _preparePlottingSections(self,rows,columns):
        """Clear currently plotted information and recreate the plot area.
        :parameter rows: Rows of Matplotlib Axes instances
        :parameter columns: Columns of Matplotlib Axes instances
        """
        self.fig.clf(keep_observers=True)
        self.gs = gridspec.GridSpec(rows,columns)

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

    def getTripleAxes(self):
        """Recreate plotting area with a three sets of axes. Arranged
        as two side by side with one below.
        :returns: List of three Matplotlib Axes instances
        """
        self._preparePlottingSections(3,2)
        self.axes = [None,None,None]
        self.axes[0] = self.fig.add_subplot(self.gs[0:2,0])
        self.axes[1] = self.fig.add_subplot(self.gs[0:2,1])
        self.axes[2] = self.fig.add_subplot(self.gs[2,:])        

        return self.axes[0],self.axes[1],self.axes[2]
        
    
    def plotMassSpectra(self):
        """Change plot area to single set of axes and plot the stacked
        mass spectra.
        """
        ax = self.getSingleAxis()
        self.settings._forPlotMassSpectra(ax,lift=self.lift)
        self.label1dPlots(ax)
        self.currentPlot = IesPlotState.MASS_SPECTRA
        self.draw()

    def plotMassSpectraWidthLimits(self):
        """Run self.plotMassSpectra() and fill in the m/z value regions which
        are to be used for extracting ATDs.
        """
        self.plotMassSpectra()
        species,z = self.settings.getSpeciesAndCharge()
        
        limits = self.settings._getMzLimits(species,z)
        ylims = self.axes[0].get_ylim()
        self.axes[0].fill_betweenx(ylims,limits[0],limits[1],color='red',alpha=0.3)
        self.axes[0].set_ylim(ylims)
        
    def refresh_plot(self):
        """Update the plotting window, given the current settings.
        """
        self.settings.checkboxStates.printStates()
        plotType = self.settings.checkboxStates.getPlotType()
        self.currentPlot = self.settings.checkboxStates.getEnum()
        
        if plotType == 1:
            ax = self.getSingleAxis()
            if self.currentPlot == IesPlotState.MASS_SPECTRA:
                if self.settings.atrOb:
                    self.plotMassSpectraWidthLimits()
                else:
                    self.plotMassSpectra()
 
            elif self.currentPlot == IesPlotState.ATD:
                self.plotCcsDistributions(ax)

            elif self.currentPlot == IesPlotState.CONTOUR:
                self.plotContourPlots(ax)

            elif self.currentPlot == IesPlotState.NONE:
                self.getSingleAxis()
                
            elif self.currentPlot == IesPlotState.CONFORMATIONS:
                self.plotConformationHeights(ax)
                
        elif plotType == 2:
            if self.currentPlot == IesPlotState.CONTOUR_ATD:
                ax1,ax2 = self.getDoubleColumnAxes()
                self.plotCcsDistributions(ax1)
                self.plotContourPlots(ax2)
            elif self.currentPlot == IesPlotState.CONTOUR_CONFORMATIONS:
                ax1,ax2 = self.getDoubleRowAxes()
                self.plotContourPlots(ax1)
                self.plotConformationHeights(ax2)
            elif self.currentPlot == IesPlotState.ATD_CONFORMATIONS:
                ax1,ax2 = self.getDoubleRowAxes()
                self.plotCcsDistributions(ax1)
                self.plotConformationHeights(ax2)
            else:
                print 'plot type = 2 but this option hasnt been implemented yet'

        elif plotType == 3:
            ax1,ax2,ax3 = self.getTripleAxes()
            self.plotCcsDistributions(ax1)
            self.plotContourPlots(ax2)
            self.plotConformationHeights(ax3)
        else:
            print 'plottype:'
            print plotType
        if self.settings.displayPeaks:
            self.drawCcsLines()
        for ax in self.axes:
            ax.set_yticks([])
        self.draw()
            
    def on_pick(self,event):
        """Clicking on the plotting area (still to be implemented).
        """
        # TODO(gns) - Still to be implemented
        if True: #self.pickingActive:
            #my attempt to make a new function
            plotState = self.settings.checkboxStates.getEnum()
            # print event.mouseevent.data
            # print dir(event.mouseevent)
            # print '==============================='
            # print dir(event)
            # copied directly
            # if self.pickedAlready:
            #     line = self.axMs.lines.pop(len(self.axMs.lines)-1)
            #     self.tempAnnotation.set_visible(False)
            # self.pickedValue = event.mouseevent.xdata
            # yval = self.ms.yvals[utils.closest(self.pickedValue,self.ms.xvals)]
            # peakId = sorted(self.ms.gPeaks.keys())[-1]+1
            # self.axMs.axvline(self.pickedValue,color='k')
            # self.tempAnnotation = self.axMs.annotate(str(peakId),[self.pickedValue,yval])
            # self.draw()
            
            # self.pickedAlready = True
    
      
    def draw(self):
        """Update plotting area.
        """
        self.canvas.draw()

    def plotContourPlots(self,ax):
        """Plot the 3D data as contour plots. Automatically switches between ATD and CCS
        depending on whether a calibration has been added.
        :parameter ax: Matplotlib Axes instance
        """
        if self.settings.calibrationOb:
            # get the data
            ccsAxis, matrices = self.settings.getCcsAxisAndGrid()

            # plot the data
            for i,matrix in enumerate(matrices):
                ax.imshow(matrix, extent=[ccsAxis[0],ccsAxis[-1],i*10,i*10+10],
                          aspect='auto',origin='lower')
                if i:
                    ax.axhline(i*10)

            # plotting parameters
            ax.set_ylabel('Relative $m/z$')
            ax.set_xlabel('CCS ($\AA^2$)')
            ax.set_yticks([])

            ax.autoscale(ax)
            ccsLims = self.autoAxesImshow(ccsAxis,matrices)
            ax.set_xlim(ccsLims)
            self.labelContourPlots(ax)
        else:
            '''see explanation in plotCcsDistributions'''
            tdAxis,matrices = self.settings.getAtdAxisAndGrid()
            for i,matrix in enumerate(matrices):
                ax.imshow(matrix,
                          extent=[tdAxis[0],tdAxis[-1],i*10,i*10+10],
                          aspect='auto',origin='lower')
                if i:
                    ax.axhline(i*10)
            # plotting parameters
            ax.set_ylabel('Relative $m/z$')
            ax.set_xlabel('t$_d$')
            ax.set_yticks([])

            ax.autoscale(ax)
            tdLims = self.autoAxesImshow(tdAxis,matrices)
            ax.set_xlim(tdLims)
            self.labelContourPlots(ax)

    def plotCcsDistributions(self,ax):
        """Plot CCS distributions.
        :parameter ax: Matplotlib Axes instance
        """
        ax.clear()
        if self.settings.calibrationOb:
            ccsAxes,ccsLines = self.settings.getCcsLines()
            # Normalising using base peak
            i = 0
            for ccsAxis,line in zip(ccsAxes,ccsLines):
                line = line/line.max()*100
                ax.plot(ccsAxis,line+(i*self.lift),color='k')
                i += 1
            self.autoAxesCcsDistributions(ax,ccsAxes,ccsLines)
            self.label1dPlots(ax)
            ax.set_ylabel('Intensity')
            ax.set_xlabel('CCS ($\AA^2$)')
        else:
            '''horribly hacked to add the ability to show ATDs
            when there hasn't been a calibration file added.
            This is so that we can compare different files where different 
            wave height and velocity values were used'''
            xaxes,lines = self.settings.getAtdLines()
            for i,(xaxis,line) in enumerate(zip(xaxes,lines)):
                line = line/line.max()*100
                ax.plot(xaxis,line+(i*self.lift),color='k')
            self.autoAxesCcsDistributions(ax,xaxes,lines)
            self.label1dPlots(ax)
            ax.set_ylabel('Intensity')
            ax.set_xlabel('t$_d$')
            
    def calculateAutoAxes(self,ccsAxes,ccsLines):
        """Create automatic x axis limits for CCS distributions.
        :parameter ccsAxes: CCS distribution x axis
        :parameter ccsLines: CCS distribution y axis
        :returns: [lowerLimit,upperLimit]
        """
        minX = min([utils.findFirstNonZeroXvalue(axis,line,zero=2) for axis,line in zip(ccsAxes,ccsLines)])
        #maxX = ax.get_xlim()[1]
        maxX = max([ utils.findLastNonZeroXvalue(axis,line)for axis,line in zip(ccsAxes,ccsLines)])
        return [minX,maxX]
        
    def autoAxesCcsDistributions(self,ax,ccsAxes,ccsLines):
        """Plot CCS distributions with automatic x axis limits.
        :parameter ax: Matplotlib Axes instance
        :parameter ccsAxes: CCS distribution x axis
        :parameter ccsLines: CCS distribution y axis
        """
        minX,maxX = self.calculateAutoAxes(ccsAxes,ccsLines)
        ax.set_xlim([minX,maxX])

    def autoAxesImshow(self,ccsAxis,matrices):
        """Create automatic CCS axis limits for 3D contour plots.
        :parameter ccsAxes: CCS axis
        :parameter matrices: Matrix of intensities
        :returns: [lowerLimit,upperLimit]
        """
        ccsLines = [ np.sum(matrix,axis=0) for matrix in matrices ]
        ccsAxes = [ ccsAxis for x in matrices ]
        minX,maxX = self.calculateAutoAxes(ccsAxes,ccsLines)
        return [minX,maxX]
        
    def label1dPlots(self,ax):
        """Label stacked mass spectra or ATDs/CCSDs (e.g. '5V', '10V' in
        voltage ramp experiments).
        :parameter ax: Matplotlib Axes instance
        """
        # needs a fully populated list of values even if just blank strings
        values = self.settings.values
        yheights = [i*self.lift + (self.lift*0.1) for i in xrange(len(values))]

        xlims = ax.get_xlim()
        x_range = xlims[1]-xlims[0]
        xposition = (x_range*0.95) + xlims[0]
        
        self._labelPlots(ax,'right',xposition,yheights,'k')


    def labelContourPlots(self,ax):
        """Draw individual contour plot labels (e.g. '5 V', '10 V' in
        volage ramps).
        :parameter ax: Matplotlib Axes instance
        """
        yheights = [5 + (10*i) for i in xrange(len(self.settings.values))]
        xlims = ax.get_xlim()
        xpos = xlims[0] + (xlims[1]-xlims[0])*0.05
        alignment = 'left'
        self._labelPlots(ax,alignment,xpos,yheights,'cyan')

    def _labelPlots(self,ax,alignment,xpos,yheights,colour):
        """Draw labels for spectra or contour plots (stacked) in a plot area.
        Function gets the values and units for labels from IesSettings() object.
        :parameter ax: Matplotlib Axes instance
        :parameter alignment: Text alignment for annotation
        :parameter xpos: Horizontal position for the labels
        :parameter yheights: Vertical positions for the labels (list)
        :parameter colour: Matplotlib colour for the labels
        """
        units = self.settings.units
        ints = True
        for i,value in enumerate(self.settings.values):
            if value:
                if value%1 != 0:
                    ints = False
        for i,value in enumerate(self.settings.values):
            if value:
                if ints:
                    s = "%.0f %s" %(value,units)
                else:
                    s = "%s %s" %(value,units)
                ax.annotate(s, (xpos,yheights[i]),
                            horizontalalignment=alignment,color=colour)

             
    #===========================================================================
    # Peak picking
    #===========================================================================

    def oneClickPicking(self,onoff):
        """Turn on one click (to add a peak) mode. This activates the
        peak picking function (self._oneClick()).
        :parameter onoff: Boolean for turning on or off peak picking
        """
        if onoff:
            self.picker = self.fig.canvas.mpl_connect('button_press_event',self._oneClick)
        else:
            self.fig.canvas.mpl_disconnect(self.picker)

    def _oneClick(self,event):
        """Function to process a click in one click per
        peak mode.
        """
        xlabel = event.inaxes.get_xlabel()
        import re
        if re.search('CCS', xlabel):
            ccs = event.xdata
            self.settings.addConformation(ccs)
            self.refresh_plot()
        else:
            # TODO(gns) - Consider removing or turning into a warning dialog
            print 'Axis label not matched'

    def togglePicking(self,onoff):
        """Turn off or on peak picking. In the form of adding a peak
        to the GUI when the toggle is turned back to off.
        :parameter onoff: Boolean 
        """
        if onoff:
            self.picker = self.fig.canvas.mpl_connect('button_press_event',self._togglePicking)
        else:
            self.fig.canvas.mpl_disconnect(self.picker)
            self.pickedValue = None

    def _togglePicking(self,event):
        """Sub function for self.togglePicking().
        """
        #print event.xdata
        xlabel = event.inaxes.get_xlabel()
        import re
        self.settings.removeConformation(self.pickedValue)
        if re.match('CCS', xlabel):
            ccs = event.xdata
            self.pickedValue = ccs
            self.settings.addConformation(ccs)
            self.refresh_plot()
        
        

    def drawCcsLines(self):
        """Draw vertical lines corresponding to the CCS value of the user defined
        conformation positions.
        """
        enum = self.settings.checkboxStates.getEnum()
        ps = IesPlotState
        justax1 = [ps.CONTOUR,ps.ATD,ps.ATD_CONFORMATIONS,ps.CONTOUR_CONFORMATIONS]
        ax1ax2 = [ps.CONTOUR_ATD,ps.CONTOUR_ATD_CONFORMATIONS]
        if enum in justax1 or enum in ax1ax2:        
            for i,ccs in enumerate(self.settings.conformations):
                ccs = round(ccs,0)
                self.axes[0].axvline(ccs,label='%d $\AA^2$'%(int(ccs)),color=utils.colourList[i])
                if enum in ax1ax2:
                    self.axes[1].axvline(ccs,label='%d $\AA^2$'%(int(ccs)),color=utils.colourList[i])
            self.axes[0].legend()
        


    def plotConformationHeights(self,ax):
        """Track the abundance of the different conformations across
        the data files and plot them.
        """
        # gather the data
        ccsAxes,ccsLines = self.settings.getCcsLines()
        ccsDic = OrderedDict()
        for ccs in self.settings.conformations:
            ccsDic[ccs] = []
            for xaxis,line in zip(ccsAxes,ccsLines):
                height = line[utils.closest(xaxis,ccs)]
                ccsDic[ccs].append(height)

        # plot it
        if not '' in self.settings.values and len(self.settings.values):
                valueList = self.settings.values
                ax.set_xlabel(self.settings.units)
        else:
            valueList = [ x for x in range(len(ccsDic.keys())) ]
            ax.set_xlabel('Arbitary progression')

        ax.set_xticks(valueList)

        for i,ccs in enumerate(ccsDic.keys()):
            ax.plot(valueList, ccsDic[ccs], label=str(ccs),color=utils.colourList[i])

            
            
        
                
        
