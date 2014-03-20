import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt
import numpy as np
import wx
import classes.MassSpectrum as MassSpectrum
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
        self.settings = settings

    def getSingleAxis(self):
        self._preparePlottingSections(1,1)
        self.axes = [None]
        self.axes[0] = self.fig.add_subplot(self.gs[0,0])
        self.axes[0].plot([],[])
        return self.axes[0]
    
        
    def _preparePlottingSections(self,rows,columns):
        self.fig.clf(keep_observers=True)
        self.gs = gridspec.GridSpec(rows,columns)

    def getDoubleColumnAxes(self):
        self._preparePlottingSections(1,2)
        self.axes = [None,None]
        self.axes[0] = self.fig.add_subplot(self.gs[0])
        self.axes[1] = self.fig.add_subplot(self.gs[1])

        return self.axes[0],self.axes[1]

    def getDoubleRowAxes(self):
        self._preparePlottingSections(2,1)
        self.axes = [None,None]
        self.axes[0] = self.fig.add_subplot(self.gs[0])
        self.axes[1] = self.fig.add_subplot(self.gs[1])

        return self.axes[0],self.axes[1]

    def getTripleAxes(self):
        self._preparePlottingSections(3,2)
        self.axes = [None,None,None]
        self.axes[0] = self.fig.add_subplot(self.gs[0:2,0])
        self.axes[1] = self.fig.add_subplot(self.gs[0:2,1])
        self.axes[2] = self.fig.add_subplot(self.gs[2,:])        

        return self.axes[0],self.axes[1],self.axes[2]
        
    
    def plotMassSpectra(self):
        ax = self.getSingleAxis()
        self.settings._forPlotMassSpectra(ax,lift=self.lift)
        self.label1dPlots(ax)
        self.currentPlot = IesPlotState.MASS_SPECTRA
        self.draw()

    def plotMassSpectraWidthLimits(self):
        self.plotMassSpectra()

        species,z = self.settings.getSpeciesAndCharge()
        
        limits = self.settings._getMzLimits(species,z)
        ylims = self.axes[0].get_ylim()
        self.axes[0].fill_betweenx(ylims,limits[0],limits[1],color='red',alpha=0.3)
        self.axes[0].set_ylim(ylims)
        
    def refresh_plot(self):
        print 'refreshing plot!'
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
                print 'plot type = 2 but this aint done yet'

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
        self.canvas.draw()

    def plotConformations(self,ax):
        ##########################
        # TODO
        ##########################
        ax.plot([1,2,3],[1,2,3])
        self.draw()

    def plotContourPlots(self,ax):
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
            '''see rant in plotCcsDistributions'''
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
        minX = min([utils.findFirstNonZeroXvalue(axis,line,zero=2) for axis,line in zip(ccsAxes,ccsLines)])
        #maxX = ax.get_xlim()[1]
        maxX = max([ utils.findLastNonZeroXvalue(axis,line)for axis,line in zip(ccsAxes,ccsLines)])
        return [minX,maxX]
        
    def autoAxesCcsDistributions(self,ax,ccsAxes,ccsLines):
        minX,maxX = self.calculateAutoAxes(ccsAxes,ccsLines)
        ax.set_xlim([minX,maxX])

    def autoAxesImshow(self,ccsAxis,matrices):
        # ccsLims = [ccsAxis[0],ccsAxis[-1]]
        # for matrix in matrices:
        #     collapsed = np.sum(matrix,axis=0)
        #     thismin = utils.findFirstNonZeroXvalue(ccsAxis,collapsed)
        #     thismax = utils.findLastNonZeroXvalue(ccsAxis,collapsed)
        #     if ccsLims[0]<thismin:
        #         ccsLims[0] = thismin
        #     if ccsLims[1]>thismax:
        #         ccsLims[1] = thismax
        # return ccsLims
        ccsLines = [ np.sum(matrix,axis=0) for matrix in matrices ]
        ccsAxes = [ ccsAxis for x in matrices ]
        minX,maxX = self.calculateAutoAxes(ccsAxes,ccsLines)
        return [minX,maxX]
        
    def label1dPlots(self,ax):
        '''needs a fully populated list of values even if just blank strings'''
        values = self.settings.values
        yheights = [i*self.lift + (self.lift*0.1) for i in xrange(len(values))]

        xlims = ax.get_xlim()
        x_range = xlims[1]-xlims[0]
        xposition = (x_range*0.95) + xlims[0]
        
        self._labelPlots(ax,'right',xposition,yheights,'k')


    def labelContourPlots(self,ax):
        yheights = [5 + (10*i) for i in xrange(len(self.settings.values))]
        xlims = ax.get_xlim()
        xpos = xlims[0] + (xlims[1]-xlims[0])*0.05
        alignment = 'left'
        self._labelPlots(ax,alignment,xpos,yheights,'cyan')

    def _labelPlots(self,ax,alignment,xpos,yheights,colour):
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
        if onoff:
            self.picker = self.fig.canvas.mpl_connect('button_press_event',self._oneClick)
        else:
            self.fig.canvas.mpl_disconnect(self.picker)

    def _oneClick(self,event):
        xlabel = event.inaxes.get_xlabel()
        import re
        if re.search('CCS', xlabel):
            ccs = event.xdata
            self.settings.addConformation(ccs)
            self.refresh_plot()
        else:
            print 'Axis label not matched'

    def togglePicking(self,onoff):
        if onoff:
            self.picker = self.fig.canvas.mpl_connect('button_press_event',self._togglePicking)
        else:
            self.fig.canvas.mpl_disconnect(self.picker)
            self.pickedValue = None

    def _togglePicking(self,event):
        print 'pick noticed'
        print event.xdata
        xlabel = event.inaxes.get_xlabel()
        import re
        self.settings.removeConformation(self.pickedValue)
        if re.match('CCS', xlabel):
            ccs = event.xdata
            self.pickedValue = ccs
            self.settings.addConformation(ccs)
            self.refresh_plot()
        
        

    def drawCcsLines(self):
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

            
            
        
                
        
