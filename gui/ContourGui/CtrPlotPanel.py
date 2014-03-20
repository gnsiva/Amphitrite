"""Class for plotting 3D data as heatmaps for ContourGui()."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import matplotlib
import matplotlib as mpl
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
from imClasses import ChargeStatePeak

import matplotlib.font_manager
prop = matplotlib.font_manager.FontProperties(size='small')

class CtrPlotPanel():
    """
    :parameter panel: wx panel to contain this plot
    """
    def __init__(self,panel):
        self.dpi = 80
        self.fig = Figure((3.0, 2.0), dpi=self.dpi)
        self.canvas = FigCanvas(panel, -1, self.fig)
        self.gs = gridspec.GridSpec(1,1)
        self.axes = [None]
                
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

        self.backgroundColour = 'white'
        self.vmax = 100
        self.cmap = mpl.colors.LinearSegmentedColormap.from_list(
                    'my_cmap',['white','b','k','k'],256)
        self.cmap2 = mpl.colors.LinearSegmentedColormap.from_list(
                    'my_cmap',['white','r','k','k'],256)
        self.differenceCmap = mpl.colors.LinearSegmentedColormap.from_list(
                    'my_cmap',['pink','r','r','white','b','b','c'],256)
        
    #====================
    # Adding other panel components        
    def setSettings(self,settings):
        """
        :parameter settings: CtrSettings() object
        """
        self.settings = settings
    def setGui(self,gui):
        """Set the GUI object to enable accessing its attributes.
        :parameter gui: ContourGui() object
        """
        self.gui = gui
    #====================
    ############################################################
    # Axes stuff
    ############################################################        
    def _preparePlottingSections(self,rows,columns):
        """Set up the plot panel using GridSpec.
        """
        self.fig.clf(keep_observers=True)
        self.gs = gridspec.GridSpec(rows,columns)

    def getSingleAxis(self):
        """Set up the plot panel with a single plot area and return it.
        :returns: Matplotlib Axes object
        """
        self._preparePlottingSections(1,1)
        self.axes = [None]
        self.axes[0] = self.fig.add_subplot(self.gs[0,0])
        self.axes[0].plot([],[])
        return self.axes[0]
        
    def getDoubleAxes(self):
        """Set up the plot panel with two plot areas and return them.
        :returns: List of two matplotlib Axes objects
        """        
        self._preparePlottingSections(3,1)
        self.axes = [None,None]
        self.axes[0] = self.fig.add_subplot(self.gs[0:2])
        self.axes[1] = self.fig.add_subplot(self.gs[2])
        return self.axes[0],self.axes[1]

    # Axes stuff
    ############################################################

    def setVmax(self,vmax):
        """The value used for the maximum color saturation value.
        """
        if utils.isNumber(vmax):
            self.vmax = float(vmax)
    
    ################################################################
    ################################################################
    # Extraction limits functions

    ################################################################
    # Mass spectra (Extraction limits)
    def plotMsExtractionLimits(self,ax,limitsD=None):
        """Fill in the area to be extracted from the m/z axis
        as arrival time or CCS.
        :parameter ax: Matplotlib Axes object to use
        :parameter limitsD: Set to True to use automatic maximum width size
        for extraction
        """
        if self.gui.checkboxAutoLimits.IsChecked():
            self.plotAtroposSelectionAutoWidth(ax,limitsD)
        else:
            self.plotAtroposSelection(ax)

    def plotAtroposSelection(self,ax):
        """Fill in the area to be extracted from the m/z axis
        as arrival time or CCS.
        :parameter ax: Matplotlib Axes object to use
        """
        ylims = ax.get_ylim()
        # if same limits all species
        for sp in self.settings.atrOb.species.keys():
            for z in self.settings.atrOb.simulatedSpecies[sp].charges:
                limits = self._getMzLimits(sp,z)
                ax.fill_betweenx(ylims,limits[0],
                                   limits[1],color='red',alpha=0.1)
                
    def plotAtroposSelectionAutoWidth(self,ax,limitsD):
        """Uses limits generated by: getAutoPeakLimits() instead of relative to 
        atropos generated fwhm values.
        :parameter ax: Matplotlib Axes object to use
        :parameter limitsD: Set to True to use automatic maximum width size
        for extraction        
        """
        ylims = ax.get_ylim()
        for sp,d in limitsD.items():
            mass = self.settings.imOb.species[sp].mass
            for z,limits in d.items():
                mz = utils.get_mz(mass,z)
                ax.fill_betweenx(ylims,limits[0],
                                   limits[1],color='red',alpha=0.1)

    def _getMzLimits(self,species,charge):
        """Get the m/z limits for a specific species charge state peak from
        the Atropos fit.
        :parameter species: Name of molecular species
        :parameter charge: Charge state to get limits for
        :returns: limits - in the form of [lower,upper]
        """
        limits = self.settings.atrOb.simulatedSpecies[species].getPeakLimits(
            charge,self.settings.widthL[species],self.settings.widthR[species])
        return limits
    ################################################################
    # Ion mobility (Extraction limits)
    def getAutoPeakLimits(self):
        # Get all the mzs of the various species charge states
        zsD,massD = self.settings.getChargeAndMassDictionaries()
        mzsD = OrderedDict()
        chargeStatePeakObs = []
        allMzs = []
        for sp, zs in zsD.items():
            mzsD[sp] = []
            for z in zs:
                mz = utils.get_mz(massD[sp],z)
                mzsD[sp].append(mz)
                allMzs.append(mz)

        # See if they overlap
        limitsD = OrderedDict()
        for sp,zs in zsD.items():
            limitsD[sp] = OrderedDict()
            for z in zs:
                csOb = ChargeStatePeak.ChargeStatePeak(massD[sp],z,zs)
                limits = csOb.getLimits(allMzs)
                limitsD[sp][z] = limits
        return limitsD
                
    ################################################################
    ################################################################

    def plotAtdForRefreshPlot(self,ax):
        '''Handles: Plain, log scale and difference atd contour plots '''
        matrix,x,y = self.settings.imOb._getMatrixWithinLimits()
        if self.gui.choiceScale.GetSelection() == 1:
            matrix = self.logScaleMatrix(matrix)

        # Difference plot
        if self.gui.checkboxDifferencePlot.IsChecked():
            matrix2,x2,y2 = self.settings.imOb2._getMatrixWithinLimits()

            if self.gui.choiceScale.GetSelection() == 1:
                matrix2 = self.logScaleMatrix(matrix2)

            x = np.concatenate((x,x2))
            y = np.concatenate((y,y2))
            matrix = matrix - matrix2
            cmap = self.differenceCmap
            vmin = self.vmax*-1
        # Single plot
        else:
            cmap = self.cmap
            vmin = 0

        ax.imshow(matrix, origin=[0,0],aspect='auto', 
                   extent=[x.min(),x.max(),y.min(),y.max()],
                  cmap=cmap,vmin=vmin,vmax=self.vmax)

        ax.set_xlabel('$m/z$')
        ax.set_ylabel('t$_d$ ($ms$)')

    def plotCcsVsMzForRefreshPlot(self,ax,limitsD):
        '''Deals with CCS vs. mass spectrum and
        charge state clemmer CCS plot'''
        ########
        # Extracting data in the imOb
        # Auto peak limits
        if self.gui.checkboxAutoLimits.IsChecked():
            self.settings.imOb.generateSpeciesSlicesExplicit(limitsD)
            if self.settings.imOb2:
                self.settings.imOb2.generateSpeciesSlicesExplicit(limitsD)
        # Relative to FWHM zstate peak limits
        else:
            for i,sp in enumerate(self.settings.species):
                if i:
                    self.settings.imOb.generateSpeciesSlicesFwhm(
                        sp,self.settings.widthL[sp],self.settings.widthR[sp])
                    if self.settings.imOb2:
                        self.settings.imOb2.generateSpeciesSlicesFwhm(
                            sp,self.settings.widthL[sp],self.settings.widthR[sp])
                        
                    if  self.gui.checkboxDifferencePlot.IsChecked():
                        self.settings.imOb.generateSpeciesSlicesFwhm(
                            sp,self.settings.widthL[sp],self.settings.widthR[sp])
                        if self.settings.imOb2:
                            self.settings.imOb2.generateSpeciesSlicesFwhm(
                                sp,self.settings.widthL[sp],self.settings.widthR[sp])

        # getting data
        dataSlices = self.settings.imOb.getDataSlices()
        dataSlices = self.scalingDataSlices(dataSlices)
        vmin = 0
        cmap = self.cmap

        # difference plot changes
        if self.gui.checkboxDifferencePlot.IsChecked():
            dataSlices2 = self.settings.imOb2.getDataSlices()
            dataSlices2 = self.scalingDataSlices(dataSlices2)
            dataSlices = self.getDataSliceDifferences(dataSlices,dataSlices2)
            vmin = self.vmax*-1
            cmap = self.differenceCmap

        #######
        # Plotting
        if not self.gui.checkboxShowChargeStatePlot.IsChecked():
            self.settings.imOb.plotCcsHeatMapFromDataSlices(
                ax,self.settings.calibrationOb,dataSlices,
                cmap=cmap,vmin=vmin,vmax=self.vmax)

            # set xlims
            xlims = self.getContourXlims()
            ax.set_xlim(xlims)
            ax.set_xlabel('$m/z$')
            species = 0

        # Clemmer charge state plot
        else:
            species = self.gui.choiceChargeStatePlotSpecies.GetStringSelection()
            self.settings.imOb.plotChargeStateContourPlotsFromDataSlices(
                ax,self.settings.calibrationOb,dataSlices,species,
                vmin=vmin,vmax=self.vmax,cmap=cmap)
            ax.set_xlabel('Charge State')
            ax.set_xticks([])


        ax.set_ylabel('CCS ($\AA^2$)')
        ax.set_axis_bgcolor(self.backgroundColour)

        if self.gui.checkboxShowPeakTops.IsChecked():
            self.plotPeakTops(ax,species)

            
    def getDataSliceDifferences(self,dataSlices1,dataSlices2):
        for sp,spSlices in dataSlices1.items():
            for z,dataSlice in spSlices.items():
                dataSlice.matrix -= dataSlices2[sp][z].matrix
        return dataSlices1
        
    def plotMsForRefreshPlot(self,limitsD):
        ax,ax2 = self.getDoubleAxes()
        colour1 = 'k'
        
        if self.gui.checkboxDifferencePlot.IsChecked():
            colour1 = 'b'
            self.settings.imOb2.massSpectrum.plot(ax2,color='r')
            
        self.settings.imOb.massSpectrum.plot(ax2,color=colour1)

        ax2.set_ylim([0,105])

        if self.settings.atrOb:
            self.plotMsExtractionLimits(ax2,limitsD)
        
        return ax,ax2

        
    def refresh_plot(self):
        # Autopeak limits (maximising charge state strip width)
        limitsD = None
        if self.gui.checkboxAutoLimits.IsChecked():
            limitsD = self.getAutoPeakLimits()

        # Mass Spectrum Panel
        if self.gui.checkboxShowMassSpectrum.IsChecked():
            ax,ax2 = self.plotMsForRefreshPlot(limitsD)
        else:
            ax = self.getSingleAxis()
            ax2 = None
        
        # ATD vs m/z
        if self.gui.radioboxPlotPanel.GetSelection() == 0:
            self.plotAtdForRefreshPlot(ax)
            if ax2:
                ax.set_xlim(ax2.get_xlim())
                
        # CCS vs m/z
        elif self.gui.radioboxPlotPanel.GetSelection() == 1:
            self.plotCcsVsMzForRefreshPlot(ax,limitsD)

        self.draw()

                
        
    def plotPeakTops(self,ax,species=0):
        '''Handles all possibilities for displaying the peak tops:
        CCS vs. m/z (not clemmer plot) -> plotCcsVsMzPeakTops'''
        #### NOT FINISHED ####

        # Get information from gui
        smths = self.gui.textCtrlSmoothes.GetValue()
        wlen = self.gui.textCtrlWlen.GetValue()
        smths, wlen = int(smths), int(wlen)     
        limit = float(self.gui.textCtrlLimit.GetValue())

        
        if not self.gui.checkboxShowChargeStatePlot.IsChecked():
            for i,sp in enumerate(self.settings.species):
                if i:
                    # CCS vs. m/z plot
                    if self.gui.radioboxPlotPanel.GetSelection() == 1:
                        # NEED AN IF STATEMENT FOR SUBSTRACTION PLOT COLOUR
                        self.plotCcsVsMzPeakTops(ax,sp)
                    # ATD <=========== still to do
        else:
            # Clemmer plot (CCS ONLY AT THIS POINT)
            self.settings.imOb.plotChargeStateContourPeaktops(
                ax,self.settings.calibrationOb,species,smths,
                wlen,limit=limit,dataType='ccs',colour='green')
            # need if loops for ATD version and subtraction plot
                        
            
    def plotCcsVsMzPeakTops(self,ax,species,colour='green'):
        """
        :parameter ax: Matplotlib Axes() object
        :parameter species: Molecular species name
        :parameter colour: Matplotlib compatible colour string
        """
        for z in self.settings.imOb.species[species].charges:
            dataSlice = self.settings.imOb.dataSlices[species][z]
            smths = self.gui.textCtrlSmoothes.GetValue()
            wlen = self.gui.textCtrlWlen.GetValue()
            smths, wlen = int(smths), int(wlen)        

            mz = dataSlice.getMsApex()
            limit = float(self.gui.textCtrlLimit.GetValue())

            ccsPeaks = dataSlice.getCcsPeaks(smths,wlen,mz,
                                             self.settings.calibrationOb,
                                             limit=limit)
            ax.scatter([mz]*len(ccsPeaks),ccsPeaks,marker='o',color=colour,alpha=0.5)
        
    def getContourXlims(self):
        """Get the m/z limits for the contour plot panel. Otherwise it
        defaults to the extremities of the extracted Charge state slices
        :returns: xlims - as [lower,upper]
        """
        if self.gui.checkboxShowMassSpectrum.IsChecked():
            xlims = self.axes[1].get_xlim()
        else:
            if self.settings.imOb:
                xlimsL1 = self.settings.imOb.xaxis[0]
                xlimsR1 = self.settings.imOb.xaxis[-1]
            # should probably add another if for imOb2
            xlims = [xlimsL1,xlimsR1]
        return xlims
                

    def refreshMsPlotOnly(self):
        """Update the mass spectrum panel of the plot alone.
        """
        if self.gui.checkboxShowMassSpectrum.IsChecked():
            self.axes[1].clear()
            self.draw()
            self.settings.imOb.massSpectrum.plot(self.axes[1])
            if self.settings.atrOb:
                self.plotAtroposSelection(self.axes[1])
            self.draw()

    def updateColourMap(self):
        """Update the colourmaps depending on what the background
        background colour is.
        """
        if self.backgroundColour == 'white':
            self.cmap = mpl.colors.LinearSegmentedColormap.from_list(
                'my_cmap',['white','b','k','k'],256)
            self.cmap2 = mpl.colors.LinearSegmentedColormap.from_list(
                'my_cmap',['white','r','k','k'],256)
        elif self.backgroundColour == 'black':
            self.cmap = mpl.colors.LinearSegmentedColormap.from_list(
                'my_cmap',['k','b','w','w'],256)
            self.cmap2 = mpl.colors.LinearSegmentedColormap.from_list(
                'my_cmap',['k','r','w','w'],256)

            
    def scalingDataSlices(self,dataSlices):
        """Change scaling between log and linear for intensity of dataSlices using
        selection in ContourGui().
        :parameter dataSlices: Dictionary of imClasses.DataSlice() objects, in the
        form of d[speciesName][chargeState]
        """
        if self.gui.choiceScale.GetSelection() == 1:
            dataSlices = self.logScaleDataSlices(dataSlices)
        # TODO(gns) - Shouldn't you do something to change back to non log scale?
        return dataSlices

    def logScaleMatrix(self,matrix):
        """Natural log matrix, do various manipulations to keep the
        values between 0-100.
        :parameter matrix: 2D numpy array of intensity values
        """
        matrix = np.log(matrix) 
        matrixMin = np.min(matrix[np.isneginf(matrix)!=True])
        matrixMax = np.max(matrix)

        # correct for values under 0
        matrix = matrix + matrixMin*-1

        # Deal with non numerical results
        matrix[np.isnan(matrix)] = 0.0
        matrix[np.isneginf(matrix)] = 0.0

        # scale the whole thing to 0-100
        matrix = (matrix/(matrixMax+matrixMin*-1))*100.

        return matrix

    def logScaleDataSlices(self,dataSlices):
        """Natural log matrices, do various manipulations to keep the
        values between 0-100.
        :parameter dataSlices: Dictionary of imClasses.DataSlice() objects, in the
        form of d[speciesName][chargeState]
        """
        loggedMin = 0
        loggedMax = 0
        for sp,spSlices in dataSlices.items():
            for z,dataSlice in spSlices.items():
                matrix = np.log(dataSlices[sp][z].matrix)

                # get minimum value for all dataSlices
                # You need to ignore -inf values returned from logging 0
                matrixMin = np.min(matrix[np.isneginf(matrix)!=True])
                if matrixMin < loggedMin:
                    loggedMin = matrixMin

                # log the maximum value for all the dataSlices for normalising to 100 later
                matrixMax = np.max(matrix)
                if matrixMax > loggedMax:
                    loggedMax = matrixMax


        for sp,spSlices in dataSlices.items():
            for z,dataSlice in spSlices.items():
                # Log the matrix
                dataSlices[sp][z].matrix = np.log(dataSlices[sp][z].matrix)
                # correct for values under 0
                dataSlices[sp][z].matrix = dataSlices[sp][z].matrix + loggedMin*-1
                
                # Deal with non numerical results
                dataSlices[sp][z].matrix[np.isnan(dataSlices[sp][z].matrix)] = 0.0
                dataSlices[sp][z].matrix[np.isneginf(dataSlices[sp][z].matrix)] = 0.0
                
                # scale the whole thing to 0-100
                dataSlices[sp][z].matrix = (dataSlices[sp][z].matrix/(loggedMax+loggedMin*-1))*100.

        return dataSlices

    def draw(self):
        """Re-draw plot area.
        """
        self.canvas.draw()



        




            
            
