"""Class for holding full ion mobility data files, containing data for
multiple charge states."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import ImData
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from lib import utils
from lib.utils import localMaxima
from imClasses import DataSlice
import collections
import lib.RawFileProcessor as RawFileProcessor
import os

class Im(ImData.ImData):

    def __init__(self):
        ImData.ImData.__init__(self)
        self.dataSlices = collections.OrderedDict()
        self.species = collections.OrderedDict()

        self.minCCS = 10000.0
        self.maxCCS = 0.0

    ###############################################################
    # Data loading functions
    ###############################################################
    def loadFolderAuto(self,folderName):
        """Automatically carries out the processing depending on what is input.
        If the folder is MassLynx raw data file, it will detect and process that
        (windows only), if the folder contains correctly named text files it will
        do the same. Finally if the folderName points to an Amphitrite data file
        ('.a') it will detect this and process accordingly.
        :parameter folderName: Absolute path to data folder or file
        """
        '''Automatically checks the type of data in the folder
        and processes it to load the ImObj'''

        processedFilenames = ['MassMobilityXaxis','MassMobilityYaxis','MassMobility']
        textFilenames = [x+'.txt' for x in processedFilenames]
        amphiFilenames = [x+'.amphi' for x in processedFilenames]

        if os.path.isdir(folderName):
            allFound = True
            processed = False

            # text files
            if not processed:
                allFound = True
                if self._checkForFiles(folderName, textFilenames):
                    print 'Processing as textfiles'
                    self.loadFolder(folderName)
                    processed = True

            # raw file process
            if not processed:
                print 'Processing from raw files'
                self.loadRawFolder(folderName)
                processed = True

        # single amphitrite data file (.a file)
        elif os.path.isfile(folderName):
            self.loadAmphiFile(folderName)
            processed = True
        else:
            print 'folder not found', folderName
            processed = False
        if not processed:
            print 'Not processed, disaster', folderName


    def loadAmphiFile(self,filename):
        """Load data from amphitrite data file ('.a').
        :parameter filename: Absolute path to data file
        """
        dataList = utils.unPickleAmphitriteProject(filename)
        if dataList:
            self.setDataFromAmphiExtract(dataList)

    def setDataFromAmphiExtract(self,dataList):
        """Set the m/z and arrival time axis as well as the
        intensity matrix from a list of those arrays.
        :parameter dataList: List of form [xaxis,yaxis,intensityArray]
        """
        self.setAxisX(dataList[0])
        self.setAxisY(dataList[1])
        self.setMatrix(dataList[2])
        self.xaxisUnaltered = self.xaxis.copy()
        self.yaxisUnaltered = self.yaxis.copy()


    def _checkForFiles(self,folder,fileNames):
        """Check if the required text files (correctly named)
        are in the given folder.
        :parameter folder: Absolute path to directory containing data
        :parameter fileNames: List of filenames to check for in the directory
        """
        allFound = True
        # TODO(gns) - this is nonsense, its probably supposed to be:
        # TODO(gns) - len(fileNames) == 0:
        if type(fileNames).__name__ == 'int':
            fileNames = [fileNames]
            
        for f in fileNames:
            if not f in os.listdir(folder):
                allFound = False
        return allFound

    def _loadX(self,filename):
        """Load m/z axis from text file.
        :parameter filename: Absolute path to data file
        """
        xaxis = np.fromfile(filename,sep=',')
        self.setAxisX(xaxis[:-2])  # NO IDEA WHERE THE -2 CAME FROM
        
    def _loadY(self,filename):
        """Load arrival time axis from text file.
        :parameter filename: Absolute path to data file
        """
        tds_orig = np.fromfile(filename,sep='\n')
        if not len(tds_orig) == 200:
            print 'Y axis conversion probably failed'
            print '%d bins detected' %len(self.tds_orig)
        self.setAxisY(tds_orig[::-1])
        
    def _loadMatrix(self,filename):
        """Load intensity matrix from text file.
        :parameter filename: Absolute path to data file
        """
        temp = np.fromfile(filename,dtype=np.float64,sep=',')
        temp = np.array_split(temp,200)
        self.setMatrix(np.flipud(temp))
        
    def loadFolder(self,folderPath):
        """Load data from text files. Files should be named 'MassMobilityXaxis.txt'
        for the m/z axis, 'MassMobilityYaxis.txt' for arrival time axis and
        'MassMobility.txt' for the matrix of intensities.
        :parameter folderPath: Absolute path to directory containing text files.
        """
        self._loadX(os.path.join(folderPath,'MassMobilityXaxis.txt'))
        self._loadY(os.path.join(folderPath,'MassMobilityYaxis.txt'))
        self._loadMatrix(os.path.join(folderPath,'MassMobility.txt'))
        self.xaxisUnaltered = self.xaxis.copy()
        self.yaxisUnaltered = self.yaxis.copy()

    def loadRawFolder(self,rawPath,grain=2):
        """Load data from MassLynx raw data set.
        :parameter rawPath: Path to raw data directory
        :parameter grain: m/z spacing for extracting the data
        (minimum 0.5).
        """
        raw = RawFileProcessor.RawFileProcessor(rawPath)
        raw.processFolder(grain)
        self._getDataFromRawProcessor(raw,rawPath)

    def _getDataFromRawProcessor(self,raw):
        """Sub function of self.loadRawFolder(). Used for processing
        MassLynx raw files.
        :parameter raw: RawFileProcessor object
        """
        self.setAxisX(raw.getAxisX())
        self.setAxisY(raw.getAxisY())
        self.xaxisUnaltered = self.xaxis.copy()
        self.yaxisUnaltered = self.yaxis.copy()
        self.setMatrix(raw.getMassMobility())

    ###############################################################
    # Preparing species
    ###############################################################
    def setSpecies(self,species):
        """Add a species object to the self.species dictionary.
        :parameter species: msClasses.Species() object
        """
        self.species[species.name] = species
        
    def _generateSlice(self,xtremes):
        """Generate a data slice (sliced in td dimension), given
        the m/z limits.
        :parameter xtremes: m/z limits in the form [lower,upper]
        :returns: imClasses.DataSlice() object
        """
        dataSlice = DataSlice()
        dataSlice.setAxisY(self.yaxisUnaltered)
        dataSlice.setAxisX(self.xaxisUnaltered[xtremes[0]:xtremes[1]])
        dataSlice.setMatrix(self.matrixUnaltered[:,xtremes[0]:xtremes[1]])
        return dataSlice

    def generateSpeciesSlicesFwhm(self,speciesName,leftMultiplier=1,rightMultiplier=1):
        """Create data slices using peak width from mass spectrum fit
        for a specific species.
        :parameter speciesName: Name of molecular species
        :parameter leftMultiplier: Value to multiply the peak fwhm below the mean
        :parameter rightMultiplier: Value to multiply the peak fwhm above the mean
        """
        self.dataSlices[speciesName] = collections.OrderedDict()
        self.generateMassSpectrum()
        left = float(leftMultiplier*self.species[speciesName].peakFwhm)/2
        right = float(rightMultiplier*self.species[speciesName].peakFwhm)/2

        for z in self.species[speciesName].charges:
            mz = self.species[speciesName].getMz(z)
            mzIndex = utils.closest(mz,self.xaxis)
            mzIndex = localMaxima(mzIndex,self.massSpectrum.xvals,self.massSpectrum.yvals, scan_range=10)
            leftI = utils.closest(mz-left,self.xaxis)
            rightI = utils.closest(mz+right,self.xaxis)
            self.dataSlices[speciesName][z] = self._generateSlice([leftI,rightI])
            self.dataSlices[speciesName][z].generateAtd()
            self.dataSlices[speciesName][z].charge = z

    def generateSpeciesSlicesExplicit(self,limitDic):
        """Create data slices using explicit m/z limits for arrival
        time data extraction, for all supplied species and charge states.
        :parameter limitDic: d[species][z] = [lowerLimit,upperLimit]
        """
        '''
        limitDic generated by ChargeStatePeak.py
        Has actual absolute mz values for the limits
        Does all species at the same time unlike self.generateSpeciesSlicesFwhm()
        which does them all together
        limitDic structure is species>z>[lowerMz,upperMz]
        '''
        self.generateMassSpectrum()

        for speciesName,subD in limitDic.items():
            self.dataSlices[speciesName] = collections.OrderedDict()
            for z,limit in subD.items():
                leftI = utils.closest(limit[0],self.xaxis)
                rightI = utils.closest(limit[-1],self.xaxis)
                self.dataSlices[speciesName][z] = self._generateSlice([leftI,rightI])
                self.dataSlices[speciesName][z].generateAtd()
                self.dataSlices[speciesName][z].charge = z



    ###############################################################
    # Save and load Atropos fits
    ###############################################################
    # TODO(gns) - There should be a better way of doing this (don't need to pickle the whole thing, definitely don't need to do the data).
    # def saveMsFit(self,filename):
    #     """Save the current state of self.massSpectrum as a pickle.
    #     :parameter filename: Absolute path for mass spectrum fit file
    #     """
    #     import cPickle as pickle
    #     ofile = open(filename, 'wb')
    #     pickle.dump(self.massSpectrum,ofile)
    #     ofile.close()

    def loadMsFit(self,filename):
        """Load and set the Amphitrite mass spectrum fit from
        file. (Erases currently set mass spectrum data).
        :parameter filename: Absolute path to pickled
        msClasses.MassSpectrum() object ('.afit' file).
        """
        # TODO(gns) - this functions execution seems inconsistent with setMsFit()
        # TODO(gns) - If you change this, update documentation as well (bit in brackets)
        import cPickle as pickle
        ifile = open(filename, 'rb')
        self.massSpectrum = pickle.load(ifile)
        self.massSpectrum.xvals = np.array([])
        self.massSpectrum.yvals = np.array([])
        self.species = self.massSpectrum.species
        ifile.close()

    def setMsFit(self,msOb):
        """Set the Amphitrite mass spectrum fit, and as a by product,
        the species objects.
        :parameter msOb: msClasses.MassSpectrum() object
        """
        xvals = self.massSpectrum.xvals.copy()
        yvals = self.massSpectrum.yvals.copy()

        self.massSpectrum = msOb
        self.species = self.massSpectrum.simulatedSpecies
        self.simulatedSpecies = self.massSpectrum.simulatedSpecies

        self.massSpectrum.xvals = xvals
        self.massSpectrum.yvals = yvals

    ###############################################################
    # Cross section calculations
    ###############################################################
    def generateCcsAxesAndGrids(self,calibrationObj,ccsInterval=1):
        """Calibrate the ion mobility data by calculating CCS equivalents
        of the arrival time data.
        :parameter calibrationObj: imClasses.Calibration() object
        :parameter ccsInterval: Interval in Angstrom**2 to use for interpolation
        of the calibrated arrival time axes.
        """
        for name,spSlices in self.dataSlices.items():
            for charge,dataSlice in spSlices.items():
                dataSlice.generateCcsAxisAndGrid(calibrationObj,ccsInterval)

    def plotAtdSlicesHeatMap(self,ax):
        """Draw 3D contour plots showing each dataslice.
        :parameter ax: Matplotlib Axes instance
        """
        # TODO(gns) - should sort/reverse-sort the charges for reproducibility
        for name,spSlices in self.dataSlices.items():
            for charge,dataSlice in spSlices.items():
                matrix =  dataSlice.matrix
                x = [dataSlice.xaxis[0],dataSlice.xaxis[-1]]
                y = [dataSlice.yaxis[0],dataSlice.yaxis[-1]]
                ax.imshow(dataSlice.matrix,extent=[x[0],x[1],y[0],y[1]],aspect='auto',origin=[0,0])
        ax.autoscale()

    def getDataSlice(self,sp,z):
        """Get the data for a given species charge state.
        :parameter sp: Species name
        :parameter z: Charge state
        :returns: imClasses.DataSlice() object
        """
        # TODO(gns) - would be nice if this generated the ds if it didn't exist yet
        ds = self.dataSlices[sp][z]
        return ds


    ###############################################################
    # Plotting functions
    ###############################################################

    def plotChargeStateAtds(self,ax,speciesName,charges=0,lift=0,colourList=0,smoothing=0,window_len=3, smoothes=1, poly_order=1,**kwargs):
        """Draw the arrival time distribution of all (or given) charge states of
        a species.
        :parameter ax: matplotlib.axes.Axes() object
        :parameter speciesName: Name of species
        :parameter charges: List of charges to use (or 0/False to use all charges)
        :parameter colour: Matplotlib.pyplot compatible colour
        :parameter lift: Absolute value to separate ATDs vertically
        :parameter smoothing: Boolean for Savitzky Golay smoothing
        :parameter window_len: Window size for smoothing
        :parameter poly_order: Polynormial to use for smoothing
        :parameter \*\*kwargs: matplotlib.pyplot.plot() arguments
        """
        if not colourList:
            colourList = utils.colourList
        colourList = colourList * 50
        if not charges:
            charges = self.species[speciesName].charges
        if lift:
            maxHeights = []
            for z in charges:
                maxHeights.append(self.dataSlices[speciesName][z].atd.yvals.max())
            lift = max(maxHeights) * float(lift)/100

        for i,z in enumerate(sorted(charges)):
            if smoothing:
                self.dataSlices[speciesName][z].atd.smoothingSG(window_len=window_len, smoothes=smoothes, poly_order=poly_order)
            x = self.dataSlices[speciesName][z].atd.xvals
            y = self.dataSlices[speciesName][z].atd.yvals
            ax.plot(x,y+i*lift,color=colourList[i],label='+'+str(z),**kwargs)

        ax.legend(loc='upper right',prop=FontProperties(size=utils.legendFontSize))
        ax.set_xlabel('t$_d$')
        ax.set_ylabel('Intensity')
        ax.set_yticks([])

    def plotChargeStateAtd(self,ax,speciesName,charge,colour='gray',normalise='bpi',label='',lift=0,smoothing=0,window_len=3, smoothes=1, poly_order=1,**kwargs):
        """Draw the arrival time distribution of a particular charge state of
        a species.
        :parameter ax: matplotlib.axes.Axes() object
        :parameter speciesName: Name of species
        :parameter charge: Charge state
        :parameter colour: Matplotlib.pyplot compatible colour
        :parameter normalise: Normalisation method; base peak intensity ('bpi') or
        per area ('area').
        :parameter label: Label to be used in legend
        :parameter lift: Absolute value to lift the trace vertically (used for stacking
        multiple ATDs
        :parameter smoothing: Boolean for Savitzky Golay smoothing
        :parameter window_len: Window size for smoothing
        :parameter poly_order: Polynormial to use for smoothing
        :parameter \*\*kwargs: matplotlib.pyplot.plot() arguments
        """
        if smoothing: self.dataSlices[speciesName][charge].atd.smoothingSG(window_len=window_len, smoothes=smoothes, poly_order=poly_order)
        x = self.dataSlices[speciesName][charge].atd.xvals
        y = self.dataSlices[speciesName][charge].atd.yvals

        if normalise == 'bpi':
            y = y/y.max()*100
        if normalise == 'area':
            y = y/sum(y)*100*np.shape(y)[0]
        ax.plot(x,y+lift,color=colour,label=str(label),**kwargs)

        # return max point
        index = y.argmax()
        av,sum_of_weights = np.average(x,weights=y,returned=True)

        # calculate weighted stdev for return
        average,stdev = utils.weightedAverageAndStd(x,y)

        return x[index],av,sum_of_weights,stdev


    def plotMsExtractionLimits(self,ax,speciesName,colourList=0,**kwargs):
        """Fill in the area between the max and min limits of FWHM around the
        peak mean.
        :parameter ax: Matplotlib Axes instance
        :parameter speciesName: Species name
        :parameter colourList: List of matplotlib colours to use. If False default colours are used
        :parameter \*\*kwargs: matplotlib.pyplot.plot() arguments
        """
        # TODO(gns) - Replaced this with MassSpectrum.plotPeakWidthLimits()
        # TODO(gns) - Check before deleting as something might still use this
        if not colourList:
            colourList = utils.colourList
        if not 'alpha' in kwargs:
            kwargs['alpha'] = 0.3
        height = ax.get_ylim()[1]
        for i,z in enumerate(sorted(self.species[speciesName].charges)):
            left = self.dataSlices[speciesName][z].xaxis.min()
            right = self.dataSlices[speciesName][z].xaxis.max()
            ax.fill_betweenx([0,height], x1=left, x2=right,color=colourList[i], **kwargs)
        ax.set_ylim([0,height])

    ################
    def getDataSlices(self):
        """Get all available data slices.
        :returns: Dictionary of imClasses.DataSlice() objects in the form
        d[speciesName][chargeState] = DataSlice()
        """
        return self.dataSlices

    ####
    def plotCcsHeatMap(self,ax,calibrationOb,**kwargs):
        """Draw a 3D m/z against CCS contour plot.
        :parameter ax: Matplotlib Axes instance
        :parameter calibrationOb: imClasses.Calibration() object
        :parameter \*\*kwargs: Matplotlib imshow compatible arguments
        """
        dataSlices = self.getDataSlices()
        self.plotCcsHeatMapFromDataSlices(ax,calibrationOb,dataSlices,**kwargs)

    def plotCcsHeatMapFromDataSlices(self,ax,calibrationOb,dataSlices,**kwargs):
        """Draw a 3D m/z against CCS contour plot given dataSlices. (Used
        when you don't want to plot all dataSlices available).
        :parameter ax: Matplotlib Axes instance
        :parameter calibrationOb: imClasses.Calibration() object
        :parameter dataSlices: Dictionary of imClasses.DataSlice() objects
        (d[chargeState] = dataSlice)
        :parameter \*\*kwargs: Matplotlib imshow compatible arguments
        """
        '''subfunction of plotCcsHeatMap().
        Allows for the gui to manipulate the data (e.g. by scaling) b4 plotting'''
        for name,spSlices in dataSlices.items():
            for charge,dataSlice in spSlices.items():
                ccsMatrix,mzs,ccsAxis = dataSlice.getCcsMatrix(calibrationOb)
                ax.imshow(ccsMatrix,extent=[mzs[0],mzs[-1],ccsAxis[0],ccsAxis[-1]],aspect='auto',origin=[0,0],**kwargs)
        ax.autoscale()
    ####

    def plotChargeStateContourPlots(self,ax,calibrationOb,species,**kwargs):
        """Draw 3D Clemmer-style contour plots given species name.
        :parameter ax: Matplotlib Axes instance
        :parameter calibrationOb: imClasses.Calibration() object
        :parameter species: Species name
        :parameter \*\*kwargs: Matplotlib imshow compatible arguments
        """
        dataSlices = self.getDataSlices()
        self.plotChargeStateContourPlotsFromDataSlices(
            ax,calibrationOb,dataSlices,species,**kwargs)



    def plotChargeStateContourPlotsFromDataSlices(self,ax,calibrationOb,dataSlices,species,**kwargs):
        """Draw 3D Clemmer-style contour plots from provided dataSlices.
        :parameter ax: Matplotlib Axes instance
        :parameter calibrationOb: imClasses.Calibration() object
        :parameter dataSlices: Dictionary of imClasses.DataSlice() objects
        (d[chargeState] = dataSlice)
        :parameter species: The molecular species the dataSlices are from
        :parameter \*\*kwargs: Matplotlib imshow compatible arguments
        """
        spSlices = dataSlices[species]
        x = 0
        for charge,dataSlice in spSlices.items():
            ccsMatrix,mzs,ccsAxis = dataSlice.getCcsMatrix(calibrationOb)
            x += 1
            ax.imshow(ccsMatrix,extent=[x,x+2,ccsAxis[0],ccsAxis[-1]],
                      aspect='auto',origin=[0,0],**kwargs)
            x += 2
        ax.autoscale()
        ax.set_xlim([0,x+1])

        self.labelChargeStateContourPlots(ax,species)

    def labelChargeStateContourPlots(self,ax,species):
        """Annotate 3D Clemmer-style contour plots with charge
        state numbering.
        :parameter ax: Matplotlib Axes instance
        :parameter species: The molecular species to be labelled
        """
        '''Compatible with ATDs and CCSs'''
        ylims = ax.get_ylim()
        height = 0.03*(ylims[1]-ylims[0]) + ylims[0]
        spSlices = self.dataSlices[species]

        x = 0
        for charge,dataSlice in spSlices.items():
            x += 2
            ax.annotate("+%d"%charge,[x,height],horizontalalignment='center')
            x += 1


    def plotChargeStateContourPeaktops(self,ax,calibrationOb,species,smths,wlen,
                                       colour='b',limit=0,dataType='ccs',**kwargs):
        """For contour plots in Clemmer-style charge state mode;
        plots the peak tops as found using derivative.
        :parameter ax: Matplotlib Axes instance
        :parameter calibrationOb: imClasses.Calibration() object
        :parameter species: Name of molecular species to use
        :parameter smths: Rounds of smoothing to do before picking peaks
        :parameter wlen: Window length for smoothing (Savitzky Golay)
        :parameter colour: Colour for peak top markers
        :parameter dataType: Either CCS ('ccs') or arrival time ('std')
        :parameter \*\*kwargs: Matplotlib scatter compatible arguments
        """
        spSlices = self.dataSlices[species]
        x = 0
        for charge,dataSlice in spSlices.items():
            x += 1
            # Get the xvalue for where to draw the peaktop marker
            perc = dataSlice.getMsApexProportion()
            xval = perc*2. + x
            mz = dataSlice.getMsApex()

            # generate peaks
            if dataType == 'ccs':
                peaks = dataSlice.getCcsPeaks(smths,wlen,mz,calibrationOb,limit)
            elif dataType == 'atd':
                peaks = dataSlice.getAtqdPeaks(smths,wlen,limit)
            else:
                print 'invalid datatype in plotChargeStateContourPeaktops()'

            # plot the peaks
            # TODO(gns) - I haven't tested the kwargs option here yet
            ax.scatter([xval]*len(peaks),peaks,marker='o',color=colour,alpha=0.5, **kwargs)
            x += 2

        ax.set_xlim([0,x+1])


    #================================================================
    # CeRamp related functions

    def calculateAtdStatistics(self,speciesName,charge):
        """Calculate the weighted mean and standard deviation for an arrival
        time distribution.
        :parameter speciesName: Name of molecular species
        :parameter charge: Charge state to use
        :returns: average, standard deviation
        """
        average,stdev = self.calculateStatistics(self.dataSlices[speciesName][charge].atd)
        return average,stdev

    def calculateCcsdStatistics(self,speciesName,charge):
        """Calculate the weighted mean and standard deviation for a CCS
        distribution.
        :parameter speciesName: Name of molecular species
        :parameter charge: Charge state to use
        :returns: average, standard deviation
        """
        average,stdev = self.calculateStatistics(self.dataSlices[speciesName][charge].atd)
        return average,stdev


    def calculateStatistics(self,twoDdata):
        """Calculate the weighted mean and standard deviation
        :parameter twoDdata: An object which inherits from msClasses.TwoDdata()
        (e.g. MassSpectrum() or Atd()).
        :returns: average, standard deviation
        """
        # works on ccsds and atds
        average,stdev = twoDdata.calculateWeightedMeanStandardDeviation()
        return average,stdev



    def getAtd(self,speciesName,charge,normalised=True):
        """Get the arrival time distribution for the specific
        charge state of a species.
        :parameter speciesName: Name of molecular species
        :parameter charge: Charge state to use
        :parameter normalised: Boolean for whether to normalise data to
        base peak intensity
        :returns: imClasses.Atd() object
        """
        self.dataSlices[speciesName][charge].atd.normalisationBpi()
        return self.dataSlices[speciesName][charge].atd

    def getDataSlice(self,speciesName,charge):
        """Get the section of the data related to this charge state
        of this species.
        :parameter speciesName: Name of molecular species
        :parameter charge: Charge state to use
        :returns: imClasses.DataSlice() object
        """
        return self.dataSlices[speciesName][charge]
