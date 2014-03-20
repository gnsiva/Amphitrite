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
        '''Automatically checks the type of data in the folder
        and processes it to load the ImObj'''

        processedFilenames = ['MassMobilityXaxis','MassMobilityYaxis','MassMobility']
        textFilenames = [x+'.txt' for x in processedFilenames]
        amphiFilenames = [x+'.amphi' for x in processedFilenames]

        if os.path.isdir(folderName):
            allFound = True
            processed = False

            # amphi
            if self._checkForFiles(folderName, amphiFilenames):
                print 'Processing as pickles'
                self.loadFromPickles(folderName)
                processed = True

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

        # single amphi folder
        elif os.path.isfile(folderName):
            self.loadAmphiFile(folderName)
            processed = True
        else:
            print 'folder not found', folderName
            processed = False
        if not processed:
            print 'Not processed, disaster', folderName


    def loadAmphiFile(self,filename):
        dataList = utils.unPickleAmphitriteProject(filename)
        if dataList:
            self.setDataFromAmphiExtract(dataList)

    def setDataFromAmphiExtract(self,dataList):
        self.setAxisX(dataList[0])
        self.setAxisY(dataList[1])
        self.setMatrix(dataList[2])
        self.xaxisUnaltered = self.xaxis.copy()
        self.yaxisUnaltered = self.yaxis.copy()


    def _checkForFiles(self,folder,fileNames):
        allFound = True
        if type(fileNames).__name__ == 'int':
            fileNames = [fileNames]
        for f in fileNames:
            if not f in os.listdir(folder):
                allFound = False
        return allFound


    def _loadX(self,filename):
        xaxis = np.fromfile(filename,sep=',')
        self.setAxisX(xaxis[:-2])  # NO IDEA WHERE THE -2 CAME FROM
    def _loadY(self,filename):
        tds_orig = np.fromfile(filename,sep='\n')
        if not len(tds_orig) == 200:
            print 'Y axis conversion probably failed'
            print '%d bins detected' %len(self.tds_orig)
        self.setAxisY(tds_orig[::-1])
    def _loadMatrix(self,filename):
        temp = np.fromfile(filename,dtype=np.float64,sep=',')
        temp = np.array_split(temp,200)
        self.setMatrix(np.flipud(temp))
    def loadFolder(self,folderPath):
        self._loadX(os.path.join(folderPath,'MassMobilityXaxis.txt'))
        self._loadY(os.path.join(folderPath,'MassMobilityYaxis.txt'))
        self._loadMatrix(os.path.join(folderPath,'MassMobility.txt'))
        self.xaxisUnaltered = self.xaxis.copy()
        self.yaxisUnaltered = self.yaxis.copy()

    def loadRawFolder(self,rawPath,grain=2):

        raw = RawFileProcessor.RawFileProcessor(rawPath)
        raw.processFolder(grain)
        self._getDataFromRawProcessor(raw,rawPath)

    def _getDataFromRawProcessor(self,raw,folder):
        self.setAxisX(raw.getAxisX())
        self.setAxisY(raw.getAxisY())
        self.xaxisUnaltered = self.xaxis.copy()
        self.yaxisUnaltered = self.yaxis.copy()
        self.setMatrix(raw.getMassMobility())

    def loadFromPickles(self,folderName):
        raw = RawFileProcessor.RawFileProcessor(folderName)
        self._getDataFromRawProcessor(raw,folderName)


    ###############################################################
    # Preparing species
    ###############################################################
    def setSpecies(self,species):
        self.species[species.name] = species
    def _generateSlice(self,xtremes):
        dataSlice = DataSlice()
        dataSlice.setAxisY(self.yaxisUnaltered)
        dataSlice.setAxisX(self.xaxisUnaltered[xtremes[0]:xtremes[1]])
        dataSlice.setMatrix(self.matrixUnaltered[:,xtremes[0]:xtremes[1]])
        return dataSlice

    def generateSpeciesSlicesFwhm(self,speciesName,leftMultiplier=1,rightMultiplier=1):

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
    def saveMsFit(self,filename):
        import cPickle as pickle
        ofile = open(filename, 'wb')
        pickle.dump(self.massSpectrum,ofile)
        ofile.close()

    def loadMsFit(self,filename):
        import cPickle as pickle
        ifile = open(filename, 'rb')
        self.massSpectrum = pickle.load(ifile)
        self.massSpectrum.xvals = np.array([])
        self.massSpectrum.yvals = np.array([])
        self.species = self.massSpectrum.species
        ifile.close()

    def setMsFit(self,msOb):
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
        for name,spSlices in self.dataSlices.items():
            for charge,dataSlice in spSlices.items():
                dataSlice.generateCcsAxisAndGrid(calibrationObj,ccsInterval)
                print np.shape(dataSlice.matrix)



    def plotAtdSlicesHeatMap(self,ax):
        for name,spSlices in self.dataSlices.items():
            for charge,dataSlice in spSlices.items():
                matrix =  dataSlice.matrix
                x = [dataSlice.xaxis[0],dataSlice.xaxis[-1]]
                y = [dataSlice.yaxis[0],dataSlice.yaxis[-1]]
                ax.imshow(dataSlice.matrix,extent=[x[0],x[1],y[0],y[1]],aspect='auto',origin=[0,0])
        ax.autoscale()

    def getDataSlice(self,sp,z):
        # TODO(gns) - would be nice if this generated the ds if it didn't exist yet
        ds = self.dataSlices[sp][z]
        return ds


    ###############################################################
    # Plotting functions
    ###############################################################

    def plotChargeStateAtds(self,ax,speciesName,charges=0,lift=0,colourList=0,smoothing=0,window_len=3, smoothes=1, poly_order=1,**kwargs):
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
        '''Lift is an absolute value as the object can't access all the max heights'''

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
        return self.dataSlices

    ####
    def plotCcsHeatMap(self,ax,calibrationOb,**kwargs):
        dataSlices = self.getDataSlices()
        self.plotCcsHeatMapFromDataSlices(ax,calibrationOb,dataSlices,**kwargs)

    def plotCcsHeatMapFromDataSlices(self,ax,calibrationOb,dataSlices,**kwargs):
        '''subfunction of plotCcsHeatMap().
        Allows for the gui to manipulate the data (e.g. by scaling) b4 plotting'''
        for name,spSlices in dataSlices.items():
            for charge,dataSlice in spSlices.items():
                ccsMatrix,mzs,ccsAxis = dataSlice.getCcsMatrix(calibrationOb)
                ax.imshow(ccsMatrix,extent=[mzs[0],mzs[-1],ccsAxis[0],ccsAxis[-1]],aspect='auto',origin=[0,0],**kwargs)
        ax.autoscale()
    ####

    def plotChargeStateContourPlots(self,ax,calibrationOb,species,**kwargs):
        '''
        Clemmer plot as contour plots
        '''
        dataSlices = self.getDataSlices()
        self.plotChargeStateContourPlotsFromDataSlices(
            ax,calibrationOb,dataSlices,species,**kwargs)



    def plotChargeStateContourPlotsFromDataSlices(self,ax,calibrationOb,dataSlices,species,**kwargs):
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
        '''
        Clemmer plot peak tops
        only compatible with plotChargeStateContourPlots()

        dataType can be 'ccs' or 'atd'
        '''
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
            ax.scatter([xval]*len(peaks),peaks,marker='o',color=colour,alpha=0.5)
            x += 2

        ax.set_xlim([0,x+1])


    #================================================================
    # CeRamp related functions

    def calculateAtdStatistics(self,speciesName,charge):
        average,stdev = self.calculateStatistics(self.dataSlices[speciesName][charge].atd)
        return average,stdev

    def calculateCcsdStatistics(self,speciesName,charge):
        average,stdev = self.calculateStatistics(self.dataSlices[speciesName][charge].atd)
        return average,stdev


    def calculateStatistics(self,twoDdata):
        # works on ccsds and atds
        average,stdev = twoDdata.calculateWeightedMeanStandardDeviation()
        return average,stdev



    def getAtd(self,speciesName,charge,normalised=True):
        self.dataSlices[speciesName][charge].atd.normalisationBpi()
        return self.dataSlices[speciesName][charge].atd

    def getDataSlice(self,speciesName,charge):
        return self.dataSlices[speciesName][charge]
