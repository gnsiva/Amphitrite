import wx,os,re
import cPickle as pickle
from collections import OrderedDict
from classes import MassSpectrum
from imClasses import Im
import gui.guiFunctions as gf
import IesCheckboxStates
import ListCtrlConformationsIes
from lib import utils
class IesSettings():

    def __init__(self,plotPanel):
        self.species = []
        self.speciesCharges = OrderedDict()
        
        self.atrOb = None
        self.calibrationOb = None
        
        self.plotPanel = plotPanel

        self.checkboxStates = IesCheckboxStates.IesCheckboxStates()
        self.listCtrlConformations = None

        self.loadedFiles = OrderedDict()
        self.massSpectra = OrderedDict() # filename as key
        self.filenames = []              # should be set by listctrlfiles

        self.widthL = 1.0
        self.widthR = 1.0

        self.speciesSelected = None
        self.chargeSelected = None

        self.units = ''
        self.values = []

        self.conformations = []

        self.displayPeaks = True

    #====================
    # Adding other components
    def setListCtrlConformations(self,listCtrl):
        self.listCtrlConformations = listCtrl

    def setPlotPanel(self,plotPanel):
        self.plotPanel = plotPanel

    #====================

    def setDisplayPeaks(self,bool):
        self.displayPeaks = bool

    def addConformation(self,ccs):
        self.listCtrlConformations.addConformation(ccs)
    def removeConformation(self,ccs):
        self.listCtrlConformations.removeConformation(ccs)
        
    def setSpeciesAndCharge(self,species,charge):
        self.speciesSelected = species
        self.chargeSelected = charge
    def getSpeciesAndCharge(self):
        return self.speciesSelected,self.chargeSelected
    
    def setUnits(self,units):
        self.units = units
        self.plotPanel.refresh_plot()
        
    def setValues(self,values):
        self.values = []
        for value in values:
            if utils.isNumber(value):
                self.values.append(float(value))
            else:
                self.values.append(None)
                                   

    def setCalibration(self,filename):
        import imClasses.Calibration
        try:
            self.calibrationOb = pickle.load(open(filename,'rb'))
        except:
            message = 'Something wrong with calibration file!'
            gf.warningDialog(message)
        
    def loadAtroposSpeciesAndCharges(self,path):
        self.atrOb = pickle.load(open(path,'rb'))

        # reset variables
        self.species = []
        self.speciesCharges = OrderedDict()
        
        for sp in sorted(self.atrOb.simulatedSpecies.keys()):
            self.species.append(sp)
            self.speciesCharges[sp] = self.atrOb.simulatedSpecies[sp].charges
        self.setImObsAtropos()
        
    def setImObsAtropos(self):
        if self.atrOb:
            for fn,imOb in self.loadedFiles.items():
                imOb.setMsFit(self.atrOb)
            
    def setFilenames(self,filenames):
        self.filenames = filenames
        
    def _forPlotMassSpectra(self,ax,lift):
        '''give the full path'''
        self._loadData(self.filenames)
        self._removeExcessLoadedData(self.filenames)
        self._plottingStackedSpectra(ax,lift)
                
    def _loadData(self,filenames):
        for fn in filenames:
            if not fn in self.loadedFiles.keys():
                imOb = Im.Im()
                imOb.loadFolderAuto(fn)              
                ms = imOb.getMassSpectrum()
                self.massSpectra[fn] = ms
                imOb.massSpectrum.normalisationBpi()
                self.loadedFiles[os.path.basename(fn)] = imOb
        self.setImObsAtropos()
                
    def _removeExcessLoadedData(self,filenames):
        fns = [os.path.basename(x) for x in filenames]
        toRemove = []
        for fn in self.loadedFiles.keys():
            if not fn in fns:
                toRemove.append(fn)
        for fn in toRemove:
            print 'Deleting: %s' %fn
            del self.loadedFiles[fn]
            
    def _plottingStackedSpectra(self,ax,lift):

        for i,msOb in enumerate(self.massSpectra.values()):
            msOb.yvals += lift*i
            msOb.plot(ax)
            msOb.yvals -= lift*i
            
    def _getMzLimits(self,species,charge):
        # limits = self.atrOb.simulatedSpecies[species].getPeakLimits(
        #          charge,self.widthL,self.widthR)
        # return limits
        mz = utils.get_mz(self.atrOb.simulatedSpecies[species].mass,charge)
        fwhm = self.atrOb.simulatedSpecies[species].peakFwhm
        left = (self.widthL*fwhm)/2.
        right = (self.widthR*fwhm)/2.
        return [mz-left,mz+right]

    def setWidthL(self,val):
        self.widthL = val
        self.plotPanel.refresh_plot()
    def setWidthR(self,val):
        self.widthR = val
        self.plotPanel.refresh_plot()

    def _getDataSlice(self,imOb):
        sp = self.speciesSelected
        z = self.chargeSelected
        return imOb.dataSlices[sp][z]
    
    def _updateExtractedSlices(self):
        for fn,imOb in self.loadedFiles.items():
            imOb.generateSpeciesSlicesFwhm(
                self.speciesSelected,self.widthL,self.widthR)
    def getCcsLines(self):
        lines = []
        ccsAxes = []
        self._updateExtractedSlices()
        for fn,imOb in self.loadedFiles.items():
            dataSlice = self._getDataSlice(imOb)
            dataSlice.generateCcsAxisAndGrid(self.calibrationOb,
                                             ccsInterval=1)
            ccsAxes.append(dataSlice.ccsAxis)
            line = dataSlice.getCcsDistribution(self.calibrationOb)
            lines.append(line/line.max()*100)
        return ccsAxes, lines

    def getAtdLines(self):
        lines = []
        xaxes = []
        self._updateExtractedSlices()
        for fn,imOb in self.loadedFiles.items():
            dataSlice = self._getDataSlice(imOb)
            xaxes.append(dataSlice.atd.xvals)
            lines.append(dataSlice.atd.yvals)
        return xaxes,lines
        

    def getCcsAxisAndGrid(self):
        matrices = []
        self._updateExtractedSlices()
        for i,(fn, imOb) in enumerate(self.loadedFiles.items()):
            ds = imOb.getDataSlice(self.speciesSelected,
                                    self.chargeSelected)
            ccsMatrix,xaxis,ccsAxis = ds.getCcsMatrix(self.calibrationOb)
            ccsMatrix = ccsMatrix.transpose()
            matrices.append(ccsMatrix)
        return ccsAxis,matrices
            
    def getAtdAxisAndGrid(self):
        matrices = []
        for i,(fn, imOb) in enumerate(self.loadedFiles.items()):
            ds = imOb.getDataSlice(self.speciesSelected,
                                    self.chargeSelected)
            matrix,xaxis,tdAxis = ds.getData()
            matrix = matrix.transpose()
            matrices.append(matrix)
        return tdAxis,matrices
