"""Class for ion mobility protein calibrant data and analysis."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
# from imClasses.Im import Im
# from msClasses.MassSpectrum import MassSpectrum
# from msClasses.Species import Species

from imClasses import Im
from msClasses import MassSpectrum, Species

im = Im()

import lib.utils as utils
import collections
import os

class Calibrant():

    def __init__(self,proteinName,folderName,peakFwhm=10):
        """
        :parameter proteinName: Calibrant name (see self._setProteinName() for options)
        :parameter folderName: Absolute path to amphitrite data file (or -1 to not load data)
        :parameter peakFwhm: m/z peak width to use when extracting arrival time data
        """
        # TODO(gns) - folderName should be changed to amphiPath or something
        self.charges = []
        self.approxMass = collections.OrderedDict()
        self.publishedCcss = collections.OrderedDict()
        self.tds = collections.OrderedDict()
        self.tdsDoublePrime = collections.OrderedDict()
        self.CcssPrime = collections.OrderedDict()
        self.leftMultiplier = 1
        self.rightMultiplier = 1
        self.peakFwhm = float(peakFwhm)
        self.imObj = None
        self.speciesObj = None

        # load protein information
        self._setProteinName(proteinName)

        if folderName != -1:
            self._loadDataAndCreateImObj(folderName)
            self._generateSpeciesObject(peakFwhm)

            # generate the matrix slices
            self._updateSpecies(self.speciesObj, self.peakFwhm, self.leftMultiplier, self.rightMultiplier)

    ###########################################################################
    # Functions for init
    ###########################################################################
    def _loadDataAndCreateImObj(self,folderName):
        """Set the amphitrite data file path, load the data and create an imClasses.Im()
        object.
        :parameter folderName: Absolute path to amphitrite data file
        """
        # TODO(gns) - rename folderName (keep consistent though, this occurs in other functions as well)
        self.folderName = folderName
        self.imObj = Im()
        self.imObj.loadFolderAuto(folderName)

    def _generateSpeciesObject(self,peakFwhm):
        """Create a species object using preset charges and masses for the given
        calibrant type.
        :parameter peakFwhm: m/z width to use for extracting arrival time data
        """
        self.speciesObj = Species(self.name)
        self.speciesObj.mass = self.approxMass
        self.speciesObj.peakFwhm = peakFwhm

    ###########################################################################
    # Calibrant set up functions
    ###########################################################################
    def setPeakFwhm(self,fwhm):
        """Set the m/z width to be used for extracting arrival time
        information and update the extraction.
        :parameter fwhm: Peak full width half maximum (m/z)
        """
        self._updateSpecies(self.speciesObj, peakFwhm=fwhm,updateApex=0)
        
    def setPeakFwhmMultipliers(self,leftMultiplier,rightMultiplier):
        """Set values to multiply the m/z value peak width for the lower
        and upper boundary for arrival time extraction. (Then update td data)
        :parameter leftMultiplier: Multiply by the peak width lower than the mean
        :parameter rightMultiplier: Multiply by the peak width above the mean
        """
        self._updateSpecies(self.speciesObj, self.peakFwhm, leftMultiplier, rightMultiplier,updateApex=0)
        
    def setCharges(self,charges):
        """Set the charges for the species. (Used for ATD extraction).
        :parameter charges: List of charge states (int)
        """
        self.speciesObj.charges = charges
        self.charges = charges
        self._updateSpecies(self.speciesObj,updateApex=0)

    def setTdValue(self,td,z):
        """Set the arrival time value for a particular charge state.
        :parameter td: Arrival time value for peak (usually apex)
        :parameter z: Charge state
        """
        self.tds[z] = td
        self._updateSpecies(self.speciesObj, self.peakFwhm, self.leftMultiplier, self.rightMultiplier, 0)

    def getChargeTdsDict(self):
        """Get a dictionary for the arrival time values for each charge
        state specified in self.charges.
        :returns: Dictionary of form d[ChargeState] = arrival time
        """
        # for the GUI (CalibrantGuiGrids)
        charges = sorted(self.charges)
        d = {z: self.tds[z] for z in charges}
        import operator
        d = collections.OrderedDict(sorted(d.items(), key=operator.itemgetter(0)))
        return d

    #===========================================================================
    # plotting functions
    #===========================================================================
    def plotMsAndExtractionLimits(self,ax,colourList=0,**kwargs):
        """Plot the calibrant's mass spectrum on the supplied axes and
        shade in the areas which is to be used for extraction of arrival
        time data.
        :parameter ax: Matplotlib Axes instances
        :parameter colourList: List of colours to use or 0/False to use defaults
        :parameter \*\*kwargs: Matplotlib.pyplot compatible commands for mass spectrum
        """
        self.imObj.massSpectrum.plot(ax, **kwargs)
        self.imObj.plotMsExtractionLimits(ax, speciesName=self.name, colourList=colourList)
        
    def plotChargeStateAtds(self,ax,lift=0,colourList=0,charges=0,smoothing=1,**kwargs):
        """Plot the arrival time distributions for the charge states of
        the calibrant.
        :parameter ax: Matplotlib Axes instances
        :parameter lift: Vertical spacing between traces
        :parameter colourList: List of colours to use or 0/False to use defaults
        :parameter charges: List of charges to display (or 0 to display all)
        :parameter smoothing: Boolean; whether or not to smooth the ATDs
        :parameter \*\*kwargs: Matplotlib.pyplot compatible commands for mass spectrum
        """
        if not charges:
            charges = self.charges
        self.imObj.plotChargeStateAtds(ax, self.name, lift=lift, charges=charges,colourList=colourList, smoothing=smoothing)
        ax.legend(loc='upper right',prop=FontProperties(size=utils.legendFontSize))

    def plotCalibrantTdPeaks(self,ax,colourList=0):
        """Plot the arrival time values which are to be used in calibration
        (usually using the ATD apex as the time value).
        :parameter ax: Matplotlib Axes instance
        :parameter colourList: List of colours to use or 0/False to use defaults
        """
        for i,z in enumerate(sorted(self.charges)):
            ax.axvline(self.tds[z],color='k',lw=0.5)
            if type(colourList).__name__ == 'int':
                ax.axvline(self.tds[z],color=utils.colourList[i],lw=3,alpha=0.3)
            else:
                ax.axvline(self.tds[z],color=colourList[i],lw=3,alpha=0.3)

    ###########################################################################
    # Generating Corrected Td and CCS values (called from Calibration() class)
    ###########################################################################
    def generateCorrectedTdsAndCcss(self,waveVelocity,gas='Nitrogen'):
        """Calculate td prime and CCS prime, using the calibration equations
        including the reduced mass equation.
        :parameter waveVelocity: Velocity in m/s
        :parameter gas: Mobililty gas used (default: 'Nitrogen')
        """
        for i,charge in enumerate(self.charges):
            mz = utils.get_mz(self.approxMass, charge)
            # td
            tdPrime = utils._calculateTdPrime(self.tds[charge], waveVelocity)
            self.tdsDoublePrime[charge] = utils._calculateTdDoublePrime(tdPrime, mz)
            #CCS
            reducedMass = utils._calculateReducedMass(mz, charge, gas)
            self.CcssPrime[charge] = utils._calculateOmegaPrime(self.publishedCcss[charge], charge, reducedMass)
            
    def getTdsDoublePrime(self):
        """Calculate td''.
        :returns: List of td'' values
        """
        return [self.tdsDoublePrime[z] for z in sorted(self.charges)]
    
    def getCcssPrime(self):
        """Calculate CCS''.
        :returns: List of CCS'' values
        """
        return [self.CcssPrime[z] for z in sorted(self.charges)]


    ###########################################################################
    # Private functions
    ###########################################################################
    def _updateSpecies(self,speciesObj=0,peakFwhm=0,leftMultiplier=0,rightMultiplier=0,updateApex=1):
        """Update the settings of the species object.
        :parameter speciesObj: msClasses.Species() object (or 0 to use self.speciesObj)
        :parameter peakFwhm: m/z width to use when extracting arrival time data
        :parameter leftMultiplier: Multiply by the peak width lower than the mean
        :parameter rightMultiplier: Multiply by the peak width above the mean
        :parameter updateApex: If True, recalculate the td apex
        """
        if peakFwhm:
            self.peakFwhm = peakFwhm
        if leftMultiplier:
            self.leftMultiplier = leftMultiplier
        if rightMultiplier:
            self.rightMultiplier = rightMultiplier
        if not speciesObj:
            speciesObj = self.speciesObj

        # check charges are within the mass spectrum limits
        zs = []
        left = self.peakFwhm * self.leftMultiplier
        right = self.peakFwhm * self.rightMultiplier
        for z in self.charges:
            mz = utils.get_mz(self.approxMass,z)
            if self.imObj.xaxis.min()+left < mz < self.imObj.xaxis.max()-right:
                zs.append(z)
            else:
                print 'charge outside m/z range'
        self.charges = zs
        self.speciesObj.charges = self.charges

        self.imObj.setSpecies(self.speciesObj)

        self.imObj.generateSpeciesSlicesFwhm(self.name, self.leftMultiplier, self.rightMultiplier)

        if updateApex:
            self.tds = collections.OrderedDict()

            for z in self.charges:
                self.imObj.dataSlices[self.name][z].atd.smoothingSG(window_len=3, smoothes=1, poly_order=1)
                td = self.imObj.dataSlices[self.name][z].getAtdApex()
                self.tds[z] = td

    def _setProteinName(self,name):
        """Set the published CCS values, their associated charge states and
        an approximate mass of the calibrant given the protein name.
        :parameter name: Calibrant name, choices: myoglobin, cytochrome c
        native,cytochrome c denatured, avidin, bsa, adh, sap5, sap10,
        concanavalin a, pyruvate kinase, blac1, blac2 and bk.
        """
        name = name.lower()
        self.name = name
        self._privatePublishedCcs = collections.OrderedDict()
        # TODO(gns) - check to make sure the myoglobin values are the helium ones
        self._privatePublishedCcs['myoglobin'] = {11:2942.,12:3044.,13:3136.,14:3143.,15:3230.,16:3313.,17:3384.,18:3489.}

        # Native proteins Bush2010 - analysed using helium as mobility gas
        self._privatePublishedCcs['cytochrome c native'] = {6:1240.,7:1280.}
        self._privatePublishedCcs['blac1'] = {7:1660.,8:1690.,9:1780.}
        self._privatePublishedCcs['blac2'] = {11:2850,12:2900.,13:2960}
        self._privatePublishedCcs['transthyretin'] = {14:3410.,15:3400.,16:3380}
        self._privatePublishedCcs['avidin'] = {15:3640.,16:3640.,17:3640.,18:3640.}
        self._privatePublishedCcs['bsa'] = {14:4090.,15:4100.,16:4060.,17:4040.}
        self._privatePublishedCcs['concanavalin a'] = {20:5550.,21:5550.,22:5480.,23:5450.}
        self._privatePublishedCcs['adh'] = {23:6940.,24:6940.,25:6830.,26:6720.}
        self._privatePublishedCcs['sap5'] = {22:7030.,23:6970.,24:6930.,25:6860.,26:6830.}
        self._privatePublishedCcs['sap10'] = {31:1040.,32:1050.,33:1060.,34:1050.,35:1070.}
        self._privatePublishedCcs['pyruvate kinase'] = {30:10300.,31:10300.,32:10300.,33:10200.,34:10000.}
        self._privatePublishedCcs['glutamate dehydrogenase'] = {39:12800.,40:12800.,41:12800.,42:12800.,43:12800.}
        self._privatePublishedCcs['groel'] = {67:20900.,68:20900.,69:20700.,70:20700.,71:20600.,72:20700.}

        # Denatured proteins Bush2010 - helium
        self._privatePublishedCcs['bk'] = {2:237.}



        # Denatured proteins - Clemmer:
        # http://www.indiana.edu/~clemmer/Research/cross%20section%20database/Proteins/protein_cs.htm
        self._privatePublishedCcs['cytochrome c denatured'] = {3:1139,
                                                       4:1153,
                                                       5:1196,5:1340,
                                                       6:1393,6:1244,6:1602,
                                                       7:1785,7:1247,7:2007,7:1620,
                                                       8:1845,8:1250,8:2061,8:1702,
                                                       9:2215,9:1964,
                                                       10:2226,
                                                       11:2303,
                                                       12:2335,
                                                       13:2391,
                                                       14:2473,
                                                       15:2579,
                                                       16:2679,
                                                       17:2723,
                                                       18:2766,
                                                       19:2800,
                                                       20:2889}

        # Mass Estimates for fitting (taken from experimental data):
        self.cal_mass = collections.OrderedDict()
        self.cal_mass['cytochrome c native'] = 12355.06
        self.cal_mass['cytochrome c denatured'] = 12323.528 #12355.06-600
        self.cal_mass['blac1'] = 18400.
        self.cal_mass['blac2'] = 36800.

        self.cal_mass['avidin'] = 64133.17
        self.cal_mass['bsa'] = 66776.04
        self.cal_mass['adh'] = 148138.02
        self.cal_mass['concanavalin a'] = 102913.46
        self.cal_mass['myoglobin'] = 16952.28
        self.cal_mass['sap5'] = 125000.
        self.cal_mass['sap10'] = self.cal_mass['sap5']*2

        # Need to update this: !!!!!!!
        self.cal_mass['pyruvate kinase'] = 227925.

        self.cal_mass['bk'] = 1060.

        ################################################
        #Function code:
        ################################################
        self.publishedCcss = self._privatePublishedCcs[name]
        self.charges = self._privatePublishedCcs[name].keys()
        self.approxMass = self.cal_mass[name]
