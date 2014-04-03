__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com>"

"""
This file isn't imported by any of the GUI apps, probably not needed in release
"""

import lib.utils as utils
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import numpy as np

class SpecialFigures():
    def __init__(self):
        pass
    
    
    def ccsDistributions(self,ax,imObj,calibrationObj,spName=0,
                         ccsInterval=1,colourList=0,lift=0,**kwargs):
        """Extract arrival time distributions for each charge state of a species,
        convert to CCS and plot them vertically stacked.

        :parameter ax: matplotlib Axes instance
        :parameter imObj: imClasses.Im() object
        :parameter calibrationObj: imClasses.Calibration() object (IM calibration)
        :parameter spName: Name of molecular species to use from imObj
        :parameter ccsInterval: Spacing of CCS values in Angstrom**2
        :parameter colourList: List of matplotlib compatible colours to use
        :parameter \*\*kwargs: matplotlib.pyplot.plot() arguments
        """
        '''use spName=0 if there is only one species in imObj'''
        if not spName:
            spName = imObj.dataSlices.keys()[0]
        if not colourist:
            colourList = utils.colourList

        ccsLines = []
        maxIntensity = 0
        
        spCharge = imObj.dataSlices[spName]
        for i,(charge,dataSlice) in enumerate(spCharge.items()):
            dataSlice.generateCcsAxisAndGrid(calibrationObj,ccsInterval)
            ccsLines.append(dataSlice.getCcsDistribution(calibrationObj))
            if ccsLines[i].max() > maxIntensity:
                maxIntensity = ccsLines[i].max()
        
        for i,(charge,dataSlice) in enumerate(spCharge.items()):
            label = '+%d' %charge
            ccsLine = ccsLines[i] + i*(maxIntensity * (float(lift)/100))
            ax.plot(dataSlice.ccsAxis,ccsLine,label=label,color=colourList[i],**kwargs)
        
        ax.legend(loc='upper left',prop=FontProperties(size=utils.legendFontSize))
        plt.xlabel('$\Omega$ ($\AA^2$)')

        
    def calibrationFigure(self,fig,calibrationObj,spToDisplay=0,autoAxesLimits=1,colourList=0,**kwargs):
        # TODO(gns) - This can probably be removed only used in:
        # 121207_test_SpecialFigureCalibration.py, 121209_test_SpecialFigureCalibration2.py
        if not spToDisplay:
            spToDisplay = calibrationObj.calibrants.keys()[0]
        if not colourList:
            colourList = utils.colourList
        
        sps = calibrationObj.calibrants.keys()          
        bandColour = colourList[sps.index(spToDisplay)]
        
        calibrant = calibrationObj.calibrants[spToDisplay]
        ax = fig.add_subplot(311) # mass spectrum
        calibrant.plotMsAndExtractionLimits(ax,bandColour,**kwargs)
        
        if autoAxesLimits:
            mzs = [utils.get_mz(calibrant.approxMass, z) for z in calibrant.charges]
            plt.xlim([min(mzs)*0.5,max(mzs)*2])
        
        
        ax = fig.add_subplot(312) # atds
        calibrant.plotChargeStateAtds(ax,colourList=colourList)
        calibrant.plotCalibrantTdPeaks(ax,colourList=colourList)
        if autoAxesLimits:
            tdPrimes = [calibrant.tdsDoublePrime[z] for z in calibrant.charges]
            plt.xlim([min(tdPrimes)*0.5,max(tdPrimes)*2.0])
        
        ax = fig.add_subplot(313) # fit
        calibrationObj.plotCalibrationCurve(ax,colourList,**kwargs)
        
        
    def atDistributions(self,ax,imObj,spName=0,colourList=0,lift=0,**kwargs):
        # TODO(gns) - Delete this function, only used in an obscure test file
        # Amphitrite_2.1/test_files/121207_test_CalibrationOnDataGrid2.0.py
        if not spName:
            spName = imObj.dataSlices.keys()[0]
        if not colourList:
            colourList = utils.colourList
        
        atdIntensities = []
        maxIntensity = 0
        
        spCharge = imObj.dataSlices[spName]
        for i,(charge,dataSlice) in enumerate(spCharge.items()):
            intensities = dataSlice.atd.yvals
            atdIntensities.append(intensities)
            if atdIntensities[i].max() > maxIntensity:
                maxIntensity = atdIntensities[i].max()
        
        for i,(charge,dataSlice) in enumerate(spCharge.items()):
            label = '+%d' %charge
            intensities = atdIntensities[i] + i*(maxIntensity * (float(lift)/100))
            ax.plot(dataSlice.atd.xvals,intensities,label=label,color=colourList[i],**kwargs)
        
        ax.legend(loc='upper left',prop=FontProperties(size=utils.legendFontSize))   
        
