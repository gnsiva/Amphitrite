"""Class for deconvoluting and plotting native mass spectra."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com>"


import TwoDdata as tdd
import Species, Gaussian
from scipy import optimize
import time
import numpy as np
import collections
import lib.utils as utils
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


class MassSpectrum(tdd.TwoDdata):
    """Class for storing and deconvoluting mass spectra

    Inherits from TwoDdata:
    def __init__(self)
    def readFile(self,filename,x_range=0,grain=1)
    def normalisationBpi(self)
    def normalisationArea(self)
    def setNormalisationType(self,type)
    def smoothingSG(self,window_len=3,smoothes=2,poly_order=1)
    def _normalisePreset(self)
    def restoreRawYvals(self)
    def getAxesWithoutNans(self)
    def calculateWeightedMeanStandardDeviation(self)
    def calculateAreaUnderCurve(self)
    def _calculateGradient(self)
    def findPeaks(self,limit=0)
    def addPeak(self,mz)
    def plot(self,ax,**kwargs)
    def plotgPeaks(self,ax,labels=0,**kwargs)
    """
    def __init__(self):
        tdd.TwoDdata.__init__(self)
        self.species = collections.OrderedDict()
        self.simulatedSpecies = collections.OrderedDict()
        self.simulatedSpectrum = []
    
    def addSpecies(self,speciesObj,allowReplace=False):
        """Add a Species() object to the MassSpectrum.species dictionary.
        This allows it to be used in the fitting procedure.
        allowReplace means adding a second species with the same name to overwrite existing
        items in the dictionary.
        """
        name = speciesObj.name
        if not name in self.species.keys():
            self.species[name] = speciesObj
        else:
            if not allowReplace:
                print "Species with this name already exists.\nOperation failed"
            else:
                print "Species with this name exists, overwriting"
                self.species[name] = speciesObj
        
        
    def leastSquaresOptimisation(self,speciesNames,oneFwhm=False):
        """Deconvolute mass spectrum using non linear least squares.
        Setting oneFwhm to true results in the same value for peak full width
        half height being used for all species.
        """
        # setting up initial parameters
        zs,p0 = [],[]
        for name in speciesNames:
            zs.append(self.species[name].charges)
            p0 += self.species[name].getp0(oneFwhm)
        
        # using scipy optimize
        def errorfunc(p,xvals,yvals,zs,oneFwhm):
            return self.forLeastSquaresOptimisation(p,xvals,zs,oneFwhm) - yvals
        
        startTime = time.time()
        p1, success = optimize.leastsq(errorfunc,p0[:],args=(self.xvals,self.yvals,zs,oneFwhm))
        print "Optimisation took:", time.time()-startTime, "s"
        
        # storing the fit
        for i,p in enumerate(self._splitPs(p1)):
            tempSp = Species.Species(speciesNames[i])
            if not oneFwhm:
                peakFwhm = p[4]
            else:
                peakFwhm = oneFwhm
            tempSp.setValsFromp1([p[0],p[1],p[2],p[3],peakFwhm])
            tempSp.setCharges(zs[i])
            self.simulatedSpecies[speciesNames[i]] = tempSp
            if oneFwhm:
                self.simulatedSpecies[speciesNames[i]].peakFwhm = p1[4]
        self.simulatedSpectrum = self._simulateSpectrum(speciesNames, self.xvals)
    
    def _simulateSpectrum(self,speciesNames,xvals):
        """Builds a simulated spectrum using the fitted parameters of given species names.
        (uses the result of deconvolution)
        """
        for i,name in enumerate(speciesNames):
            if not i: 
                combined = self.simulatedSpecies[name].simulateSpecies(xvals)
            else: 
                combined += self.simulatedSpecies[name].simulateSpecies(xvals)
        return combined
    
    def errorCalculations(self,xvals):
        """Calculates and prints to console error statistics for comparing the simulated
        spectrum to the original data.
        + average error per data point
        + average error per m/z unit
        + RMSD (root mean square deviation)
        """
        dataBpi = self.yvals/self.yvals.max()*100
        simBpi = self.simulatedSpectrum/self.simulatedSpectrum.max()*100
        difference = np.abs(dataBpi-simBpi)
        
        avErrorPerDataPoint = np.sum(difference)/len(difference)
        print "Average error per datapoint: %.2f %%" %avErrorPerDataPoint
        
        mzGap = np.average(np.diff(xvals))
        avErrorPerMz = avErrorPerDataPoint/mzGap
        print "Average error per $m/z$ unit: %.2f %%" %avErrorPerMz
        
        rmsd = np.sqrt(np.sum(dataBpi-simBpi)**2/len(dataBpi))
        print "RMSD: %.2f" %rmsd
    
    def forLeastSquaresOptimisation(self,p0,xvals,zs,oneFwhm):
        """Creates the simulated mass spectrum for MassSpectrum.leastSquaresOptimisation().
        p0 - parameters to be used in simulation
        zs - charges to be simulated
        oneFwhm - if True use same FWHM for each species.
        """
        ps = self._splitPs(p0)
        combined = np.zeros(len(xvals))
        zGauss = Gaussian.Gaussian()
        
        for i,p in enumerate(ps):
            zGauss.setParameters(p[1], p[2], p[3])
            for j,z in enumerate(zs[i]):
                centre = utils.get_mz(p[0], z)
                amplitude = zGauss.calculateAmplitude(centre)
                if oneFwhm:
                    combined += utils.draw_peaks['hybrid'](xvals,amplitude,centre,p0[4])
                else:
                    combined += utils.draw_peaks['hybrid'](xvals,amplitude,centre,p[4])
        return combined
    
    def _splitPs(self,pLong):
        """Splits 1D list created by leastSquaresOptimisation into a multidimensional list.
        One dimension per species.
        """
        #5 params per species
        return [pLong[pos:pos + 5] for pos in xrange(0, len(pLong),5)]
    
            
    def getgPeaksFromIds(self,ids):
        """Generates a new gPeaks dictionary containing only the gPeaks
        with the supplied ids.
        gPeak format is gPeak[id] = [mz,intensity]
        """
        gPeaksOut = {}
        for id in ids:
            gPeaksOut[id] = self.gPeaks[id]
        return gPeaksOut

    # ================================================================
    # Plotting functions

    
    def plotSimulatedSpectrum(self,ax,xvals,speciesNames=0,**kwargs):
        """Plots just the simulated mass spectrum as generated by self.leastSquaresOptimisation(),
        with a default color of red.
        All matplotlib.pyplot.plot() arguments are allowed in **kwargs.
        """
        # TODO(gns) speciesNames is an uneeded argument, remove and bug test
        # TODO(gns)
        #ax = utils.checkAx(ax) # to add in later
        if not 'color' in kwargs:
            kwargs['color'] = 'r'
        ax.plot(xvals,self.simulatedSpectrum,**kwargs)   
    
    def plotTheoreticalMzs(self,ax,name,absolute=False,mass=0,**kwargs):
        """Used to draw the positions of calculated charge state peaks.
        Can be used with a Species() object stored in self.species by entering a 'name'
        or can be calculated directly from 'mass'
        If a value for absolute is given, it is used as the y axis position for charge
        state labels.
        All matplotlib.pyplot.plot() arguments are allowed in **kwargs.
        """
        self.species[name].plotTheoreticalMzs(ax,self.xvals,self.yvals,
                                              absolute,mass,**kwargs)

                    
    def plotAllComponents(self,ax,xvals,yvals,liftPercentage=10,speciesNames=0,colourList=0,traceLabels='mass',labelSize='small',**kwargs):
        """Plot mass spectrum, each of the deconvoluted elements and the simulated spectrum
        traceLabels can be 'mass' or 'name', with name the Species name is printed on the
        right tail of the trace.
        Species names allow you to change the order of the list
        All matplotlib.pyplot.plot() arguments are allowed in **kwargs.
        """
        # TODO - if you work out why you want to be able to supply the axes, write it down here!
        # colours
        if not colourList:
            colourList = utils.colourList 
        colourList = colourList * 50
        
        # heights and widths
        max = yvals.max()
        lift = liftPercentage/100. *max
        xleft = xvals.max() - (xvals[0].max()-xvals.min())*0.005 
        
        # plotting species
        if type(speciesNames).__name__ == 'int':
            speciesNames = self.simulatedSpecies.keys()
        j = 0
        for i,name in enumerate(speciesNames):
            tempys = self.simulatedSpecies[name].simulateSpecies(xvals)
            ax.plot(xvals,tempys+(j*lift+ lift/3.),color=colourList[j])
            if traceLabels == 'mass':
                ax.annotate('%.2f Da' %self.simulatedSpecies[name].mass,xy=(xleft,(j*lift+ lift/2.5)),size=labelSize,horizontalalignment='right')
            elif traceLabels == 'name':
                ax.annotate('%s' %name,xy=(xleft,(j*lift+ lift/2.5)),size=labelSize,horizontalalignment='right')
            j += 1
            
        # plotting combined
        ax.plot(xvals,self.simulatedSpectrum+(j*lift+ lift/3.),color='r',label='Simulated')
            
        j += 1
        
        # data
        ax.plot(xvals,yvals+(j*lift+ lift/3.),color='k',label='Experimental')
            
        # plotting parameters
        ax.set_yticks([])
        ax.legend(loc='upper right',prop=FontProperties(size=utils.legendFontSize))
        ax.set_xlabel('$m/z$')
        ax.set_ylabel('Intensity')  
        
        
    def writeCoordinatesToTxt(self,filename):
        """Write the base mass spectrum data to file
        Each line is in the format of mz\tintensity (tab separated).
        """
        file = open(filename,'w')
        for x,y in zip(self.xvals,self.yvals):
            print>>file, "%.3f\t%.3f" %(x,y)
        file.close()

