"""Class to hold the information for molecular species which can used with MassSpectrum()."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com>"


import numpy as np
import lib.utils as utils
import Gaussian
import collections
from scipy import integrate

class Species():
    '''Holds the data for a single species.
    Separate objects should be created for before and
    after (auto) simulation'''
    
    def __init__(self,name):
        self.name = name
        self.zGauss = Gaussian.Gaussian()
        self.charges = []
        self.mass = 0.0
        self.calcMassError = 0.0
        self.peakFwhm = 0.0


    def plotTheoreticalMzs(self,ax,xvals,yvals,absolute=0,mass=0,**kwargs):
        """Used to draw the positions of calculated charge state peaks.
        If a value for absolute is given, it is used as the y axis position for charge
        state labels.
        All matplotlib.pyplot.plot() arguments are allowed in **kwargs.
        """
        if not mass:
            mass = self.mass
        oneDone = False
        for i,z in enumerate(np.arange(100)+1):
            mz = utils.get_mz(mass, z)
            if xvals[0] < mz < xvals[-1]:
                # the line
                if not oneDone: ax.axvline(mz, label="%s - %.1f Da +/- %.1f" %(self.name,self.mass,self.calcMassError), **kwargs)
                else: ax.axvline(mz, **kwargs)
                oneDone = True
                # annotation
                if z/2 == float(z)/2: lift = .1 * yvals.max()
                else: lift = .2 * yvals.max()    
                if not 'color' in kwargs:
                    color = 'blue'
                else:
                    color = kwargs['color']      
                if not absolute:
                    ax.text(mz,lift,'+%d' %z, color=color)
                else:
                    ax.text(mz,absolute,'+%d' %z, color=color)
    
    def getp0(self,oneFwhm=0):
        """Return the list of initial parameters for leastSquaresOptimisation()
        (in MassSpectrum()) from this Species() object.
        """
        if oneFwhm:    
            return [self.mass,self.zGauss.amplitude,self.zGauss.centre,self.zGauss.fwhm,oneFwhm]
        else:
            return [self.mass,self.zGauss.amplitude,self.zGauss.centre,self.zGauss.fwhm,self.peakFwhm]

   

    def setValsFromp1(self,p1):
        """After leastSquaresOptimisation() set the fitted parameters (p1)
        to this Species() object.
        """
        self.mass = p1[0]
        self.zGauss.amplitude = p1[1]
        self.zGauss.centre = p1[2]
        self.zGauss.fwhm = p1[3]
        self.peakFwhm = p1[4]     


    def setVals(self,mass,z_mu,z_amp,z_fwhh,p_fwhh=0):
        """Set parameters for this object.
        """
        # TODO(gns) - This function doesn't seem necessary, try to remove
        # Not even sure if it would work, setValsFromp1 should break if p1 only has 4 elements
        if p_fwhh:
            self.setValsFromp1(mass,z_mu,z_amp,z_fwhh,p_fwhh)
        else:
            self.setValsFromp1(mass,z_mu,z_amp,z_fwhh)
    
    def simulateSpecies(self,xvals,peakShape='hybrid'):
        """Simulates a mass spectrum using the object's attributes
        Valid peak shapes are: 'hybrid', 'gaussian' & 'lorentzian
        if one_fwhh is to be used be sure to set self.peakFwhm before
        calling this function.
        """
        combined = np.zeros(len(xvals),dtype='float')  
        for z in self.charges:
            centre = utils.get_mz(self.mass, z)
            amplitude = self.zGauss.calculateAmplitude(centre)
            combined += utils.draw_peaks[peakShape](xvals,amplitude,centre,self.peakFwhm)

        return combined
    
    def getMz(self,z):
        """Calculate the m/z value for a given charge state using the mass 
        from self.mass
        """
        return utils.get_mz(self.mass, z)
    
    def calculateMassAndCharges(self,mzs):
        """Calculate the mass of a molecular species using the given m/z
        values.
        Primarily used as a subfunction of self.calculateMass().
        """
        iarray = mzs[:]   
        iarray.sort()
        iarray.reverse()
        charges = collections.OrderedDict()      
        lowest = 10000000
        lowest_z = 0
        zs = xrange(1,101)
        for z in zs:
            charges[z] = []
            for i,mz in enumerate(iarray):
                charges[z].append(utils.get_mass(mz,z+i))
        
        for z in charges.keys():
            sd = np.std(charges[z])
            if sd < lowest:
                lowest = sd
                lowest_z = z
        
        # calculating error
        total_error = []
        for mass in charges[lowest_z]:
            total_error.append(abs(np.average(charges[lowest_z])-mass))
        average_error = np.average(total_error)
        
        return np.average(charges[lowest_z]), average_error, [lowest_z+i for i in xrange(len(mzs))]
        
    def calculateMass(self,mzs):
        """Calculate the mass of a molecular species using the given m/z
        values.
        Returns mass and error, where error is the average absolute
        difference between the given m/z values and the m/z from the
        calculated mass.
        """
        mass,error,charges = self.calculateMassAndCharges(mzs)
        return mass,error


    def setSpecies(self,gPeaks,peakFwhm=10):
        """Setup this Species() object using gPeaks (d[id]=[mz,intensity])
        Calculates the mass and estimates parameters for the other
        Species() parameters for later use as initial values for
        leastSquaresOptimisation() (MassSpectrum()).
        """
        xvals,yvals = [],[]
        for id in gPeaks:
            xvals.append(gPeaks[id][0])
            yvals.append(gPeaks[id][1])
        self.mass, self.calcMassError = self.calculateMass(xvals)
        
        if len(xvals) > 2:
            self.zGauss.optimiseParameters(xvals, yvals, setValues=1)
        else:
            self.zGauss.estimateParameters(xvals, yvals, setValues=1)
        if self.charges == []: self.charges = self._estimateCharges()
        self.peakFwhm = peakFwhm 
        
        
    def setSpeciesGivenMass(self,mass,xvals,yvals,peakFwhm,zs):
        """Similar to Species.setSpecies(). Except that instead of using
        gPeaks, the Species() object is setup automatically when given
        values for mass and the charges to be analysed.
        """
        self.mass = mass
        pseudogPeaks = {}

        for z in zs:
            xval = utils.get_mz(self.mass, z)
            i = utils.closest(xval, xvals)
            yval = yvals[i]
            pseudogPeaks[z] = [xval,yval]
        self.setSpecies(pseudogPeaks, peakFwhm)
            
    
    def setCharges(self,charges):
        """For this molecular species, set the charges which
        are detected in the mass spectrometer.
        """
        self.charges = charges
    
    def _estimateCharges(self,limit=1):
        """Estimate charges to be simulated by fitting the charge state
        Gaussian distribution.
        Used as a subfunction for self.setSpecies()
        Limit is given as a percentage of the total height of the Gaussian
        and is used as the cutoff point for whether a charge state is to
        be included or discarded.
        """
        # TODO(gns) - perhaps lower the limit for atropos at least
        # probably the default value as well.
        zs = np.arange(1,151)
        charges = []
        for z in zs:
            xval = utils.get_mz(self.mass, z)
            height = self.zGauss.calculateAmplitude(xval)
            if height > self.zGauss.amplitude*(float(limit)/100):
                charges.append(z)
        print charges
        return charges
    
    def getTotalArea(self,xvals):
        """Calculate the area under the deconvoluted mass spectrum for just this
        species using the trapezium method of integration.
        """
        return integrate.trapz(self.simulateSpecies(xvals), xvals)#, dx, axis)
    
    def getTotalIntensity(self,xvals,yvals):
        """Using the supplied data, sum the intensity of the highest
        point of the data at the m/z value for each charge state in
        self.charges using self.mass as the mass.
        """
        '''Gets experimental intensity'''
        intensity = 0
        for charge in self.charges:
            mz = utils.get_mz(self.mass, charge)
            i = utils.closest(mz, xvals)
            intensity += yvals[i]
        return intensity
        
        
    
    def getPeakLimits(self,z,leftMultiplier=1.0,rightMultiplier=1.0):
        """Returns lower and upper m/z limits for charge state peaks
        using self.peakFwhm as the original value, and using multipliers
        to alter the left and right limits independently.
        """
        mz = self.getMz(z)
        left = mz - (self.peakFwhm/2)*leftMultiplier
        right = mz + (self.peakFwhm/2)*rightMultiplier
        return [left,right]
