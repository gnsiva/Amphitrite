"""Class for fitting and plotting 3 parameter Gaussian distributions."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com>"

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize


class Gaussian():

    def __init__(self):
        self.centre = 0.0
        self.fwhm = 0.0
        self.amplitude = 0.0

    def _gaussian(self,xvals,amplitude,centre,fwhm):
        """Generate a Gaussian given the input parameters.
        Returns the associated array of y values.
        xvals must be a numpy array (or a scalar value)
        """
        return amplitude*np.exp((-(xvals-centre)**2)/(2*(fwhm/2.3548200450309493)**2))

    def getData(self,xvals):
        """Generates the associated Gaussian y values for the supplied xvals
        (which should be a numpy array or scalar).
        self.amplitude, self.centre and self.fwhm should already have been set
        """
        yvals = self._gaussian(xvals,self.amplitude,self.centre,self.fwhm)
        return yvals

    def calculateAmplitude(self,xval):
        """Return Gaussian... think this is an unfinished function
        """
        # TODO(gns) - Check that this function isn't used anywhere and delete it
        # even if it is used somewhere else getData() does the same thing
        return self._gaussian(xval,self.amplitude,self.centre,self.fwhm)


    def setParameters(self,amplitude,centre,fwhm):
        """Set the values for a three parameter Gaussian"""
        self.amplitude = amplitude
        self.centre = centre
        self.fwhm = fwhm

    def estimateParameters(self,xvals,yvals,setValues=True):
        """Estimate values for the 3 parameter Gaussian
        The estimation is crude and is usually used to determine initial values for
        optimise parameters, which is much more robust.
        If setValues is True, use the estimated parameters as the values
        for self.amplitude, self.centre and self.fwhm.
        Otherwise values are returned as a list [amplitude,centre,fwhm]
        """
        fwhm = max(xvals) - min(xvals)
        amplitude = max(yvals)
        centre = np.average(xvals)
        if setValues:
            self.amplitude = amplitude
            self.centre = centre
            self.fwhm = fwhm
        else:
            return amplitude,centre,fwhm

    def optimiseParameters(self,xvals,yvals,setValues=1):
        """Use non linear least squares to fit the parameters of the
        Gaussian.
        setValues means the optimised parameters are
        set to this object, else the parameters are returned
        as a dictionary
        """
        fitfunc = lambda p,x: self._gaussian(x,p[0],p[1],p[2])
        errorfunc = lambda p,x,y: fitfunc(p,x)-y

        h,c,f = self.estimateParameters(xvals,yvals,setValues=0)
        p0 = [h,c,f]
        p1,success = optimize.leastsq(errorfunc,p0[:], args=(xvals,yvals))

        if not success:
            print 'Gaussian charge state distribution estimation failed'
            # TODO - dangerous quit statement
            quit()

        if setValues:
            self.amplitude = p1[0]
            self.centre = p1[1]
            self.fwhm = p1[2]
        else:
            d = {}
            d['amplitude'] = p1[0]
            d['centre'] = p1[1]
            d['fwhm'] = p1[2]
            return d

    def plot(self,xvals,ax,**kwargs):
        """Plot the Gaussian using the objects Gaussian parameters
        All matplotlib.pyplot.plot() arguments are allowed in **kwargs.
        """
        yvals = self._gaussian(xvals,self.amplitude,self.centre,self.fwhm)
        ax.plot(xvals,yvals,**kwargs)
