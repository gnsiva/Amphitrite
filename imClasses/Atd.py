"""Class for arrival time distributions (ATD) - plotting, peak identification and calculations."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com>"

import msClasses.TwoDdata as tdd
from lib import utils

class Atd(tdd.TwoDdata):
    
    def __init__(self):
        tdd.TwoDdata.__init__(self)
    
    def setAtdYvals(self,yvals):
        """Set the y axis values for the object.
        
        :parameter yvals: Intensity values for the Atd() object
        :type numpy array:
        :returns: None
        """
        self.rawyvals = yvals.copy()
        self.yvals = self.rawyvals.copy()
    
    def setAtdXvals(self,xvals):
        """Set the x axis values for the object.
        
        :parameter yvals: Arrival time values for the Atd() object
        :type numpy array:
        :returns: None
        """
        # TODO(gns) - there should be a function for this and the above in TwoDdata()
        self.xvals =xvals.copy() 

    def plotCcs(self,ax,ccsAxis,**kwargs):
        """Plot a CCS distribution, handles figure labelling.
        
        :parameter ax: matplotlib axes instance or False
        :parameter ccsAxis: CCS axis which corresponds to the arrival time axis
        :type numpy array:
        :parameter \*\*kwargs: matplotlib.pyplot.plot() arguments
        :returns: None
        """
        if not 'color' in kwargs:
            kwargs['color'] = 'black'
        if not 'lw' in kwargs:
            kwargs['lw'] = 0.8
        ax.plot(ccsAxis,self.yvals,**kwargs)
        ax.set_ylabel('Intensity')
        ax.set_xlabel('$\Omega$ ($\AA^2$)')
        
    
    def autoAxisX(self,ax,limit=0.1,percExtra=10,ccsAxisValues=False):
        """Automatically crop the x axis of the axes instance.
        
        :parameter ax: matplotlib.axes() object
        :parameter limit: highest value allowed to be excluded (1 = 1% of BPI)
        :parameter percExtra: how far past the limit values to set axis limits
        :parameter ccsAxisValues: Provide to use CCS axis instead of td.
        :returns: None
        """
        if not ccsAxisValues:
            xvals = self.xvals
        else:
            xvals = ccsAxisValues
        
        ys = self.yvals/self.yvals.max()*100        
        xsAboveLimit = xvals[ys>=limit]
        left = xsAboveLimit[0] - (float(percExtra)/100*xvals.max())
        right = xsAboveLimit[-1] + (float(percExtra)/100*xvals.max())
        ax.set_xlim([left,right])
    
    def autoAxisXoff(self,ax,ccsAxisValues=False):
        """Revert axis limits to original values.

        :parameter ax: matplotlib.axes() object
        :parameter ccsAxisValues: Provide to use CCS axis instead of td.
        :returns: None
        """
        if not ccsAxisValues:
            left = self.xvals.min()
            right = self.xvals.max()
            ax.set_xlim([left,right])
        else:
            left = ccsAxisValues.min()
            right = ccsAxisValues.max()


