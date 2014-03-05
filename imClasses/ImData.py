"""Class for describing ion mobility data. This is used for both
small sections of the full set of mobility data (DataSlice) and whole
datasets (its the superclass of Im()).
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com>"

import Atd
import classes.MassSpectrum as MassSpectrum
import numpy as np
import matplotlib.pyplot as plt
import lib.utils as utils

class ImData():

    def __init__(self):
        # TODO(gns) rename axes, like tdAxis and mzAxis or something
        # xaxis is mz axis
        # yaxis is the arrival time axis
        self.xaxis = np.array([])
        self.yaxis = np.array([])
        self.matrix = np.array([])
        self.matrixUnaltered = np.array([])
        self.xlims = []
        self.ylims = []

        self.xaxisUnaltered = np.array([])
        self.yaxisUnaltered = np.array([])

        self.atd = Atd.Atd()
        self.massSpectrum = MassSpectrum.MassSpectrum()


    # Setters
    def setAxisX(self,x):
        """Set x (m/z) axis

        :parameter x: numpy array
        """
        self.xaxis = x
        self.xlims = [x.min(),x.max()]
        if len(self.xaxisUnaltered) == 0:
            self.xaxisUnaltered = self.xaxis.copy()

    def setAxisY(self,y):
        """Set y (arrival time) axis

        :parameter y: numpy array
        """
        self.yaxis = y
        self.ylims = [y.min(),y.max()]
        if len(self.yaxisUnaltered) == 0:
            self.yaxisUnaltered = self.yaxis.copy()

    def setMatrix(self,matrix):
        """Set the matrix of intensity values (which aligns with the x
        and y axes). 

        :parameter matrix: numpy array (2D)
        """
        self.matrix = matrix
        self.matrixUnaltered = self.matrix.copy()

    # limiting axes
    def limitAxisX(self,lims):
        """Set limits for the m/z axis.

        :parameter lims: list in format of [lower,upper] limits
        """
        self.xlims = [lims[0],lims[1]]
    def limitAxisY(self,lims):
        """Set limits for the arrival time axis.

        :parameter lims: list in format of [lower,upper] limits
        """
        self.ylims = [lims[0],lims[1]]

    # normalisationimport matplotlib.pyplot as plt
    def normalisationBpi(self):
        """Normalise intensity matrix to the peak height at its highest
        point.
        """
        self.matrix = (self.matrixUnaltered/self.matrixUnaltered.max())*100
    def normalisationArea(self):
        """Normalise intensity matrix to the area under the matrix.

        Calculation is done by dividing by the sum of the intensity matrix and
        multipling by the number of entries.
        """
        shape = np.shape(self.matrix)
        self.matrix = self.matrixUnaltered/np.sum(self.matrixUnaltered)*(shape[0]*shape[1])
    def normalisationBpiByMs(self):
        """Not implemented yet."""
        # TODO(gns) - probably remove this, don't think it will work
        # TODO(gns) - its mentioned in the below, so you need to fix that
        # SpectralAveragingGui/SaSettings.py:
        '''
        if not len(self.massSpectrum.xvals):
            self.generateMassSpectrum()
        yvals = self.massSpectrum.yvals
        bpis = yvals/yvals.max()
        # the whole column should sum to the corresponding value in bpi
        # i.e. column = column/np.sum(column)*bpi
        self.matrix = self.matrixUnaltered
        for i,bpi in enumerate(bpis):
            self.matrix[i] = self.matrix[i]/np.sum(self.matrix[i])
            self.matrix[i] = self.matrix[i]*bpi
        '''

    ################################################################
    # generate 2D data
    def generateAtd(self):
        """Sum the intensity matrix to create an arrival time distribution."""
        # Y axis of grid is x axis of atd
        self.atd.setAtdXvals(self.yaxis.copy())
        intensity = np.sum(self.matrix,axis=1)
        self.atd.setAtdYvals(intensity)

    def generateMassSpectrum(self,keepExisting=True):
        """Sum the intensity matrix to create a mass spectrum.

        :parameter keepExisting: if True don't overwrite exisiting mass spectrum for this object
        """
        # Check if there is already a mass spectrum for this object
        if not len(self.massSpectrum.xvals):
            self._definitelyGenerateMassSpectrum()
        elif len(self.massSpectrum.xvals):
            if keepExisting == False:
                self._definitelyGenerateMassSpectrum()
            else:
                pass

    def _definitelyGenerateMassSpectrum(self):
        """Sum the intensity matrix to create a mass spectrum. If there is already
        a mass spectrum for this object, overwrite it.
        """
        self.massSpectrum.xvals = self.xaxis.copy()
        self.massSpectrum.yvals = np.sum(self.matrix,axis=0)
        if not len(self.massSpectrum.xvals) == len(self.massSpectrum.yvals):
            print 'Unequal x and y axis lengths (truncating)'
            self.massSpectrum.xvals = self.massSpectrum.xvals[len(self.massSpectrum.xvals)-len(self.massSpectrum.yvals):]
            self.massSpectrum.yvals = self.massSpectrum.yvals[:len(self.massSpectrum.xvals)]
        self.massSpectrum.rawyvals = self.massSpectrum.yvals.copy()

    def getMassSpectrum(self,keepExisting=True):
        """Sum the intensity matrix to create a mass spectrum and return it.

        :parameter keepExisting: If False recreate the mass spectrum before returning it
        :returns: msClasses.MassSpectrum() object
        """
        # TODO(gns) - make the same function for ATDs
        self.generateMassSpectrum()
        return self.massSpectrum


    ###############################################################
    # Plotting functions
    ###############################################################

    def _getMatrixWithinLimits(self):
        """Applies the self.xlims and self.ylims limits and returns the reduced
        matrix with x and y axes.

        :returns: matrix, x axis, yaxis
        """
        lowIndexX = utils.closest(self.xlims[0], self.xaxis)
        highIndexX = utils.closest(self.xlims[1], self.xaxis)
        x = self.xaxis[lowIndexX:highIndexX]

        lowIndexY = utils.closest(self.ylims[0], self.yaxis)
        highIndexY = utils.closest(self.ylims[1], self.yaxis)
        y = self.yaxis[lowIndexY:highIndexY]

        matrix = self.matrix[lowIndexY:highIndexY,lowIndexX:highIndexX]
        return matrix, x, y

    def plotHeatmap(self,ax=0,**kwargs):
        """Plot a heatmap of this object's data.

        :parameter ax: matplotlib Axes instance or boolean
        :parameter \*\*kwargs: matplotlib.pyplot.plot() arguments
        """
        matrix,x,y = self._getMatrixWithinLimits()
        self.plotHeatmapFromMatrix(matrix,x,y,ax,**kwargs)

    def plotHeatmapFromMatrix(self,matrix,x,y,ax=False,**kwargs):
        """Create a heatmap from provided information instead of using
        this object's properties.

        :parameter matrix: numpy array (intensity data)
        :parameter x: numpy array (m/z axis)
        :parameter y: numpy array (arrival time axis)
        :parameter ax: matplotlib Axes instance or boolean
        :parameter \*\*kwargs: matplotlib.pyplot.plot() arguments
        """
        if not ax:
            ax = plt
        ax.imshow(matrix, origin=[0,0],aspect='auto',
                   extent=[x.min(),x.max(),y.min(),y.max()],**kwargs)
        
    def show(self):
        """Convenience function with calls matplotlib.pyplot.show()."""
        plt.show()
