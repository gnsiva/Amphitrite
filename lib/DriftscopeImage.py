"""Driftscope image class:
Draws a simple driftscope style contour plot with the
mass spectrum below it.
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import re
import cPickle as pickle

class DriftscopeImage():
    def __init__(self,filename=False):
        """
        :parameter filename: Absolute path to Amphitrite data file - draw figure
        and save in same folder as data file. OR False, then just instantiate
        class.
        """
        self.mzs = None
        self.tds = None
        self.mobility = None
    
    def loadFile(self,filename):
        """Load an amphitrite data file ('.a').
        :parameter filename: Absolute path to file
        """
        # load data
        npObj = np.load(filename)
        self.mzs = npObj[0]
        self.tds = npObj[1]
        self.mobility = npObj[2]
        self.its = np.sum(self.mobility,axis=0)
        # process data
        self.normaliseMobility()
        
    #===========================================================================
    # Data preparation
    #===========================================================================
    def normaliseMobility(self):
        """Normalise the intensity of the intensity matrix
        to its highest point. (Necessary for proper colouring of
        heatmap).
        """
        self.mobility = (self.mobility/self.mobility.max())*1000
        self.mobility = np.log2(self.mobility)    
        infinite = np.isinf(self.mobility)
        self.mobility[infinite] = 0
        self.mobility[self.mobility<1] = 0
        
    #===========================================================================
    # THE GRAPH
    #===========================================================================
    def driftscope(self,mz_range=0,a_range=0):
        """Draw the driftscope image. Top region is the contour
        plot, bottom section is mass spectrum.
        """
        #================
        # Set up colour scheme for driftscope image
        # TODO(gns) - Change the colour map used 
        grayNums = ['0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.5']
        c_array = ['black']+grayNums+['r','magenta','purple','b','c','w']
        cmap = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',c_array,256)

        #================
        # Setup ranges to plot
        if not mz_range:
            mz_range = [self.mzs[0],self.mzs[-1]]
        if not a_range:
            a_range = [self.tds[0],self.tds[-1]]

        mzLims = [mz_range[0],mz_range[1]]
        tdLims = [a_range[0],a_range[1]]

        #================
        # Do plotting
        
        f = plt.figure(figsize=[6,8]) #FIGSIZE
        gs = gridspec.GridSpec(2,1)

        # Contour plot
        ax1 = plt.subplot(gs[0],xlim=mzLims,ylim=tdLims,ylabel='Arrival Time ($ms$)',xlabel='$m/z$')
        ax1.imshow(self.mobility,aspect='auto',origin='lower',cmap=cmap,extent=mzLims+tdLims)

        # Mass Spectrum
        ax2 = plt.subplot(gs[1], xlabel="$m/z$",xlim=mzLims,yticks=[],ylabel='Intensity')
        ax2.plot(self.mzs,self.its,color='black')

        try:
            plt.tight_layout()
        except:
            pass
    
    #===========================================================================
    # The Function
    #===========================================================================
    def drawFigure(self,filename):
        """Run the program. Open Amphitrite data file, draw a
        driftscope-like diagram and save the figure in the same
        directory as the data file.
        """
        outputDir = os.path.dirname(filename)
        self.loadFile(filename)
        self.driftscope()

        searcher = re.search('(.*)\.a$',filename)
        if searcher:
            fn = searcher.group(1) + "_ATDandMS.png"
        else:
            fn = filename + "_ATDandMS.png"
        plt.savefig(fn)

# Example usage        
'''
fn = "/home/gns/Downloads/BetaLactoglobulin.a"
ob = DriftscopeImage(fn)
'''

