"""Class for holding additional functions and settings for ContourGui()."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import wx,os,re,copy
import cPickle as pickle
from collections import OrderedDict
import numpy as np

from msClasses import MassSpectrum
from imClasses import Im
import gui.guiFunctions as gf
from lib import utils
from AmphitriteEnums import *

class CtrSettings():

    def __init__(self):
        
        # Atropos stuff
        self.atrOb = None
        self.species = []
        self.speciesSelected = None
        
        self.calibrationOb = None

        # Boundaries
        self.widthL = OrderedDict()
        self.widthR = OrderedDict()
        self.widthL['Apply to all'] = 1.0
        self.widthR['Apply to all'] = 1.0

        # components
        self.plotPanel = None
        self.gui = None
        
        # Data
        self.imOb = None
        self.imOb2 = None

        self.defaultDirectory = '/home/ganesh/Dropbox/workspaces/Amphitrite_2.0/gui/tempData'

        # Misc
        self.limit = "0"

    #====================
    # Adding other panel components
    def setPlotPanel(self,plotPanel):
        """Set the plot panel object to enable accessing its attributes.
        :parameter gui: CtrPlotting() object
        """
        self.plotPanel = plotPanel
    def setGui(self,yourself):
        """Set the GUI object to enable accessing its attributes.
        :parameter gui: ContourGui() object
        """
        self.gui = yourself
    #====================

    ############################################################
    # Width functions
    ############################################################

    def setImOb(self,imOb):
        """Set the primary data object.
        :parameter imOb: imClasses.Im() object
        """
        self.imOb = imOb
        self.setImObsAtropos()
        self.plotPanel.refresh_plot()

    def setImOb2(self,imOb):
        """Set the secondary data object for difference plots.
        :parameter imOb: imClasses.Im() object
        """
        self.imOb2 = imOb
        self.setImObsAtropos()
        self.plotPanel.refresh_plot()        

    def setWidthL(self,val,species):
        """Set the peak width multiplier for the lower m/z value limit.
        :parameter val: Peak width multiplier
        :parameter species: Moleculae species to apply the value
        """
        self.widthL[species] = val

    def setWidthR(self,val,species):
        """Set the peak width multiplier for the upper m/z value limit.
        :parameter val: Peak width multiplier
        :parameter species: Molecular species to apply the value
        """
        self.widthR[species] = val

    def initialiseSpeciesWidths(self,species):
        """Set the peak width multipliers to 1.
        :parameter species: Molecular species
        """
        self.widthL[species] = 1.0
        self.widthR[species] = 1.0
    # Width functions
    ############################################################



    ############################################################
    # Atropos related functions
    ############################################################

    def loadAtroposSpecies(self,path):
        """Load pickled mass spectrum fit.
        :parameter path: Path to pickled Atropos() fit object
        """
        self.atrOb = pickle.load(open(path,'rb'))

        # reset variables
        self.species = ['Apply to all']
        for sp in sorted(self.atrOb.simulatedSpecies.keys()):
            self.species.append(sp)

        for sp in self.species:
            if not sp in self.widthL.keys():
                self.initialiseSpeciesWidths(sp)
            
                # update choices 
                self.gui.choiceSpecies.Append(sp)

        self.setImObsAtropos()


    def setImObsAtropos(self):
        """Set the atropos MS fit to the primary and secondary
        imClasses.Im() objects.
        """
        if self.atrOb:
            if self.imOb:
                #self.imOb.setMsFit(self.atrOb)
                self.imOb.setMsFit(copy.deepcopy(self.atrOb))
            if self.imOb2:
                #self.imOb2.setMsFit(self.atrOb)
                self.imOb2.setMsFit(copy.deepcopy(self.atrOb))

    # Atropos related functions
    ############################################################


    def setCalibration(self,path):
        """Load pickled calibration object.
        :parameter path: Absolute path to imClasses.Calibration() object
        """
        import imClasses.Calibration
        try:
            self.calibrationOb = pickle.load(open(path,'rb'))
        except:
            message = 'Something wrong with calibration file!'
            gf.warningDialog(message)

        
    def getChargeAndMassDictionaries(self):
        """Return the charge and mass dictionaries for each species.
        :returns: zsD and massD in the form of d[speciesName] = value
        """
        zsD = OrderedDict()
        massD = OrderedDict()
        for i,sp in enumerate(self.species):
            if i:
                zsD[sp] = self.atrOb.species[sp].charges
                massD[sp] = self.atrOb.species[sp].mass
        return zsD,massD
        
