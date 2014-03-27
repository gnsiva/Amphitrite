"""Class for processing the state of display checkboxes in IesGui."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

from AmphitriteEnums import *

class IesCheckboxStates():

    def __init__(self):
        self.atds = False
        self.contourPlots = False
        self.conformations = False
        
        self.massSpectra = True

        self.checkBoxes = [self.atds,self.contourPlots,self.conformations,
                           self.massSpectra]
        
        self.singleCheckBoxEnums = [IesPlotState.ATD,IesPlotState.CONTOUR,
                                    IesPlotState.CONFORMATIONS,
                                    IesPlotState.MASS_SPECTRA]

    def printStates(self):
        """Print to console the states of each checkbox.
        """
        print '================================'
        print '====== Check box states'
        print '%d - atds' %self.atds
        print '%d - contour' %self.contourPlots
        print '%d - conformations' %self.conformations
        print '%d - mass spectra' %self.massSpectra
        print '================================'
        
    def getEnum(self):
        """Return the enum for the current state shown in checkboxes.
        :parameter enum: Enum from AmphitriteEnums.py
        """
        self.updateCheckboxList()
        # Mass Spectrum
        if self.massSpectra:
            return IesPlotState.MASS_SPECTRA
        
        # Single plot (including None)
        if sum(self.checkBoxes[:-1]) == 1:
            i = self.checkBoxes[:-1].index(True)
            return self.singleCheckBoxEnums[i]
        elif sum(self.checkBoxes) == 0:
            return IesPlotState.NONE

        # Double plot
        if sum(self.checkBoxes[:-1]) == 2:
            if self.atds:
                if self.contourPlots:
                    return IesPlotState.CONTOUR_ATD
                elif self.conformations:
                    return IesPlotState.ATD_CONFORMATIONS
                
            if self.contourPlots:
                if self.conformations:
                    return IesPlotState.CONTOUR_CONFORMATIONS

        # Triple plot
        if sum(self.checkBoxes[:-1]) == 3:
            return IesPlotState.CONTOUR_ATD_CONFORMATIONS

    def getPlotType(self):
        """Return the plot type which should be displayed according to
        checkbox states.
        :returns: (int) 1: single plot are, 2: Double plot and 3: Triple plot
        """
        self.updateCheckboxList()
        print self.checkBoxes
        # Mass Spectrum
        if self.massSpectra:
            return 1
        
        # Single plot (including None)
        if sum(self.checkBoxes[:-1]) == 1:
            return 1
        elif sum(self.checkBoxes[:-1]) == 0:
            return 1
        
        # Double plot
        if sum(self.checkBoxes[:-1]) == 2:
            return 2

        # Triple plot
        if sum(self.checkBoxes[:-1]) == 3:
            return 3

        def updateCheckboxList(self):
            """Update the self.checkBoxes property after changing the state
            of one of the checkboxes in this object.
            """
            self.checkBoxes = [self.atds,self.contourPlots,self.conformations,
                               self.massSpectra]
