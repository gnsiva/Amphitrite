"""Class stores parameters for the GUI as well as carries out the calculations."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

from lib import utils

class AtroposGuiSettings():
    
    def __init__(self):
        self.openFilePath = ''
        self.msPanel = None
        self.listCtrlSmoothing = None
        self.listCtrlAssigningSpecies = None
        self.listCtrlOptimise = None
        
        self.atroposPath = ''
        
        self.speciesColumn = None
        self.speciesSimulated = []
        self.addedPeakMzs = []
        
    def setMsPanel(self,msPanel):
        """Store the plotting panel in this object."""
        self.msPanel = msPanel
    def setListCtrlSmoothing(self,listCtrlSmoothing):
        """Store the smoothing parameter list control in this object."""
        self.listCtrlSmoothing = listCtrlSmoothing
    def setListCtrlAssigningSpecies(self,listCtrlAssigningSpecies):
        """Store the assigning species parameters list control in this object."""
        self.listCtrlAssigningSpecies = listCtrlAssigningSpecies
    def setListCtrlOptimise(self,listCtrlOptimise):
        """Store the optimisation results list control in this object."""
        self.listCtrlOptimise = listCtrlOptimise

    def setSpeciesColumn(self,col):
        """Set the currently selected species in ListCtrlAssigningSpecies
        by the column index.

        :parameter col: column index
        """
        self.assignSpeciesColumn = col
    def getSpeciesColumn(self):
        """Return the last accessed species column in ListCtrlAssigningSpecies

        :returns: column index
        """
        return self.assignSpeciesColumn
    
    def getLimit(self):
        """Get the lower limit for peak finding. Value is given as a
        percentage of base peak intensity.
        """
        limit = self.listCtrlSmoothing.listCtrl.GetItem(4,1).GetText()
        if utils.isNumber(limit):
            return float(limit)
        else:
            print 'Limit value is not a number (using default - 0.0)'    
            return 0

    
    def addPeakMz(self,mz):
        """When the user adds a peak manually use this function to keep
        the peak list up to date for this object.

        :parameter mz: m/z value for peak being added
        """
        self.addedPeakMzs.append(mz)
        
    def leastSquaresOptimisation(self,oneFwhm):
        """Carry out deconvolution of the mass spectrum.

        :parameter oneFwhm: boolean (use single FWHM value for all species)
        """
        speciesNames = self.listCtrlAssigningSpecies.getSpeciesToSimulate()
        self.speciesToSimulate = speciesNames
        self.msPanel.ms.leastSquaresOptimisation(speciesNames,oneFwhm)
       
        self.msPanel.plotLeastSquaresSimulation(speciesNames)
        self.listCtrlOptimise.setSpeciesColumns(speciesNames)
        self.speciesSimulated = speciesNames

    def plotLeastSquaresSimulationForToggle(self):
        """Tell the plot panel to display the deconvolution. Only plot
        the species in self.speciesSimulated.
        """
        self.msPanel.plotLeastSquaresSimulation(self.speciesSimulated)

    def getChargesToSimulate():
        """Get charge list from ListCtrlAssigningSpecies.
        Does not automatically check for valid input.
        """
        zs = self.listCtrlAssigningSpecies.listCtrl.GetItem(4,1).GetText()
        return utils.getHyphenCommaList(zs)
        
            
            
