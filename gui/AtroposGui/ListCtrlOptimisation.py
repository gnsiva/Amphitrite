"""ListCtrl for holding the results of the deconvolution."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import wx
from EditableListCtrl import EditableListCtrl
from lib import utils
import collections

class ListCtrlOptimisation(wx.Panel):
    
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,-1)
        
        self.listCtrl = EditableListCtrl(self)
        self.listCtrl.setColumnToIgnore(0)
        
        self.line = collections.OrderedDict()
        self.line["Mass"] = 0
        self.line["Total Area"] = 1
        self.line["% Area"] = 2 
        self.line["Mu"] = 3
        self.line["Amplitude"] = 4 
        self.line["FWHM"] = 5
        self.line["Peak FWHM"] = 6
        self.line["Total Intensity"] = 7
        self.line["% Intensity"] = 8
        self.setParameterColumn()
        
        sizer_1 = wx.FlexGridSizer(1, 1, 0, 0) 
        sizer_1.Add(self.listCtrl, 1, wx.EXPAND)  
        sizer_1.AddGrowableRow(0, 1)
        sizer_1.AddGrowableCol(0, 1)
        self.SetSizer(sizer_1)          
     

        
    def setParameterColumn(self):
        """Fill in the first column of the ListCtrl, which is the
        names for all the parameters.
        """
        self.listCtrl.InsertColumn(0,"Parameter",width=100)
        
        for name, row in self.line.items():
            self.listCtrl.InsertStringItem(row, "" )
            self.listCtrl.SetStringItem(row, 0, name)
        '''
        self.listCtrl.SetStringItem(0, 0, "Mass")
        self.listCtrl.SetStringItem(1, 0, "Total Area")
        self.listCtrl.SetStringItem(2, 0, "Mu")
        self.listCtrl.SetStringItem(3, 0, "Amplitude")
        self.listCtrl.SetStringItem(4, 0, "FWHM")
        self.listCtrl.SetStringItem(5, 0, "Peak FWHM")
        self.listCtrl.SetStringItem(6, 0, "Total Intensity")
        '''
            
            
    def setSettings(self,settings):
        """Set the AtroposGuiSettings object, so its functions and
        attributes can be accessed more readily.

        :parameter settings: AtroposGuiSettings() object
        """
        self.settings = settings

    def addSpeciesColumn(self,speciesName):
        """Add a new column to the ListCtrl for a molecular species.

        :parameter species: Name of species for this column
        """
        cols = self.listCtrl.GetColumnCount()
        self.listCtrl.InsertColumn(cols,speciesName,width=100)
    
    def removeSpeciesColumns(self):
        """Clear the ListCtrl. Is run before refilling with updated
        information.
        """
        self.listCtrl.DeleteAllColumns()
        self.setParameterColumn()        
    
    def setSpeciesColumns(self,simulatedSpecies):
        """Fill in ListCtrl using the result of deconvolution.

        :parameter simulatedSpecies: list of the species names which were deconvoluted
        """
        self.removeSpeciesColumns()
        for i,speciesName in enumerate(simulatedSpecies):
            # get species name and column
            speciesOb = self.settings.msPanel.ms.simulatedSpecies[speciesName]
            column = i+1
            
            self.addSpeciesColumn(speciesName)
            
            # Fill in the data
            self.listCtrl.SetStringItem(self.line["Mass"],column,"%.2f" %speciesOb.mass)
            
            # % Total Area - To Calculate
            totArea = speciesOb.getTotalArea(self.settings.msPanel.ms.xvals)
            self.listCtrl.SetStringItem(self.line["Total Area"],column,"%.1f" %totArea)
            self.listCtrl.SetStringItem(self.line["Mu"],column,"%.1f" %speciesOb.zGauss.centre)
            self.listCtrl.SetStringItem(self.line["Amplitude"],column,"%.1f" %speciesOb.zGauss.amplitude)
            self.listCtrl.SetStringItem(self.line["FWHM"],column,"%.1f" %speciesOb.zGauss.fwhm)
            self.listCtrl.SetStringItem(self.line["Peak FWHM"],column,"%.1f" %speciesOb.peakFwhm)
            
            # Total Intensity - To Calculate
            totIts = speciesOb.getTotalIntensity(self.settings.msPanel.ms.xvals,
                                           self.settings.msPanel.ms.yvals)
            self.listCtrl.SetStringItem(self.line["Total Intensity"],column,"%.1f" %totIts)
        
    def calculatePercentageArea(self):
        """Calculate the percentage abundance (by area) of each molecular species
        deconvolution. Fill in the ListCtrl with the results.
        """
        cols = self.listCtrl.GetColumnCount()
        areas = []
        for column in xrange(cols):
            if column:
                areas.append(float(self.listCtrl.GetItem(self.line["Total Area"],column).GetText()))
            else:
                areas.append(0)
        for column in xrange(cols):
            if column:
                percArea = (areas[column]/sum(areas))*100
                self.listCtrl.SetStringItem(self.line["% Area"],column,"%.1f" %percArea)   

    def calculatePercentageIntensity(self):
        """Calculate the percentage abundance (by peak intensity) of each molecular
        species deconvolution. Fill in the ListCtrl with the results.
        """
        cols = self.listCtrl.GetColumnCount()
        peaks = []
        for column in xrange(cols):
            if column:
                peaks.append(float(self.listCtrl.GetItem(self.line["Total Intensity"],column).GetText()))
            else:
                peaks.append(0)
        for column in xrange(cols):
            if column:
                percIts = (peaks[column]/sum(peaks))*100
                self.listCtrl.SetStringItem(self.line["% Intensity"],column,"%.1f" %percIts)        
                     
         
        
          
