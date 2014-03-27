#!/usr/bin/env python
# -*- coding: utf-8 -*-
# generated by wxGlade 0.6.4 on Wed Apr 10 11:18:42 2013

import wx, os
import ListCtrlFilesIes,ListCtrlConformationsIes
import wx.lib.agw.multidirdialog as MDD
import gui.guiFunctions as gf
import cPickle as pickle
import IesSettings,IesPlotPanel

# begin wxGlade: extracode
# end wxGlade

"""Program for analysing IM-MS data for multiple conditions e.g. protein
unfolding experiments."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

class IesGui(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: IesGui.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.window_1 = wx.SplitterWindow(self, -1, style=wx.SP_3D | wx.SP_BORDER)
        self.window_1_pane_1 = wx.Panel(self.window_1, -1)
        self.label_7 = wx.StaticText(self.window_1_pane_1, -1, "Open files")
        self.buttonOpenFiles = wx.Button(self.window_1_pane_1, -1, "Open")
        self.buttonCloseFile = wx.Button(self.window_1_pane_1, -1, "Close")
        self.panelListCtrlFiles = wx.Panel(self.window_1_pane_1, -1)
        self.label_8 = wx.StaticText(self.window_1_pane_1, -1, "Value Units")
        self.textCtrlUnits = wx.TextCtrl(self.window_1_pane_1, -1, "")
        self.sizer_2_staticbox = wx.StaticBox(self.window_1_pane_1, -1, "Files")
        self.label_2 = wx.StaticText(self.window_1_pane_1, -1, "Atropos file ")
        self.textCtrlAtropos = wx.TextCtrl(self.window_1_pane_1, -1, "")
        self.buttonAtropos = wx.Button(self.window_1_pane_1, -1, "Open")
        self.label_5 = wx.StaticText(self.window_1_pane_1, -1, "Species")
        self.choiceSpecies = wx.Choice(self.window_1_pane_1, -1, choices=[])
        self.label_6 = wx.StaticText(self.window_1_pane_1, -1, "Charge State")
        self.choiceChargeState = wx.Choice(self.window_1_pane_1, -1, choices=[])
        self.label_3 = wx.StaticText(self.window_1_pane_1, -1, "Width multiplier: L  ", style=wx.ALIGN_RIGHT)
        self.textCtrlWidthL = wx.TextCtrl(self.window_1_pane_1, -1, "1.0")
        self.label_4 = wx.StaticText(self.window_1_pane_1, -1, "R  ")
        self.textCtrlWidthR = wx.TextCtrl(self.window_1_pane_1, -1, "1.0")
        self.static_line_1 = wx.StaticLine(self.window_1_pane_1, -1)
        self.label_1 = wx.StaticText(self.window_1_pane_1, -1, "Calibration file ")
        self.textCtrlCalibration = wx.TextCtrl(self.window_1_pane_1, -1, "")
        self.buttonCalibration = wx.Button(self.window_1_pane_1, -1, "Open")
        self.sizer_3_staticbox = wx.StaticBox(self.window_1_pane_1, -1, "Input")
        self.buttonProcess = wx.Button(self.window_1_pane_1, -1, "Process")
        self.sizer_4_staticbox = wx.StaticBox(self.window_1_pane_1, -1, "")
        self.panelPlotPanel = wx.Panel(self.window_1_pane_1, -1)
        self.window_1_pane_2 = wx.Panel(self.window_1, -1)
        self.checkboxContour = wx.CheckBox(self.window_1_pane_2, -1, "Contour plots")
        self.checkboxAtds = wx.CheckBox(self.window_1_pane_2, -1, "ATDs / CCSDs")
        self.checkboxConformations = wx.CheckBox(self.window_1_pane_2, -1, "Conformation Tracking")
        self.checkboxMassSpectra = wx.CheckBox(self.window_1_pane_2, -1, "Mass Spectra")
        self.sizer_5_staticbox = wx.StaticBox(self.window_1_pane_2, -1, "Display")
        self.checkboxPeakMarkers = wx.CheckBox(self.window_1_pane_2, -1, "Peak Markers")
        self.buttonOneClick = wx.ToggleButton(self.window_1_pane_2, -1, "One Click")
        self.buttonToggle = wx.ToggleButton(self.window_1_pane_2, -1, "Toggle")
        self.sizer_6_staticbox = wx.StaticBox(self.window_1_pane_2, -1, "Pick Peak")
        self.panelListCtrlConformations = wx.Panel(self.window_1_pane_2, -1)
        self.buttonEdit = wx.Button(self.window_1_pane_2, -1, "Edit")
        self.buttonDelete = wx.Button(self.window_1_pane_2, -1, "Delete")
        self.sizer_7_staticbox = wx.StaticBox(self.window_1_pane_2, -1, "Conformations")

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_BUTTON, self.eventButtonOpenFiles, self.buttonOpenFiles)
        self.Bind(wx.EVT_BUTTON, self.eventButtonCloseFile, self.buttonCloseFile)
        self.Bind(wx.EVT_TEXT, self.eventTextCtrlUnits, self.textCtrlUnits)
        self.Bind(wx.EVT_BUTTON, self.eventButtonAtropos, self.buttonAtropos)
        self.Bind(wx.EVT_CHOICE, self.eventChoiceSpecies, self.choiceSpecies)
        self.Bind(wx.EVT_CHOICE, self.eventChoiceChargeState, self.choiceChargeState)
        self.Bind(wx.EVT_TEXT, self.eventTextCtrlWidthL, self.textCtrlWidthL)
        self.Bind(wx.EVT_TEXT, self.eventTextCtrlWidthR, self.textCtrlWidthR)
        self.Bind(wx.EVT_BUTTON, self.eventButtonCalibration, self.buttonCalibration)
        self.Bind(wx.EVT_BUTTON, self.eventButtonProcess, self.buttonProcess)
        self.Bind(wx.EVT_CHECKBOX, self.eventCheckboxContour, self.checkboxContour)
        self.Bind(wx.EVT_CHECKBOX, self.eventCheckboxAtds, self.checkboxAtds)
        self.Bind(wx.EVT_CHECKBOX, self.eventCheckboxConformations, self.checkboxConformations)
        self.Bind(wx.EVT_CHECKBOX, self.eventCheckboxMassSpectra, self.checkboxMassSpectra)
        self.Bind(wx.EVT_CHECKBOX, self.eventCheckboxPeakMarkers, self.checkboxPeakMarkers)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.eventButtonOneClick, self.buttonOneClick)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.eventButtonToggle, self.buttonToggle)
        self.Bind(wx.EVT_BUTTON, self.eventButtonEdit, self.buttonEdit)
        self.Bind(wx.EVT_BUTTON, self.eventButtonDelete, self.buttonDelete)
        self.Bind(wx.EVT_SPLITTER_UNSPLIT, self.eventWindow_1Unsplit, self.window_1)
        # end wxGlade

        ## Importing each other
        self.settings = IesSettings.IesSettings(self.plotPanel)
        self.settings.setListCtrlConformations(self.listCtrlConformations)
        self.settings.setPlotPanel(self.plotPanel)

        self.plotPanel.setSettings(self.settings)
        
        self.listCtrlFiles.setSettings(self.settings)
        self.listCtrlFiles.setPlotPanel(self.plotPanel)

        self.listCtrlConformations.setSettings(self.settings)
        self.listCtrlConformations.setPlotPanel(self.plotPanel)

        self.window_1.SetSashGravity(1.0)

    def __set_properties(self):
        # begin wxGlade: IesGui.__set_properties
        self.SetTitle("Injection Energy Study GUI Beta")
        self.SetSize((944, 612))
        self.buttonCloseFile.Enable(False)
        self.choiceSpecies.Enable(False)
        self.choiceChargeState.Enable(False)
        self.checkboxContour.Enable(False)
        self.checkboxAtds.Enable(False)
        self.checkboxConformations.Enable(False)
        self.checkboxMassSpectra.SetValue(1)
        self.checkboxPeakMarkers.SetValue(1)
        # end wxGlade

        #self.listCtrlFiles = ListCtrlFilesIes.ListCtrlFiles(self.panelListCtrlFiles)
        self.panelListCtrlFiles.Hide()
        self.panelListCtrlFiles = ListCtrlFilesIes.ListCtrlFiles(self.window_1_pane_1)
        self.listCtrlFiles = self.panelListCtrlFiles
        self.plotPanel = IesPlotPanel.IesPlotPanel(self.panelPlotPanel,self)
        self.panelListCtrlConformations.Hide()
        self.panelListCtrlConformations = ListCtrlConformationsIes.ListCtrlConformations(self.window_1_pane_2)
        self.listCtrlConformations = self.panelListCtrlConformations
        
    def __do_layout(self):
        # begin wxGlade: IesGui.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_10 = wx.FlexGridSizer(3, 1, 0, 0)
        self.sizer_7_staticbox.Lower()
        sizer_7 = wx.StaticBoxSizer(self.sizer_7_staticbox, wx.HORIZONTAL)
        grid_sizer_13 = wx.FlexGridSizer(2, 1, 0, 0)
        grid_sizer_14 = wx.GridSizer(1, 2, 0, 0)
        self.sizer_6_staticbox.Lower()
        sizer_6 = wx.StaticBoxSizer(self.sizer_6_staticbox, wx.HORIZONTAL)
        grid_sizer_12 = wx.FlexGridSizer(2, 1, 0, 0)
        grid_sizer_15 = wx.FlexGridSizer(1, 2, 0, 0)
        self.sizer_5_staticbox.Lower()
        sizer_5 = wx.StaticBoxSizer(self.sizer_5_staticbox, wx.HORIZONTAL)
        grid_sizer_11 = wx.FlexGridSizer(5, 1, 0, 0)
        grid_sizer_1 = wx.FlexGridSizer(1, 2, 0, 0)
        grid_sizer_2 = wx.FlexGridSizer(3, 1, 0, 0)
        self.sizer_4_staticbox.Lower()
        sizer_4 = wx.StaticBoxSizer(self.sizer_4_staticbox, wx.HORIZONTAL)
        grid_sizer_9 = wx.FlexGridSizer(1, 1, 0, 0)
        self.sizer_3_staticbox.Lower()
        sizer_3 = wx.StaticBoxSizer(self.sizer_3_staticbox, wx.HORIZONTAL)
        grid_sizer_3 = wx.FlexGridSizer(4, 1, 0, 0)
        grid_sizer_16 = wx.FlexGridSizer(1, 3, 0, 0)
        grid_sizer_5 = wx.FlexGridSizer(6, 2, 0, 0)
        grid_sizer_4 = wx.FlexGridSizer(2, 3, 0, 0)
        self.sizer_2_staticbox.Lower()
        sizer_2 = wx.StaticBoxSizer(self.sizer_2_staticbox, wx.HORIZONTAL)
        grid_sizer_6 = wx.FlexGridSizer(3, 1, 0, 0)
        grid_sizer_8 = wx.FlexGridSizer(1, 2, 0, 0)
        grid_sizer_7 = wx.FlexGridSizer(1, 3, 0, 0)
        grid_sizer_7.Add(self.label_7, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_7.Add(self.buttonOpenFiles, 0, 0, 0)
        grid_sizer_7.Add(self.buttonCloseFile, 0, 0, 0)
        grid_sizer_7.AddGrowableCol(0)
        grid_sizer_6.Add(grid_sizer_7, 1, wx.EXPAND, 0)
        grid_sizer_6.Add(self.panelListCtrlFiles, 1, wx.EXPAND, 0)
        grid_sizer_8.Add(self.label_8, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_8.Add(self.textCtrlUnits, 0, 0, 0)
        grid_sizer_8.AddGrowableCol(0)
        grid_sizer_6.Add(grid_sizer_8, 1, wx.EXPAND, 0)
        grid_sizer_6.AddGrowableRow(1)
        grid_sizer_6.AddGrowableCol(0)
        sizer_2.Add(grid_sizer_6, 1, wx.EXPAND, 0)
        grid_sizer_2.Add(sizer_2, 1, wx.EXPAND, 0)
        grid_sizer_4.Add(self.label_2, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_4.Add(self.textCtrlAtropos, 0, wx.EXPAND, 0)
        grid_sizer_4.Add(self.buttonAtropos, 0, 0, 0)
        grid_sizer_4.AddGrowableCol(1)
        grid_sizer_3.Add(grid_sizer_4, 1, wx.EXPAND, 0)
        grid_sizer_5.Add(self.label_5, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_5.Add(self.choiceSpecies, 0, wx.EXPAND, 0)
        grid_sizer_5.Add(self.label_6, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_5.Add(self.choiceChargeState, 0, wx.EXPAND, 0)
        grid_sizer_5.Add(self.label_3, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_5.Add(self.textCtrlWidthL, 0, wx.EXPAND, 0)
        grid_sizer_5.Add(self.label_4, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_5.Add(self.textCtrlWidthR, 0, wx.EXPAND, 0)
        grid_sizer_5.AddGrowableCol(1)
        grid_sizer_3.Add(grid_sizer_5, 1, wx.EXPAND, 0)
        grid_sizer_3.Add(self.static_line_1, 0, wx.EXPAND, 0)
        grid_sizer_16.Add(self.label_1, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_16.Add(self.textCtrlCalibration, 0, wx.EXPAND, 0)
        grid_sizer_16.Add(self.buttonCalibration, 0, 0, 0)
        grid_sizer_3.Add(grid_sizer_16, 1, wx.EXPAND, 0)
        sizer_3.Add(grid_sizer_3, 1, wx.EXPAND, 0)
        grid_sizer_2.Add(sizer_3, 1, wx.EXPAND, 0)
        grid_sizer_9.Add(self.buttonProcess, 0, wx.EXPAND, 0)
        grid_sizer_9.AddGrowableCol(0)
        sizer_4.Add(grid_sizer_9, 1, wx.EXPAND, 0)
        grid_sizer_2.Add(sizer_4, 1, wx.EXPAND, 0)
        grid_sizer_2.AddGrowableRow(0)
        grid_sizer_1.Add(grid_sizer_2, 1, wx.EXPAND, 0)
        grid_sizer_1.Add(self.panelPlotPanel, 1, wx.EXPAND, 0)
        self.window_1_pane_1.SetSizer(grid_sizer_1)
        grid_sizer_1.AddGrowableRow(0)
        grid_sizer_1.AddGrowableCol(1)
        grid_sizer_11.Add(self.checkboxContour, 0, 0, 0)
        grid_sizer_11.Add(self.checkboxAtds, 0, 0, 0)
        grid_sizer_11.Add(self.checkboxConformations, 0, 0, 0)
        grid_sizer_11.Add(self.checkboxMassSpectra, 0, 0, 0)
        sizer_5.Add(grid_sizer_11, 1, wx.EXPAND, 0)
        grid_sizer_10.Add(sizer_5, 1, wx.EXPAND, 0)
        grid_sizer_12.Add(self.checkboxPeakMarkers, 0, 0, 0)
        grid_sizer_15.Add(self.buttonOneClick, 0, wx.EXPAND, 0)
        grid_sizer_15.Add(self.buttonToggle, 0, wx.EXPAND, 0)
        grid_sizer_15.AddGrowableCol(0)
        grid_sizer_15.AddGrowableCol(1)
        grid_sizer_12.Add(grid_sizer_15, 1, wx.EXPAND, 0)
        grid_sizer_12.AddGrowableCol(0)
        sizer_6.Add(grid_sizer_12, 1, wx.EXPAND, 0)
        grid_sizer_10.Add(sizer_6, 1, wx.EXPAND, 0)
        grid_sizer_13.Add(self.panelListCtrlConformations, 1, wx.EXPAND, 0)
        grid_sizer_14.Add(self.buttonEdit, 0, wx.EXPAND, 0)
        grid_sizer_14.Add(self.buttonDelete, 0, wx.EXPAND, 0)
        grid_sizer_13.Add(grid_sizer_14, 1, wx.EXPAND, 0)
        grid_sizer_13.AddGrowableRow(0)
        grid_sizer_13.AddGrowableCol(0)
        sizer_7.Add(grid_sizer_13, 1, wx.EXPAND, 0)
        grid_sizer_10.Add(sizer_7, 1, wx.EXPAND, 0)
        self.window_1_pane_2.SetSizer(grid_sizer_10)
        grid_sizer_10.AddGrowableRow(2)
        grid_sizer_10.AddGrowableCol(0)
        self.window_1.SplitVertically(self.window_1_pane_1, self.window_1_pane_2, 733)
        sizer_1.Add(self.window_1, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        self.Layout()
        # end wxGlade

    def eventButtonOpenFiles(self, event):  # wxGlade: IesGui.<event_handler>
        """Open FileDialog which is capable of opening multiple files.
        """
        home = os.path.expanduser('~')
        dlg = wx.FileDialog(self,message="Choose Amphitrite file(s)",
                            wildcard="Amphitrite IM file (*.a)|*.a",
                            style=wx.MULTIPLE)

        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            self.listCtrlFiles.addFiles(paths)
        else:
            print 'multidir dialog error'
        dlg.Destroy()
        

    def eventChoiceSpecies(self, event):  # wxGlade: IesGui.<event_handler>
        """Triggered when the species from a mass spectrum fit selected.
        """
        self._populateCharges()
        
    def eventChoiceChargeState(self, event):  # wxGlade: IesGui.<event_handler>
        """Triggered when user selects the charge state of a species from a mass
        spectrum fit. Draws extraction range onto mass spectrum.
        """
        sp = self._getChoiceSpecies()
        z = self._getChoiceChargeState()
        self.settings.setSpeciesAndCharge(sp,z)
        self.plotPanel.plotMassSpectraWidthLimits()
        self.plotPanel.draw()

    def eventButtonProcess(self, event):  # wxGlade: IesGui.<event_handler>
        """Extract the ATDs from the mass spectrum, and display them in the plot window.
        """
        # Turn on appropriate toggles
        self.checkboxAtds.SetValue(True)
        self.settings.checkboxStates.atds = True
        
        self.checkboxMassSpectra.SetValue(False)
        self.setCheckboxEnableDisable(False)
        self.settings.checkboxStates.massSpectra = False
        
        self.plotPanel.refresh_plot()



    def eventButtonOneClick(self, event):  # wxGlade: IesGui.<event_handler>
        """Turns on peak picking (by clicking on an ATD or CCSD). Each click
        adds another peak.
        """
        # deal with other toggle picking button
        if self.buttonToggle.GetValue():
            self.buttonToggle.SetValue(False)
            self.plotPanel.togglePicking(False)
        
        self.plotPanel.oneClickPicking(self.buttonOneClick.GetValue())


    def eventButtonToggle(self, event):  # wxGlade: IesGui.<event_handler>
        """Turns on peak picking (by clicking on an ATD or CCSD). A peak isn't
        registered until the button is pressed again (this allows for more accurate
        peak positioning).
        """
        # deal with other toggle picking button
        if self.buttonOneClick.GetValue():
            self.buttonOneClick.SetValue(False)
            self.plotPanel.oneClickPicking(False)
        
        self.plotPanel.togglePicking(self.buttonToggle.GetValue())


    def eventButtonEdit(self, event):  # wxGlade: IesGui.<event_handler>
        # TODO(gns) - unfinished function
        print "Event handler `eventButtonEdit' not implemented!"
        event.Skip()

    def eventButtonDelete(self, event):  # wxGlade: IesGui.<event_handler>
        # TODO(gns) - unfinished funciton
        print "Event handler `eventButtonDelete' not implemented!"
        event.Skip()

    def eventButtonCalibration(self, event):  # wxGlade: IesGui.<event_handler>
        """Open an amphitrite IM calibration file.
        """
        path = gf.openCalibrationFile(self)
        if path:
            try:
                self.settings.setCalibration(path)
                self.textCtrlCalibration.SetValue(path)
            except:
                gf.warningDialog('Problem with calibration file\nTry again')

    def eventButtonAtropos(self, event):  # wxGlade: IesGui.<event_handler>
        """Open an amphitrite mass spectrum fit file.
        """
        path = gf.openAtroposFile(self)
        if path:
            try:
                self.settings.loadAtroposSpeciesAndCharges(path)
                self.textCtrlAtropos.SetValue(path)
                self._populateSpeciesChoices()
            except:
                gf.warningDialog('Problem with afit file\nTry again')

        #else:
        #    print 'Atropos path problem: %s' %path

    def _populateSpeciesChoices(self):
        """Update dropdown list showing the different molecular species in
        the mass spectrum fit.
        """
        self.choiceSpecies.Enable(True)
        self.choiceSpecies.Clear()
        
        spChoices = self.settings.speciesCharges.keys()
        for sp in spChoices:
            self.choiceSpecies.Append(item=sp)

    def _populateCharges(self):
        """Update dropdown list showing available charges (from the selected
        species (dropdown) in the mass spectrum fit).
        """        
        sp = self._getChoiceSpecies()
        zs = self.settings.speciesCharges[sp]

        self.choiceChargeState.Enable(True)
        self.choiceChargeState.Clear()
        for z in zs:
            self.choiceChargeState.Append(str(z))


    def _getChoiceSpecies(self):
        """Get the currently selected species from the dropdown.
        :returns: Species name (string)
        """
        i = self.choiceSpecies.GetCurrentSelection()
        sp = self.settings.species[i]
        self.setBackToMs()
        return sp
    
    def _getChoiceChargeState(self):
        """Get the currently selected charge state from dropdown.
        :returns: Charge state (int)
        """
        sp = self._getChoiceSpecies()
        i = self.choiceChargeState.GetCurrentSelection()
        z = self.settings.speciesCharges[sp][i]
        self.setBackToMs()
        return z

    def setBackToMs(self):
        """Change the display to show mass spectra.
        """
        self.checkboxMassSpectra.SetValue(True)
        self.settings.checkboxStates.massSpectra = True
        self.setCheckboxEnableDisable()
        try:
            self.plotPanel.plotMassSpectraWidthLimits()
        except:
            self.plotPanel.plotMassSpectra()
        
            
    def eventButtonCloseFile(self, event):  # wxGlade: IesGui.<event_handler>
        # TODO(gns) - unfinished function (should probably be renamed remove or something)
        print "Event handler `eventButtonCloseFile' not implemented"
        event.Skip()

    def eventTextCtrlWidthL(self, event):  # wxGlade: IesGui.<event_handler>
        """Set the left peak FWHM multiplier for extracting arrival times.
        """
        val = gf.checkIfNumberTextCtrl(self.textCtrlWidthL)
        if type(val).__name__ != 'str':
            self.settings.setWidthL(val)
        else:
            message = 'Please only enter numbers in this box!'
            gf.warningDialog(message)
            self.textCtrlWidthL.SetValue(str(self.settings.widthL))         

    def eventTextCtrlWidthR(self, event):  # wxGlade: IesGui.<event_handler>
        """Set the right peak FWHM multiplier for extracting arrival times.
        """        
        val = gf.checkIfNumberTextCtrl(self.textCtrlWidthR)
        if type(val).__name__ != 'str':
            self.settings.setWidthR(val)
        else:
            message = 'Please only enter numbers in this box!'
            gf.warningDialog(message)
            self.textCtrlWidthR.SetValue(str(self.settings.widthR))

    def setCheckboxEnableDisable(self,MassSpectra=True):
        """Update display setting checkboxes.
        see eventCheckboxMassSpectra
        """
        if MassSpectra:
            others = False
            self.settings.checkboxStates.massSpectrum = True    
        else:
            others = True
            self.settings.checkboxStates.massSpectrum = False
            
        self.checkboxAtds.Enable(others)
        self.checkboxContour.Enable(others)
        self.checkboxConformations.Enable(others)
        self.checkboxPeakMarkers.Enable(others)
        
            
    def eventCheckboxMassSpectra(self, event):  # wxGlade: IesGui.<event_handler>
        """Checkbox to display mass spectra.
        """
        state = event.IsChecked()
        self.setCheckboxEnableDisable(state) #settings state is updated here too
        if state:
            self.settings.checkboxStates.massSpectra = True
        else:
            self.settings.checkboxStates.massSpectra = False
        
        self.plotPanel.refresh_plot()

    def eventCheckboxContour(self, event):  # wxGlade: IesGui.<event_handler>
        """Checkbox to toggle the display of contour plots of m/z against
        arrival time or CCS.
        """
        state = event.IsChecked()
        if state:
            self.settings.checkboxStates.contourPlots = True
        else:
            self.settings.checkboxStates.contourPlots = False
        self.plotPanel.refresh_plot()

    def eventCheckboxAtds(self, event):  # wxGlade: IesGui.<event_handler>
        """Checkbox to toggle the display of arrival time distributions.
        """
        state = event.IsChecked()
        if state:
            self.settings.checkboxStates.atds = True
        else:
            self.settings.checkboxStates.atds = False
        self.plotPanel.refresh_plot()

    def eventCheckboxConformations(self, event):  # wxGlade: IesGui.<event_handler>
        """Checkbox to toggle the display of conformation tracking of
        peaks in the arrival time/CCS distributions.
        """
        state = event.IsChecked()
        if state:
            self.settings.checkboxStates.conformations = True
        else:
            self.settings.checkboxStates.conformations = False
        self.plotPanel.refresh_plot()

    def eventCheckboxPeakMarkers(self, event):  # wxGlade: IesGui.<event_handler>
        """Checkbox to toggle the display of user selected peak positions in
        contour plots and ATDs/CCSDs.
        """
        self.settings.setDisplayPeaks(event.IsChecked())
        self.plotPanel.refresh_plot()
        
    def eventTextCtrlUnits(self, event):  # wxGlade: IesGui.<event_handler>
        """Units to use when displaying annotations for individual spectra.
        e.g. 'eV' for unfolding experiments.
        """
        self.settings.setUnits(self.textCtrlUnits.GetValue())

    def eventWindow_1Unsplit(self, event):  # wxGlade: IesGui.<event_handler>
        # TODO(gns) - unfinished, maybe remove
        event.Veto()

# end of class IesGui
if __name__ == "__main__":
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    frame_1 = IesGui(None, -1, "")
    app.SetTopWindow(frame_1)
    frame_1.Show()
    app.MainLoop()