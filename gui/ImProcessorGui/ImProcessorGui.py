#!/usr/bin/env python
# -*- coding: utf-8 -*-
# generated by wxGlade 0.6.4 on Fri Apr  5 14:39:32 2013

import wx
import ImProcessorSettings
import gui.guiFunctions as gf
import os
# begin wxGlade: extracode
# end wxGlade

"""Program for converting Waters Synapt ion mobility data files into
Amphitrite compatible format."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

class ImProcessorGui(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: ImProcessorGui.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.textCtrlInputDirectory = wx.TextCtrl(self, -1, "")
        self.buttonProjectDirectory = wx.Button(self, -1, "Open")
        self.sizer_6_staticbox = wx.StaticBox(self, -1, "Project Directory")
        self.listCtrlFiles = wx.ListCtrl(self, -1, style=wx.LC_REPORT | wx.SUNKEN_BORDER)
        self.sizer_5_staticbox = wx.StaticBox(self, -1, "Files")
        self.label_1 = wx.StaticText(self, -1, "Grain (m/z)")
        self.textCtrlGrain = wx.TextCtrl(self, -1, "1")
        self.sizer_2_staticbox = wx.StaticBox(self, -1, "Options")
        self.checkboxAlternateDirectory = wx.CheckBox(self, -1, "Alternate output directory")
        self.textCtrlAlternateDirectory = wx.TextCtrl(self, -1, "")
        self.buttonAlternateDirectory = wx.Button(self, -1, "Open")
        self.sizer_3_staticbox = wx.StaticBox(self, -1, "")
        self.buttonProcess = wx.Button(self, -1, "Process")

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_TEXT, self.eventTextCtrlInputDirectory, self.textCtrlInputDirectory)
        self.Bind(wx.EVT_BUTTON, self.eventButtonProjectDirectory, self.buttonProjectDirectory)
        self.Bind(wx.EVT_TEXT, self.eventTextCtrlGrain, self.textCtrlGrain)
        self.Bind(wx.EVT_CHECKBOX, self.eventToggleAlternateDirectory, self.checkboxAlternateDirectory)
        self.Bind(wx.EVT_TEXT, self.eventTextCtrlAlternateDirectory, self.textCtrlAlternateDirectory)
        self.Bind(wx.EVT_BUTTON, self.eventButtonAlternateDirectory, self.buttonAlternateDirectory)
        self.Bind(wx.EVT_BUTTON, self.eventProcess, self.buttonProcess)
        # end wxGlade
        self.settings = ImProcessorSettings.Settings()
        self.settings.setGui(self)

    def __set_properties(self):
        # begin wxGlade: ImProcessorGui.__set_properties
        self.SetTitle("ImProcessor")
        self.SetSize((350, 500))
        self.textCtrlAlternateDirectory.Enable(False)
        self.buttonAlternateDirectory.Enable(False)
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: ImProcessorGui.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_1 = wx.FlexGridSizer(1, 1, 0, 0)
        grid_sizer_2 = wx.FlexGridSizer(3, 1, 0, 0)
        self.sizer_3_staticbox.Lower()
        sizer_3 = wx.StaticBoxSizer(self.sizer_3_staticbox, wx.HORIZONTAL)
        sizer_4 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_3 = wx.FlexGridSizer(1, 2, 0, 0)
        self.sizer_2_staticbox.Lower()
        sizer_2 = wx.StaticBoxSizer(self.sizer_2_staticbox, wx.HORIZONTAL)
        grid_sizer_4 = wx.FlexGridSizer(1, 2, 0, 0)
        self.sizer_5_staticbox.Lower()
        sizer_5 = wx.StaticBoxSizer(self.sizer_5_staticbox, wx.HORIZONTAL)
        self.sizer_6_staticbox.Lower()
        sizer_6 = wx.StaticBoxSizer(self.sizer_6_staticbox, wx.HORIZONTAL)
        grid_sizer_5 = wx.FlexGridSizer(1, 2, 0, 0)
        grid_sizer_5.Add(self.textCtrlInputDirectory, 0, wx.EXPAND, 0)
        grid_sizer_5.Add(self.buttonProjectDirectory, 0, 0, 0)
        grid_sizer_5.AddGrowableCol(0)
        sizer_6.Add(grid_sizer_5, 1, wx.EXPAND, 0)
        grid_sizer_1.Add(sizer_6, 1, wx.EXPAND, 0)
        sizer_5.Add(self.listCtrlFiles, 1, wx.EXPAND, 0)
        grid_sizer_2.Add(sizer_5, 1, wx.EXPAND, 0)
        grid_sizer_4.Add(self.label_1, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_4.Add(self.textCtrlGrain, 0, wx.EXPAND, 0)
        grid_sizer_4.AddGrowableCol(0)
        grid_sizer_4.AddGrowableCol(1)
        sizer_2.Add(grid_sizer_4, 1, wx.EXPAND, 0)
        grid_sizer_2.Add(sizer_2, 1, wx.EXPAND, 0)
        sizer_4.Add(self.checkboxAlternateDirectory, 0, 0, 0)
        grid_sizer_3.Add(self.textCtrlAlternateDirectory, 0, wx.EXPAND, 0)
        grid_sizer_3.Add(self.buttonAlternateDirectory, 0, 0, 0)
        grid_sizer_3.AddGrowableCol(0)
        sizer_4.Add(grid_sizer_3, 1, wx.EXPAND, 0)
        sizer_3.Add(sizer_4, 1, wx.EXPAND, 0)
        grid_sizer_2.Add(sizer_3, 1, wx.EXPAND, 0)
        grid_sizer_2.Add(self.buttonProcess, 0, wx.EXPAND, 0)
        grid_sizer_2.AddGrowableRow(0)
        grid_sizer_2.AddGrowableCol(0)
        grid_sizer_1.Add(grid_sizer_2, 1, wx.EXPAND, 0)
        grid_sizer_1.AddGrowableRow(1)
        grid_sizer_1.AddGrowableCol(0)
        sizer_1.Add(grid_sizer_1, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        self.Layout()
        # end wxGlade

        self.listCtrlFiles.InsertColumn(0,'Filename',width=300)
        self.grainValue = 1
        
        
    def eventButtonProjectDirectory(self, event):  # wxGlade: ImProcessorGui.<event_handler>
        """Open directory dialog for selecting the directory where the files
        are.
        """
        dlg = wx.DirDialog(self, message="Open Data Directory",
                            style=wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            if dlg.GetPath() != '':
                path = self._checkIfRawFileReturnParent(dlg.GetPath())
                self.settings.setInDirPath(path)
                self.textCtrlAlternateDirectory.SetValue(dlg.GetPath())
                self.textCtrlInputDirectory.SetValue(self.settings.inDirPath)
                # set list ctrl up
                rawfiles = self.settings.listRawFiles(self.settings.inDirPath)
        dlg.Destroy()

    def eventTextCtrlGrain(self, event):  # wxGlade: ImProcessorGui.<event_handler>
        """Set the grain of the m/z binning when converting the data files.
        """
        val = event.GetString()
        if val != '':
            try:
                val = float(val)
                if val % 1 == 0:
                    val = int(val)
                self.grainValue = val
            except:
                message = 'Please only enter numbers!'
                gf.warningDialog(message)
                self.textCtrlGrain.SetValue(str(self.grainValue))
        
    def _checkIfRawFileReturnParent(self,path):
        """Check directory to see if it ends in .raw (is a Synapt data file).
        :parameter path: Absolute path to raw file directory
        """
        if path.rstrip('/')[-4:] == '.raw':
            path = os.path.dirname(path)
        return path

    def eventButtonAlternateDirectory(self, event):  # wxGlade: ImProcessorGui.<event_handler>
        """Specify a different directory to put converted Amphitrite
        data files. 
        """
        dlg = wx.DirDialog(self, message="Open Data Output Directory",
                            style=wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            #self.settings.setOutDirPath(dlg.GetPath())
            self.textCtrlAlternateDirectory.SetValue(dlg.GetPath())
        dlg.Destroy()
        
                  
    def eventToggleAlternateDirectory(self, event):  # wxGlade: ImProcessorGui.<event_handler>
        """Under normal operation the Amphitrite data files are stored in
        the same directory as the Synapt files. This toggle indicates that an
        alternate directory should be used.
        """
        val = self.checkboxAlternateDirectory.GetValue()
        if val == True:
            self.textCtrlAlternateDirectory.Enable()
            self.buttonAlternateDirectory.Enable()
        elif val == False:
            self.textCtrlAlternateDirectory.Enable(False)
            self.buttonAlternateDirectory.Enable(False)
            
                    
    def eventTextCtrlInputDirectory(self, event):  # wxGlade: ImProcessorGui.<event_handler>
        """Processing input directory when the path is entered into the
        textCtrl instead of using the DirDialog.
        """
        path = self.textCtrlInputDirectory.GetValue()
        if os.path.isdir(path):
            rawfiles = self.settings.listRawFiles(path)
            self._updateListCtrl()
            self.settings.setInDirPath(path)
        
        
    def eventTextCtrlAlternateDirectory(self, event):  # wxGlade: ImProcessorGui.<event_handler>
        """This function is called when the path is entered into the alternate
        directory textCtrl. Currently does nothing though.
        """
        # TODO(gns) - Remove this function
        #self.settings.setOutDirPath(self.textCtrlAlternateDirectory.GetValue())
        pass

    def _updateListCtrl(self):
        """Delete all the items in the listCtrl and repopulate.
        """
        self.listCtrlFiles.DeleteAllItems()
        for  k,v in self.settings.rawFilePaths.items():
            self.listCtrlFiles.InsertStringItem(index=0,label=k)
    
    def _getSelectedItemsListCtrl(self):
        """Get the data files the user has selected to be converted.
        :returns: List of the selected filenames (not full paths).
        """
        selection = []
        if self.listCtrlFiles.GetSelectedItemCount():
            index = self.listCtrlFiles.GetFirstSelected()
            selection.append(self.listCtrlFiles.GetItemText(index))
            while len(selection) != self.listCtrlFiles.GetSelectedItemCount():
                index = self.listCtrlFiles.GetNextSelected(index)
                selection.append(self.listCtrlFiles.GetItemText(index))
            selection = [str(s) for s in selection]
        return selection

    def eventProcess(self, event):  # wxGlade: ImProcessorGui.<event_handler>
        """Convert the selected Synapt data files into Amphitrite
        files. (Process button event.)
        """
        # TODO(gns) - important sounding comment below here, check what its about
        ## NEED TO CHECK GRAIN AND OTHER TEXT BOXES
        selected = self._getSelectedItemsListCtrl()
        if not len(selected):
            message = 'Please select one or more files and try again'
            gf.warningDialog(message)
        else:
            grain = self._processGrain()
            self.settings.extract_rawfiles(selected,grain,self)
        
    
    def _processGrain(self):
        """Check that the value for m/z binning widths (grain) is valid.
        :returns: grain (float)
        """
        grain = self.textCtrlGrain.GetValue()
        try:
            grain = float(grain)
            if grain < 0.5:
                # warning message
                message = 'Grain spacing too low: %s\n \
Using minimum value (0.5 m/z units)' %grain
                dial = wx.MessageDialog(None, message, 'Warning',
                                        wx.OK | wx.ICON_EXCLAMATION)
                dial.ShowModal()                
                # set grain value
                grain = 0.5
        except:
            # warning message
            message = 'Non-number used for grain: %s\n \
Using default value (2.0 m/z units)' %grain
            gf.warningDialog(message)
            # set grain value
            grain = 2.0
        return grain



# end of class ImProcessorGui
if __name__ == "__main__":
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    frame_1 = ImProcessorGui(None, -1, "")
    app.SetTopWindow(frame_1)
    frame_1.Show()
    app.MainLoop()
