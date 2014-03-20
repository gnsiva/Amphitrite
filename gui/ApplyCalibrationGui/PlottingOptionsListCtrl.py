import wx, ApplyCalibrationGuiSettings
from EditableListCtrl import EditableListCtrl
from lib import utils

class PlottingOptionsListCtrl(wx.Panel):
    
    def __init__(self,parent):
        wx.Panel.__init__(self,parent)
        
        self.listCtrl = EditableListCtrl(self)
        self.listCtrl.setColumnToIgnore(0)
        
        self.listCtrl.InsertColumn(0,"Option",width=80)
        self.listCtrl.InsertColumn(1,"Value",width=80)
        
        sizer_1 = wx.FlexGridSizer(1, 1, 0, 0) 
        sizer_1.Add(self.listCtrl, 1, wx.EXPAND)  
        sizer_1.AddGrowableRow(0, 1)
        sizer_1.AddGrowableCol(0, 1)
        self.SetSizer(sizer_1)          
        
        self.listCtrl.InsertStringItem(0, "" )
        self.listCtrl.SetStringItem(0, 0, "Smoothing")
        self.listCtrl.SetStringItem(0, 1, "False")
        
        self.listCtrl.InsertStringItem(1, "" )
        self.listCtrl.SetStringItem(1,0, "Line Width")      
        self.listCtrl.SetStringItem(1, 1, "1")
              
        self.listCtrl.InsertStringItem(2, "" )
        self.listCtrl.SetStringItem(2,0, "Auto Axis")   
        self.listCtrl.SetStringItem(2, 1, "False")
            
        self.listCtrl.InsertStringItem(3, "" )
        self.listCtrl.SetStringItem(3,0, "Colour")   
        self.listCtrl.SetStringItem(3, 1, "Black")
        
        self.settings = ApplyCalibrationGuiSettings.ApplyCalibrationGuiSettings()
    
    def setSettings(self,settings):
        """Set the Gui Settings object.
        :parameter settings: ApplyCalibrationGuiSettings() object
        """
        self.settings = settings
    
    def getSettingValue(self,setting):
        """Get one of the plotting setting values.
        :parameter setting: String name of attribute ('Smoothing', 'Line Width',
        'Auto Axis', 'Colour')
        """
        row = None
        if setting == 'Smoothing':
            row = 0
        elif setting == 'Line Width':
            row = 1
        elif setting == 'Auto Axis':
            row = 2
        elif setting == 'Colour':
            row = 3
        else:
            print 'Unrecognised setting: %s (PlottingOptionsListCtrl.getSettingValue)' %setting
            return None

        return self.listCtrl.GetItem(row,1).GetText()


    def checkFieldsValid(self):
        """Check if each of the plotting options contains a valid input.
        :returns: valid,errors - Boolean, list of the names of values with errors
        """
        valid = True
        errors = []
        smoothing = self.getSettingValue('Smoothing')
        if not utils.isBinaryResponse(smoothing):
            valid = False
            errors.append('Smoothing')
        aaxis = self.getSettingValue('Auto Axis')
        if not utils.isBinaryResponse(aaxis):
            valid = False
            errors.append('Auto Axis') 
        colour = self.getSettingValue('Colour')
        if not utils.isMplColour(colour):
            valid = False
            errors.append('Colour')
        lw = self.getSettingValue('Line Width')
        if not utils.isNumber(lw):
            valid = False
            errors.append('Line Width')
        return valid,errors
    
    def processSmoothing(self):
        """Carry out the smoothing on the arrival time distribution.
        """
        smoothing = self.getSettingValue('Smoothing')
        smoothing = utils.getBinaryReponse(smoothing)
        if smoothing:
            self.settings.atdPanel.atd.smoothingSG()
        else:
            self.settings.atdPanel.atd.restoreRawYvals()
        
    def processLineWidth(self):
        """Get the float value for line width from ListCtrl.
        :returns: Line width
        """
        lw = self.getSettingValue('Line Width')
        return float(lw)
    def processAutoAxis(self):
        """Update the plot to auto axis on or off depending on the data
        in the ListCtrl.
        """
        autoaxis = utils.getBinaryReponse(self.getSettingValue('Auto Axis'))
        if autoaxis:
            self.settings.atdPanel.autoAxisX()
        else:
            self.settings.atdPanel.autoAxisXoff()
    def processColour(self):
        """Get the value for colour from the ListCtrl.
        :returns: Colour
        """
        return self.getSettingValue('Colour').lower()
        
    def updateSettings(self):
        """Check if the values in the ListCtrl are valid, then update
        the plot window.
        """
        valid,errors = self.checkFieldsValid()
        if not valid:
            self.ShowMessage(errors)
        else:
            self.processSmoothing()
            lw = self.processLineWidth()
            color = self.processColour()
            self.settings.atdPanel.clear()
            if self.settings.displayedView == 'td':
                self.settings.atdPanel.plotAtd(color=color,lw=lw)
            else:
                self.settings.atdPanel.plotCcs(color=color,lw=lw)
            autoaxis = self.processAutoAxis()
    
    def ShowMessage(self,errors):
        """Create a pop window (MessageBox) displaying errors.
        :parameter errors: Errors as to be shown in window
        """
        wx.MessageBox('Fails:\n%s' %str(errors)[1:-1], 'Info', 
        wx.OK | wx.ICON_INFORMATION)


     
        
        
        
