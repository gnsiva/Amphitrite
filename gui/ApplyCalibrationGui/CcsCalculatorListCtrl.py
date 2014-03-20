"""ListCtrl (table) for the conversion of individual arrival time
values to CCS.
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import wx, ApplyCalibrationGuiSettings
from EditableListCtrl import EditableListCtrl
    
class CcsCalculatorListCtrl(wx.Panel):

    """
    :parameter parent: The wx area to place the ListCtrl
    """
    def __init__(self,parent):
        wx.Panel.__init__(self,parent)
        
        self.listCtrl = EditableListCtrl(self)
        self.listCtrl.setColumnToIgnore(3)
        
        self.listCtrl.InsertColumn(0,"td",width=40)
        self.listCtrl.InsertColumn(1,"m/z",width=60)
        self.listCtrl.InsertColumn(2,"z",width=30)
        self.listCtrl.InsertColumn(3,"CCS")
        
        sizer_1 = wx.FlexGridSizer(1, 1, 0, 0) 
        sizer_1.Add(self.listCtrl, 1, wx.EXPAND)  
        sizer_1.AddGrowableRow(0, 1)
        sizer_1.AddGrowableCol(0, 1)
        self.SetSizer(sizer_1)        
        
        for i in xrange(10):
            self.listCtrl.InsertStringItem(i,"")
    
        self.settings = ApplyCalibrationGuiSettings.ApplyCalibrationGuiSettings()
    
    def setSettings(self,settings):
        """Set the Gui Settings object.
        :parameter settings: ApplyCalibrationGuiSettings() object
        """
        self.settings = settings
    
    def calculate(self):
        """Calculate the CCS values and update the table.
        """
        rows = self.listCtrl.GetItemCount()
        
        for row in xrange(rows):
            td = self.isNumber(self.listCtrl.GetItem(row, 0).GetText())
            mz = self.isNumber(self.listCtrl.GetItem(row, 1).GetText())
            z = self.isNumber(self.listCtrl.GetItem(row, 2).GetText())
            if not 0 in [td,mz,z]:
                ccsValue = self.calcCcs(td, mz, z)
                ccs = self.listCtrl.GetItem(row, 3)
                ccs.SetText(str(ccsValue))
                self.listCtrl.SetItem(ccs)
            
            
    def isNumber(self,s):
        """Check if value is a number.
        :parameter s: String value for proposed number
        :returns: Boolean False if NaN, or float(s)
        """
        # TODO(gns) - This should be handled by a lib.utils function
        try:
            return float(s)
        except:
            return 0
    
    def calcCcs(self,td,mz,z):
        """Calculate a CCS value.
        :parameter td: Arrival time value (ms)
        :parameter mz: Mass to charge ratio
        :parameter z: Charge state
        :returns: ccs (float)
        """
        ccs = self.settings.calibration.apply1dCalibration(mz, td, z)
        return ccs












