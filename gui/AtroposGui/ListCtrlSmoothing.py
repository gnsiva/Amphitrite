"""ListCtrl for holding the information about smoothing and the lower
limit of base peak intensity to be included in first derivative peak finding.
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import wx
from EditableListCtrl import EditableListCtrl
from lib import utils

class ListCtrlSmoothing(wx.Panel):
    
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,-1)
        
        self.listCtrl = EditableListCtrl(self)
        self.listCtrl.setColumnToIgnore(0)
        
        self.listCtrl.InsertColumn(0,"Option",width=100)
        self.listCtrl.InsertColumn(1,"Value",width=80)
        
        sizer_1 = wx.FlexGridSizer(1, 1, 0, 0) 
        sizer_1.Add(self.listCtrl, 1, wx.EXPAND)  
        sizer_1.AddGrowableRow(0, 1)
        sizer_1.AddGrowableCol(0, 1)
        self.SetSizer(sizer_1)          
        
        self.listCtrl.InsertStringItem(0, "" )
        self.listCtrl.SetStringItem(0, 0, "Smoothing")
        self.listCtrl.SetStringItem(0, 1, "True")
        
        self.listCtrl.InsertStringItem(1, "" )
        self.listCtrl.SetStringItem(1, 0, "Window length")
        self.listCtrl.SetStringItem(1, 1, "5")

        self.listCtrl.InsertStringItem(2, "" )
        self.listCtrl.SetStringItem(2, 0, "Smoothes")
        self.listCtrl.SetStringItem(2, 1, "2")

        self.listCtrl.InsertStringItem(3, "" )
        self.listCtrl.SetStringItem(3, 0, "Polynomial")
        self.listCtrl.SetStringItem(3, 1, "1")  
        
        self.listCtrl.InsertStringItem(4, "" )
        self.listCtrl.SetStringItem(4, 0, "Limit")
        self.listCtrl.SetStringItem(4, 1, "1") 

    def setSettings(self,settings):
        """Set the AtroposGuiSettings object, so its functions and
        attributes can be accessed more readily.

        :parameter settings: AtroposGuiSettings() object
        """
        self.settings = settings
        
    def getValues(self):
        """Get the smoothing parameters (not limit).
        Checks for valid input.
        """
        # TODO(gns) - make alert boxes instead of printing to terminal
        strings = self.listCtrl.getColumnStrings(1)
        
        smoothing = utils.getBinaryReponse(strings[0])
        if smoothing == 'Error':
            smoothing = True
            print 'Smoothing entry is not a binary value (using default - True)'
        
        if utils.isNumber(strings[1]):
            wlen = float(strings[1])
        else:
            wlen = 5
            print 'Wlen entry is not a number (using default - %d)' %wlen
            
        if utils.isNumber(strings[2]):
            smoothes = int(strings[2])
        else:
            smoothes = 2
            print 'Wlen entry is not a number (using default - %d)' %smoothes
    
        if utils.isNumber(strings[3]):
            poly = float(strings[3])
        else:
            poly = 1
            print 'Polynomial entry is not a number (using default - %d)' %poly        

        
        return smoothing,wlen,smoothes,poly
        

