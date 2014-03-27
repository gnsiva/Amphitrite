"""Table for holding conformations used in conformation tracking in IesGui()."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import wx,os,re
from gui.EditableListCtrl import EditableListCtrl
from lib import utils
from collections import OrderedDict
import gui.guiFunctions as gf
import IesSettings,IesPlotPanel

class ListCtrlConformations(wx.Panel):

    def __init__(self,parent):
        wx.Panel.__init__(self,parent,-1)
        
        self.listCtrl = EditableListCtrl(self)
        #self.listCtrl.setColumnToIgnore(0)
        
        self.listCtrl.InsertColumn(0,"CCS",width=100)
        for i in xrange(20):
            self.listCtrl.InsertStringItem(i,'')
        
        sizer_1 = wx.FlexGridSizer(1, 1, 0, 0) 
        sizer_1.Add(self.listCtrl, 1, wx.EXPAND)  
        sizer_1.AddGrowableRow(0, 1)
        sizer_1.AddGrowableCol(0, 1)
        self.SetSizer(sizer_1)          


        self.Bind(wx.EVT_LIST_END_LABEL_EDIT, self.OnEndLabelEdit, self.listCtrl)

        self.settings = IesSettings.IesSettings('')
        self.plotPanel = None

        self.conformations = []

    def setSettings(self,settings):
        """Set the Gui settings class.
        :parameter settings: IesSettings() object
        """
        self.settings = settings
        
    def setPlotPanel(self,plotPanel):
        """Add the Matplotlib plot panel object to this class.
        :parameter plotPanel: IesPlotPanel() object
        """
        self.plotPanel = plotPanel        

    def addConformation(self,ccs):
        """Add a conformation to the ListCtrl.
        :parameter ccs: CCS value of the conformation
        """
        # TODO(gns) - Check for number first and use a warning dialog
        self.conformations.append(float(ccs))
        self.updateListCtrl()

    def removeConformation(self,ccs):
        """Remove conformation from the ListCtrl and conformation tracking.
        :parameter ccs: CCS value of the conformation
        """
        # TODO(gns) - Need to integrate this with the button on the GUI
        if ccs in self.conformations:
            self.conformations.remove(ccs)
        else:
            print 'CCS value not found in conformations'
        self.updateListCtrl()
        
    def OnEndLabelEdit(self,event):
        """Editing CCS value of a selected conformation entry in the ListCtrl.
        """
        newVal = event.GetText()
        i = event.GetIndex()

        if utils.isNumber(newVal):
            newVal = float(newVal)
            self.conformations[i] = newVal
            self.updateListCtrl()
            self.plotPanel.refresh_plot()
        else:
            if not newVal == '':
                message = 'Only enter a numerical value here!'
                gf.warningDialog(message)
                event.Veto()
            else:
                del self.conformations[i]
                self.updateListCtrl()

    def updateListCtrl(self):
        """Refresh the data in the ListCtrl.
        """
        self.conformations = sorted(self.conformations)
        self.settings.conformations = self.conformations
        for i in xrange(self.listCtrl.GetItemCount()):
            self.listCtrl.SetItemText(i,'')  # UNTESTED

        for i in xrange(len(self.conformations)):
            self.listCtrl.InsertStringItem(i,"%.0f" %self.conformations[i])


