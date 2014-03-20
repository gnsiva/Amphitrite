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
        self.settings = settings
    def setPlotPanel(self,plotPanel):
        self.plotPanel = plotPanel        

    def addConformation(self,ccs):
        '''check if it is a number beforehand'''
        self.conformations.append(float(ccs))
        self.updateListCtrl()

    def removeConformation(self,ccs):
        if ccs in self.conformations:
            self.conformations.remove(ccs)
        else:
            print 'CCS value not found in conformations'
        self.updateListCtrl()
        
    def OnEndLabelEdit(self,event):
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
        self.conformations = sorted(self.conformations)
        self.settings.conformations = self.conformations
        for i in xrange(self.listCtrl.GetItemCount()):
            self.listCtrl.SetItemText(i,'')  # UNTESTED

        for i in xrange(len(self.conformations)):
            self.listCtrl.InsertStringItem(i,"%.0f" %self.conformations[i])


