"""Table for holding the data files, as well as associated values and units,
for IesGui()."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import wx,os,re
from gui.EditableListCtrl import EditableListCtrl
from lib import utils
from collections import OrderedDict
import gui.guiFunctions as gf
import IesSettings,IesPlotPanel


class ListCtrlFiles(wx.Panel):
    
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,-1)
        #self.settings = AtroposGuiSettings.AtroposGuiSettings()
        
        self.listCtrl = EditableListCtrl(self)
        self.listCtrl.setColumnToIgnore(1)
        
        self.listCtrl.InsertColumn(0,"Value",width=50)
        self.listCtrl.InsertColumn(1,"Filename",width=200)
        
        sizer_1 = wx.FlexGridSizer(1, 1, 0, 0) 
        sizer_1.Add(self.listCtrl, 1, wx.EXPAND)  
        sizer_1.AddGrowableRow(0, 1)
        sizer_1.AddGrowableCol(0, 1)
        self.SetSizer(sizer_1)          


        self.Bind(wx.EVT_LIST_END_LABEL_EDIT, self.OnEndLabelEdit, self.listCtrl)

        self.files = OrderedDict() #[basename] = fullpath
        self.unitValues = OrderedDict() #[basename] = value
        self.settings = IesSettings.IesSettings('')
        self.plotPanel = None


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
        
    def OnEndLabelEdit(self,event):
        """Function for editing the value labels for each of the
        data files.
        """
        val = event.GetText()
        message = 'Only enter numerical values here'
        if len(val):
            if utils.isNumber(val):
                self.listCtrl.SetItemText(event.GetIndex(),val)
                self.updateUnitValues()
                self.plotPanel.refresh_plot()
            else:
                gf.warningDialog(message)
                event.Veto()

    def addFiles(self,l):
        """Add data files to the ListCtrl.
        """
        lfixed = [self.fixMultiDirPath(f) for f in l]
        for f in lfixed:
            if gf.checkIfAmphiProject:
                basename = os.path.basename(f)
                self.files[basename] = f
                self._addFileToListCtrl(basename)

        self.settings.setFilenames(lfixed)
        self.plotPanel.plotMassSpectra()
        self.updateUnitValues()

    def updateUnitValues(self):
        """Update the values entered next to the filenames to the
        IesSettings() object.
        """
        self.unitValues = OrderedDict()
        for i,fn in enumerate(self.files.keys()):
            s = self.listCtrl.GetItem(i,0).GetText()
            self.unitValues[fn] = s
        self.settings.setValues(self.unitValues.values())

    def _addFileToListCtrl(self,basename):
        """Subfunction of self.addFiles
        """
        count = self.listCtrl.GetItemCount()
        self.listCtrl.InsertStringItem(count,"")
        self.listCtrl.SetStringItem(count,1,basename)


    def fixMultiDirPath(self,path):
        """Fix problems with the output from MultiDirDialog.
        """
        # TODO(gns) - DANGER THIS NEEDS TO BE FIXED
        # TODO(gns) - DANGER THIS NEEDS TO BE FIXED
        # TODO(gns) - DANGER THIS NEEDS TO BE FIXED
        path = re.sub("Home directory","/home/ganesh",path)
        path = re.sub("Macintosh HD","",path)
        if path[-1] == '/' or path[0] == '\\':
            path = path[:-1]

        return path
    


