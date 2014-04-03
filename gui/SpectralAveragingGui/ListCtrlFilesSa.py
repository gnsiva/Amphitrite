"""Table containing files to be used in SpectralAveragingGui()."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import wx,os,re
from EditableListCtrl import EditableListCtrl
from lib import utils
from collections import OrderedDict
import gui.guiFunctions as gf
import SaSettings,SaPlotPanel


class ListCtrlFiles(wx.Panel):
    
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,-1)
        self.listCtrl = EditableListCtrl(self)
        self.listCtrl.setColumnToIgnore(0)
        
        self.listCtrl.InsertColumn(1,"Filename",width=300)
        
        sizer_1 = wx.FlexGridSizer(1, 1, 0, 0) 
        sizer_1.Add(self.listCtrl, 1, wx.EXPAND)  
        sizer_1.AddGrowableRow(0, 1)
        sizer_1.AddGrowableCol(0, 1)
        self.SetSizer(sizer_1)          


        #self.Bind(wx.EVT_LIST_END_LABEL_EDIT, self.OnEndLabelEdit, self.listCtrl)

        self.files = OrderedDict() #[basename] = fullpath
        self.unitValues = OrderedDict() #[basename] = value
        self.settings = None
        self.plotPanel = None

    #====================
    # Adding other panel components
    def setSettings(self,settings):
        """Set the Gui settings object.
        :parameter settings: SaSettings() object
        """
        self.settings = settings
        
    def setPlotPanel(self,plotPanel):
        """Set the plotting area object.
        :parameter plotPanel: SaPlotPanel() object
        """
        self.plotPanel = plotPanel
    #====================        
        

    def addFiles(self,l):
        """Add Amphitrite data files ('.a') for spectral averaging.
        :parameter l: List of paths
        """
        lfixed = [self.fixMultiDirPath(f) for f in l]
        for f in lfixed:
            if gf.checkIfAmphiProject:
                basename = os.path.basename(f)
                self._addFileToListCtrl(basename)

        self.settings.setImFilenames(lfixed)
        self.plotPanel.refresh_plot()

    def addTxtFiles(self,l):
        """Add spectrum list text files. 
        :parameter l: List of paths
        """
        self.settings.setTxtFilenames(l)
        self.plotPanel.refresh_plot()

        fns = self.settings.massSpectra.keys()
        for s in dir(self.listCtrl):
            print s
        self.listCtrl.DeleteAllItems()
        for fn in fns:
            self._addFileToListCtrl(os.path.basename(fn))

        
    def _addFileToListCtrl(self,basename,column=0):
        """Subfunction of self.addFiles.
        :parameter basename: Filename without path
        :parameter column: Column to add the file to
        """
        # TODO(gns) - pretty sure the column will always be 0, should remove the option and fix the comment.
        count = self.listCtrl.GetItemCount()
        self.listCtrl.InsertStringItem(count,"")
        self.listCtrl.SetStringItem(count,column,basename)


    def fixMultiDirPath(self,path):
        """Soon to be removed...
        """
        # TODO(gns) - Get rid of this and switch to using Amphitrite pickled data.
        # TODO(gns) - DANGER THIS NEEDS TO BE FIXED
        # TODO(gns) - DANGER THIS NEEDS TO BE FIXED
        # TODO(gns) - DANGER THIS NEEDS TO BE FIXED        
        path = re.sub("Home directory","/home/ganesh",path)
        path = re.sub("Macintosh HD","",path)
        if path[-1] == '/' or path[0] == '\\':
            path = path[:-1]
        return path



