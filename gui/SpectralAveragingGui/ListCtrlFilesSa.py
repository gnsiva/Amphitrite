import wx,os,re
from gui.EditableListCtrl import EditableListCtrl
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
        self.settings = settings
    def setPlotPanel(self,plotPanel):
        self.plotPanel = plotPanel
    #====================        
        

    def addFiles(self,l):
        lfixed = [self.fixMultiDirPath(f) for f in l]
        for f in lfixed:
            if gf.checkIfAmphiProject:
                basename = os.path.basename(f)
                self._addFileToListCtrl(basename)

        self.settings.setImFilenames(lfixed)
        self.plotPanel.refresh_plot()

    def addTxtFiles(self,l):
        self.settings.setTxtFilenames(l)
        self.plotPanel.refresh_plot()

        fns = self.settings.massSpectra.keys()
        for s in dir(self.listCtrl):
            print s
        self.listCtrl.DeleteAllItems()
        for fn in fns:
            self._addFileToListCtrl(os.path.basename(fn))

        
    def _addFileToListCtrl(self,basename,column=0):
        '''subfunction of self.addFiles'''
        count = self.listCtrl.GetItemCount()
        self.listCtrl.InsertStringItem(count,"")
        self.listCtrl.SetStringItem(count,column,basename)


    def fixMultiDirPath(self,path):
        path = re.sub("Home directory","/home/ganesh",path)
        path = re.sub("Macintosh HD","",path)
        if path[-1] == '/' or path[0] == '\\':
            path = path[:-1]
        return path



