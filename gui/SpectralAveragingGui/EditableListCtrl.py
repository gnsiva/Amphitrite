"""Super class for creating an list control (table), where
the cells can be edited by the user.
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import wx
from wx.lib.mixins.listctrl import TextEditMixin, ListCtrlAutoWidthMixin

#class EditableListCtrl(wx.ListCtrl, TextEditMixin, ListCtrlAutoWidthMixin):
class EditableListCtrl(wx.ListCtrl, TextEditMixin):
    def __init__(self, parent):
        
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT | wx.SUNKEN_BORDER)#, ID, pos, size, style)
        TextEditMixin.__init__(self) 
        #ListCtrlAutoWidthMixin.__init__(self)
        
        self.Bind(wx.EVT_LIST_BEGIN_LABEL_EDIT, self.OnBeginLabelEdit)
        self.columnToIgnore = None
        
    
    def OnBeginLabelEdit(self, event):
        """Test whether this cell should be editable, and
        ignore the input if it shouldn't.
        """
        if event.m_col == self.columnToIgnore:
            event.Veto()

  
    def setColumnToIgnore(self,col):
        """Set which column should not be editable.
        :parameter col: Column index (int)
        """
        self.columnToIgnore = col

    def getColumnStrings(self,col):
        """Get a list of strings containing the data in each
        cell in a column.
        :parameter col: Column index (int)
        :returns: List of strings
        """
        rows = self.GetItemCount()
        strings = []
        for i in xrange(rows):
            strings.append(self.GetItem(i,col).GetText())
        
        return strings
