"""Class for creating a ListCtrl (~table) where the user can enter text."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

# TODO(gns) - the atropos version doesn't use ListCtrlAutoWidthMixin
# TODO(gns) - consider removing if that's not important and inherit from the same class.

import wx
from wx.lib.mixins.listctrl import TextEditMixin, ListCtrlAutoWidthMixin

class EditableListCtrl(wx.ListCtrl, TextEditMixin, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT | wx.SUNKEN_BORDER)#, ID, pos, size, style)
        TextEditMixin.__init__(self) 
        ListCtrlAutoWidthMixin.__init__(self)
        
        self.Bind(wx.EVT_LIST_BEGIN_LABEL_EDIT, self.OnBeginLabelEdit)
        self.columnToIgnore = None
        
    
    def OnBeginLabelEdit(self, event):
        """Check if this cell is editable otherwise,
        do not take the text input.
        """
        if event.m_col == self.columnToIgnore:
            event.Veto()

  
    def setColumnToIgnore(self,col):
        """Do not process events from given column index.
        (Mainly used to make cells uneditable)

        :parameter col: column index
        """
        self.columnToIgnore = col

    def getColumnStrings(self,col):
        """Get the entries in each data cell for the given column.

        :parameter col: column index
        :returns: list of strings with one entry per data cell
        """
        rows = self.GetItemCount()
        strings = []
        for i in xrange(rows):
            strings.append(self.GetItem(i,col).GetText())
        
        return strings
