"""Super class for creating a wx panel containing a
Matplotlib plotting area.
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.pyplot as plt
import numpy as np
import wx

matplotlib.rc('xtick', labelsize=8) 
matplotlib.rc('ytick', labelsize=8) 

class PlotPanel():
    
    def __init__(self,panel):
        self.dpi = 80
        self.fig = Figure((3.0, 2.0), dpi=self.dpi)
        self.canvas = FigCanvas(panel, -1, self.fig)
        self.ax = self.fig.add_subplot(111)       

        self.toolbar = NavigationToolbar(self.canvas)  
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)
        panel.SetSizer(self.vbox)
        self.canvas.draw() 
        
        self.ax.plot([1,2,3],[1,2,3])


    def draw(self):
        """Update the figure.
        """
        self.canvas.draw()
