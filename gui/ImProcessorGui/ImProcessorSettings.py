"""Class to hold functions and attributes related to ImProcessorGui()."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import os, wx
#from lib import RawFileProcessor
import lib.RawFileProcessor_v2 as RawFileProcessor
from collections import OrderedDict

import gui.guiFunctions as gf

class Settings():
    
    def __init__(self):
        self.inDirPath = ''
        #self.outDirPath = ''
        self.rawFilePaths = OrderedDict()
        self.gui = None

    def setGui(self,yourself):
        """Set the Gui object so the attributes can be accessed from here.
        :parameter yourself: self from the ImProcessorGui() object
        """
        self.gui = yourself
        
    def setInDirPath(self,path):
        """Set the input directory, containing the Synapt data files.
        :parameter path: Absolute path
        """
        # if self.inDirPath == self.outDirPath:
        #     self.outDirPath = path
        self.inDirPath = path

    # def setOutDirPath(self,path):
    #     self.outDirPath = path

    def getOutDirPath(self):
        """Get output directory path. If alternate directory checkbox is selected
        get the value from the output directory textCtrl, otherwise return the input
        directory.
        :returns: Absolute output path
        """
        if self.gui.checkboxAlternateDirectory.IsChecked():
            dirPath = self.gui.textCtrlAlternateDirectory.GetValue()
            if os.path.isdir(dirPath):
                return dirPath
            else:
                message = \
'''
Output directory not found:
%s 
Using default directory:
%s
''' %(dirPath,self.inDirPath)
                gf.warningDialog(message)
                return self.inDirPath
        else:
            return self.inDirPath
                
    def listRawFiles(self,dirPath):
        """Get list of raw files in a directory, also store them in self.rawFilePaths
        dictionary.
        :parameter dirPath: Absolute path to directory containing raw files
        :returns: rawFiles - list of base names e.g. ['bsa.raw','cona.raw']
        """
        files = os.listdir(dirPath)
        rawFiles = []
        self.rawFilePaths = OrderedDict()
        
        for fn in files:
            if str(fn)[-4:] == '.raw':
                rawFiles.append(str(fn))
                self.rawFilePaths[str(fn)] = os.path.join(dirPath,str(fn))
        
        self.rawFilePaths = OrderedDict(sorted(self.rawFilePaths.items(), key=lambda x: x[0],reverse=True))
        
        return rawFiles
    
    def extract_rawfiles(self,filePaths,grain,yourself):
        """Run the conversion process to convert Synapt data files to
        Amphitrite files.
        :parameter filePaths: List of absolute paths to each Synapt file
        :parameter grain: m/z binning width
        :parameter yourself: self of ImProcessorGui() object
        """
        outDirPath = self.getOutDirPath()
        for path in filePaths:
            # please wait dialog
            bi = wx.BusyInfo("Processing %s" %path, yourself)
            wx.Yield() 
            
            # Processing
            fullPath = str(self.rawFilePaths[path])
            rp = RawFileProcessor.RawFileProcessor(fullPath)
            rp.setOutputFolder(outDirPath)
            rp.processFolder(grain)
            #rp.makePreview()
              
            # close please wait dialog
            bi.Destroy()
