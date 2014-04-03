"""File for running cppapplication.exe to create Amphitrite
data files from MassLynx raw files.
Deprecated - Use RawFileProcessor_v2.py
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import os
import shutil
import re
import subprocess
import cPickle as pickle
import numpy as np
import time
import utils

class RawFileProcessor():
    
    def __init__(self,rawPath):
        self.path = rawPath
        self.rawfolder = os.path.basename(self.path.rstrip('/'))
        if not self.path.rstrip('/')[-2:] == '.a':
            self.outputfolder = self.path.rstrip('.raw') + '.a'
        else:
            self.outputfolder = self.path
    
    def setOutputFolder(self,outputFolderPath):
        rawFileName = os.path.basename(self.rawfolder)
        rawFileName = rawFileName.rstrip('.raw/') + '.a'

        self.outputfolder = os.path.join(outputFolderPath,rawFileName)
        
    def processFolder(self,grain=2):
        '''1. Copy raw file to working directory
        2. Run CppApplication
        3. Read text files and make imObj
        4. Delete text files, remove raw file
        5. Make new folder for processed data
        6. Dump pickles there'''
        
        if not self._checkIfProcessed():
            # 1
            if not os.path.isdir(os.path.join('.',self.rawfolder)):
                shutil.copytree(self.path,os.path.join('.',self.rawfolder))
            # 2
            #print 'raw folder', self.rawfolder
            #print 'path', self.path      
            print \
'''=================================
Arguments passed
================================='''
            print ['cppapplication.exe',self.rawfolder,"0","1",str(grain),"0"]
            print \
'''================================='''
            p = subprocess.call(['cppapplication.exe',self.rawfolder,"0","1",str(grain),"0"])
            #print p
            #print 'cwd', os.getcwd()
            
            for file in ['MassMobility.txt','MassMobilityXaxis.txt','MassMobilityYaxis.txt']:
                try:
                    os.rename(file, os.path.join(self.rawfolder,file))
                except:
                    print 'waiting for cppapplication'
                    subprocess.call(['cppapplication.exe',str(self.rawfolder),"0",str(grain),"0"])
                    time.sleep(5)
                    try:
                        os.rename(file, os.path.join(self.rawfolder,file))
                    except:
                        print 'still waiting'
                        time.sleep(10)
                        try:
                            os.rename(file, os.path.join(self.rawfolder,file))
                        except:
                            print 'Couldnt open file: %s' %self.rawfolder
                            shutil.rmtree(self.rawfolder)
            if not os.path.isdir(self.outputfolder):
                os.mkdir(self.outputfolder)
            self._processAxisX()
            self._processAxisY()
            self._processMassMobililty()
            shutil.rmtree(path=self.rawfolder)
            
            print 'File processed: %s' %self.rawfolder
    
    def _checkIfProcessed(self):
        processed = False
        amphiFns = ['MassMobilityXaxis.amphi','MassMobilityYaxis.amphi','MassMobility.amphi']
        if os.path.isdir(self.outputfolder):
            if utils.isInDir(self.outputfolder, amphiFns):
                processed = True
        
        # Legacy support for text files (TO BE REMOVED)
        textFns = ['MassMobilityXaxis.txt','MassMobilityYaxis.txt','MassMobility.txt']
        if utils.isInDir(self.path,textFns):
            processed = True
        return processed
        
    def getAxisX(self):
        return self._unPickle(os.path.join(self.outputfolder,'MassMobilityXaxis.amphi'))
    def getAxisY(self):
        return self._unPickle(os.path.join(self.outputfolder,'MassMobilityYaxis.amphi'))
    def getMassMobility(self):
        return self._unPickle(os.path.join(self.outputfolder,'MassMobility.amphi'))


    def _processMassMobililty(self,removeTxt=1):
        path = os.path.join(self.rawfolder, 'MassMobility.txt')
        text = open(path,'r').readlines()
        
        if removeTxt:
            os.remove(path)
        
        lines = len(text)
        file = open('temp.xsg','w')
        for i,line in enumerate(text):
            if i != (lines-1):
                print>> file, line.rstrip('\n')
            else:
                print>> file, line.rstrip(',\n')
        file.close()
        
        ifile = open('temp.xsg','r')
        temp = np.fromfile(ifile,dtype=np.float64,sep=',')
        ifile.close()
        os.remove('temp.xsg')
        
        temp = np.array_split(temp,200)
        massMobility = np.flipud(temp)     
        
        self._pickle(massMobility, os.path.join(self.outputfolder,'MassMobility.amphi'))
        
    
    def _processAxisX(self,removeTxt=1):
        path = os.path.join(self.rawfolder,'MassMobilityXaxis.txt')
        ifile = open(path,'r')
        xAxis = np.fromfile(ifile,dtype='float64',sep=',')
        ifile.close()
        if removeTxt:
            os.remove(path)        
        self._pickle(xAxis[:-2], os.path.join(self.outputfolder,'MassMobilityXaxis.amphi'))
         
    def _processAxisY(self,removeTxt=1):
        path = os.path.join(self.rawfolder,'MassMobilityYaxis.txt')
        ifile = open(path,'r')
        yAxis = np.fromfile(path,sep='\n')
        ifile.close()
        if removeTxt:
            os.remove(path)        
        yAxis = yAxis[::-1]
        self._pickle(yAxis, os.path.join(self.outputfolder,'MassMobilityYaxis.amphi'))    
    

    def _pickle(self,obj,filename):
        obj.dump(filename)
    
    def _unPickle(self,filename):
        ifile = open(os.path.join(filename),'rb')
        obj = pickle.load(ifile)
        ifile.close()
        return obj        
                
    def makePreview(self):
        import fast_driftscope_image as fdi
        image = fdi.Driftscope_image()
        image.load_folder(self.outputfolder)
        image.normalise_mobility()
        image.driftscope()
        imagename = self.rawfolder.rstrip('/')[:-4] + '_preview.png'
        image.savefig(os.path.dirname(self.outputfolder),imagename)
        
        


                
                
                
