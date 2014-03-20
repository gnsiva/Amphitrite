import os
import shutil
import re
import subprocess
import cPickle as pickle
import numpy as np
import time
import utils
import gui.guiFunctions as gf

class RawFileProcessor():
    
    def __init__(self,rawFilePath):
        self.matrix = np.array([])
        self.xAxis = np.array([])
        self.yAxis = np.array([])
        
        self.rawFilePath = rawFilePath.rstrip('/')
        #self.path = rawFilePath
        self.rawFileBasename = os.path.basename(self.rawFilePath)
        if self.rawFileBasename[-2:] != '.a':
            self.outputFileBasename = self.rawFileBasename.rstrip('.raw') + '.a'
        else:
            self.outputFileBasename = self.rawFileBasename

        self.outputFile = os.path.join(
                            os.path.dirname(self.rawFilePath),self.outputFileBasename)
    
    def setOutputFolder(self,outputFolderPath):
        self.outputFile = os.path.join(outputFolderPath,self.outputFileBasename)
        
        
    def processFolder(self,grain=2):
        '''1. Copy raw file to working directory
        2. Run CppApplication
        3. Read text files and make imObj
        4. Delete text files, remove raw file
        5. Make new folder for processed data
        6. Dump pickles there'''
        
        if not self._checkIfProcessed():
            # 1
            if not os.path.isdir(os.path.join('.',self.rawFileBasename)):
                shutil.copytree(self.rawFilePath,os.path.join('.',self.rawFileBasename))
            # 2
            print \
'''=================================
Arguments passed
================================='''
            print ['cppapplication.exe',self.rawFileBasename,"0","1",str(grain),"0"]
            print \
'''================================='''
            p = subprocess.call(['cppapplication.exe',self.rawFileBasename,"0","1",str(grain),"0"])
            converted = True
            for file in ['MassMobility.txt','MassMobilityXaxis.txt','MassMobilityYaxis.txt']:
                try:
                    os.rename(file, os.path.join(self.rawFileBasename,file))
                except:
                    print 'waiting for cppapplication'
                    #subprocess.call(['cppapplication.exe',str(self.rawFileBasename),"0",str(grain),"0"])
                    time.sleep(5)
                    try:
                        os.rename(file, os.path.join(self.rawFileBasename,file))
                    except:
                        print 'still waiting'
                        time.sleep(10)
                        try:
                            os.rename(file, os.path.join(self.rawFileBasename,file))
                        except:
                            print 'Couldnt open file: %s' %self.rawFileBasename
                            #shutil.rmtree(self.rawFileBasename)
                            converted = False

            if converted:
                self._processAxisX()
                self._processAxisY()
                self._processMassMobililty()
                self.pickleAmphitriteProject()
            
            shutil.rmtree(path=self.rawFileBasename)
            print 'File processed: %s' %self.rawFileBasename
    
    def _checkIfProcessed(self):
        processed = False
        if self.outputFileBasename in os.listdir(os.path.dirname(self.outputFile)):
            processed = True
        return processed
        
    def _processMassMobililty(self,removeTxt=1):
        path = os.path.join(self.rawFileBasename, 'MassMobility.txt')
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
        
        self.matrix = massMobility
        
    
    def _processAxisX(self,removeTxt=1):
        path = os.path.join(self.rawFileBasename,'MassMobilityXaxis.txt')
        ifile = open(path,'r')
        xAxis = np.fromfile(ifile,dtype='float64',sep=',')
        ifile.close()
        if removeTxt:
            os.remove(path)        

        self.xAxis = xAxis[:-2]
         
    def _processAxisY(self,removeTxt=1):
        path = os.path.join(self.rawFileBasename,'MassMobilityYaxis.txt')
        ifile = open(path,'r')
        yAxis = np.fromfile(path,sep='\n')
        ifile.close()
        if removeTxt:
            os.remove(path)        
        yAxis = yAxis[::-1]

        self.yAxis = yAxis

    def pickleAmphitriteProject(self):
        utils.pickleAmphitriteProject(self.outputFile,
                                   self.xAxis,self.yAxis,self.matrix)
    
    def unPickleAmphitriteProject(self,filename):
        [xAxis,yAxis,matrix] = utils.unPickleAmphitriteProject(filename)
        return [xAxis,yAxis,matrix]
                
    def makePreview(self):
        import fast_driftscope_image as fdi
        image = fdi.Driftscope_image()
        image.load_file(self.outputFile)
        image.normalise_mobility()
        image.driftscope()
        imagename = self.outputFile.rstrip('.a') + '_preview.png'
        image.savefig(os.path.dirname(self.outputFile),imagename)
        
        


                
                
                
