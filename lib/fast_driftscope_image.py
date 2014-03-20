import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import cPickle as pickle

class Driftscope_image():
    def __init__(self):
        self.atds = {}
        self.height_normed_atds = {}
        self.area_normed_atds = {}
        self.Z = {}
        self.zs4CCS = {}
    
    ############################################################################
    # For single full amphitrite data files (as opposed to using folders
    ############################################################################
    def load_file(self,filename):
        # load data
        npObj = np.load(filename)
        self.mzs_orig = npObj[0]
        self.tds_orig = npObj[1]
        self.mobility_orig = npObj[2]
        self.its_orig = np.sum(self.mobility_orig,axis=0)
        # process data
        self.limit_axes()
        self.normalise_mobility()
        
    ############################################################################
    #===========================================================================
    # File loading (Using folder of pickled data)
    #===========================================================================
    def load_x(self,filename):
        self.mzs_orig = np.load(filename)
        
    def load_y(self,filename):
        self.tds_orig = np.load(filename)
        
    def load_mobility(self,filename):
        self.mobility_orig = np.load(filename) 
        self.its_orig = np.sum(self.mobility_orig,axis=0)
        self.limit_axes()
    
    def load_folder(self,folderPath):
        self.load_x(os.path.join(folderPath,'MassMobilityXaxis.amphi'))
        self.load_y(os.path.join(folderPath,'MassMobilityYaxis.amphi'))
        self.load_mobility(os.path.join(folderPath,'MassMobility.amphi'))
        self.normalise_mobility()
        
    def limit_axes(self):
        '''still to be implemented, 
        ensures that original data is preserved'''
        self.mzs = self.mzs_orig
        self.tds = self.tds_orig
        self.mobility = self.mobility_orig
        self.its = self.its_orig
    #===========================================================================


    # The rest of the functions
    #===========================================================================
    # Data preparation
    #===========================================================================
    def normalise_mobility(self):
        self.mobility = (self.mobility/self.mobility.max())*1000
        self.mobility = np.log2(self.mobility)    
        infinite = np.isinf(self.mobility)
        self.mobility[infinite] = 0
        self.mobility[self.mobility<1] = 0
        
        
    #===========================================================================
    # THE GRAPH
    #===========================================================================
    def driftscope(self,mz_range=0,a_range=0):
        c_array = ['black','0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.5','r','magenta','purple','b','c','w']
        cmap = mpl.colors.LinearSegmentedColormap.from_list('my_cmap',c_array,256)
        
        # Limit m/z range        
        if mz_range:
            pass
        else:
            mz_range = [self.mzs[0],self.mzs[-1]]
        
        # atd range
        if a_range:
            pass
        else:
            a_range = [self.tds[0],self.tds[-1]]
        
        import matplotlib.gridspec as gridspec
        
        f = plt.figure(figsize=[6,8]) #FIGSIZE
        gs = gridspec.GridSpec(2,1)
        ax1 = plt.subplot(gs[0],xlim=[mz_range[0],mz_range[1]],ylim=[a_range[0],a_range[-1]],ylabel='Arrival Time ($ms$)',xlabel='$m/z$')
        ax1.imshow(self.mobility,aspect='auto',origin='lower',cmap=cmap,extent=[self.mzs[0],self.mzs[-1],self.tds[0],self.tds[-1]])
        
        ax2 = plt.subplot(gs[1], xlabel="$m/z$",xlim=[mz_range[0],mz_range[1]],yticks=[],ylabel='Intensity')
        ax2.plot(self.mzs,self.its,color='black')

        plt.tight_layout()
    
    def savefig(self,outputfolder,rawfilename):
        plt.savefig(os.path.join(outputfolder,rawfilename))
    #===========================================================================
    # The Function
    #===========================================================================
    def draw_figure(self,projectFolder,dataFolder,dataLocation):
        '''Not used by GUI'''
        self.load_folder(folderPath=dataLocation+'/')
        self.driftscope()
        plt.savefig(os.path.join(projectFolder,(dataFolder+'_ATDandMS.png')))


'''
import sys

projectFolder = sys.argv[1]    #full path
dataFolder = sys.argv[2]    #just folder name (no trailing /)
dataLocation = sys.argv[3]     #full path

#pic_folder = '/home/ganesh/Dropbox/PhD/Amphitrite/CppApplication/JUN_120503_HeatingExp/'
#name = '2_ADH_60_3'
#folder = p.p['heating03']+name + '/'

im_ob = IM()
im_ob.load_folder(folderPath=dataLocation+'/')
im_ob.driftscope()
plt.savefig(projectFolder+dataFolder+'_ATDandMS.png')
#plt.savefig(projectFolder+dataFolder+'_ATDandMS.eps')
'''
        
