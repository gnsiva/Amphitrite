import os,sys
import cPickle as pickle
import numpy as np
import utils
import RawFileProcessor_v2 as RawFileProcessor

filename = sys.argv[1]
#lib = os.path.join('E:/','workspaces','Amphitrite_2.1','lib')
lib = os.path.join('/home','ganesh','Dropbox','workspaces','Amphitrite_2.1','lib')

rfp = RawFileProcessor.RawFileProcessor(os.path.join(lib,filename))
rfp.rawFileBasename = lib

# move to the directory of the rawfile
os.chdir(lib)
print os.getcwd()

rfp._processAxisX(removeTxt=1)
rfp._processAxisY(removeTxt=1)
rfp._processMassMobililty(removeTxt=1)

ofile = filename[:-4] + '.a'
rfp.outputFile = os.path.join(lib,ofile)#ofile)
rfp.pickleAmphitriteProject()

rfp.makePreview()
