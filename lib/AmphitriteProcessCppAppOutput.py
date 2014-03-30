"""Small utility for processing MassLynx raw data files.
(Just produced for my personal use).

Run as:
python AmphitriteProcessCppAppOutput.py <filename>

Must be run in the Amphitrite/lib folder.
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import os,sys
import cPickle as pickle
import numpy as np
import utils
import RawFileProcessor_v2 as RawFileProcessor

filename = sys.argv[1]
#lib = os.path.join('/home','ganesh','Dropbox','workspaces','Amphitrite_2.1','lib')
lib = os.getcwd()

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
