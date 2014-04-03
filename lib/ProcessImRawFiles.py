"""Program for processing raw files. Still uses old RawFileProcessor."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

import RawFileProcessor
import paths as p
import os,sys

if len(sys.argv) == 1:
    print 'Enter parent directory as command line argument'
    quit()

rootFolder = sys.argv[1]
if not os.path.isdir(rootFolder):
    print 'Directory not found:', rootFolder

#print rootFolder
#quit()

# get raw file names
folders = []
for (path,dirs,files) in os.walk(rootFolder):
    folders = dirs
    print 'Folders to process:'
    for folder in folders:
        print folder
    break


# iterate over raw files
for i,folder in enumerate(folders):
    r = RawFileProcessor.RawFileProcessor(os.path.join(rootFolder,folder))
    r.processFolder()
    print '''############################
#
# %s done
#
############################
''' %(folder)

print 'done'

'''
Paste Bin:
Tekmor:
F:
cd F:\workspaces\Amphitrite_2.0\lib

Phoebe:
cd C:\Users\ganesh\Dropbox\workspaces\Amphitrite_2.0\lib
python ProcessImRawFiles.py 
'''
