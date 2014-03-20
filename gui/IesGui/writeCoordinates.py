from classes import TwoDdata, MassSpectrum
from imClasses import Im

fn = '/home/ganesh/MyData/2013/04_April/130404_harpal/130404_bsa_001.a'
imOb = Im.Im()
imOb.loadFolderAuto(fn)
msOb = imOb.getMassSpectrum()

import matplotlib.pyplot as plt
f = plt.figure()
ax = f.add_subplot(111)

#msOb.plot(ax)
#plt.show()

fout = 'bsa_001.txt'
msOb.writeCoordinatesToTxt(fout)
