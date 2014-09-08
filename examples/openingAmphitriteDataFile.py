"""
This is an example script for how to extract an ATD from an Amphitrite datafile.
For use in CIU applications, I recommend a for loop for the file name, as all the rest
of the settings should be the same.

n.b. the data is held as numpy arrays, so don't be surprised if it looks a bit weird 
when you print the vals
"""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

from imClasses import Im

dataFn = 'data/datafile.a'

imob = Im()
imob.loadAmphiFile(dataFn)

"""
To get x and y axes and the intensity matrix
(this is for the whole data file)
"""
print "Mass spectrum axis"
print imob.xaxis
print "Arrival time axis"
print imob.yaxis
print "Intensity matrix"
print imob.matrix

"""
To only analyse a certain peak, use an Atropos fit. 
You can do this programmatically or using the GUI
"""

atroposFn = 'data/atroposfit.afit'
imob.loadMsFit(atroposFn)

"""
Here you want to access a particular molecular species and charge state.
These would have to be the values used when you created the atropos fit
"""
speciesName = 'Apo'
charge = 13

imob.generateSpeciesSlicesFwhm(speciesName)
ds = imob.getDataSlice(speciesName,charge)

"""
Now you can access the ATD data for this charge state
"""

print "ATD x axis (arrival time)"
print ds.atd.xvals 
print "ATD intensity values (y values)"
print ds.atd.yvals 




