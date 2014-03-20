import msClasses.Species as Species
import Calibration, Calibrant, Im
import matplotlib.pyplot as plt
from collections import OrderedDict
import numpy as np


class MsMsCalibrant():
    """This class allows you to build up a single Calibrant() object,
    made up from  individual MSMS files of different charge states
    Using getCalibrantOb returns a fully functional Calibrant() object
    which should be fully compatible with other classes
    """
    def __init__(self,proName):
        self.proFiles = []
        self.proName = proName.lower()
        self.charges = []
        self.d = None  # key = charge, value = filename

        # needed to check validity of proName
        calibrantOb = Calibrant.Calibrant(proName,-1)
        if not self.proName in calibrantOb.cal_mass.keys():
            print 'Calibrant name not recognised, exiting...'
            print 'See imClasses/Calibrant.py for options'
            quit()

    def setUpProteinFiles(self,d):
        """Supply dictionary of charges (key) and filenames (value)
        for a single protein
        """
        self.d = d
        for z,fn in d.items():
            pFile = MsMsChargeState()
            pFile.setProName(self.proName)
            pFile.setCharge(z)
            pFile.setFilename(fn)
            pFile.generateCalibrationOb()
            self.proFiles.append(pFile)

    def getCalibrantOb(self,peakFwhm=20):
        """This creates a single Calibrant() object
        from all the protein files
        """
        self.createBlankCalibrantOb(peakFwhm)
        self.checkIfAxesMatch()
        self.combineChargeStateData()

        # The rest of the __init__ steps in Calibrant()
        cOb = self.calibrantOb
        cOb._generateSpeciesObject(peakFwhm)
        cOb._updateSpecies(cOb.speciesObj,cOb.peakFwhm,
                           cOb.leftMultiplier,cOb.rightMultiplier)
        cOb.setCharges(self.d.keys())
        cOb._updateSpecies()


        return self.calibrantOb

    def setCharges(self,charges=0):
        if charges:
            self.calibrantOb.setCharges(charges)
        else:
            if self.d:
                self.calibrantOb.setCharges(self.d.keys())
            else:
                print 'Dictionary of charges and files not given'
                print 'Using default charges'


    def createBlankCalibrantOb(self,peakFwhm):
        """Setup calibrationOb without any data
        data can then be added with combineChargeStateData()
        """
        self.calibrantOb = Calibrant.Calibrant(self.proName,-1)
        self.calibrantOb.imObj = Im.Im()
        spOb = Species.Species(self.proName)
        spOb.mass = self.calibrantOb.approxMass
        spOb.peakFwhm = peakFwhm


    def checkIfAxesMatch(self):
        """Check that the x and y axes from each charge state
        data file is the same
        This allows the data in the matrix to be summed directly
        """
        # TODO(gns) Change this so that it can use the x and y axes
        # as alignment for summing the matrices
        # This still won't work if the increments rather than the
        # extremities are different
        # You could then technically use np.interp if it is a real problem
        match = True
        imObj = self.proFiles[0].calibrantOb.imObj
        for proFile in self.proFiles:
            if not np.array_equal(imObj.xaxis,proFile.calibrantOb.imObj.xaxis):
                match = False
            if not np.array_equal(imObj.yaxis,proFile.calibrantOb.imObj.yaxis):
                match = False
            if match == False:
                print "Axes do not match between input files"
                print "Unable to combine, exiting"
                quit()

    def combineChargeStateData(self):
        """Sum the data matrices from each proteinFile
        also copy the x and y axes to the self.calibrantOb
        """
        imObj = self.proFiles[0].calibrantOb.imObj
        self.calibrantOb.imObj.setAxisX(imObj.xaxis)
        self.calibrantOb.imObj.setAxisY(imObj.yaxis)

        matrix = np.zeros(np.shape(imObj.matrix),dtype='float')
        for proFile in self.proFiles:
            tempMatrix = proFile.calibrantOb.imObj.matrix
            matrix += tempMatrix/tempMatrix.max() * 100
        self.calibrantOb.imObj.setMatrix(matrix)


class MsMsChargeState():
    """This file generates a Calibrant() object
    for a single charge state (isolated using MSMS)
    Usage: Use all the setters, then run generateCalibrationOb()
    """
    def __init__(self):
        self.proName = ''
        self.charge = None
        self.filename = ''
        self.calibrantOb = None

    def setProName(self,name):
        self.proName = name
    def setCharge(self,charge):
        self.charge = charge
    def setFilename(self,filename):
        self.filename = filename

    def generateCalibrationOb(self):
        self.calibrantOb = Calibrant.Calibrant(self.proName,self.filename)
    def getCalibrantOb(self):
        return self.calibrantOb


if __name__ == "__main__":

    pInfo = MsMsCalibrant('concanavalin a')

    d = OrderedDict()
    folder = '/home/ganesh/MyData/2013/11_November/131106_just_a_files/'
    d[19] = folder+'131106_CAL4130904_CONA_19Z.a' # bush doesn't see 19...
    d[20] = folder+'131106_CAL4130904_CONA_20Z.a'

    pInfo.setUpProteinFiles(d)
    calibrantOb = pInfo.getCalibrantOb()

    calibrantOb.imObj.plotHeatmap()
    calibrantOb.setCharges(d.keys())
    calibrantOb._updateSpecies()

    f = plt.figure()
    axAtds = f.add_subplot(211)
    calibrantOb.plotChargeStateAtds(axAtds)

    print calibrantOb.imObj.dataSlices['concanavalin a'][20].atd.xvals

    axMs = f.add_subplot(212)
    calibrantOb.plotMsAndExtractionLimits(axMs)

    plt.show()
