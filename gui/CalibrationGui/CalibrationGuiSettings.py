"""Class to hold the settings for the ion mobility calibration GUI."""

__author__ = "Ganesh N. Sivalingam <g.n.sivalingam@gmail.com"

class CalibrationGuiSettings():
    def __init__(self):
        self.gas = 'nitrogen'
        self.waveVelocity = '350'
        self.calPath = ''
        self.calNameCurrent = ''
        self.xlimsMs = {}
        self.ylimsMs = {}
        self.xlimsTds = {}
        self.ylimsTds = {}
        self.list_ctrl = None
        
        self.saveFilePath = ''
        self.openFilePath = ''
        
    # gas
    def setGas(self,gas):
        """Set the mobility gas used for mobility separation. (Nitrogen for
        Waters Synapt G1.)
        """
        self.gas = gas
    def getGas(self):
        """Return mobility gas."""
        return self.gas
    
    # wave velocity
    def setWaveVelocity(self,wv):
        "Set the ion mobility cell travelling wave velocity."
        try:
            float(wv)
            self.waveVelocity = wv
        except:
            print'Error: non-number wave velocity value'
    def getWaveVelocity(self):
        "Get the ion mobility cell travelling wave velocity."
        return float(self.waveVelocity)
    def getWaveVelocityString(self):
        "Set the ion mobility cell travelling wave velocity as a string."
        return self.waveVelocity
    
    # cal path
    def setCalPath(self,calPath):
        """Set the path to a data file containing a '.a' file for
        a calibrant.
        :parameter calPath: Absolute path (string)
        """
        self.calPath = calPath
    def getCalPath(self):
        """Get the path to the currently selected calibrant '.a' file.
        :return: Absolute path (string)
        """
        return self.calPath
    
    def setCalNameCurrent(self,name):
        """Set the type of calibrant associated with calPath.
        :parameter name: Protein name (e.g. bsa, see self.getCalNameCurrent()
        for full list of options)
        """
        self.calNameCurrent = name
    def getCalNameCurrent(self):
        """Return the name of the currently selected calibrant, in the
        format which is compatible with backend (imClasses.Calibrant()
        and imClasses.Calibration())
        :returns: Protein name (string)
        """
        name = self.calNameCurrent
        if name == "Myoglobin (denatured)":
            name = 'myoglobin'
        elif name == "Serum Amyloid P (5-mer)":
            name = "SAP5"
        elif name == "Serum Amyloid P (10-mer)":
            name = "SAP10"
        elif name == "Cytochrome c (denatured)":
            name = "cytochrome c denatured"
        elif name == "Cytochrome c (native)":
            name = "cytochrome c native"
        elif name == "Beta-lactoglobulin Monomer":
            name = "blac1"
        elif name == "Beta-lactoglobulin Dimer":
            name = "blac2"
        elif name == "Bovine Serum Albumin":
            name = "BSA"
        elif name == "Alcohol Dehydrogenase":
            name = "ADH"
        elif name == 'Bradykinin':
            name = 'bk'
        return name.lower()
    
    
