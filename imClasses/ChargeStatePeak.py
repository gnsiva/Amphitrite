from lib import utils
import numpy as np

class ChargeStatePeak():

    def __init__(self,mass,charge,speciesAllCharges):
        # self.mass = None
        # self.charge = None
        # self.species = None

        self.mass = mass
        self.charge = charge

        self.speciesAllCharges = sorted(speciesAllCharges,reverse=True)
        self.mz = utils.get_mz(self.mass,self.charge)

        self.limits = [None,None]


    def calcLimits(self,allMzs):
        allMzs = sorted(allMzs)
        allMzs = np.array(allMzs)

        diffs = allMzs - self.mz
        zeroDiffIs = [i for i,v in enumerate(diffs == 0) if v]
        lowerMz,upperMz = None, None
        
        # down
        if zeroDiffIs[0] != 0:
            downI = zeroDiffIs[0] - 1
            lowerMz = self.mz - ((self.mz - allMzs[downI])/2)
            
        # up
        if zeroDiffIs[-1] < len(allMzs)-1:
            upI = zeroDiffIs[-1] + 1
            upperMz = self.mz + ((allMzs[upI] - self.mz)/2)

        # correct end charges
        lowerMz,upperMz = self.correctForFirstLastCharge(lowerMz,upperMz)
        self.limits = [lowerMz,upperMz]


    def getLimits(self,allMzs):
        self.calcLimits(allMzs)
        return self.limits

    
    def correctForFirstLastCharge(self,lowerMz,upperMz):

        if len(self.speciesAllCharges) != 1:
            # fix for leftmost charge
            if self.charge == self.speciesAllCharges[0]:
                diff = upperMz - self.mz
                lowerMz = self.mz - diff
            # fix for right
            if self.charge == self.speciesAllCharges[-1]:
                diff = self.mz - lowerMz
                upperMz = self.mz + diff
        else:
            # NEED BUG WARNING
            # I am fudging it here and giving arbitary width
            # (corner case -> only one charge for the species)
            lowerMz = self.mz - 20
            upperMz = self.mz + 20
        return lowerMz,upperMz
