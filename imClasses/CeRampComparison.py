from imClasses import Im,CeRamp
import matplotlib.pyplot as plt
import numpy as np
from msClasses import *
import os
from lib import utils
import lib.SG as sg
from collections import OrderedDict
legendFontSize = 11

class CeRampComparison():
    def __init__(self):
        self.ceRamps = []
        self.proteinNames = []
        self.fnRoots = []
        self.smoothDistributions = False # if True then 3,2 sav golay smoothing will be performed on incoming yaxes
        self.smoothingSettings = [2,2,1] # wlen, smoothes, poly order

        self.proteinChargeDic = {} #k=protein name,v=charge, for calc'ing eV
        self.useEvs = False


    def useElectronVolts(self,proteinChargeDictionary):
        '''k=protein name,v=charge, for calc-ing eV'''
        self.proteinChargeDic = proteinChargeDictionary
        self.useEvs = True

    def addCeRamp(self,ceRamp):
        proteinName = ceRamp.proteinName
        fnRoot = ceRamp.fnRoot

        if not proteinName in self.proteinNames:
            self.proteinNames.append(proteinName)
        if not fnRoot in self.fnRoots:
            self.fnRoots.append(fnRoot)

        self.ceRamps.append(ceRamp)


    def createProteinRamps(self,name,dire,fnRoots,noOfFiles,startNo,
                           afit,sp,z,lm,rm,
                           vStart,vEnd,vStep,cal=0):

        for fnRoot in fnRoots:
            rmp = CeRamp.CeRamp(name)
            rmp.setupProteinRamp(dire,fnRoot,noOfFiles,startNo,
                                 afit,sp,z,lm,rm,
                                 vStart,vEnd,vStep,cal)
            rmp.calculateStatistics()
            if self.smoothDistributions:
                s = self.smoothingSettings
                rmp.smoothDistributions(s[0],s[1],s[2])
            self.addCeRamp(rmp)


    def getProteinInformation(self,protein):
        # TODO - hopefully retire this
        avs,stdevs = [],[]
        cbCheckList = [] #checking if each ramp was calibrated

        for ceRamp in protein:
            name = ceRamp.proteinName
            avData,stdevData,calibrated = ceRamp.getStatistics()
            avs.append(avData)
            stdevs.append(stdevData)
            cbCheckList.append(calibrated)
            vs = ceRamp.getVoltages()
        return name,vs,avs,stdevs,cbCheckList



    def plotDeltaUnfolding(self,ax,proteinName,spreadType='min/max',colour=0):
        protein = self.getProteinCeRamps(proteinName)

        name,vs,avs,stdevs,cbCheckList = self.getProteinInformation(protein)
        if self.useEvs:
            charge = self.proteinChargeDic[name]
            vs = np.array(vs)*charge
        if not colour:
            col = 'k' #utils.colourList2[i]
        else:
            col = colour

        if len(protein) > 1:
            # if there is more than one replicate
            # =================================== ADD MORE SPREADTYPES HERE
            if spreadType == 'min/max':
                mnU,mxU,avU = self.getMinMaxAverage(avs)

            startPoint = avU[0]
            for val in [mnU,mxU,avU]:
                val = val - startPoint
            # plotting
            self.plotShadedBoundaryAndAverage(ax,name,vs,mnU,mxU,avU,col)
        # if there's only one replicate
        else:
            values = avs[0] - avs[0][0]
            ax.plot(vs,values,color=col,label=name)
            ax.scatter(vs,values,color=col)

        # Axes labels (TURN THIS INTO A FUNCTION)
        #if calibrated:
        if cbCheckList[0]:
            ax.set_ylabel('$\Delta$CCS ($\AA^2$)')
        else:
            ax.set_ylabel('$\Delta$ Arrival Time (ms)')
        # TODO(gns) - commented out for unfolding report fig 1
        #ax.legend(loc='best',prop={'size':legendFontSize})
        if self.useEvs:
            ax.set_xlabel('Collision Energy (eV)')
        else:
            ax.set_xlabel('Voltage')

    def compareProteinUnfoldingCurves(self,axU,proteinNames,spreadType='min/max',colour = 0):
        # TODO(gns) - split unfoldingAndStdev() and this one to reduce redundancy
        listOlistOramps = []
        for proName in proteinNames:
            listOlistOramps.append(self.getProteinCeRamps(proName))

        for i,protein in enumerate(listOlistOramps):
            name,vs,avs,stdevs,cbCheckList = self.getProteinInformation(protein)
            if self.useEvs:
                charge = self.proteinChargeDic[name]
                vs = np.array(vs)*charge
            if not colour:
                col = utils.colourList2[i]
            else:
                col = colour

            # if there is more than one replicate
            if len(protein) > 1:
                # =================================== ADD MORE SPREADTYPES HERE
                if spreadType == 'min/max':
                    mnU,mxU,avU = self.getMinMaxAverage(avs)
                # plotting
                self.plotShadedBoundaryAndAverage(axU,name,vs,mnU,mxU,avU,col)
            else:
                axU.plot(vs,avs[0],color=col,label=name)
                axU.scatter(vs,avs[0],color=col)

            # Axes labels (TURN THIS INTO A FUNCTION)
            #if calibrated:
            if cbCheckList[0]:
                axU.set_ylabel('CCS ($\AA^2$)')
            else:
                axU.set_ylabel('Arrival Time (ms)')
            axU.legend(loc='best',prop={'size':legendFontSize})
            if self.useEvs:
                axU.set_xlabel('Collision Energy (eV)')
            else:
                axU.set_xlabel('Voltage')


    def compareProteinsUnfoldingAndStdev(self,axU,axS,proteinNames,spreadType='min/max',colour=0):
        listOlistOramps = []
        for proName in proteinNames:
            listOlistOramps.append(self.getProteinCeRamps(proName))

        calibratedCheckList = []
        for i,protein in enumerate(listOlistOramps):
            name,vs,avs,stdevs,cbCheckList = self.getProteinInformation(protein)
            if self.useEvs:
                charge = self.proteinChargeDic[name]
                vs = np.array(vs)*charge
            if not colour:
                col = utils.colourList2[i]
            else:
                col = colour

            # if there is more than one replicate
            if len(protein) > 1:
                # =================================== ADD MORE SPREADTYPES HERE
                if spreadType == 'min/max':
                    mnU,mxU,avU = self.getMinMaxAverage(avs)
                    mnS,mxS,avS = self.getMinMaxAverage(stdevs)

                # plotting
                self.plotShadedBoundaryAndAverage(axU,name,vs,mnU,mxU,avU,col)
                self.plotShadedBoundaryAndAverage(axS,name,vs,mnS,mxS,avS,col)
            else:
                axU.plot(vs,avs[0],color=col,label=name)
                axS.plot(vs,stdevs[0],color=col,label=name)

                axS.scatter(vs,stdevs[0],color=col)
                axU.scatter(vs,avs[0],color=col)

            # Axes labels (TURN THIS INTO A FUNCTION)
            #if calibrated:
            if cbCheckList[0]:
                axU.set_ylabel('CCS ($\AA^2$)')
            else:
                axU.set_ylabel('Arrival Time (ms)')
            axU.legend(loc='best',prop={'size':legendFontSize})
            if self.useEvs:
                axS.set_xlabel('Collision Energy (eV)')
                axU.set_xlabel('Collision Energy (eV)')
            else:
                axS.set_xlabel('Voltage')
                axU.set_xlabel('Voltage')
            axS.set_ylabel('Standard Deviation')

    def averageDistributions(self,ceRamps):

        xs = [ceRamp.xaxes for ceRamp in ceRamps]
        ys = [ceRamp.yaxes for ceRamp in ceRamps]

        # SHOULD CHECK IF THEY ARE ALL THE SAME
        xaxis = xs[0][0]
        # Averaging yvals
        ys = np.array(ys)
        avYs = np.average(ys,axis=0)
        print len(ys[0][0]), len(avYs),len(avYs[0]),len(xaxis)

        return xaxis,avYs

    def plotAllDistributions(self,ax,proteinName,average=False,colour='k',single=0,lift=5,oddonly=False):
        ceRamps = self.getProteinCeRamps(proteinName)
        # Get data
        vs = ceRamps[0].voltages
        if average:
            xaxis,yaxes = self.averageDistributions(ceRamps)
        else:
            xaxis = ceRamps[single].xaxes[0]
            yaxes = ceRamps[single].yaxes

        # Plotting
        tempVs = []
        for i,yaxis in enumerate(yaxes):
            if oddonly:
                if not i%2: # if not even
                    y = yaxis+i/2.*lift
                    if not i:
                        ax.plot(xaxis,y,color=colour,label=proteinName)
                    else:
                        ax.plot(xaxis,y,color=colour)
                    tempVs.append(vs[i])
            else:
                y = yaxis+i*lift
                if not i:
                    ax.plot(xaxis,y,color=colour,label=proteinName)
                else:
                    ax.plot(xaxis,y,color=colour)



        # Labels
        if ceRamps[0].calibrateDistributions:
            ax.set_xlabel('CCS ($\AA^2$)')
        else:
            ax.set_xlabel('Arrival Time (ms)')
        ax.set_yticks([])
        ax.set_ylabel('Intensity')

        if oddonly:
            utils.label1dPlots(ax,lift,tempVs,'V')
        else:
            utils.label1dPlots(ax,lift,vs,'V')
        ax.legend(loc='upper right',prop={'size':legendFontSize})


    def plotShadedBoundaryAndAverage(self,ax,name,x,low,high,average,colour='b'):
        ax.plot(x,average,color=colour,label=name)
        ax.fill_between(x,low,high,color=colour,alpha=0.2)


    def getMinMaxAverage(self,listOfSignal):
        arrayOfSignal = np.array(listOfSignal)
        minLine = np.min(arrayOfSignal,axis=0)
        maxLine = np.max(arrayOfSignal,axis=0)
        averageLine = np.average(arrayOfSignal,axis=0)
        return minLine,maxLine,averageLine


    def getProteinCeRamps(self,proteinName):
        out = []
        for rmp in self.ceRamps:
            if rmp.proteinName == proteinName:
                out.append(rmp)
        return out


    #================================================================
    # Single ramp plotting
    #================================================================
    # TODO - you need to update the multi ramp ones above
    # they are pretty terrible.
    # Base the code of these functions rather than the existing code
    # TODO - update... these functions could even be directly on CeRamp
    # they don't need the min/max stuff that this class was created for
    def plotCcsMeanCurve(self,ax,proteinName,colour='k'):
        """Only for use with single rmps per protein"""
        rmp = self.getProteinCeRamps(proteinName)[0]
        vs = self.getEvsOrVs(proteinName,rmp)
        ax.plot(vs,rmp.meansCcsds,colour+"o-",label=proteinName)
        self.labelAxisEvOrVs(ax)
        ax.set_ylabel("$\overline{\Omega}$ ($\AA^2$)")

    def plotCcsDeltaCurve(self,ax,proteinName,colour='k'):
        rmp = self.getProteinCeRamps(proteinName)[0]
        vs = self.getEvsOrVs(proteinName,rmp)
        ax.plot(vs,rmp.deltasCcsds,colour+"o-")
        self.labelAxisEvOrVs(ax)
        ax.set_ylabel("$\Delta\overline{\Omega}$ ($\AA^2$)")

    def plotCcsStdevCurve(self,ax,proteinName,colour='k'):
        rmp = self.getProteinCeRamps(proteinName)[0]
        vs = self.getEvsOrVs(proteinName,rmp)
        ax.plot(vs,rmp.stdevsCcsds,colour+"o-")
        self.labelAxisEvOrVs(ax)
        ax.set_ylabel("$\overline{\sigma}$ ($A^2$)")

    def plotCcsMeanAndStdevCurve(self,ax,proteinName,colour='k'):
        rmp = self.getProteinCeRamps(proteinName)[0]
        vs = self.getEvsOrVs(proteinName,rmp)
        upper = rmp.meansCcsds + rmp.stdevsCcsds
        lower = rmp.meansCcsds - rmp.stdevsCcsds

        ax.plot(vs,rmp.meansCcsds,colour+"o-",label=proteinName)
        ax.fill_between(vs,lower,upper,color=colour,alpha=0.3)
        self.labelAxisEvOrVs(ax)
        ax.set_ylabel("$\Omega$ ($\AA^2$)")


    def plotCcsVarianceCurve(self,ax,proteinName,colour='k'):
        rmp = self.getProteinCeRamps(proteinName)[0]
        vs = self.getEvsOrVs(proteinName,rmp)
        ax.plot(vs,rmp.variancesCcsds,colour+"o-")
        self.labelAxisEvOrVs(ax)
        ax.set_ylabel("$\overline{\sigma}^2$")

    def plotCcsAucCurve(self,ax,proteinName,colour='k'):
        rmp = self.getProteinCeRamps(proteinName)[0]
        vs = self.getEvsOrVs(proteinName,rmp)
        ax.plot(vs,rmp.aucsCcsds,colour+"o-")
        self.labelAxisEvOrVs(ax)
        ax.set_ylabel('Area under curve')

    def plotAllDistributions(self,ax,proteinName,colour='k',lift=5):
        rmp = self.getProteinCeRamps(proteinName)[0]
        vs = self.getEvsOrVs(proteinName,rmp)
        for i,(x,y) in enumerate(zip(rmp.xaxes,rmp.yaxes)):
            ax.plot(x,y+i*lift,color=colour)
        if self.useEvs:
            utils.label1dPlots(ax,lift,vs,'eV')
        else:
            utils.label1dPlots(ax,lift,vs,'V')
        rmp.labelDistributionAxes(ax)


    #================
    # Plotting support functions
    def labelAxisEvOrVs(self,ax):
        if self.useEvs:
            ax.set_xlabel('Collision Energy (eV)')
        else:
            ax.set_xlabel('Collision Energy (V)')

    def getEvsOrVs(self,proteinName,rmp):
        if self.useEvs:
            return self.proteinChargeDic[proteinName] * rmp.voltages
        else:
            return rmp.voltages.copy()
