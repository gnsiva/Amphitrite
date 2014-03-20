# TODO(gns) - You shouldn't import paths.py in the class file as its only relevant to you
import paths as p
import imClasses.Im as Im
import matplotlib.pyplot as plt
import numpy as np
from msClasses import *
import os
from lib import utils
import lib.SG as sg

legendFontSize = 12

class CeRamp():
    def __init__(self,proteinName):
        self.proteinName = proteinName
        self.fnRoot = None
        self.filenames = []

        self.afit = None
        self.speciesName = None
        self.charge = None
        self.imOb = None

        self.leftMultiplier = 1
        self.rightMultiplier = 1

        self.calibrationOb = None
        self.calibrateDistributions = False
        #________________
        # Outputs
        self.dataSlices = None
        self.xaxes = None
        self.yaxes = None

        # Run self.calculateStatistics() to fill in these values
        self.meansAtds = None
        self.stdevsAtds = None
        self.deltasAtds = None
        self.variancesAtds = None
        self.aucsAtds = None

        self.meansCcsds = None
        self.stdevsCcsds = None
        self.deltasCcsds = None
        self.variancesCcsds = None
        self.aucsCcsds = None

    #================================================================
    # Setup (auto)
    def setupProteinRamp(self,dire,fnRoot,noOfFiles,startNo,
                           afit,sp,z,lm,rm,
                           vStart,vEnd,vStep,cal=0):

        self.fileNameGenerator(dire+fnRoot,noOfFiles,startNo)
        self.setAtroposSettings(afit,sp,z)
        self.voltagesGenerator(vStart,vEnd,vStep)
        self.setMultiplers(lm,rm)
        self.extractDataFromFiles()

        if cal:
            self.setCalibration(cal)
        #self.calculateStatistics()

    #================================================================
    # Filenames and voltages
    def fileNameGenerator(self,rootName,noOfFiles,start=0):
        self.fnRoot = os.path.basename(rootName)
        if start == 0:
            numbers = np.arange(noOfFiles)
        else:
            numbers = np.arange(start,noOfFiles+start)
        for number in numbers:
            self.filenames.append("%s%03d.a" %(rootName,number))

    def setFilenames(self,filenames,rootName):
        self.filenames = filenames
        self.fnRoot = os.path.basename(rootName)

    def getFilenames(self,filenames):
        return self.filenames

    def voltagesGenerator(self,start,stop,step,autocrop=True):
        self.voltages = np.arange(start,stop+0.001,step)
        # If an error points here its probably because you
        # havent set self.charge yet
        self.eVs = self.voltages * self.charge
        if autocrop:
            self.voltages = self.voltages[:len(self.filenames)]
            self.eVs = self.eVs[:len(self.filenames)]

    def setVoltages(self,voltages):
        self.voltages = voltages
        self.eVs = self.voltages *self.eVs

    def getVoltages(self):
        return self.voltages
    def getEvs(self):
        return self.eVs
    #================================================================
    # aFits, Calibrations and Multipliers
    def setAtroposSettings(self,afit,speciesName,charge):
        self.afit = afit
        self.imOb = Im.Im()
        self.imOb.loadMsFit(self.afit)

        self.speciesName = speciesName
        self.charge = charge

    def setMultiplers(self,leftMultiplier,rightMultiplier):
        self.leftMultiplier = leftMultiplier
        self.rightMultiplier = rightMultiplier

    def setCalibration(self,calibration,changeToCcsDs=True):
        '''Takes filename for calibration'''
        calibrationOb = np.load(calibration)
        self.calibrationOb = calibrationOb
        self.generateCcsDistributions()
        if changeToCcsDs:
            self.useCcsDs()
            self.updateAxes()

    def generateCcsDistributions(self):
        for ds in self.dataSlices:
            ds.generateCcsDistribution(self.calibrationOb)

    def useAtds(self):
        if self.calibrateDistributions == True:
            self.calibrateDistributions = False
            self.updateAxes()
        else:
            self.calibrateDistributions = False

    def useCcsDs(self):
        if self.calibrationOb:
            if self.calibrateDistributions == False:
                self.calibrateDistributions = True
                self.updateAxes()
            else:
                self.calibrateDistributions = True
        else:
            print 'No calibration object found!'

    #================================================================
    # Processing
    def extractDataFromFiles(self):
        self.means,self.stdevs = [],[]
        #self.atds = []
        self.dataSlices = []

        for name in self.filenames:
            self.imOb.loadAmphiFile(name)
            self.imOb.generateSpeciesSlicesFwhm(speciesName=self.speciesName,
                                                leftMultiplier=self.leftMultiplier,
                                                rightMultiplier=self.rightMultiplier)
            dataSlice = self.imOb.getDataSlice(self.speciesName,self.charge)

            # Storage
            dataSlice.atd.normalisationBpi()
            self.dataSlices.append(dataSlice)

            self.updateAxes()

    def updateAxes(self):
        self.xaxes,self.yaxes = [],[]

        for ds in self.dataSlices:
            if self.calibrateDistributions:
                self.xaxes.append(ds.ccsd.xvals)
                self.yaxes.append(ds.ccsd.yvals)
            else:
                self.xaxes.append(ds.atd.xvals)
                self.yaxes.append(ds.atd.yvals)

    def removeNans(self,xaxis,yaxis):
        out_x = xaxis[np.invert(np.isnan(xaxis))]
        out_y = yaxis[np.invert(np.isnan(xaxis))]
        return out_x,out_y


    #================================================================
    # Calculating statistics
    def calculateStatistics(self):
        # For CCSs
        if self.calibrationOb:
            if self.dataSlices[0].ccsd:
                dataSlices = [ob.ccsd for ob in self.dataSlices]
                self._calculateCcsdsStatistics(dataSlices)
            else:
                print 'No CcsD is found in dataSlice'
        else:
            print 'No calibration object found'
        # For Atds
        dataSlices = [ob.atd for ob in self.dataSlices]
        self._calculateAtdsStatistics(dataSlices)

    def _calculateCcsdsStatistics(self,dataSlices):
        avs,stdevs = self._statisticsCalcAvStd(dataSlices)
        self.meansCcsds = avs
        self.stdevsCcsds = stdevs
        self.deltasCcsds = self._statisticsCalcDelta(avs)
        self.variancesCcsds = self._statisticsCalcVariance(stdevs)
        self.aucsCcsds = self._statisticsCalcAuc(dataSlices)

    def _calculateAtdsStatistics(self,dataSlices):
        avs,stdevs = self._statisticsCalcAvStd(dataSlices)
        self.meansAtds = avs
        self.stdevsAtds = stdevs
        self.deltasAtds = self._statisticsCalcDelta(avs)
        self.variancesAtds = self._statisticsCalcVariance(stdevs)
        self.aucsAtds = self._statisticsCalcAuc(dataSlices)


    def _statisticsCalcAvStd(self,twoDdataList):
        avs = np.zeros(len(twoDdataList),dtype='float')
        stdevs = np.zeros(len(twoDdataList),dtype='float')

        for i,tD in enumerate(twoDdataList):
            av,stdev = tD.calculateWeightedMeanStandardDeviation()
            avs[i] = av
            stdevs[i] = stdev
        return avs,stdevs
    def _statisticsCalcDelta(self,avs):
        deltas = np.abs(avs - avs[0])
        return deltas
    def _statisticsCalcVariance(self,stdevs):
        variances = stdevs**2
        return variances
    def _statisticsCalcAuc(self,twoDdataList):
        """Function for calculating the area under the curve
        Intended for use on Atds and CCSDs
        Second variable is asking for xaxis
        Either self.eVs or self.voltages
        """
        aucs = np.zeros(len(twoDdataList),dtype='float')
        for i,dist in enumerate(twoDdataList):
            aucs[i] = dist.calculateAreaUnderCurve()
        return aucs


    def getAtdStatistics(self):
        return self.meansAtds,self.stdevsAtds
    def getCcsStatistics(self):
        return self.meansCcsds,self.stdevsCcsds

    def getStatistics(self):
        '''returns stats of currently being used form (e.g. atd or ccs)
        final return is boolean for if calibrated results are sent'''
        if self.calibrateDistributions:
            return self.meansCcsds,self.stdevsCcsds,True
        else:
            return self.meansAtds,self.stdevsAtds,False
        return self.meansCcsds,self.stdevsCcsds


    #================================================================
    # Smoothing distributions
    def smoothDistributions(self,window_len=2,smoothes=1,poly_order=1):
        for ds in self.dataSlices:
            if self.calibrationOb:
                ds.ccsd.smoothingSG(window_len,smoothes,poly_order)
            ds.atd.smoothingSG(window_len,smoothes,poly_order)
        self.updateAxes()


    #================================================================
    # Plotting
    def plotAllDistributions(self,ax,colour='k',lift=5):
        # TODO - fix this
        # this has basically been rewritten in CeRampComparison.py
        for i,(x,y) in enumerate(zip(self.xaxes,self.yaxes)):
            ax.plot(x,y+i*lift,color=colour)
        utils.label1dPlots(ax,lift,self.voltages,'V')

        self.labelDistributionAxes(ax)


    def plotSingleUnfoldingCurve(self,ax,colour='k',marker='.'):
        if self.calibrateDistributions:
            ax.plot(self.voltages,self.meansCcsds,'-'+colour+marker)
            ax.set_ylabel('CCS ($\AA^2$)')
        else:
            ax.plot(self.voltages,self.meansAtds,'-'+colour+marker)
            ax.set_ylabel('Arrival Time (ms)')
        ax.set_xlabel('Trap Voltage')


    def plotSingleStdevCurve(self,ax,colour='k',marker='.'):
        if self.calibrateDistributions:
            ax.plot(self.voltages,self.stdevsCcsds,'-'+colour+marker)
            ax.set_ylabel('Standard Deviation')
        else:
            ax.plot(self.voltages,self.stdevsAtds,'-'+colour+marker)
            ax.set_ylabel('Standard Deviation')
        ax.set_xlabel('Trap Voltage')

    def drawTriplePlot(self,figsize=0,lift=5):
        # setting up figure
        if not figsize:
            figsize = [7,10]
        f = plt.figure(figsize=figsize)
        import matplotlib as mpl
        gs = mpl.gridspec.GridSpec(4,1)

        # plotting atds
        ax = f.add_subplot(gs[0:2])
        self.plotAllDistributions(ax,lift=lift)
        # plotting unfolding
        ax2 = f.add_subplot(gs[2])
        self.plotSingleUnfoldingCurve(ax2)
        # plotting stdevs
        ax3 = f.add_subplot(gs[3])
        self.plotSingleStdevCurve(ax3)

        plt.tight_layout()

    #================================================================
    # Fitting

    def leastSquaresStaticMeans(self,xaxis,yaxis,peakTops,amp=60):
        fwhm = self.getStartingWidth()
        amp = float(amp)
        p0 = []
        for mean in peakTops:
            p0.append(amp)
            p0.append(fwhm)

        xaxis,yaxis = self.removeNans(xaxis,yaxis)

        # Running fitting
        import time
        from scipy import optimize
        startTime = time.time()
        p1, success = optimize.leastsq(self.fitFuncStaticMeans,p0[:],
                                       args=(xaxis,yaxis,peakTops))#,maxfev=100000)
        print "Optimisation took:", time.time()-startTime, "s"
        p1 = np.abs(p1)
        sim = self.distributionFromP(p1,xaxis,peakTops)
        return xaxis,sim,p1

    def fitFuncStaticMeans(self,p,xaxis,yaxis,peakTops):
        p = np.abs(p)
        simulation = self.distributionFromP(p,xaxis,peakTops)
        return np.abs(yaxis - simulation)


    def distributionFromP(self,p,xaxis,peakTops):
        '''returns simulated distribution'''
        simulation = np.zeros(len(xaxis),dtype=np.float)
        gGen = self.gaussParameterGeneratorStaticMeans(p,peakTops)
        for mean in peakTops:
            gP = gGen.next()
            simulation += utils.gaussian(xaxis,gP[1],gP[0],gP[2])
        return simulation

    def gaussParameterGeneratorStaticMeans(self,p,peakTops):
        '''return the 3 gauss parameters from a distributions'''
        for i,mean in enumerate(peakTops):
            gaussP = np.zeros(3,dtype=np.float)
            gaussP[0] = peakTops[i]         # mean
            gaussP[1] = p[i*2]      # amp
            gaussP[2] = p[i*2+1]    # fwhm
            yield gaussP

    def getStartingWidth(self,factor=100):
        '''returns a quarter of the width of the xaxis'''
        xaxis = self.xaxes[0]
        xaxis = xaxis[np.invert(np.isnan(xaxis))]  # remove nans
        return (xaxis.max()-xaxis.min())/factor

    #================================
    # plotting (fitting only)
    def fitAndPlotAll(self,ax,peakTops,lift=20):
        self.plotAllDistributions(ax,lift=lift)
        for i in range(len(self.xaxes)):

            xaxis = self.xaxes[i]
            yaxis = self.yaxes[i]

            sim = self.leastSquaresStaticMeans(xaxis,yaxis,peakTops)

            liftInd = i*lift
            self.plotStaticMeansGaussians(ax,sim[2],xaxis,peakTops,liftInd)
            self.plotStaticMeansSimulation(ax,sim[2],xaxis,peakTops,liftInd)

    def plotStaticMeansGaussians(self,ax,p1,xaxis,peakTops,lift):
        '''plots the individual gaussians determined by the fit'''
        gGen = self.gaussParameterGeneratorStaticMeans(p1,peakTops)
        for j,mean in enumerate(peakTops):
            gP = gGen.next()
            y = utils.gaussian(xaxis,gP[1],gP[0],gP[2])
            ax.plot(xaxis,y+lift,color=utils.colourList[j])

        self.labelDistributionAxes(ax)

    def plotStaticMeansSimulation(self,ax,p1,xaxis,peakTops,lift):
        ySim = self.distributionFromP(p1,xaxis,peakTops)
        ax.plot(xaxis,ySim+lift,color='r')
        self.labelDistributionAxes(ax)


    def labelDistributionAxes(self,ax):
        if self.calibrateDistributions:
            ax.set_xlabel('CCS ($\AA^2$)')
        else:
            ax.set_xlabel('Arrival Time (ms)')
        ax.set_ylabel('Intensity')
        ax.set_yticks([])

    #================================================================
    # Genetic algorthim fitting
    def geneticAlgorithmStaticMeans(self,xaxis,yaxis,peakTops,maxFwhm,popsize,
                                    generations,elitism=False,mutationRate=False):
        '''i is the index for the voltage you want to fit for
        peakTops are the means of the gaussians you want to fit
        if not False elitism should be set to a number
        if not False mutationRate should be set to a number (0->1)'''

        import StaticMeansChromosome as Smc
        from pyevolve import G1DList, GSimpleGA, Selectors, Scaling, Crossovers
        from pyevolve import Initializators, Mutators, Consts, DBAdapters

        chromosome = Smc.StaticMeansChromosome(peakTops,maxFwhm)
        cLen = chromosome.getChromosomeLength()
        cLimits = chromosome.getChromosomeLimits()

        #================
        # Genome initialisation
        genome = G1DList.G1DList(cLen)
        genome.mutator.set(Mutators.G1DListMutatorGish)
        genome.initializator.set(Initializators.G1DListInitializatorGish)
        genome.evaluator.set(self.fitFuncStaticMeansGeneticAlgorithm)
        #================
        # Passing parameters
        xaxis,yaxis = self.removeNans(xaxis,yaxis)
        params = {'xaxis':xaxis,'yaxis':yaxis,'peakTops':peakTops}
        genome.setParams(**params)
        #================
        # Genetic algorithm initialisation
        ga = GSimpleGA.GSimpleGA(genome)
        ga.setPopulationSize(popsize)

        ga.setParams(limits=cLimits)
        ga.selector.set(Selectors.GRouletteWheel)

        pop = ga.getPopulation()
        pop.scaleMethod.set(Scaling.SigmaTruncScaling)

        ga.setGenerations(generations)
        ga.minimax = Consts.minimaxType["minimize"]
        #ga.setMultiProcessing(True)
        #================
        # Elitism
        if not elitism:
            ga.setElitism(False)
        else:
            ga.setElitism(True)
            ga.setElitismReplacement(elitism)
        #================
        # Mutation rate
        if mutationRate:
            ga.setMutationRate(mutationRate)
        #================
        # Crossover
        # NOT DONE YET (SEE 130226_GA_dummy_data_tests_2.0.4_cm_testing.py)

        #================
        # Run it
        ga.evolve(freq_stats=100,limits=cLimits)

        #================
        # Return output
        best = ga.bestIndividual()
        print best
        simulation = self.distributionFromC(list(best),xaxis,peakTops)
        return xaxis,simulation,list(best)

    def fitFuncStaticMeansGeneticAlgorithm(self,chromosome):
        error = 0.0
        xaxis = chromosome.getParam('xaxis')
        yaxis = chromosome.getParam('yaxis')
        peakTops = chromosome.getParam('peakTops')

        simulation = self.distributionFromC(chromosome,xaxis,peakTops)
        error = np.sum(np.abs(yaxis-simulation))
        return error

    def splitChromosomeToGaussians(self,chromosome,peakTops):
        '''this function shouldnt really be here
        it already exists in StaticMeansChromosome.
        (only does one distribution at a time)
        per gaussian [mu,amp,fwhm]'''
        cGaussians = [[] for g in peakTops]
        counter = 0
        for i,peakTop in enumerate(peakTops):
            cGaussians[i].append(peakTop)
            cGaussians[i].append(chromosome[counter])
            counter += 1
            cGaussians[i].append(chromosome[counter])
            counter += 1
        return cGaussians

    def distributionFromC(self,chromosome,xaxis,peakTops):
        '''chromosome output from a static means genetic algorithm'''
        simulation = np.zeros(len(xaxis),dtype=np.float)
        cGaussians = self.splitChromosomeToGaussians(chromosome,peakTops)
        for gC in cGaussians:
            simulation += utils.gaussian(xaxis,gC[1],gC[0],gC[2])
        return simulation

    #================================
    # Plotting for genetic algorthim fitting
    def plotStaticMeansGeneticAlgorithm(self,ax,chromosome,xaxis,peakTops,lift=0):
        simulation = self.distributionFromC(chromosome,xaxis,peakTops)
        xaxis,simulation = self.removeNans(xaxis,simulation)
        ax.plot(xaxis,simulation+lift,color='r')
        self.labelDistributionAxes(ax)

    def plotStaticMeansGaGaussians(self,ax,chromosome,xaxis,peakTops,lift=0):
        cGaussians = self.splitChromosomeToGaussians(chromosome,peakTops)
        for i,gC in enumerate(cGaussians):
            y = utils.gaussian(xaxis,gC[1],gC[0],gC[2])
            ax.plot(xaxis,y+lift,color=utils.colourList[i])

    def plotStaticMeanConformationTracking(self,ax,bests,peakTops,xaxes,absolute=True):
        stackOfDistributions = []
        for i,best in enumerate(bests):
            conformations = []
            cGaussians = self.splitChromosomeToGaussians(best,peakTops)
            # [mu,amp,fwhm]
            xaxis = xaxes[i][np.invert(np.isnan(xaxes[i]))]

            total = 0.0
            for j,g in enumerate(cGaussians):
                y = utils.gaussian(xaxis,g[1],g[0],g[2])
                integral = np.trapz(y,x=xaxis)
                conformations.append(integral)
                total += integral

            if not absolute:
                for j in xrange(len(conformations)):
                    conformations[j] = conformations[j]/total*100

            stackOfDistributions.append(conformations)
        stackOfDistributions = np.array(stackOfDistributions)

        for i in xrange(len(stackOfDistributions[0])):
            ax.plot(self.voltages,stackOfDistributions[:,i],color=utils.colourList[i],
                    label=str(peakTops[i])+' $\AA^2$')
            ax.scatter(self.voltages,stackOfDistributions[:,i],color=utils.colourList[i])
        #ax.legend(loc='upper right',prop={'size':legendFontSize})
        ax.set_xlabel('Voltage')
        if not absolute:
            ax.set_ylabel('Conformational Area (%)')
        else:
            ax.set_ylabel('Conformational Area')

        # stop the y axis going well negative for no reason
        ylims = ax.get_ylim()
        lower = (ylims[1]*0.01)*-1
        ax.set_ylim(lower,ylims[1])


    def _convertAbsoluteConformationAbundanceToPerc(self,conformations):
        con = np.array(conformations)
        conformations_out = np.zeros(shape=np.shape(con),dtype='float')
        for i in xrange(len(conformations[0])):
            # going round the voltages
            total = 0.0

            for j in xrange(len(conformations)):
                total += conformations[j][i]

            abundances = [conformations[k][i]/total*100 for k in xrange(len(conformations))]
            for j in xrange(len(conformations)):
                conformations_out[j][i] = abundances[j]

        return conformations_out

    def plotPeakHeightConformationTracking(self,ax,peakTops,absolute=True):
        '''Absolute flag gives absolution intensities whereas
        False gives percentage of total'''
        conformations = [[] for p in peakTops]

        for i,peakTop in enumerate(peakTops):
            for xaxis,yaxis in zip(self.xaxes,self.yaxes):
                xaxis,yaxis = self.removeNans(xaxis,yaxis)
                index = utils.closest(peakTop,xaxis)
                conformations[i].append(yaxis[index])

        if not absolute:
            conformations = self._convertAbsoluteConformationAbundanceToPerc(conformations)

        for i,conformation in enumerate(conformations):
            ax.plot(self.voltages,conformation,color=utils.colourList[i],
                    label=str(peakTops[i])+' $\AA^2$')
            ax.scatter(self.voltages,conformation,color=utils.colourList[i])

        if not absolute:
            ax.set_ylabel('Peak Height Abundance (%)')
        else:
            ax.set_ylabel('Peak Height')
        ax.set_xlabel('Voltage')

        # stop the y axis going well negative for no reason
        ylims = ax.get_ylim()
        lower = (ylims[1]*0.01)*-1
        ax.set_ylim(lower,ylims[1])

        #ax.legend(loc='best',prop={'size':legendFontSize})

    def plotPeakTopPositions(self,ax,peakTops,lift=20):
        for i,(xaxis,yaxis) in enumerate(zip(self.xaxes,self.yaxes)):
            xaxis,yaxis = self.removeNans(xaxis,yaxis)
            ax.plot(xaxis,yaxis+lift*i,color='k')

        for j,(xaxis,yaxis) in enumerate(zip(self.xaxes,self.yaxes)):
            for i,peakTop in enumerate(peakTops):
                if not j:
                    plt.axvline(peakTop,color=utils.colourList[i],lw=0.8,
                                label=str(peakTops[i])+' $\AA^2$')
                else:
                    plt.axvline(peakTop,color=utils.colourList[i],lw=0.8)

        ax.set_xlabel('CCS $\AA^2$')
        ax.set_ylabel('Intensity')
        ax.set_yticks([])
        ax.legend(loc='upper left',prop={'size':legendFontSize})

    def plotCiuContourPlot(self,ax=0):
        # TODO make a function for this
        # encapsulate it out and add it everywhere else
        if ax == 0:
            f = plt.figure()
            ax = f.add_subplot(111)

        # Get rid of the NaNs
        yaxes = []
        x = self.xaxes[0]
        for y in self.yaxes:
            xaxis,yaxis = self.removeNans(x,y)
            yaxes.append(yaxis)

        # Need to interpolate the damn Xaxis
        new_x = np.arange(xaxis.min(),xaxis.max(),dtype='float')

        yaxes2 = []
        for yaxis in yaxes:
            yaxes2.append(np.interp(new_x,xaxis,yaxis))

        ccsAxis = new_x

        # Make the matrix and get the orientation right
        matrix = np.array(yaxes2)

        matrix = np.transpose(matrix)

        vs = self.voltages

        ax.imshow(matrix,aspect='auto',origin='lower',extent=[vs[0],vs[-1],ccsAxis[0],ccsAxis[-1]])
        ax.set_xlabel('Voltage')
        ax.set_ylabel('CCS ($\AA^2$)')
