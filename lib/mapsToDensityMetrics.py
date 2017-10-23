from vxlsPerAtmAnalysisPlots import plotVxlsPerAtm, plotDensForAtm
from perAtomClusterAnalysis import perAtomClusterAnalysis
from densityAnalysisPlots import edens_scatter
from savevariables import save_objectlist
from itertools import izip as zip
from PDBFileManipulation import PDBtoList
from readMap import readMap
import matplotlib.pyplot as plt
from errors import error
import numpy as np
import sys
import time
import os

# check if seaborn library is present and include if so
from checkDependencies import checkDependencies
c = checkDependencies()
if c.checkSeaborn():
    import seaborn as sns


class maps2DensMetrics():

    # assign values within a density map to specific atoms, using
    # the an atom-tagged map to determine which regions of space
    # are to be assigned to each atom

    def __init__(self,
                 filesIn='', filesOut='', pdbName='', atomTagMap='',
                 densityMap='', FCmap='',  plotScatter=False, plotHist=False,
                 logFile='./untitled.log', calcFCmap=True):

        # the input directory
        self.filesIn = filesIn

        # the output directory
        self.filesOut = filesOut

        # the pdb file name
        self.pdbName = pdbName

        # atom-tagged map name
        self.map1 = atomTagMap

        # density map name (typically Fo-Fo)
        self.map2 = densityMap

        # FC map name
        self.map3 = FCmap

        # (bool) plot scatter plots or not
        self.plotScatter = plotScatter

        # (bool) plot histogram plots or not
        self.plotHist = plotHist

        # log file name
        self.log = logFile

        # whether FC map should be generated
        self.calcFCmap = calcFCmap

    def maps2atmdensity(self):

        # the map run method for this class. Will read in an atom-tagged map
        # and density map and assign density values for each individual atom
        # (as specified within the atom-tagged map). From these summary metrics
        # describing the density map behaviour in the vicinity of each refined
        # atom can be calculated

        self.readPDBfile()
        self.readAtomMap()
        self.readDensityMap()
        self.reportDensMapInfo()
        self.checkMapCompatibility()

        if self.calcFCmap:
            self.readFCMap()

        self.createVoxelList()

        if self.plotHist:
            self.plotDensHistPlots()

        self.calcDensMetrics()

        if self.plotScatter:
            self.plotDensScatterPlots()

        self.pickleAtomList()

    def readPDBfile(self):

        # read in pdb file info here. A list of atom objects
        # is created, to which density metric information
        # will be added as additional attributes in the
        # methods included below

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln='Reading pdb file: {}'.format(self.pdbName))

        # read in the pdb file to fill list of atom objects
        self.PDBarray = PDBtoList('{}{}.pdb'.format(
            self.filesIn, self.pdbName))
        self.stopTimer()

        # make sure array of atoms ordered by atom number
        self.PDBarray.sort(key=lambda x: x.atomnum)

    def readAtomMap(self):

        # read in the atom-tagged map

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln='Reading atom-tagged map file...')
        self.lgwrite(ln='Atom map name: {}'.format(self.map1))

        self.atmmap, self.atomIndices = readMap(
            dirIn=self.filesIn, dirOut=self.filesOut, mapName=self.map1,
            mapType='atom_map', log=self.log)
        self.stopTimer()

        # find number of atoms in structure
        numAtms = len(self.PDBarray)

        # find atom numbers present in list (repeated atom numbers removed)
        seen = set()
        seenAdd = seen.add
        uniqAtms = [x for x in self.atmmap.vxls_val if not
                    (x in seen or seenAdd(x))]

        # find set of atoms numbers not present
        # (i.e atoms not assigned to voxels)
        AtmsNotPres = set(range(1, numAtms+1)) - set(uniqAtms)
        self.lgwrite(
            ln='Number of atoms not assigned to voxels: ' +
               '{}'.format(len(AtmsNotPres)))

    def readDensityMap(self):

        # read in the density map

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln='Reading density map file...')
        self.lgwrite(ln='Density map name: {}'.format(self.map2))

        self.densmap = readMap(dirIn=self.filesIn, dirOut=self.filesOut,
                               mapName=self.map2, mapType='density_map',
                               atomInds=self.atomIndices, log=self.log)
        self.stopTimer()

    def readFCMap(self):

        # read in the FC (calculated structure factor) density map.
        # This method should not be called if no FC density map
        # has been provided in the current run.

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln='Reading Fcalc density map file...')
        self.lgwrite(ln='Density map name: {}'.format(self.map3))

        self.FCmap = readMap(dirIn=self.filesIn, dirOut=self.filesOut,
                             mapName=self.map3, mapType='density_map',
                             atomInds=self.atomIndices, log=self.log)

        self.stopTimer()

    def reportDensMapInfo(self,
                          numSfs=4):

        # report the density map summary information to a log file

        totalNumVxls = np.product(self.atmmap.nxyz.values())
        structureNumVxls = len(self.densmap.vxls_val)
        totalMean = self.densmap.density['mean']
        structureMean = np.mean(self.densmap.vxls_val)
        solvNumVxls = totalNumVxls - structureNumVxls
        solvMean = (totalNumVxls*totalMean -
                    structureNumVxls*structureMean)/solvNumVxls

        self.lgwrite(
            ln='\nFor voxels assigned to structure:\n' +
               '\tmean structure density : {}\n'.format(
                round(structureMean, numSfs)) +
               '\tmax structure density : {}\n'.format(
                round(max(self.densmap.vxls_val), numSfs)) +
               '\tmin structure density : {}\n'.format(
                round(min(self.densmap.vxls_val), numSfs)) +
               '\tstd structure density : {}\n'.format(
                round(np.std(self.densmap.vxls_val), numSfs)) +
               '\t# voxels included : {}\n'.format(structureNumVxls) +
               '\nFor voxels assigned to solvent:\n' +
               '\tmean solvent-region density : {}\n'.format(
                round(solvMean), numSfs) +
               '\t# voxels included : {}'.format(solvNumVxls))

    def checkMapCompatibility(self):

        # check that atom-tagged and density map
        # can be combined successfully. This
        # requirement is met if the maps have the
        # the same map header information. Grid
        # dimensions are permitted to deviate
        # between the two maps, however this is
        # flagged at run time

        self.printStepNumber()
        self.lgwrite(
            ln='Checking that maps have same dimensions and sampling...')

        self.startTimer()
        # Check that the maps have the same dimensions, grid sampling,..
        if (self.atmmap.axis != self.densmap.axis or
            self.atmmap.gridsamp != self.densmap.gridsamp or
            self.atmmap.start != self.densmap.start or
            self.atmmap.nxyz != self.densmap.nxyz or
                self.atmmap.type != self.densmap.type):

            error(text='Incompatible map properties',
                  log=self.log, type='error')

        elif self.atmmap.celldims != self.densmap.celldims:
            self.lgwrite(ln='Not exact same map grid dimensions..')
            # now check if grid dims same to a
            # specific dp and consider continuing
            stop = True
            for i in list(reversed(range(7))):
                count = 0
                for key in self.atmmap.celldims.keys():
                    roundedAtmmapDim = np.round(self.atmmap.celldims[key], i)
                    roundedDensmapDim = np.round(self.densmap.celldims[key], i)
                    if roundedAtmmapDim == roundedDensmapDim:
                        count += 1
                if count == 6:
                    self.lgwrite(
                        ln='Map grid dimensions same to {}dp\n'.format(i) +
                           '--> continuing with processing anyway')
                    stop = False
                    break
            if stop:
                    error(text='Map grid dimensions still not same to 0dp',
                          log=self.log, type='error')

        else:
            self.success()
            self.lgwrite(
                ln='The atom and density map are of compatible format!')
        self.stopTimer()

        self.lgwrite(
            ln='Total number of voxels assigned to atoms: {}'.format(
                len(self.atmmap.vxls_val)))

    def createVoxelList(self,
                        plotVoxelsAssignedToAtoms=False,
                        getVoxelXYZs=False):

        # create dictionary of voxels with atom numbers as keys

        self.startTimer()
        self.printStepNumber()
        self.lgwrite(ln='Combining voxel density and atom values...')
        self.success()
        vxlDic = {atm: [] for atm in self.atmmap.vxls_val}
        xyzDic = {atm: [] for atm in self.atmmap.vxls_val}

        self.densmap.reshape1dTo3d()
        self.densmap.abs2xyz_params()
        for atm, dens in zip(self.atmmap.vxls_val, self.densmap.vxls_val):
            vxlDic[atm].append(dens)

        # The following is not essential for run (TODO maybe omit)
        if getVoxelXYZs:
            xyz_list = self.densmap.getVoxXYZ(
                self.atomIndices, coordType='fractional')
            for atm, xyz in zip(self.atmmap.vxls_val, xyz_list):
                xyzDic[atm].append(xyz)

            # can visualise voxel positions if required
            if plotVoxelsAssignedToAtoms:
                self.plot3DvoxelPositions(xyzDic)
            self.xyzsPerAtom = xyzDic

        self.vxlsPerAtom = vxlDic

        if self.calcFCmap:
            vxlDic2 = {atm: [] for atm in self.atmmap.vxls_val}
            for atm, dens in zip(self.atmmap.vxls_val, self.FCmap.vxls_val):
                vxlDic2[atm].append(dens)
            self.FCperAtom = vxlDic2

        self.deleteMapsAttributes()
        self.stopTimer()

    def plot3DvoxelPositions(self,
                             xyzDic={}, xMin=0.48, xMax=1, yMin=0,
                             yMax=0.34, zMin=0.75, zMax=0.8):

        # Allow the plotting of the 3D locations of all
        # voxels assigned to atoms within the grid limits
        # specified by (xMin,xMax), (yMin,yMax) and (zMin,zMax).
        # This should not be selected in a standard run of
        # the code, however can be used for testing that
        # voxels have been suitably assigned to atoms

        import matplotlib.pyplot as plt
        xSubset, ySubset, zSubset = [], [], []

        count = 0
        colors = []
        for atm in xyzDic.keys():
            xyz = xyzDic[atm]
            found = False
            for i in range(len(xyz)):
                if xyz[i][0] > yMin and xyz[i][0] < yMax:
                    if xyz[i][1] > yMin and xyz[i][1] < yMax:
                        if xyz[i][2] > zMin and xyz[i][2] < zMax:
                            if not found:
                                count += 1
                            xSubset.append(xyz[i][0])
                            ySubset.append(xyz[i][1])
                            zSubset.append(xyz[i][2])
                            colors.append(count)
                            found = True

        colorsFinal = []
        for col in colors:
            if col % 5 == 0:
                colorsFinal.append('r')
            elif col % 5 == 1:
                colorsFinal.append('b')
            elif col % 5 == 2:
                colorsFinal.append('b')
            elif col % 5 == 3:
                colorsFinal.append('b')
            elif col % 5 == 4:
                colorsFinal.append('b')
            else:
                colorsFinal.append('b')

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xSubset, ySubset, zSubset, s=10, c=colorsFinal)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()

    def deleteMapsAttributes(self):

        # Provide the option to delete atmmap and
        # densmap attributes to save memory, if
        # they are no longer needed during a run

        del self.atmmap
        if self.calcFCmap:
            del self.FCmap
        del self.densmap.vxls_val

    def plotDensHistPlots(self,
                          getVoxelStats=False, perAtmDensHist=False):

        # Create histogram or kde plots of number of voxels per atom

        self.startTimer()
        self.printStepNumber()
        self.lgwrite(ln='Plotting histogram plots of voxels per atom...')
        self.lgwrite(ln='Plots written to "{}plots"'.format(self.filesOut))

        stats = plotVxlsPerAtm(pdbName=self.pdbName, where=self.filesOut,
                               vxlsPerAtom=self.vxlsPerAtom, plotType='both',
                               returnStats=getVoxelStats)

        if stats != '':
            print 'mean: {}\nstd: {}\nmax: {}\nmin: {}'.format(*stats)

        if perAtmDensHist:
            plotDensForAtm(pdbName=self.pdbName, where=self.filesOut,
                           vxlsPerAtom=self.vxlsPerAtom, plotType='both',
                           PDBarray=self.PDBarray)

        self.stopTimer()

    def calcDensMetricsForAtom(self,
                               atom=[], plotDistn=False, clustAnalys=False):

        # calculate density metrics for a particular atom.
        # This method includes the option to perform
        # cluster analysis on the voxel values assigned
        # to this atom, however this should not be selected
        # for a standard run of the code

        try:
            atomVxls = self.vxlsPerAtom[atom.atomnum]
        except KeyError:
            error(
                text='Warning!: No voxels assigned to an atom. Consider ' +
                     'increasing per-atom search radius parameter in RIDL ' +
                     'input .txt file.',
                log=self.log, type='warning')
            atomVxls = [np.nan]

        if len(atomVxls) != 0:
            atom.meandensity = np.mean(atomVxls)
            atom.mediandensity = np.median(atomVxls)
            atom.mindensity = min(atomVxls)
            atom.maxdensity = max(atomVxls)
            atom.stddensity = np.std(atomVxls)
            atom.min90tile = np.percentile(atomVxls, 10)
            atom.max90tile = np.percentile(atomVxls, 90)
            atom.min95tile = np.percentile(atomVxls, 5)
            atom.max95tile = np.percentile(atomVxls, 95)
            atom.numvoxels = len(atomVxls)

            posVals = [w for w in atomVxls if w > 0]
            if posVals != []:
                atom.meanPosOnly = np.mean(posVals)
            else:
                atom.meanPosOnly = 0

            negVals = [w for w in atomVxls if w < 0]
            if negVals != []:
                atom.meanNegOnly = np.mean(negVals)
            else:
                atom.meanNegOnly = 0

            if self.calcFCmap:
                # if the user has opted to calculate an Fcalc
                # map in addition to the difference map, then
                # additional metrics can be derived using this
                # map. These metrics typically use the Fcalc
                # map density at each voxel to weight the
                # contribution that each voxel's difference map
                # value should play when calculating damage metrics.
                # Effectively, a voxel far from an atom (but still
                # included in the search radius around that atom)
                # should not contribute to a damage indicator
                # as much as a voxel close to the atomic centre

                atomFCvals = self.FCperAtom[atom.atomnum]
                # currently set all negative values to zero.
                # This has the effect of ignoring Fcalc density that
                # is less than the map mean. This is implemented
                # such that all per-voxel weights (see below) are
                # positive and so therefore sensible weighted-means
                # can be calculated. This may need to be reconsidered
                # for future use!
                atomFCvals = [v if v > 0 else 0 for v in atomFCvals]

                atomFCvalsMaxNormed = np.array(atomFCvals)/max(atomFCvals)

                minIndex = np.array(atomVxls).argmin()
                weightedVxls = np.multiply(atomVxls, atomFCvalsMaxNormed)

                atom.densityWeightedMean = np.mean(weightedVxls)
                atom.densityWeightedMin = np.min(weightedVxls)
                atom.densityWeightedMax = np.max(weightedVxls)

                # the following attribute provides an indication of the
                # fraction of the local maximum Fcalc map density around
                # the current atom at the point where the minimum difference
                # map value has been located to be. A higher value (closer to
                # 1) indicates that the min density value is found at an
                # electron density-rich region of space, whereas a lower
                # value (closer to 0) indicates that the min density value is
                # located away from where the majority of the electron density
                # assigned to the atom is predicted to be.
                atom.fracOfMaxAtomDensAtMin = atomFCvalsMaxNormed[minIndex]

                posVals = [w for w in weightedVxls if w > 0]
                negVals = [w for w in weightedVxls if w < 0]
                posValsSum = np.sum(posVals)
                negValsSum = np.sum(negVals)

                posWeights = [v for v, w in zip(
                    atomFCvalsMaxNormed, weightedVxls) if w > 0]
                negWeights = [v for v, w in zip(
                    atomFCvalsMaxNormed, weightedVxls) if w < 0]
                posWeightsSum = np.sum(posWeights)
                negWeightsSum = np.sum(negWeights)

                if posVals != []:
                    atom.densityWeightedMeanPosOnly = posValsSum/posWeightsSum
                else:
                    atom.densityWeightedMeanPosOnly = 0

                if negVals != []:
                    atom.densityWeightedMeanNegOnly = negValsSum/negWeightsSum
                else:
                    atom.densityWeightedMeanNegOnly = 0

                if plotDistn:
                    self.plotFCdistnPlot(
                        atomsToPlot=['CYS-SG'], atomOfInterest=atom,
                        atomFCvals=atomFCvals, FCatMin=atomFCvals[minIndex],
                        atomFCvalsMaxNorm=atomFCvalsMaxNormed,
                        reliability=atom.reliability)

            if clustAnalys:
                # provides the user with the option to also run
                # per-atom cluster analysis on the spatial
                # distribution of voxels assigned to a single atom.
                # This would be useful to distinguish 'clumps' of
                # positive or negative difference density, in order
                # to decide whether an atom may have shifted
                # position upon irradiation.
                # It should be noted that this option takes a
                # significant time to run, and should be deselected
                # in a standard run of the code

                self.clustDoneOnAtm.append(atom.getAtomID())
                clustAnalysis = perAtomClusterAnalysis(
                    atmNum=atom.atomnum, atmId=atom.getAtomID(),
                    vxlsPerAtom=self.vxlsPerAtom, xyzsPerAtom=self.xyzsPerAtom,
                    densMapObj=self.densmap, prevAtmMidPt=self.atomMidPts[-1])

                self.atomMidPts.append(clustAnalysis.midPt)

                atom.negClusterVal = clustAnalysis.output[0]
                atom.totDensShift = clustAnalysis.output[-1]

                self.densByRegion.append(clustAnalysis.densByRegion)

    def calcDensMetrics(self,
                        plotDistn=False, clustAnalys=False,
                        showProgress=True, parallel=False,
                        makeTrainSet=False, inclOnlyGluAsp=False):

        # determine density summary metrics per atom.
        # 'includeOnlyGluAsp' allows calculations to
        # be performed only for Glu/asp carboxylates
        # (this is not typically suitable and will
        # cause later analysis to break), however allows
        # quicker generation of per-atom training sets
        # for glu/asp groups over a structure.
        # Training sets for supervised learning
        # classification can be created by setting the
        # 'makeTrainSet' input to True

        if makeTrainSet:
            clustAnalys = True
            inclOnlyGluAsp = True

        self.startTimer()
        self.printStepNumber()
        self.lgwrite(ln='Calculating electron density statistics per atom...')

        total = len(self.PDBarray)

        if parallel:
            # TODO: this would be great to implement at some point
            print 'Parallel processing not currently implemented!'
            pass
        else:

            self.densByRegion = []
            self.clustDoneOnAtm = []
            self.atomMidPts = [[np.nan]*3]

            for i, atom in enumerate(self.PDBarray):

                if inclOnlyGluAsp:

                    atmTypes = ['GLU-CD', 'GLU-OE1', 'GLU-OE2',
                                'ASP-OD1', 'ASP-OD2', 'ASP-CG']

                    tag = '-'.join(atom.getAtomID().split('-')[2:])
                    if tag not in atmTypes:
                        continue

                if showProgress:
                    sys.stdout.write('\r')
                    sys.stdout.write(
                        '{}%'.format(round(100*float(i)/total, 3)))
                    sys.stdout.flush()

                self.calcDensMetricsForAtom(atom=atom, plotDistn=plotDistn,
                                            clustAnalys=clustAnalys)

            if makeTrainSet:
                self.makeTrainingSet()

        self.success()
        self.stopTimer()

        # delete vxlsPerAtom since no longer needed
        del self.vxlsPerAtom

        for atom in self.PDBarray:
            atom.getAdditionalMetrics()

    def makeTrainingSet(self,
                        killNow=True, standardise=False):

        # make a training set of per-atom density values on
        # which a supervised-learning classifier could be trained.
        # This should not be included in a standard run of the code

        print 'Preparing classifier training dataset'

        if standardise:
            from sklearn.preprocessing import StandardScaler
            X = StandardScaler().fit_transform(self.densByRegion)
        else:
            X = self.densByRegion

        # get bfactors for atoms on which densByRegion is known
        bfactors = []
        for atmID in self.clustDoneOnAtm:
            for atm in self.PDBarray:
                if atm.getAtomID() == atmID:
                    bfactors.append(atm.Bfactor)
                    break

        # Write classification features to output file here
        f = lambda x: '{}clusterTrainingSet-{}.trset'.format(self.filesOut, x)
        i = 1
        while os.path.isfile(f(i)):
            i += 1
        print 'Writing calculated features to file: "{}"'.format(f(i))
        csvIn = open(f(i), 'w')
        for i, (atmID, dens) in enumerate(zip(self.clustDoneOnAtm, X)):
            csvIn.write(atmID+',' +
                        ','.join([str(np.round(d, 3)) for d in dens]) +
                        ',{}\n'.format(bfactors[i]))
        csvIn.close()

        if killNow:
            import sys
            sys.exit()

    def plotFCdistnPlot(self,
                        plot=True, atomOfInterest='',
                        atomsToPlot=['GLU-CD', 'CYS-SG'], atomFCvals=[],
                        atomFCvalsMaxNorm=[], FCatMin=[], reliability=[],
                        plotType='.png', axesFont=18):

        # plot a kde & histrogram distribution plot for the FCalc values for an
        # atom, both raw, and after being divided by the maximum FCalc value
        # attained for that atom (normalised-FCalc). The plot will also include
        # vertical lines indicating the FCalc and normalised-FCalc values
        # attained for the voxel where the most negative density map (not FC
        # map) voxel within the local region around the atom (this is the
        # voxel corresponding to the DLoss metric value).

        for tag in atomsToPlot:
            if tag in atomOfInterest.getAtomID():
                sns.set_style("dark")
                sns.set_context(rc={"figure.figsize": (10, 6)})
                fig = plt.figure()
                ax = plt.subplot(111)
                sns.distplot(np.array(atomFCvals), label='Fcalc')
                sns.distplot(np.array(atomFCvalsMaxNorm),
                             label='Fcalc/max(Fcalc)')
                ylims = ax.get_ylim()

                plt.plot((FCatMin, FCatMin),
                         (ylims[0], ylims[1]),
                         label='Fcalc, at position of min diff density')
                plt.plot((reliability, reliability),
                         (ylims[0], ylims[1]),
                         label='Fcalc/max(Fcalc), at ' +
                               'position of min diff density')
                leg = plt.legend(frameon=1)
                frame = leg.get_frame()
                frame.set_color('white')
                plt.xlabel('Per-voxel density map values', fontsize=axesFont)
                plt.ylabel('Normed-frequency', fontsize=axesFont)
                plt.title('Distribution of Fcalc density values: {}'.format(
                    atomOfInterest.getAtomID()))
                fig.savefig('{}testDistnPlot-{}{}'.format(
                    self.filesOut, atomOfInterest.getAtomID(), plotType))

    def plotDensScatterPlots(self,
                             printText=False, clustAnalys=False):

        # plot scatter plots for density metrics for
        # quick assessment of whether per-atom metrics
        # are behaving as expecting

        self.startTimer()
        self.fillerLine(style='line')
        self.lgwrite(
            ln='Plotting scatter plots for electron density statistics...',
            forcePrint=printText)

        plotVars = [['meandensity', 'maxdensity'],
                    ['meandensity', 'mediandensity'],
                    ['meandensity', 'mindensity'],
                    ['mindensity', 'maxdensity'],
                    ['meandensity', 'stddensity'],
                    ['mindensity', 'min90tile'],
                    ['maxdensity', 'max90tile'],
                    ['min90tile', 'min95tile'],
                    ['max90tile', 'max95tile'],
                    ['meandensity', 'meanPosOnly'],
                    ['meandensity', 'meanNegOnly'],
                    ['mindensity', 'meanNegOnly'],
                    ['maxdensity', 'meanPosOnly']]

        # # only include below if per-atom clusters are
        # # calculated - currently very slow
        if clustAnalys:
            plotVars += [['negClusterVal', 'meandensity'],
                         ['negClusterVal', 'mindensity'],
                         ['totDensShift', 'meandensity'],
                         ['totDensShift', 'mindensity']]

        if self.calcFCmap:
            plotVars.append(['meandensity', 'densityWeightedMean'])
            plotVars.append(['mindensity', 'densityWeightedMin'])
            plotVars.append(['maxdensity', 'densityWeightedMax'])
            plotVars.append(['maxdensity', 'densityWeightedMeanPosOnly'])
            plotVars.append(['mindensity', 'densityWeightedMeanNegOnly'])
            plotVars.append(['meanNegOnly', 'densityWeightedMeanNegOnly'])
            plotVars.append(['meanPosOnly', 'densityWeightedMeanPosOnly'])
            plotVars.append(
                ['densityWeightedMean', 'densityWeightedMeanPosOnly'])
            plotVars.append(
                ['densityWeightedMean', 'densityWeightedMeanNegOnly'])

        for pVars in plotVars:
            logStr = edens_scatter(outputDir=self.filesOut, metrics=pVars,
                                   PDBarray=self.PDBarray,
                                   pdbName=self.pdbName)
            self.lgwrite(ln=logStr)

    def pickleAtomList(self):

        # save list of atom objects to a .pkl file

        self.pklFileName = save_objectlist(self.PDBarray, self.pdbName)

    def startTimer(self):

        # start a timer

        self.timeStart = time.time()

    def stopTimer(self,
                  includeInLog=False):

        # stop a timer (must run startTimer before)

        elapsedTime = time.time() - self.timeStart
        if includeInLog:
            self.lgwrite(
                ln='section time: {}s\n'.format(round(elapsedTime, 3)))
        sys.stdout.flush()

    def success(self):

        # report success to log file

        self.lgwrite(ln='---> success')

    def fillerLine(self,
                   style='blank'):

        # print a filler line (several styles)
        # to command line

        if style == 'stars':
            ln = '\n***'
        elif style == 'line':
            ln = '\n'+'-'*30
        elif style == 'blank':
            ln = '\n'
        self.lgwrite(ln=ln)

    def lgwrite(self,
                ln='', strip=True, forcePrint=False):

        # write line to log file

        self.log.writeToLog(str=ln, strip=strip, forcePrint=forcePrint)

    def printStepNumber(self):

        # print a string indicating the current pipeline
        # step number directory to the command line

        try:
            self.stepNumber
        except AttributeError:
            self.stepNumber = 1
        self.lgwrite(ln='\n_______' +
                        '\nSTEP {})'.format(self.stepNumber))

        self.stepNumber += 1
