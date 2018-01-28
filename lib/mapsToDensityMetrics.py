from vxlsPerAtmAnalysisPlots import plotVxlsPerAtm, plotDensForAtm
from perAtomClusterAnalysis import perAtomXYZAnalysis
from densityAnalysisPlots import edens_scatter
from PDBFileManipulation import PDBtoList
from readMap import readMap
import matplotlib.pyplot as plt
from errors import error
import numpy as np
import sys
import time
import os
# try to import seaborn here. Don't worry if not (not essential for main run)
try:
    import seaborn as sns
except ImportError:
    pass

if sys.version_info[0] < 3:
    from itertools import izip as zip


class maps2DensMetrics():

    # assign values within a density map to specific atoms, using
    # the an atom-tagged map to determine which regions of space
    # are to be assigned to each atom

    def __init__(self,
                 filesIn='', filesOut='', pdbName='', atomTagMap='',
                 densityMap='', FCmap='',  plotScatter=False, plotHist=False,
                 logFile='./untitled.log', calcFCmap=True,
                 doXYZanalysis=False):

        # the input directory
        self.filesIn = filesIn

        # the output directory
        self.filesOut = filesOut

        # the pdb file name
        self.pdbName = pdbName

        # atom-tagged map name
        self.atomMapIn = atomTagMap

        # density map name (typically Fo-Fo)
        self.densMapIn = densityMap

        # FC map name
        self.FCmapIn = FCmap

        # (bool) plot scatter plots or not
        self.plotScatter = plotScatter

        # (bool) plot histogram plots or not
        self.plotHist = plotHist

        # log file name
        self.log = logFile

        # whether FC map should be generated
        self.calcFCmap = calcFCmap

        # whether to do analysis based on xyz of each voxel
        self.doXYZanalysis = doXYZanalysis

    def maps2atmdensity(self,
                        mapsAlreadyRead=False):

        # the map run method for this class. Will read in an atom-tagged map
        # and density map and assign density values for each individual atom
        # (as specified within the atom-tagged map). From these summary metrics
        # describing the density map behaviour in the vicinity of each refined
        # atom can be calculated

        if not mapsAlreadyRead:
            self.readPDBfile()
            self.readAtomMap()
        self.readDensityMap()
        self.reportDensMapInfo()
        self.checkMapCompatibility()

        if self.calcFCmap and not mapsAlreadyRead:
            self.readFCMap()
            self.reportDensMapInfo(mapType='calc')

        self.createVoxelList()

        if self.plotHist:
            self.plotDensHistPlots()

        self.calcDensMetrics(showProgress=False)

        if self.plotScatter:
            self.plotDensScatterPlots()

    def readPDBfile(self):

        # read in pdb file info here. A list of atom objects
        # is created, to which density metric information
        # will be added as additional attributes in the
        # methods included below

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln='Reading pdb file: {}'.format(self.pdbName))

        # read in the pdb file to fill list of atom objects
        self.PDBarray = PDBtoList('{}{}'.format(
            self.filesIn, self.pdbName))
        self.stopTimer()

        # make sure array of atoms ordered by atom number
        self.PDBarray.sort(key=lambda x: x.atomnum)

    def readAtomMap(self):

        # read in the atom-tagged map

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln='Reading atom-tagged map file...\n' +
                        'Atom map name: {}'.format(self.atomMapIn))

        self.atmmap, self.atomIndices = readMap(
            dirIn=self.filesIn, dirOut=self.filesOut, mapName=self.atomMapIn,
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
        self.lgwrite(ln='Reading density map file...\n' +
                        'Density map name: {}'.format(self.densMapIn))

        self.densmap = readMap(dirIn=self.filesIn, dirOut=self.filesOut,
                               mapName=self.densMapIn, mapType='density_map',
                               atomInds=self.atomIndices, log=self.log)
        self.stopTimer()

    def readFCMap(self):

        # read in the FC (calculated structure factor) density map.
        # This method should not be called if no FC density map
        # has been provided in the current run.

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln='Reading Fcalc density map file...\n' +
                        'Density map name: {}'.format(self.FCmapIn))

        self.FCmap = readMap(dirIn=self.filesIn, dirOut=self.filesOut,
                             mapName=self.FCmapIn, mapType='density_map',
                             atomInds=self.atomIndices, log=self.log)

        self.stopTimer()

    def reportDensMapInfo(self,
                          numSfs=4, mapType='density'):

        # report the density map summary information to a log file

        if mapType == 'density':
            mp = self.densmap
        elif mapType == 'calc':
            mp = self.FCmap

        totalNumVxls = np.product(self.atmmap.nxyz.values())
        structureNumVxls = len(mp.vxls_val)
        totalMean = mp.density['mean']
        structureMean = np.mean(mp.vxls_val)
        solvNumVxls = totalNumVxls - structureNumVxls
        solvMean = (totalNumVxls*totalMean -
                    structureNumVxls*structureMean)/solvNumVxls

        self.lgwrite(
            ln='\nFor voxels assigned to structure:\n' +
               '\tmean structure density : {}\n'.format(
                round(structureMean, numSfs)) +
               '\tmax structure density : {}\n'.format(
                round(max(mp.vxls_val), numSfs)) +
               '\tmin structure density : {}\n'.format(
                round(min(mp.vxls_val), numSfs)) +
               '\tstd structure density : {}\n'.format(
                round(np.std(mp.vxls_val), numSfs)) +
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
                        inclOnlyGluAsp=False):

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

        self.vxlsPerAtom = vxlDic

        # The following is not essential for run (TODO maybe omit)
        if self.doXYZanalysis:
            xyz_list = self.densmap.getVoxXYZ(
                self.atomIndices, coordType='fractional')

            for atm, xyz in zip(self.atmmap.vxls_val, xyz_list):
                xyzDic[atm].append(xyz)

            # get the mid points for each atom from the set of voxels
            # per atom, whilst accounting for symmetry (the asym unit
            # may not contain 1 single whole molecule, but split up)
            xyzDic2 = {}
            for atom in self.PDBarray:

                # this is more of testing reasons that any clear use
                if inclOnlyGluAsp:
                    atmTypes = ['GLU-CD', 'GLU-OE1', 'GLU-OE2',
                                'ASP-OD1', 'ASP-OD2', 'ASP-CG',
                                'CYS-SG', 'CYS-CB', 'CYS-CA',
                                'MET-SD', 'MET-CE', 'MET-CG']
                    tag = '-'.join(atom.getAtomID().split('-')[2:])
                    if tag not in atmTypes:
                        continue

                xyzAnalysis = perAtomXYZAnalysis(
                    atomObj=atom, vxlRefPoint=np.mean(xyz_list, 0),
                    densPerVxl=np.round(np.array(vxlDic[atom.atomnum]), 3),
                    xyzsPerAtom=xyzDic[atom.atomnum], densMapObj=self.densmap)
                xyzAnalysis.getxyzPerAtom()
                atom.vxlMidPt = xyzAnalysis.findVoxelMidPt()
                xyzDic2[atom.getAtomID()] = xyzAnalysis.keptPts
                self.xyzsPerAtom = xyzDic2

        if self.calcFCmap:
            vxlDic2 = {atm: [] for atm in self.atmmap.vxls_val}
            for atm, dens in zip(self.atmmap.vxls_val, self.FCmap.vxls_val):
                vxlDic2[atm].append(dens)
            self.FCperAtom = vxlDic2

        self.deleteMapsAttributes()
        self.stopTimer()

    def deleteMapsAttributes(self):

        # Provide the option to delete atmmap and
        # densmap attributes to save memory, if
        # they are no longer needed during a run

        # del self.atmmap
        # if self.calcFCmap:
        #     del self.FCmap
        del self.densmap.vxls_val

    def plotDensHistPlots(self,
                          getVoxelStats=False, perAtmDensHist=False):

        # Create histogram or kde plots of number of voxels per atom

        self.startTimer()
        self.printStepNumber()
        self.lgwrite(ln='Plotting histogram plots of voxels per atom...\n' +
                        'Plots written to "{}plots"'.format(self.filesOut))

        stats = plotVxlsPerAtm(pdbName=self.pdbName, where=self.filesOut,
                               vxlsPerAtom=self.vxlsPerAtom, plotType='both',
                               returnStats=getVoxelStats)

        if stats != '':
            print('mean: {}\nstd: {}\nmax: {}\nmin: {}'.format(*stats))

        if perAtmDensHist:
            plotDensForAtm(pdbName=self.pdbName, where=self.filesOut,
                           vxlsPerAtom=self.vxlsPerAtom, plotType='both',
                           PDBarray=self.PDBarray)

        self.stopTimer()

    def calcDensMetricsForAtom(self,
                               atom=[], plotDistn=False):

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
                # if the user has opted to calculate an Fcalc map in addition
                # to the difference map, then additional metrics can be
                # derived using this map. These metrics typically use the Fcalc
                # map density at each voxel to weight the contribution that
                # each voxel's difference map value should play when
                # calculating damage metrics. Effectively, a voxel far from an
                # atom (but still included in the search radius around that
                # atom) should not contribute to a damage indicator as much as
                # a voxel close to the atomic centre

                atomFCvals = self.FCperAtom[atom.atomnum]
                # NOTE: currently set all negative values to zero. This has
                # effect of ignoring Fcalc density that is less than the map
                # mean. This is implemented such that all per-voxel weights
                # (see below) are positive and so therefore sensible
                # weighted-means can be calculated. This may need to be
                # reconsidered for future use!
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
                    # typically only to be used for testing purposes
                    self.plotFCdistnPlot(
                        atomsToPlot=['GLU-CD', 'CYS-SG'], atomOfInterest=atom,
                        atomFCvals=atomFCvals, FCatMin=atomFCvals[minIndex],
                        atomFCvalsMaxNorm=atomFCvalsMaxNormed)

            if self.doXYZanalysis:
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

                clustAnalysis = perAtomXYZAnalysis(
                    atomObj=atom, vxlMidPt=atom.vxlMidPt,
                    knownRefPoint=self.knownRefPt1,
                    knownRefPoint2=self.knownRefPt2)

                clustAnalysis.keptPts = self.xyzsPerAtom[atom.getAtomID()]
                clustAnalysis.partitionPtsByVec()

                # atom.negClusterVal = clustAnalysis.topNegClustMean
                # atom.totDensShift = clustAnalysis.netDensShift

                self.densByRegion.append(clustAnalysis.densByRegion)

    def calcDensMetrics(self,
                        plotDistn=False, showProgress=True, parallel=False,
                        makeTrainSet=False, inclOnlyGluAsp=False,
                        doRandomSubset=False):

        # determine density summary metrics per atom. 'includeOnlyGluAsp'
        # allows calculations to be performed only for Glu/asp carboxylates
        # (this is not typically suitable and will cause later analysis to
        # break), however allows quicker generation of per-atom training sets
        # for glu/asp groups over a structure. Training sets for supervised
        # learning classification can be created by setting the 'makeTrainSet'
        # input to True

        if makeTrainSet:
            self.doXYZanalysis = True
            inclOnlyGluAsp = True

        self.startTimer()
        self.printStepNumber()
        self.lgwrite(ln='Calculating electron density statistics per atom...')

        total = len(self.PDBarray)

        if parallel:
            # TODO: this would be great to implement at some point
            print('Parallel processing not currently implemented!')
            pass
        else:

            self.densByRegion = []
            self.clustDoneOnAtm = []

            for i, atom in enumerate(self.PDBarray):

                # this is more of testing reasons that any clear use
                if inclOnlyGluAsp:

                    atmTypes = ['GLU-CD', 'GLU-OE1', 'GLU-OE2',
                                'ASP-OD1', 'ASP-OD2', 'ASP-CG',
                                'CYS-SG', 'MET-SD']

                    tag = '-'.join(atom.getAtomID().split('-')[2:])
                    if tag not in atmTypes:
                        continue

                    if self.doXYZanalysis:

                        num = '-'.join(atom.getAtomID().split('-')[:2])

                        if tag == 'GLU-CD':
                            lookFor1 = num+'-GLU-OE1'
                            lookFor2 = num+'-GLU-OE2'
                        elif tag == 'GLU-OE1':
                            lookFor1 = num+'-GLU-CD'
                            lookFor2 = num+'-GLU-OE2'
                        elif tag == 'GLU-OE2':
                            lookFor1 = num+'-GLU-CD'
                            lookFor2 = num+'-GLU-OE1'
                        elif tag == 'ASP-CG':
                            lookFor1 = num+'-ASP-OD1'
                            lookFor2 = num+'-ASP-OD2'
                        elif tag == 'ASP-OD1':
                            lookFor1 = num+'-ASP-CG'
                            lookFor2 = num+'-ASP-OD2'
                        elif tag == 'ASP-OD2':
                            lookFor1 = num+'-ASP-CG'
                            lookFor2 = num+'-ASP-OD1'
                        elif tag == 'CYS-SG':
                            lookFor1 = num+'-CYS-CB'
                            lookFor2 = num+'-CYS-CA'
                        elif tag == 'MET-SD':
                            lookFor1 = num+'-MET-CE'
                            lookFor2 = num+'-MET-CG'

                        if i < 10:
                            srt = 0
                        else:
                            srt = i-10
                        if i < len(self.PDBarray)-10:
                            stp = i + 10
                        else:
                            stp = len(self.PDBarray)

                        for atm2 in self.PDBarray[srt:stp]:
                            if atm2.getAtomID() == lookFor1:
                                self.knownRefPt1 = atm2.vxlMidPt
                            elif atm2.getAtomID() == lookFor2:
                                self.knownRefPt2 = atm2.vxlMidPt

                # only calculate metrics for a random subset of atoms
                # - for testing purposes
                if doRandomSubset:
                    import random
                    if random.uniform(0, 1) > 0.01:
                        continue
                    else:
                        print('Random atom used: ' + atom.getAtomID())

                if showProgress:
                    sys.stdout.write('\r')
                    sys.stdout.write(
                        '{}%'.format(round(100*float(i)/total, 3)))
                    sys.stdout.flush()

                self.calcDensMetricsForAtom(atom=atom, plotDistn=plotDistn)
                atom.getAdditionalMetrics()

            if makeTrainSet:
                self.makeTrainingSet()

        self.success()
        self.stopTimer()

        # delete vxlsPerAtom since no longer needed
        del self.vxlsPerAtom

        # ############################################################################
        # # TEST: cluster the density values per atom based off xyz.
        # # KEEP THIS COMMENTED WHEN USING THE CODE
        # from sklearn.cluster import KMeans
        # from sklearn.decomposition import PCA

        # d = self.densByRegion

        # numClusts = 5
        # reduced_data = PCA(n_components=2).fit_transform(d)
        # kmeans = KMeans(init='k-means++', n_clusters=numClusts)
        # kmeans.fit(reduced_data)

        # # Step size of the mesh. Decrease to increase the quality of the VQ.
        # h = .02     # point in the mesh [x_min, x_max]x[y_min, y_max].

        # # Plot decision boundary. For that, we will assign a color to each
        # x_min = reduced_data[:, 0].min()-0.5*np.abs(reduced_data[:, 0].min())
        # x_max = reduced_data[:, 0].max()+0.5*np.abs(reduced_data[:, 0].max())
        # y_min = reduced_data[:, 1].min()-0.5*np.abs(reduced_data[:, 1].min())
        # y_max = reduced_data[:, 1].max()+0.5*np.abs(reduced_data[:, 1].max())

        # xx, yy = np.meshgrid(
        #     np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

        # # Obtain labels for each point in mesh. Use last trained model.
        # Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

        # Zatoms = kmeans.predict(reduced_data)
        # atmNames = [x for _, x in sorted(zip(Zatoms, self.clustDoneOnAtm))]
        # Zatoms.sort()
        # for Za, atom in zip(Zatoms, atmNames):
        #     print('{} --> {}'.format(atom, Za))

        # # Put the result into a color plot
        # Z = Z.reshape(xx.shape)
        # plt.figure(1)
        # plt.clf()
        # plt.imshow(Z, interpolation='nearest',
        #            extent=(xx.min(), xx.max(), yy.min(), yy.max()),
        #            cmap=plt.cm.Paired,
        #            aspect='auto', origin='lower')

        # plt.plot(reduced_data[:, 0], reduced_data[:, 1], 'k.', markersize=2)
        # # Plot the centroids as a white X
        # centroids = kmeans.cluster_centers_
        # plt.scatter(centroids[:, 0], centroids[:, 1],
        #             marker='x', s=169, linewidths=3,
        #             color='w', zorder=10)

        # import pylab as pl
        # for i in range(numClusts):
        #     pl.text(centroids[i, 0], centroids[i, 1],
        #             str(i), color="white", fontsize=20)

        # plt.title('K-means clustering on per-atom density (PCA-reduced data)\n'
        #           'Centroids are marked with white cross')
        # plt.xlim(x_min, x_max)
        # plt.ylim(y_min, y_max)
        # plt.xticks(())
        # plt.yticks(())
        # plt.show()
        # import sys
        # sys.exit()
        # ############################################################################

    def makeTrainingSet(self,
                        killNow=True, standardise=False):

        # make a training set of per-atom density values on
        # which a supervised-learning classifier could be trained.
        # NOTE: This should NOT be included in a standard run

        print('Preparing classifier training dataset')

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
        print('Writing calculated features to file: "{}"'.format(f(i)))
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
                        atomFCvalsMaxNorm=[], FCatMin=[],
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
