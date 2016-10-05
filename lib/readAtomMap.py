from residueFormatter import densper_resatom_NOresidueclass,densper_res
from vxlsPerAtmAnalysisPlots import plotVxlsPerAtm,plotDensForAtm
from densityAnalysisPlots import edens_scatter
from deleteListIndices import multi_delete 
from savevariables import save_objectlist
from itertools import izip as zip, count
from PDBFileManipulation import PDBtoList
from map2VoxelClassList import readMap
import matplotlib.pyplot as plt 
from progbar import progress
from logFile import logFile
import numpy as np
import sys  
import time
import math

# check if seaborn library is present and include if so
from checkSeabornPresent import checkSeabornPresent as checkSns
seabornFound = checkSns()
if seabornFound is True:
    import seaborn as sns

class voxel_density():

    # A class for .map file voxel

    def __init__(self,
                 density = 0,
                 atmnum  = 0):

        self.density = density
        self.atmnum  = atmnum

class vxls_of_atm():
    
    # A class for collecting voxels per atom

    def __init__(self,
                 vxls   = [],
                 atmnum = 0):

        self.vxls   = vxls
        self.atmnum = atmnum

class maps2DensMetrics():

    # assign values within a density map to specific atoms, using
    # the an atom-tagged map to determine which regions of space
    # are to be assigned to each atom

    def __init__(self,
                 filesIn     = '',
                 filesOut    = '',
                 pdbName     = '',
                 atomTagMap  = '',
                 densityMap  = '',
                 FCmap       = '',
                 plotScatter = False,
                 plotHist    = False,
                 plotBar     = False,
                 logFile     = logFile,
                 calcFCmap   = True):

        self.filesIn     = filesIn     # the input directory
        self.filesOut    = filesOut    # the output directory
        self.pdbName     = pdbName     # the pdb file name
        self.map1        = atomTagMap  # atom-tagged map
        self.map2        = densityMap  # density map (typically Fo-Fo)
        self.map3        = FCmap       # FC map
        self.plotScatter = plotScatter # (bool) plot scatter plots or not
        self.plotHist    = plotHist    # (bool) plot histogram plots or not
        self.plotBar     = plotBar     # (bool) plot bar plots or not
        self.log         = logFile
        self.calcFCmap   = calcFCmap   # whether FC map should be generated

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

        if self.calcFCmap is True:
            self.readFCMap()
        self.createVoxelList(useFCmap = self.calcFCmap)

        if self.plotHist is True:
            self.plotDensHistPlots()
        self.calculateDensMetrics(FCweighted = self.calcFCmap)

        if self.plotScatter is True:
            self.plotDensScatterPlots()

        if self.plotBar is True:
            self.plotPerResidueBoxPlots()
        self.pickleAtomList()

    def readPDBfile(self):

        # read in pdb file info here

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln = 'Reading pdb file: {}'.format(self.pdbName))

        # read in the pdb file to fill list of atom objects
        self.PDBarray = PDBtoList('{}{}.pdb'.format(self.filesIn,self.pdbName))
        self.stopTimer()

        # make sure array of atoms ordered by atom number
        self.PDBarray.sort(key = lambda x: x.atomnum)
           
        # need to get VDW radius for each atom:
        for atom in self.PDBarray:
            atom.VDW_get()  
        
    def readAtomMap(self):

        # read in the atom map

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln = 'Reading atom-tagged map file...')
        self.lgwrite(ln = 'Atom map name: {}'.format(self.map1))

        self.atmmap,self.atomIndices = readMap(dirIn   = self.filesIn,
                                               dirOut  = self.filesOut,
                                               mapName = self.map1,
                                               mapType = 'atom_map',
                                               log     = self.log)  
        self.stopTimer()

        # find number of atoms in structure
        num_atoms = len(self.PDBarray)

        # find atom numbers present in list (repeated atom numbers removed)
        seen = set()
        seen_add = seen.add
        uniq_atms = [x for x in self.atmmap.vxls_val if not (x in seen or seen_add(x))] 
        
        # find set of atoms numbers not present (i.e atoms not assigned to voxels)
        Atms_notpres = set(range(1,num_atoms+1)) - set(uniq_atms)
        self.lgwrite(ln = 'Number of atoms not assigned to voxels: {}'.format(len(Atms_notpres)))

    def readDensityMap(self):

        # read in the density map

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln = 'Reading density map file...')
        self.lgwrite(ln = 'Density map name: {}'.format(self.map2))
        
        self.densmap = readMap(dirIn    = self.filesIn,
                               dirOut   = self.filesOut,
                               mapName  = self.map2,
                               mapType  = 'density_map',
                               atomInds = self.atomIndices,
                               log      = self.log)  
        self.stopTimer()

    def readFCMap(self):

        # read in the FC (calculated structure factor) density map

        self.printStepNumber()
        self.startTimer()
        self.lgwrite(ln = 'Reading Fcalc density map file...')
        self.lgwrite(ln = 'Density map name: {}'.format(self.map3))

        self.FCmap = readMap(dirIn    = self.filesIn,
                             dirOut   = self.filesOut,
                             mapName  = self.map3,
                             mapType  = 'density_map',
                             atomInds = self.atomIndices,
                             log      = self.log)  

        self.stopTimer()

    def reportDensMapInfo(self):

        # print density map summary information to command line

        totalNumVxls     = np.product(self.atmmap.nxyz.values())
        structureNumVxls = len(self.densmap.vxls_val)
        totalMean        = self.densmap.density['mean']
        structureMean    = np.mean(self.densmap.vxls_val)
        solvNumVxls      = totalNumVxls - structureNumVxls
        solvMean         = (totalNumVxls*totalMean - structureNumVxls*structureMean)/solvNumVxls

        txt = '\nFor voxels assigned to structure:\n'+\
              '\tmean structure density : {}\n'.format(round(structureMean,4))+\
              '\tmax structure density : {}\n'.format(round(max(self.densmap.vxls_val),4))+\
              '\tmin structure density : {}\n'.format(round(min(self.densmap.vxls_val),4))+\
              '\tstd structure density : {}\n'.format(round(np.std(self.densmap.vxls_val),4))+\
              '\t# voxels included : {}\n'.format(structureNumVxls)+\
              '\nFor voxels assigned to solvent:\n'+\
              '\tmean solvent-region density : {}\n'.format(round(solvMean),4)+\
              '\t# voxels included : {}'.format(solvNumVxls)
        self.lgwrite(ln = txt)

    def checkMapCompatibility(self):

        # check that atom-tagged and density map 
        # can be combined successfully

        self.printStepNumber()
        self.lgwrite(ln = 'Checking that maps have same dimensions and sampling properties...' )

        self.startTimer()
        # Check that the maps have the same dimensions, grid sampling,..
        if (self.atmmap.axis != self.densmap.axis or
            self.atmmap.gridsamp != self.densmap.gridsamp or 
            self.atmmap.start != self.densmap.start or
            self.atmmap.nxyz != self.densmap.nxyz or
            self.atmmap.type != self.densmap.type):

            self.lgwrite(ln = 'Incompatible map properties --> terminating script')
            sys.exit()

        elif self.atmmap.celldims != self.densmap.celldims:
            self.lgwrite(ln = 'Not exact same map grid dimensions..')

            # now check if grid dims same to a specific dp and consider continuing
            stop = True
            for i in list(reversed(range(7))):
                count = 0
                for key in self.atmmap.celldims.keys():
                    if np.round(self.atmmap.celldims[key],i) == np.round(self.densmap.celldims[key],i):
                        count += 1
                if count == 6:
                    str = 'Map grid dimensions same to {}dp\n'.format(i)+\
                          '--> continuing with processing anyway'
                    self.lgwrite(ln = str)
                    stop = False
                    break
            if stop == True:
                err = 'Map grid dimensions still not same to 0dp\n'+\
                      ' --> terminating script'
                self.lgwrite(ln = err)
                sys.exit()

        else:
            self.success()
            self.lgwrite(ln = 'The atom and density map are of compatible format!')
        self.stopTimer()

        str = 'Total number of voxels assigned to atoms: {}'.format(len(self.atmmap.vxls_val))
        self.lgwrite(ln = str)

    def createVoxelList(self,
                        useFCmap = True):

        # create dictionary of voxels with atom numbers as keys 

        self.startTimer()
        self.printStepNumber()
        self.lgwrite(ln = 'Combining voxel density and atom values...')
        self.success()
        vxlDic = {atm:[] for atm in self.atmmap.vxls_val}
        xyzDic = {atm:[] for atm in self.atmmap.vxls_val}

        self.densmap.reshape1dTo3d()
        self.densmap.abs2xyz_params()
        for atm,dens in zip(self.atmmap.vxls_val,self.densmap.vxls_val):
            vxlDic[atm].append(dens)

        xyz_list = self.densmap.getVoxXYZ(self.atomIndices,coordType = 'fractional')
        for atm,xyz in zip(self.atmmap.vxls_val,xyz_list):
            xyzDic[atm].append(xyz)

        self.vxlsPerAtom = vxlDic
        self.xyzsPerAtom = xyzDic # not essential for run

        if useFCmap:
            vxlDic2 = {atm:[] for atm in self.atmmap.vxls_val}
            for atm,dens in zip(self.atmmap.vxls_val,self.FCmap.vxls_val):
                vxlDic2[atm].append(dens)
            self.FCperAtom = vxlDic2

        # delete atmmap and densmap now to save memory
        del self.atmmap
        del self.FCmap
        del self.densmap.vxls_val
        self.stopTimer()

    def plotDensHistPlots(self,
                          getVoxelStats  = False,
                          perAtmDensHist = False):

        # histogram & kde plots of number of voxels per atom

        self.startTimer()
        self.printStepNumber()
        self.lgwrite(ln = 'Plotting histogram plots of voxels per atom...')
        self.lgwrite(ln = 'Plots written to "{}plots"'.format(self.filesOut))

        stats = plotVxlsPerAtm(pdbName     = self.pdbName,
                               where       = self.filesOut,
                               vxlsPerAtom = self.vxlsPerAtom,
                               plotType    = 'both',
                               returnStats = getVoxelStats)

        if stats != '':
            print 'mean: {}\nstd: {}\nmax: {}\nmin: {}'.format(*stats)

        if perAtmDensHist:
            plotDensForAtm(pdbName     = self.pdbName,
                           where       = self.filesOut,
                           vxlsPerAtom = self.vxlsPerAtom,
                           plotType    = 'both',
                           PDBarray    = self.PDBarray)

        self.stopTimer()

    def calculateDensMetrics(self,
                             FCweighted   = True,
                             plotDistn    = False,
                             plotxyzScat  = False,
                             showProgress = True):

        # determine density summary metrics per atom

        self.startTimer()
        self.printStepNumber()
        self.lgwrite(ln = 'Calculating electron density statistics per atom...')

        total = len(self.PDBarray)
        for i,atom in enumerate(self.PDBarray):

            if showProgress:
                sys.stdout.write('\r')
                sys.stdout.write('{}%'.format(round(100*float(i)/total,3)))
                sys.stdout.flush()

            try:
                atomVxls = self.vxlsPerAtom[atom.atomnum]
            except KeyError:
                print 'Warning!: No voxels assigned to an atom. Consider increasing '+\
                      'per-atom search radius parameter in RIDL input .txt file.'
                atomVxls = [np.nan]

            if FCweighted is True:
                # calculate reliability measures based on electron 
                # density probability at position of min density 
                atomFCvals = self.FCperAtom[atom.atomnum]
                atomFCvalsMaxNormalised = np.array(atomFCvals)/max(atomFCvals)

                minIndex     = np.array(atomVxls).argmin()
                reliability  = atomFCvalsMaxNormalised[minIndex]
                FCatMin      = atomFCvals[minIndex]
                weightedVxls = np.multiply(atomVxls,atomFCvalsMaxNormalised)

            if len(atomVxls) != 0:
                atom.meandensity   = np.mean(atomVxls)
                atom.mediandensity = np.median(atomVxls)
                atom.mindensity    = min(atomVxls)
                atom.maxdensity    = max(atomVxls)
                atom.stddensity    = np.std(atomVxls)
                atom.min90tile     = np.percentile(atomVxls,10)
                atom.max90tile     = np.percentile(atomVxls,90)
                atom.min95tile     = np.percentile(atomVxls,5)
                atom.max95tile     = np.percentile(atomVxls,95)
                atom.numvoxels     = len(atomVxls)

                if FCweighted is True:
                    atom.reliability  = reliability
                    atom.wMean        = np.mean(weightedVxls)

                    if plotDistn is True:
                        self.plotFCdistnPlot(atomOfInterest    = atom,
                                             atomFCvals        = atomFCvals,
                                             atomFCvalsMaxNorm = atomFCvalsMaxNormalised,
                                             FCatMin           = FCatMin,
                                             reliability       = reliability)

                if plotxyzScat:
                    atmId = atom.getAtomID()
                    # if 'B-17-GLU-N' in atmId or 'GLU-OE1' in atmId or 'GLU-OE2' in atmId or 'GLU-CD' in atmId:
                    if 'CYS-SG' in atmId:
                        clustAnalysis = self.plotPerAtomVoxelsScatter(atmNum = atom.atomnum,
                                                                      atmId  = atmId,
                                                                      atmXYZ = atom.getXYZ())
                        atom.negCluster  = clustAnalysis[0]

                # # test here - remove after use...
                # fig = plt.figure()
                # ax = plt.subplot(111)
                # # plt.plot(atomFCvals,'r',label = 'FC')
                # # plt.plot(atomVxls,'b',label = 'Fo-Fo')
                # # plt.xlabel('Order')
                # # plt.ylabel('Voxel value')
                # # plt.legend()
                # plt.scatter(atomFCvalsMaxNormalised,atomVxls)
                # xRange = [min(atomFCvalsMaxNormalised),max(atomFCvalsMaxNormalised)]
                # plt.plot(xRange,[atom.wMean]*2,'r')
                # plt.plot(xRange,[atom.meandensity]*2,'g')
                # plt.xlabel('normed-FC')
                # plt.ylabel('Fo-Fo')
                # plt.title('Compare Fc and Fo-Fo per atom: {}'.format(atom.getAtomID()))
                # fig.savefig('{}FcAndFoFoPlot-{}.png'.format(self.filesOut,atom.getAtomID()))

        self.success()
        self.stopTimer()

        # delete the vxlsPerAtom list now to save memory
        del self.vxlsPerAtom

        # get additional metrics per atom
        for atom in self.PDBarray:
            atom.getAdditionalMetrics()

    def plotPerAtomVoxelsScatter(self,
                                 atmNum        = 0,
                                 atmId         = '',
                                 atmXYZ        = [0.0,0.0,0.0],
                                 showDens      = True,
                                 showSymm      = False,
                                 test          = 'insulin',
                                 makePlot      = False,
                                 report        = True,
                                 findClusters  = True,
                                 plotCentroids = False,
                                 numPosClusts  = 10,
                                 mainVecOnly   = False):

        # plot a per-atom scatter plot of xyz positions of each
        # voxel assigned to a specfied atom in 3d space        
        # BIT BELOW IS STILL UNFINISHED AND SHOULD STAY COMMENTED

        if atmNum not in self.xyzsPerAtom.keys():
            return

        makePlot = True
        # if 'CYS-SG' in atmId or 'B-17-GLU-N' in atmId or 'GLU-OE1' in atmId or 'GLU-OE2' in atmId or 'GLU-CD' in atmId:
        # # if 'CYS-SG' in atmId:
        #     makePlot = False
        # else: 
        #     return

        k = atmNum
        X = np.array([xyz[0] for xyz in self.xyzsPerAtom[k]])
        Y = np.array([xyz[1] for xyz in self.xyzsPerAtom[k]])
        Z = np.array([xyz[2] for xyz in self.xyzsPerAtom[k]])

        # keep densities here to 3dp since sym copies may have 
        # differing densities to higher accuracy than this
        density = np.round(np.array(self.vxlsPerAtom[k]),3) 

        symOps = self.densmap.getSymOps(X,Y,Z)

        # for any sym op that takes a point outside of unit cell,
        # move translate this back into the unit cell of interest
        symOps = np.mod(symOps,1).tolist()

        # for case where original xyz points are close to a unit
        # cell boundary, include additional adjacent unit cells
        translateUnitCell = []
        threshold = 0.01
        for i in range(len(symOps[0][0])):
            x = symOps[0][0][i]
            y = symOps[0][1][i]
            z = symOps[0][2][i]
            trans = []
            if x < threshold:
                trans = [[-1],[0],[0]]
            if x > 1-threshold:
                trans = [[1],[0],[0]]
            if y < threshold:
                trans = [[0],[-1],[0]]
            if y > 1-threshold:
                trans = [[0],[1],[0]]
            if z < threshold:
                trans = [[0],[0],[-1]]
            if z > 1-threshold:
                trans = [[0],[0],[1]]
            if trans != []:
                if trans not in translateUnitCell:
                    translateUnitCell.append(trans)

        l = len(translateUnitCell)
        if l != 0:
            for t in translateUnitCell:
                for sym in np.array(symOps) + t:
                    symOps.append(sym)
            if report:
                print 'Necessary to consider {} adjacent unit cell(s)'.format(l)

        # find suitable reference point near centre of asym unit
        centreOfInterest = np.mean(symOps[0], axis = 1)
        originalXYZs = np.array([[symOps[0][j][ii] for j in range(3)] for ii in range(len(X))])
        dists = list(np.linalg.norm(originalXYZs - np.array(centreOfInterest),axis = 1))
        del originalXYZs
        minInd = dists.index(min(dists))
        refPoint = [X[minInd],Y[minInd],Z[minInd]] + [density[minInd]]

        keptPoints = [refPoint]
        symOrOriginal = [0]

        samexyzs = 0
        dupcount = 0
        for ii in range(len(X)):

            symxyzs   = [[symOp[j][ii] for j in range(3)] for symOp in symOps]
            symPoints = [sym + [density[ii]] for sym in symxyzs]
            copy = False
            for i,symPt in enumerate(symPoints):
                if symPt in keptPoints:
                    copy = True
                    dupcount += 1

                    if i == 0:
                        symOrOriginal[keptPoints.index(symPt)] = 0
                    break
            if copy:
                continue

            distToFound = []
            for pt in keptPoints:
                distToFound.append(np.linalg.norm(np.array(symxyzs) - np.array(pt[:3]),axis = 1))
            
            minDist    = list(np.min(distToFound, axis = 0))
            closestInd = minDist.index(min(minDist))
            closestSym = symPoints[closestInd]

            if min(minDist) < 1e-10:
                samexyzs += 1 
                # makePlot = True

                # if new point is part of original asym unit, use this point instead
                if closestInd == 0:
                    dists = [dist[0] for dist in distToFound]
                    nearestPtId = dists.index(min(minDist))
                    keptPoints.pop(nearestPtId)
                    symOrOriginal.pop(nearestPtId)
                else:
                    continue

            keptPoints.append(closestSym)

            if closestInd != 0:
                symOrOriginal.append(1)
            else:
                symOrOriginal.append(0)

        if report:
            print 'Atom Identity: {}'.format(atmId)
            print 'Number points found: {}'.format(len(keptPoints))
            print '\tNumber that are original: {}'.format(len(symOrOriginal)-np.sum(symOrOriginal))
            print '\tNumber that are sym copies: {}'.format(np.sum(symOrOriginal))
            print '\tNumber of duplicate positions: {}'.format(dupcount+samexyzs)
            print '\t\tidentical: {}'.format(dupcount)
            print '\t\tsame xyz: {}'.format(samexyzs)

        xs = [xyz[0] for xyz in keptPoints]
        ys = [xyz[1] for xyz in keptPoints]
        zs = [xyz[2] for xyz in keptPoints]
        ds = [xyz[3] for xyz in keptPoints]

        # cluster analysis begins here
        if findClusters:
            from sklearn.cluster import KMeans

            posVals = [pt for pt in keptPoints if pt[3] >= 0]
            negVals = [pt for pt in keptPoints if pt[3] < 0]

            # find negative point clusters
            if len(negVals) > 1:
                negClustering = True
                negKmeans = KMeans(n_clusters = 2)
                negClusters = negKmeans.fit_predict(negVals)
                negCluster1 = [val for i,val in enumerate(negVals) if negClusters[i] == 0]
                negCluster2 = [val for i,val in enumerate(negVals) if negClusters[i] == 1]
                negMean1    = np.mean([c[3] for c in negCluster1])
                negMean2    = np.mean([c[3] for c in negCluster2])
                negCentroids = negKmeans.cluster_centers_
                if negMean2 > negMean1:
                    negClustOfInterest = negCluster1
                    mainNegCentroid = negCentroids[0, :]
                    negClusters = np.mod(np.array(negClusters)-1,2)
                else:
                    negClustOfInterest = negCluster2
                    mainNegCentroid = negCentroids[1, :]   
            else:
                negClustering = False
                negClusters = [1]

            # find positive point clusters
            if len(posVals) > 1:

                # reduce number of clusters if too many specified
                if numPosClusts > len(posVals):
                    numPosClusts = len(posVals) 

                posClustering = True
                posKmeans = KMeans(n_clusters = numPosClusts)
                posClusters = posKmeans.fit_predict([pt[:3] for pt in posVals])
                posClustList = []
                for j in range(numPosClusts):
                    posClustList.append([val for i,val in enumerate(posVals) if posClusters[i] == j])      
                posClustMeans = [np.mean([c[3] for c in clust]) for clust in posClustList]
                posCentroids = posKmeans.cluster_centers_
                maxPosInd = posClustMeans.index(max(posClustMeans))
                posClusterOfInterest = posClustList[maxPosInd]
            else:
                posClustering = False
                posClusters = [1]

            if posClustering and negClustering:
                # find dominant direction of density shift from neg to pos density

                if mainVecOnly:
                    clusterLoop = [posClusterOfInterest]
                else:
                    clusterLoop = posClustList

                resultantVectors = []
                for clust in clusterLoop:
                    quiverPts  = { i : [] for i in ['x','y','z','u','v','w'] }
                    for negPt in negClustOfInterest:
                        densDiffs  = []
                        directions = []
                        for posPt in clust:
                            densDiffs.append(posPt[3]-negPt[3])
                            vectorChange = np.array(posPt[:3])-np.array(negPt[:3])
                            directions.append(vectorChange/np.linalg.norm(vectorChange))
                        totalComponent = np.dot(densDiffs,directions)
                        for n,xyz in enumerate(['x','y','z']):
                            quiverPts[xyz].append(negPt[n])
                        for n,uvw in enumerate(['u','v','w']):
                            quiverPts[uvw].append(totalComponent[n])
                    resultantVectors.append(np.mean([quiverPts[uvw] for uvw in ['u','v','w']], axis = 1))
                vectorMagnitudes = [np.linalg.norm(vec) for vec in resultantVectors]
                normedVectors    = [vec/length for vec,length in zip(resultantVectors,vectorMagnitudes)]
                # weightedVectors  = [0.03*vec/max(vectorMagnitudes) for vec in resultantVectors]
                weightedVectors  = [0.03*vec/20 for vec in resultantVectors]

                print 'Vector magnitudes: {}'.format(vectorMagnitudes)
                print 'Normed vectors: {}'.format(normedVectors)
                print 'Magnitude-weighted vectors: {}'.format(weightedVectors)
                del resultantVectors


                vecStart = self.densmap.abc2xyz(asymIndices = mainNegCentroid[:3],
                                                fracInput   = True)

                # append vector positions to a pdb file 
                fname = 'Dshift-vectors.pdb'
                try:
                    with open(fname) as f:
                        for i, l in enumerate(f):
                            pass
                    lineCount = i + 1
                except IOError:
                    lineCount = 0
                vectorPDBfile = open(fname,'a')

                from classHolder import singlePDB
                from PDBFileManipulation import writePDBline
                vecPosAsAtom = singlePDB(atomnum    = lineCount,
                                         residuenum = atmId.split('-')[1],
                                         atomtype   = atmId.split('-')[3],
                                         basetype   = 'V0',
                                         chaintype  = atmId.split('-')[0],
                                         X_coord    = vecStart[0],
                                         Y_coord    = vecStart[1],
                                         Z_coord    = vecStart[2])
                pdbLine = writePDBline(vecPosAsAtom,1)
                vectorPDBfile.write(pdbLine+'\n')

                for i,vec in enumerate(weightedVectors):
                    vecEnd   = self.densmap.abc2xyz(asymIndices = np.array(mainNegCentroid[:3]) + vec,
                                                    fracInput   = True)
                    # print 'cgo_arrow [{}],[{}],hlength=0.25,radius=0.05,hradius=0.075'.format(','.join([str(round(v,2)) for v in vecStart]),','.join([str(round(v,2)) for v in vecEnd]))

                    vecPosAsAtom = singlePDB(atomnum    = lineCount+i+1,
                                             residuenum = atmId.split('-')[1],
                                             atomtype   = atmId.split('-')[3],
                                             basetype   = 'V'+str(i+1),
                                             chaintype  = atmId.split('-')[0],
                                             X_coord    = vecEnd[0],
                                             Y_coord    = vecEnd[1],
                                             Z_coord    = vecEnd[2])
                    pdbLine = writePDBline(vecPosAsAtom,1)
                    vectorPDBfile.write(pdbLine+'\n')

        # scatter plots begin here
        if makePlot:
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = '3d')

            colors = []
            cn,cp = 0,0
            for d in ds:
                if d < 0:
                    cmap = plt.cm.Reds
                    colors.append(cmap((negClusters[cn]+1)/float(4)))
                    cn += 1
                else:
                    cmap = plt.cm.Greens
                    colors.append(cmap((posClusters[cp]+1)/float(numPosClusts+2)))
                    cp += 1

            if showSymm:
                ax.scatter(xs, ys, zs, c = symOrOriginal, s = 10, alpha = 0.6)

            elif showDens:
                ax.scatter(xs, ys, zs, c = colors, s = 1000*np.abs(ds), alpha = 0.6)

                if findClusters:
                    if plotCentroids:
                        if negClustering:
                            ax.scatter(negCentroids[:3, 0],
                                       negCentroids[:3, 1],
                                       negCentroids[:3, 2],
                                       marker     = 'x',
                                       s          = 100,
                                       linewidths = 3,
                                       color      = 'k')

                        if posClustering:
                            ax.scatter(posCentroids[:, 0],
                                       posCentroids[:, 1],
                                       posCentroids[:, 2],
                                       marker     = 'd',
                                       s          = 100,
                                       linewidths = 3,
                                       color      = 'b')

                    if negClustering and posClustering:
                        for vec,length in zip(normedVectors,vectorMagnitudes):
                            ax.quiver(mainNegCentroid[0],
                                      mainNegCentroid[1],
                                      mainNegCentroid[2],
                                      vec[0],
                                      vec[1],
                                      vec[2],
                                      length=0.001*length,
                                      pivot = 'tail',
                                      color = 'k',
                                      arrow_length_ratio = 0.1)
            plt.title(atmId)
            plt.show()

        output = []
        if findClusters:
            if negClustering:
                if negMean1 > negMean2:
                    output.append(negMean2)
                else:
                    output.append(negMean1)
                if posClustering:
                    output.append(normedVectors)
                    output.append(vectorMagnitudes)
        return output

    def plotFCdistnPlot(self,
                        plot              = True,
                        atomOfInterest    = '',
                        atomsToPlot       = ['GLU-CD','CYS-SG'],
                        atomFCvals        = [],
                        atomFCvalsMaxNorm = [],
                        FCatMin           = [],
                        reliability       = [],
                        plotType          = '.png',
                        axesFont          = 18):

        # plot a kde & histrogram distribution plot for the FCalc values for an 
        # atom, both raw, and after being divided by the maximum FCalc value
        # attained for that atom (normalised-FCalc). The plot will also include 
        # vertical lines indicating the FCalc and normalised-FCalc values attained
        # for the voxel where the most negative density map (not FC map) voxel 
        # within the local region around the atom (this is the voxel corresponding
        # to the DLoss metric value).

        for tag in atomsToPlot:
            if tag in atomOfInterest.getAtomID():
                sns.set_style("dark")
                sns.set_context(rc = {"figure.figsize": (10, 6)})
                fig = plt.figure()
                ax = plt.subplot(111)
                sns.distplot(np.array(atomFCvals), label = 'Fcalc')
                sns.distplot(np.array(atomFCvalsMaxNorm), label = 'Fcalc/max(Fcalc)')
                ylims = ax.get_ylim()

                plt.plot((FCatMin,FCatMin),
                         (ylims[0],ylims[1]),
                         label='Fcalc, at position of min diff density')
                plt.plot((reliability,reliability),
                         (ylims[0],ylims[1]),
                         label = 'Fcalc/max(Fcalc), at position of min diff density')
                leg = plt.legend(frameon = 1)
                frame = leg.get_frame()
                frame.set_color('white')
                plt.xlabel('Per-voxel density map values', fontsize = axesFont)
                plt.ylabel('Normed-frequency', fontsize = axesFont)
                plt.title('Distribution of Fcalc density values: {}'.format(atomOfInterest.getAtomID()))
                fig.savefig('{}testDistnPlot-{}{}'.format(self.filesOut,atomOfInterest.getAtomID(),plotType))

    def plotDensScatterPlots(self,
                             printText = False):

        # plot scatter plots for density metrics 

        self.startTimer()
        self.fillerLine(style = 'line')
        str = 'Plotting scatter plots for electron density statistics...'
        self.lgwrite(ln = str,forcePrint = printText)

        plotVars = [['meandensity','maxdensity'],
                    ['meandensity','mediandensity'],
                    ['meandensity','mindensity'],
                    ['mindensity','maxdensity'],
                    ['meandensity','stddensity'],
                    ['mindensity','min90tile'],
                    ['maxdensity','max90tile'],
                    ['min90tile','min95tile'],
                    ['max90tile','max95tile']]

        # # only include below if per-atom clusters are
        # # calculated - currently very slow
        # plotVars += [['negCluster','meandensity'],
        #              ['negCluster','mindensity']]

        if self.calcFCmap:
            plotVars.append(['meandensity','wMean'])
            plotVars.append(['mindensity','wMean'])

        for pVars in plotVars:
            logStr = edens_scatter(outputDir = self.filesOut,
                                   metrics   = pVars,
                                   PDBarray  = self.PDBarray,
                                   pdbName   = self.pdbName)
            self.lgwrite(ln = logStr)

    def plotPerResidueBoxPlots(self,
                               lineStyle = 'line'):

        # perform residue analysis for datatset, outputting boxplots
        # for each atom specific to each residue, and also a combined 
        # boxplot across all residues in structures.

        self.fillerLine(style = lineStyle)
        for densMet in ('mean','min','max'):
            resArray = densper_resatom_NOresidueclass(where    = self.filesOut,
                                                      PDBarray = self.PDBarray,
                                                      densMet  = densMet,
                                                      pdbName  = self.pdbName)
        densper_res(where    = self.filesOut,
                    resArray = resArray,
                    pdbName  = self.pdbName)

        resArray = []
        self.stopTimer()

    def pickleAtomList(self):

        # save list of atom objects to a .pkl file

        self.pklFileName = save_objectlist(self.PDBarray,self.pdbName)

    def startTimer(self):

        # start a timer

        self.timeStart = time.time()

    def stopTimer(self,
                  includeInLog = False):

        # stop a timer (must run startTimer before)

        elapsedTime = time.time() - self.timeStart
        if includeInLog is True:
            self.lgwrite(ln = 'section time: {}s\n'.format(round(elapsedTime,3)))
        sys.stdout.flush()

    def success(self):

        # report success to log file

        self.lgwrite(ln = '---> success')

    def fillerLine(self,
                   style = 'blank'):

        # print a filler line (several styles)
        # to command line

        if style == 'stars':
            ln = '\n***'
        elif style == 'line':
            ln = '\n'+'-'*30
        elif style == 'blank':
            ln = '\n'
        self.lgwrite(ln = ln)

    def lgwrite(self,
                ln         = '',
                strip      = True,
                forcePrint = False):

        # write line to log file

        self.log.writeToLog(str        = ln, 
                            strip      = strip,
                            forcePrint = forcePrint)

    def printStepNumber(self):

        # print a string indicating the current pipeline 
        # step number directory to the command line

        try:
            self.stepNumber
        except AttributeError:
            self.stepNumber = 1
        ln =  '\n_______'+\
              '\nSTEP {})'.format(self.stepNumber)
        self.lgwrite(ln = ln)

        self.stepNumber += 1




