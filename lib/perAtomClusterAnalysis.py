import numpy as np
import os
from sklearn.cluster import KMeans
from scipy.interpolate import griddata


class perAtomXYZAnalysis():

        # class to determine the full set of voxels around an atom (including
        # symmetry related voxels). Then can find clusters of positive and
        # negative dnesity values. Can determine 3d density shift vectors
        # describing the shift in electron density from negative to positive
        # regions (based on these particular defined clusters of voxels)

    def __init__(self,
                 atomObj=0, densMapObj=[], densPerVxl=[],
                 xyzsPerAtom=[], knownRefPoint=[0, 0, 0], showDens=True,
                 knownRefPoint2=[], showSymm=False, makePlot=False,
                 report=False, plotCentroids=False, numPosClusts=5,
                 writeVecsToFile=False,  mainVecOnly=True,
                 doInterpolation=False, vxlMidPt=[],
                 autoRun=False, vxlRefPoint=[]):

        self.atmId = atomObj.getAtomID()
        self.knownRefPoint = knownRefPoint
        self.knownRefPoint2 = knownRefPoint2
        self.showDens = showDens
        self.showSymm = showSymm
        self.makePlot = makePlot
        self.report = report
        self.plotCentroids = plotCentroids
        self.numPosClusts = numPosClusts
        self.mainVecOnly = mainVecOnly
        self.writeVecsToFile = writeVecsToFile
        self.doInterpolation = doInterpolation
        self.vxlRefPoint = vxlRefPoint
        self.densMapObj = densMapObj
        self.xyzsPerAtom = xyzsPerAtom
        self.densPerVxl = densPerVxl
        self.vxlMidPt = vxlMidPt

        if autoRun:
            self.run()

    def run(self):

        # run main function of class

        self.getxyzPerAtom()

        self.findVoxelMidPt()

        props = self.doMetricClustAnalysis()

        if self.makePlot:
            self.make3dScatterPlot(
                negClustInfo=props[0], posClustInfo=props[1],
                vectorDic=props[2], plotVectors=True)

    def getxyzPerAtom(self):

        # get list of xyz positions of voxels for atom

        XYZ = np.transpose(self.xyzsPerAtom)
        symOps = self.getExtendedSymRelatedPoints(XYZ=XYZ)
        refPoint = self.findVoxRefPoint(XYZ)
        props = self.decideWhichVoxToKeep(
            refPoint=refPoint, symOps=symOps, X=XYZ[0])

        if self.doInterpolation:
            keptPointsNew = self.interpolateGrid(props[0])
            self.interpolatedPts = np.ndarray.flatten(np.array(keptPointsNew))

        self.keptPts = props[0]
        self.symOrOrig = props[1]

    def doMetricClustAnalysis(self):

        (negClustInfo, posClustInfo) = self.findPosNegClusters()

        if self.foundNegClust and self.foundPosClust:
            vecDic = self.findNegToPosDensShift(negClustInfo, posClustInfo)

            if self.writeVecsToFile:
                self.writeVecStartStopToPDBFile(
                    mainNegCentroid=negClustInfo['main centroid'],
                    weightedVectors=vecDic['weighted vectors'])
        else:
            vecDic = {}
            self.netDensShift = np.nan

        return negClustInfo, posClustInfo, vecDic

    def partitionPtsByVec(self,
                          method='known ref point'):

        # split voxels around an atom up by one of several methods
        # and then determine the behaviour of voxels in each region

        self.splitDensPointsBasedOnVector(method=method)
        self.calcSumDensForPartition()

    def getbasicSymRelatedPoints(self,
                                 XYZ=[]):

        # retrieve the symmetry operations for the map file.
        # For any sym op that takes a point outside of unit cell,
        # move translate this back into the unit cell of interest

        symOps = self.densMapObj.getSymOps(XYZ[0], XYZ[1], XYZ[2])
        symOps = np.mod(symOps, 1).tolist()

        return symOps

    def getExtendedSymRelatedPoints(self,
                                    XYZ=[], threshold=0.01, doNull=True):

        # retrieve the symmetry operations for the map file.
        # for case where original xyz points are close to a unit
        # cell boundary, include additional adjacent unit cells

        if doNull:
            symOps = [XYZ]
            return symOps

        symOps = self.getbasicSymRelatedPoints(XYZ)

        translateUnitCell = []

        for i in range(len(symOps[0][0])):
            x = symOps[0][0][i]
            y = symOps[0][1][i]
            z = symOps[0][2][i]
            trans = []
            if x < threshold:
                trans = [[-1], [0], [0]]
            if x > 1-threshold:
                trans = [[1], [0], [0]]
            if y < threshold:
                trans = [[0], [-1], [0]]
            if y > 1-threshold:
                trans = [[0], [1], [0]]
            if z < threshold:
                trans = [[0], [0], [-1]]
            if z > 1-threshold:
                trans = [[0], [0], [1]]
            if trans != []:
                if trans not in translateUnitCell:
                    translateUnitCell.append(trans)

        l = len(translateUnitCell)
        if l != 0:
            for t in translateUnitCell:
                for sym in np.array(symOps) + t:
                    symOps.append(sym)
            if self.report:
                print 'Must consider {} adjacent unit cell(s)'.format(l)

        return symOps

    def findVoxRefPoint(self,
                        XYZ=[]):

        # find suitable reference voxel point (out of those assigned to atom)
        # that is nearest to centre of asym unit. This will be the reference
        # point from to which the distance to over voxels will be determined

        if self.vxlRefPoint == []:
            centOfInterest = np.mean(XYZ, axis=1)
        else:
            centOfInterest = self.vxlRefPoint

        xyzs = [[XYZ[j][ii] for j in range(3)] for ii in range(len(XYZ[0]))]
        dists = list(np.linalg.norm(np.array(xyzs) - np.array(centOfInterest),
                     axis=1))
        minInd = dists.index(min(dists))
        refPoint = [XYZ[0][minInd], XYZ[1][minInd],
                    XYZ[2][minInd]] + [self.densPerVxl[minInd]]

        return refPoint

    def decideWhichVoxToKeepAlt(self,
                                refPoint=[], XYZ=[], dens=[], refThres=2,
                                performEndCheck=True, doNothing=False):

        # this version of the decision function can only work when the map
        # is centred on the molecule. In this case, a single well formed
        # cluster must be present for the atom, however symmetry-related
        # additional voxels may be present, which should be removed

        xyzs = [[XYZ[j][ii] for j in range(3)] for ii in range(len(XYZ[0]))]
        xyzds = [xyz+[d] for xyz, d in zip(xyzs, dens)]

        if doNothing:
            return xyzds, [0]*len(xyzds)

        A = np.array(xyzs) - np.array(refPoint[:3])
        distToRef = np.linalg.norm(A, axis=-1)
        xyzds = [xyzd for _, xyzd in sorted(zip(distToRef, xyzds))]
        distToRef.sort()
        keptPts = []
        for xyzd, dist in zip(xyzds, distToRef):
            if dist < refThres:
                keptPts.append(xyzd)
                continue
            A = np.array(keptPts)[:, None, :3] - np.array(xyzd[0:3])
            distToFound = np.linalg.norm(A, axis=2)
            minDist = np.min(distToFound, axis=0)[0]
            if minDist < refThres:
                keptPts.append(xyzd)

        return keptPts, [0]*len(keptPts)

    def decideWhichVoxToKeep(self,
                             refPoint=[], symOps=[], X=[],
                             threshold=1e-10, refThres=0.01,
                             performEndCheck=True, doNothing=False,
                             returnAll=False):

        # loop through all the located voxels assigned to the atom
        # and decide which ones form a single well defined cluster
        # in space. This will then be what is treated as the
        # voxel set of this particular atom

        # provide the null choice to return original points
        # (but not sym-related points)
        if doNothing:
            xyz = [[symOps[0][j][ii] for j in range(3)] + [self.densPerVxl[ii]] for ii in range(len(X))]
            symOrOrig = [0 for i in range(len(xyz))]
            return xyz, symOrOrig

        # provide the choice to not exclude any sym-related points
        if returnAll:
            symPointsAll = []
            symOrOrig = []
            for i in range(len(X)):
                symxyzs = [[symOp[j][i] for j in range(3)] for symOp in symOps]
                symPoints = [sym + [self.densPerVxl[i]] for sym in symxyzs]
                symPointsAll += symPoints
                symOrOrig += [0]+[1]*(len(symxyzs)-1)
            return symPointsAll, symOrOrig

        keptPoints = np.array([refPoint])
        keptPointsL = [refPoint]
        symOrOrig = [0]
        samexyzs = 0
        dupcount = 0

        xyz = [[symOps[0][j][ii] for j in range(3)] for ii in range(len(X))]
        A = np.array(xyz) - np.array(refPoint[:3])
        distToRef = np.linalg.norm(A, axis=-1)

        for ii, refDist in enumerate(distToRef):

            if refDist > refThres:

                symxyzs = [[symOp[j][ii] for j in range(3)] for symOp in symOps]
                symPoints = [sym + [self.densPerVxl[ii]] for sym in symxyzs]

                copy = False
                for i, symPt in enumerate(symPoints):
                    if symPt in keptPointsL:
                        copy = True
                        dupcount += 1

                        if i == 0:
                            symOrOrig[keptPointsL.index(symPt)] = 0
                        break
                if copy:
                    continue

                A = np.array(symxyzs) - keptPoints[:, None, :3]
                distToFound = np.linalg.norm(A, axis=2)

                minDist = np.min(distToFound, axis=0).tolist()
                minOfMinDist = min(minDist)
                closestInd = minDist.index(minOfMinDist)
                closestSym = symPoints[closestInd]

                if minOfMinDist < threshold:
                    samexyzs += 1

                    # if new point is part of original asym unit,
                    # use this point instead
                    if closestInd == 0:
                        dists = [dist[0] for dist in distToFound]
                        nearestPtId = dists.index(min(minDist))
                        keptPoints = np.delete(keptPoints, nearestPtId, axis=0)
                        keptPointsL = keptPointsL[:nearestPtId]+keptPointsL[nearestPtId+1:]
                        symOrOrig = symOrOrig[:nearestPtId]+symOrOrig[nearestPtId+1:]
                    else:
                        continue

            else:
                closestInd = 0
                closestSym = xyz[ii] + [self.densPerVxl[ii]]

                distsToPoint = np.linalg.norm(keptPoints[:, :3]-np.array(xyz[ii]), axis=-1)
                minDist = min(distsToPoint)
                if minDist < threshold:
                    nearestPtId = distsToPoint.tolist().index(minDist)
                    keptPoints = np.delete(keptPoints, nearestPtId, axis=0)
                    keptPointsL = keptPointsL[:nearestPtId]+keptPointsL[nearestPtId+1:]
                    symOrOrig = symOrOrig[:nearestPtId]+symOrOrig[nearestPtId+1:]

            # add the newly found point to the list of found points
            keptPoints = np.append(keptPoints, [closestSym], axis=0)
            keptPointsL.append(closestSym)

            if closestInd != 0:
                symOrOrig.append(1)
            else:
                symOrOrig.append(0)

        # check for included points with same xyz coordinates
        # - this should not happen
        if performEndCheck:
            for i, keptPoint in enumerate(keptPoints):
                keptPoints2 = np.delete(keptPoints, i, axis=0)
                dists = np.linalg.norm(
                    keptPoints2[:, :3]-keptPoint[:3], axis=-1)
                if min(dists) < threshold:
                    print 'Error! 2 density points included with identical xyz'
                    import sys
                    sys.exit()

        if self.report:
            self.decideWhichVoxToKeepReport(
                keptPoints, symOrOrig, dupcount, samexyzs)

        return keptPoints.tolist(), symOrOrig

    def decideWhichVoxToKeepReport(self,
                                   keptPoints=[], symOrOrig=[], dupcount=[],
                                   samexyzs=[]):

        # print to command line how the previous method did

        print 'Atom Identity: {}'.format(self.atmId)
        print 'Number points found: {}'.format(len(keptPoints))
        print '\tNumber that are original: {}'.format(
            len(symOrOrig) - np.sum(symOrOrig))
        print '\tNumber that are sym copies: {}'.format(np.sum(symOrOrig))
        print '\tNumber of duplicate positions: {}'.format(dupcount + samexyzs)
        print '\t\tidentical: {}'.format(dupcount)
        print '\t\tsame xyz: {}'.format(samexyzs)

    def findPosNegClusters(self):

        # now that a complete set of voxels have been found around an atom
        # determine the location of clusters of positive and negative density

        posVals = [pt for pt in self.keptPts if pt[3] >= 0]
        negVals = [pt for pt in self.keptPts if pt[3] < 0]

        negClustInfo = self.findNegClusters(negVals=negVals)
        posClustInfo = self.findPosClusters(posVals=posVals)

        return negClustInfo, posClustInfo

    def findNegClusters(self,
                        n_clusters=2, negVals=[]):

        # find negative point clusters

        negClustInfo = {}

        if len(negVals) > 1:
            self.foundNegClust = True

            negKmeans = KMeans(n_clusters=n_clusters)
            negClusters = negKmeans.fit_predict(negVals)
            negCluster1 = [val for i, val in enumerate(negVals) if negClusters[i] == 0]
            negCluster2 = [val for i, val in enumerate(negVals) if negClusters[i] == 1]
            mean1 = np.mean([c[3] for c in negCluster1])
            mean2 = np.mean([c[3] for c in negCluster2])
            centroids = negKmeans.cluster_centers_

            if mean2 > mean1:
                mainClust = negCluster1
                mainCentrID = 0
                clusters = np.mod(np.array(negClusters)-1, 2)
                self.topNegClustMean = mean1

            else:
                mainClust = negCluster2
                mainCentrID = 1
                clusters = negClusters
                self.topNegClustMean = mean2

            negClustInfo['clusters'] = clusters
            negClustInfo['main cluster'] = mainClust
            negClustInfo['centroids'] = centroids
            negClustInfo['main centroid'] = centroids[mainCentrID, :]
            negClustInfo['mean 1'] = mean1
            negClustInfo['mean 2'] = mean2

        else:
            self.foundNegClust = False
            negClustInfo['clusters'] = [1]
            self.topNegClustMean = np.nan

        return negClustInfo

    def findPosClusters(self,
                        numPosClusts=10, posVals=[]):

        # find positive point clusters

        posClustInfo = {}

        if len(posVals) > 1:
            # reduce number of clusters if too many specified
            if numPosClusts > len(posVals):
                numPosClusts = len(posVals)

            self.foundPosClust = True
            posKmeans = KMeans(n_clusters=numPosClusts)
            clusters = posKmeans.fit_predict([pt[:3] for pt in posVals])
            posClustList = []

            for j in range(numPosClusts):
                posClustList.append([val for i, val in enumerate(posVals) if clusters[i] == j])

            posClustMeans = [np.mean([c[3] for c in clust]) for clust in posClustList]
            maxPosInd = posClustMeans.index(max(posClustMeans))

            posClustInfo['clusters'] = clusters
            posClustInfo['centroids'] = posKmeans.cluster_centers_
            posClustInfo['main cluster'] = maxPosInd
            posClustInfo['cluster list'] = posClustList
            posClustInfo['means'] = posClustMeans

        else:
            self.foundPosClust = False
            posClustInfo['clusters'] = [1]
            posClustInfo['cluster list'] = [1]

        return posClustInfo

    def findNegToPosDensShift(self,
                              negClustInfo={}, posClustInfo={},
                              printText=False, performCheck=True):

        # now that clusters have been located, find the shift vectors
        # which determine the magnitude and direction in which density
        # appears to have shifted. A choice can be made between only
        # considering one main positive cluster of interest or all of them

        if self.mainVecOnly:
            ind = posClustInfo['main cluster']
            clusterLoop = [posClustInfo['cluster list'][ind]]
        else:
            clusterLoop = posClustInfo['cluster list']

        resultantVectors = []

        for clust in clusterLoop:
            vecChange = np.array(clust)[:, :3] - np.array(
                negClustInfo['main cluster'])[:, None, :3]
            densDiffs = np.array(clust)[:, 3] - np.array(
                negClustInfo['main cluster'])[:, None, 3]

            norms = np.linalg.norm(vecChange, axis=-1)
            directions = vecChange/(norms[:, :, None])

            if performCheck:
                for vecs in vecChange:
                    for vec in vecs:
                        if np.linalg.norm(vec) == 0.0:
                            print 'Error! zero length density shift vector'
                            import sys
                            sys.exit()

            totalComponent = (densDiffs[:, :, None]*directions).mean(1)
            resultantVectors.append(np.mean(np.array(totalComponent), axis=0))

        vectorDic = {}
        vMagnitudes = np.linalg.norm(resultantVectors, axis=1)
        vectorDic['normed vectors'] = np.array(resultantVectors)/(vMagnitudes[:, None])
        vectorDic['weighted vectors'] = np.array(resultantVectors)*(0.03/20)
        vectorDic['magnitudes'] = vMagnitudes

        self.netDensShiftVec = np.sum(resultantVectors, 0)
        self.netDensShift = np.linalg.norm(self.netDensShiftVec)

        if printText:
            for k in vectorDic.keys():
                print '{}: {}'.format(k, np.around(vectorDic[k], decimals=2))

        return vectorDic

    def findMinToMaxVector(self):

        # as an alternative to determining positive and negative clusters
        # and then the density shift vectors between them, find location
        # of max and min points and calculate direction vector between them

        minPt = self.findVoxelMinPt()
        maxPt = self.findVoxelMaxPt()

        dirVector = np.array(maxPt) - np.array(minPt)

        self.minToMaxVector = dirVector

    def findVoxelMaxPt(self):

        # find the max point assigned to the atom

        densities = self.getDensities()
        maxPtId = densities.index(max(densities))
        maxX = self.getXcoords()[maxPtId]
        maxY = self.getYcoords()[maxPtId]
        maxZ = self.getZcoords()[maxPtId]

        return [maxX, maxY, maxZ]

    def findVoxelMinPt(self):

        # find the min point assigned to the atom

        densities = self.getDensities()
        minPtId = densities.index(min(densities))
        minX = self.getXcoords()[minPtId]
        minY = self.getYcoords()[minPtId]
        minZ = self.getZcoords()[minPtId]

        return [minX, minY, minZ]

    def findVoxelMidPt(self,
                       points=[]):

        # find the mid point assigned to the atom

        midX = np.mean(self.getXcoords())
        midY = np.mean(self.getYcoords())
        midZ = np.mean(self.getZcoords())

        self.vxlMidPt = [midX, midY, midZ]

        return [midX, midY, midZ]

    def findRefPointToAtomVector(self,
                                 refPoint=[]):

        # define a direction vector from a known reference point
        # to the current atom (were the centres are defined by
        # the mid voxel position for that atom)

        dirVector = np.array(self.vxlMidPt) - np.array(refPoint)

        return dirVector/np.linalg.norm(dirVector)

    def findVoxelsAlongVector(self,
                              dirVec=[1, 0, 0],
                              posVec=[0, 0, 0], thres=5e-3,
                              method='interpolation'):

        # supplied with a 3d vector, the voxels lying on the
        # line defined by direction vector 'dirVec' and passing
        # through the point posVec are returned (or within a
        # specified distance threshold 'thres')

        x1 = np.array(posVec)
        x2 = x1 + np.array(dirVec)
        distsToMid = []
        densKept = []

        if method == 'projection':
            for pt in self.keptPts:
                x0minusx1 = np.array(pt[:3]) - x1
                x0minusx2 = np.array(pt[:3]) - x2
                x2minusx1 = x2 - x1
                top = np.linalg.norm(np.cross(x0minusx1, x0minusx2))
                dist = float(top)/np.linalg.norm(x2minusx1)

                if dist < thres:
                    proj = x1 + np.dot(x0minusx1, x2minusx1)/np.dot(
                        x2minusx1, x2minusx1) * x2minusx1
                    distToMid = np.linalg.norm(
                        proj-np.array(posVec))/np.linalg.norm(np.array(dirVec))
                    for i in range(3):
                        if np.abs(proj[i]-x1[i]-distToMid*dirVec[i]) > 1e-3:
                            distToMid = -distToMid
                            break
                    distsToMid.append(distToMid)
                    densKept.append(pt[3])

        elif method == 'interpolation':
            m = self.maxVecDist
            pointsOnVector = [x1+t*np.array(dirVec) for t in np.linspace(-m, m)]
            xs = [[pt[0]] for pt in pointsOnVector]
            ys = [[pt[1]] for pt in pointsOnVector]
            zs = [[pt[2]] for pt in pointsOnVector]
            newPts = self.interpolateGrid(keptPoints=self.keptPts,
                                          targetPts=[xs, ys, zs])
            for t, pt in zip(np.linspace(-m, m), newPts):
                distsToMid.append(t)
                densKept.append(pt[3])

        self.voxelsAlongVector = {'dists': distsToMid,
                                  'densities': densKept}

    def splitDensPointsBasedOnVector(self,
                                     partition=2, method='density shift'):

        # with a density shift vector defined, split the set of x,y,z
        # density voxel points into sections, based on where they lie
        # with respect to the density shift vector

        if method == 'density shift':
            self.findMinToMaxVector()
            try:
                shiftVec = self.netDensShiftVec/np.linalg.norm(
                    self.netDensShiftVec)
            except AttributeError:
                shiftVec = np.array([1, 0, 0])

        elif method == 'min to max':
            self.findMinToMaxVector()
            shiftVec = self.minToMaxVector

        elif method == 'known ref point':
            shiftVec = self.findRefPointToAtomVector(self.knownRefPoint)
            # this may fail if the ref point is badly chosen
            if np.isnan(reduce(lambda x, y: x*y, shiftVec)):
                print 'Error: failed to calculate vector to ref point. ' +\
                      'Using a random coordinate system for atom instead.'
                shiftVec = np.random.randn(3)
            if self.knownRefPoint2 != []:
                shiftVec2 = self.findRefPointToAtomVector(self.knownRefPoint2)
                if np.isnan(reduce(lambda x, y: x*y, shiftVec2)):
                    print 'Error: failed to calculate vector to ref point. ' +\
                          'Using a random coordinate system for atom instead.'
                    shiftVec2 = np.random.randn(3)
        else:
            print 'Error in "method" parameter assignment'
            return

        # get a (non-unique) orthonormal basis including shiftVec
        if shiftVec2 == []:
            # if only 1 vector known
            orth1 = np.random.randn(3)
            orth1 -= orth1.dot(shiftVec) * shiftVec
            orth1 /= np.linalg.norm(orth1)
            orth2 = np.cross(shiftVec, orth1)
        else:
            # if 2 (possibly non orthogonal) vectors known
            orth1 = np.cross(shiftVec, shiftVec2)
            orth1 /= np.linalg.norm(orth1)
            orth2 = np.cross(shiftVec, orth1)
            orth2 /= np.linalg.norm(orth2)

        # find angles between vector basis and each voxel point
        vecs = np.array(self.keptPts)[:, :3] - np.array(self.vxlMidPt)
        dists = np.linalg.norm(vecs, axis=-1)[:, None]
        normVecs = vecs / dists

        angles1 = np.arccos(np.clip(np.dot(normVecs, shiftVec), -1.0, 1.0))
        angles2 = np.arccos(np.clip(np.dot(normVecs, orth1), -1.0, 1.0))
        angles3 = np.arccos(np.clip(np.dot(normVecs, orth2), -1.0, 1.0))

        if partition == 1:
            order = [[i, j, k] for i in range(2) for j in range(2) for k in range(2)]
        elif partition == 2:
            order = [[0, j, k, d] for j in range(1, 3) for k in range(1, 3) for d in range(2)] + [[i, 0, 0, d] for i in [1, 2] for d in range(2)]
            self.medVecDist = np.median(dists)
            self.maxVecDist = np.max(dists)

        # self.findVoxelsAlongVector(dirVec=shiftVec, posVec=self.vxlMidPt)

        positions = []
        for i, angs in enumerate([angles1, angles2, angles3]):
            # note that one point will have angles 'nan' - it is the mid point
            # by default this will enter the [0,0,0] group
            if partition == 1:
                positions.append([1 if ang < np.pi/2 else 0 for ang in angs])
            elif partition == 2:
                if i == 0:
                    positions.append([1 if ang < np.pi/4 else 2 if ang > 3*np.pi/4 else 0 for ang in angs])
                elif i == 1:
                    positions.append(map(int, [0 if positions[0][j] in (1, 2) else 1 if np.isnan(ang) else np.ceil(2*ang/np.pi) for j, ang in enumerate(angs)]))
                else:
                    positions.append(map(int, [0 if positions[0][j] in (1, 2) else 1 if np.isnan(ang) else np.ceil(2*ang/np.pi) for j, ang in enumerate(angs)]))
                    positions.append([0 if d <= self.medVecDist else 1 for d in dists])
        positions = np.transpose(positions)

        groups = [order.index(list(pos)) for pos in positions]

        self.ptPartition = {'groups': groups, 'number': len(order)}

    def calcSumDensForPartition(self,
                                metric='mean'):

        # for a set of x,y,z voxels 'points' and a chosen classification
        # scheme self.ptPartition, calculate the min of the density
        # within each partition.

        densByRegion = []

        for i in range(self.ptPartition['number']):
            subset = [pt[3] for pt, prt in zip(
                self.keptPts, self.ptPartition['groups']) if prt == i]
            if subset != []:
                if metric == 'min':
                    sumDens = np.min(subset)
                elif metric == 'mean':
                    sumDens = np.mean(subset)
                elif metric == 'minmax':
                    sumDens = [np.min(subset), np.max(subset)]
            else:
                if metric == 'minmax':
                    sumDens = [np.nan]*2
                else:
                    sumDens = np.nan
            densByRegion.append(np.round(sumDens, 3))

        meanNoNan = np.nanmean(densByRegion, 0)
        for i in range(self.ptPartition['number']):
            if metric == 'minmax':
                val = densByRegion[i][0]
            else:
                val = densByRegion[i]
            if np.isnan(val):
                densByRegion[i] = meanNoNan

        # flatten list to 1d if needed
        if metric == 'minmax':
            densByRegion = [a for b in densByRegion for a in b]

        # densByRegion = [min(densByRegion)]

        self.densByRegion = densByRegion

    def writeVecStartStopToPDBFile(self,
                                   mainNegCentroid=[],
                                   outputPDBname='Dshift-vectors.pdb',
                                   weightedVectors=[]):

        # write located density shift vector start and stop
        # xyz positions to a pdb file so that they can then be
        # visualised in pymol (or otherwise). 'densMapObj' is
        # taken from readAtomMap.py (self.densmap).

        from classHolder import singlePDB
        from PDBFileManipulation import writePDBline

        vecStart = self.densMapObj.abc2xyz(
            asymIndices=mainNegCentroid[:3], fracInput=True)

        # append vector positions to a pdb file if it exists
        lineCount = 0
        if os.path.isfile(outputPDBname):
            try:
                with open(outputPDBname) as f:
                    i = -1
                    for i, l in enumerate(f):
                        pass
                lineCount = i + 1
            except IOError:
                pass

        vectorPDBfile = open(outputPDBname, 'a')

        atmIDparts = self.atmId.split('-')

        vecPosAsAtom = singlePDB(
            atomnum=lineCount, residuenum=atmIDparts[1],
            atomtype=atmIDparts[3], basetype='V0',
            chaintype=atmIDparts[0], X_coord=vecStart[0],
            Y_coord=vecStart[1], Z_coord=vecStart[2])

        pdbLine = writePDBline(vecPosAsAtom, 1)
        vectorPDBfile.write(pdbLine + '\n')

        for i, vec in enumerate(weightedVectors):
            vecEnd = self.densMapObj.abc2xyz(
                asymIndices=np.array(mainNegCentroid[:3]) + vec,
                fracInput=True)

            vecPosAsAtom = singlePDB(
                atomnum=lineCount + i + 1, residuenum=atmIDparts[1],
                atomtype=atmIDparts[3], basetype='V' + str(i+1),
                chaintype=atmIDparts[0], X_coord=vecEnd[0],
                Y_coord=vecEnd[1], Z_coord=vecEnd[2])

            pdbLine = writePDBline(vecPosAsAtom, 1)
            vectorPDBfile.write(pdbLine + '\n')

        vectorPDBfile.close()

    def interpolateGrid(self,
                        method='linear', targetPts=[]):

        kpArray = np.array(self.keptPts)
        XYZ = kpArray[:, :3]
        dens = kpArray[:, 3]
        minX = min(XYZ[:, 0])
        maxX = max(XYZ[:, 0])
        minY = min(XYZ[:, 1])
        maxY = max(XYZ[:, 1])
        minZ = min(XYZ[:, 2])
        maxZ = max(XYZ[:, 2])

        self.deltaX = maxX - minX
        self.deltaY = maxY - minY
        self.deltaZ = maxZ - minZ

        if targetPts == []:
            # if no specific target points specified,
            # interpolate on whole 3d grid
            gridx, gridy, gridz = np.mgrid[
                minX:maxX:10j, minY:maxY:10j, minZ:maxZ:10j]
        else:
            gridx = np.array(targetPts[0])
            gridy = np.array(targetPts[1])
            gridz = np.array(targetPts[2])

        interpDens = griddata(XYZ, dens, (gridx, gridy, gridz), method=method)

        gridxFlat = np.ndarray.flatten(gridx)
        gridyFlat = np.ndarray.flatten(gridy)
        gridzFlat = np.ndarray.flatten(gridz)
        interpDensFlat = np.ndarray.flatten(interpDens)

        interpolatedPts = [[gridxFlat[i], gridyFlat[i], gridzFlat[i], interpDensFlat[i]] for i in range(len(interpDensFlat))]

        return interpolatedPts

    def getXcoords(self):

        # get x coordinates of density voxel points

        return [xyz[0] for xyz in self.keptPts]

    def getYcoords(self):

        # get y coordinates of density voxel points

        return [xyz[1] for xyz in self.keptPts]

    def getZcoords(self):

        # get z coordinates of density voxel points

        return [xyz[2] for xyz in self.keptPts]

    def getDensities(self):

        # get list of densities for voxel points

        return [xyz[3] for xyz in self.keptPts]

    def make3dScatterPlot(self,
                          negClustInfo={}, posClustInfo={},
                          vectorDic={}, transparency=0.6,
                          plotVectors=True, colByCluster=False,
                          plotCartesian=False):

        # make a 3d scatter plot of the located voxels around an atom
        # in xyz cartesian space. The colour scheme can take multiple
        # options (see below for details)

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        xs = self.getXcoords()
        ys = self.getYcoords()
        zs = self.getZcoords()
        ds = self.getDensities()

        if plotCartesian:
            xyzCart = self.densMapObj.abc2xyz(
                fracInput=True, asymIndices=[xs, ys, zs])
            xs = xyzCart[0]
            ys = xyzCart[1]
            zs = xyzCart[2]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        colors = []
        cn, cp = 0, 0
        for d in ds:
            if d < 0:
                if colByCluster:
                    cmap = plt.cm.Reds
                    colors.append(
                        cmap((negClustInfo['clusters'][cn]+1)/float(4)))
                    cn += 1
                else:
                    colors.append('#ec685b')
            else:
                if colByCluster:
                    cmap = plt.cm.Greens
                    colors.append(
                        cmap((posClustInfo['clusters'][cp]+1)/float(
                            len(posClustInfo['cluster list'])+2)))
                    cp += 1
                else:
                    colors.append('#b8f082')

        if self.partitionPtsByVec:
            colors = []
            cmap = plt.cm.rainbow
            for g in self.ptPartition['groups']:
                colors.append(cmap(float(g)/self.ptPartition['number']))

        if self.showSymm:
            colors, sizes = [], []
            for s in self.symOrOrig:
                if s == 0:
                    colors.append('r')
                    sizes.append(50)
                else:
                    colors.append('y')
                    sizes.append(10)
            ax.scatter(xs, ys, zs, c=colors,
                       s=sizes, alpha=transparency)

        elif self.showDens:
            ax.scatter(xs, ys, zs, c=colors,
                       s=1000*np.abs(ds), alpha=transparency)

            if colByCluster:
                if self.plotCentroids:
                    if self.foundNegClust:
                        ax.scatter(negClustInfo['centroids'][:3, 0],
                                   negClustInfo['centroids'][:3, 1],
                                   negClustInfo['centroids'][:3, 2],
                                   marker='x', s=100, linewidths=3, color='k')

                    if self.foundPosClust:
                        ax.scatter(posClustInfo['centroids'][:, 0],
                                   posClustInfo['centroids'][:, 1],
                                   posClustInfo['centroids'][:, 2],
                                   marker='d', s=100, linewidths=3, color='b')

                if plotVectors:
                    if self.foundNegClust and self.foundPosClust:
                        for vec, length in zip(vectorDic['normed vectors'],
                                               vectorDic['magnitudes']):
                            ax.quiver(negClustInfo['main centroid'][0],
                                      negClustInfo['main centroid'][1],
                                      negClustInfo['main centroid'][2],
                                      vec[0], vec[1], vec[2],
                                      length=0.1*length, pivot='tail',
                                      color='k',  arrow_length_ratio=0.1)

            elif plotVectors:
                vec = self.netDensShiftVec
                ax.quiver(self.vxlMidPt[0], self.vxlMidPt[1],  self.vxlMidPt[2],
                          vec[0], vec[1], vec[2], length=self.medVecDist,
                          pivot='tail', color='k', arrow_length_ratio=0.1)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.title(self.atmId)
        plt.show()
