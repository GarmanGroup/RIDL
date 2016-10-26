import numpy as np
import os
from sklearn.cluster import KMeans

class perAtomClusterAnalysis():

        # class to determine the full set of voxels around an atom (including symmetry 
        # related voxels). Then can find clusters of positive and negative dnesity 
        # values. Can determine 3d density shift vectors describing the shift in 
        # electron density from negative to positive regions (based on these
        # particular defined clusters of voxels)  

    def __init__(self,
                 atmNum          = 0,
                 atmId           = '',
                 densMapObj      = [],
                 vxlsPerAtom     = [],
                 xyzsPerAtom     = [],
                 showDens        = True,
                 showSymm        = False,
                 makePlot        = False,
                 report          = False,
                 findClusters    = True,
                 plotCentroids   = False,
                 numPosClusts    = 10,
                 writeVecsToFile = False,
                 mainVecOnly     = False,
                 autoRun         = True):

        self.atmNum          = atmNum
        self.atmId           = atmId
        self.showDens        = showDens
        self.showSymm        = showSymm
        self.makePlot        = makePlot
        self.report          = report
        self.findClusters    = findClusters
        self.plotCentroids   = plotCentroids
        self.numPosClusts    = numPosClusts
        self.mainVecOnly     = mainVecOnly
        self.writeVecsToFile = writeVecsToFile

        self.densMapObj      = densMapObj
        self.xyzsPerAtomDic  = xyzsPerAtom
        self.vxlsPerAtomDic  = vxlsPerAtom

        self.output          = []

        if autoRun:
            self.run()

    def run(self):

        # run main function of class

        (X,Y,Z)  = self.getxyzPerAtom(xyzsPerAtomDic = self.xyzsPerAtomDic)
        density  = self.getVxlDensitiesPerAtom(vxlsPerAtomDic = self.vxlsPerAtomDic)

        symOps   = self.getExtendedSymRelatedPoints(densMapObj = self.densMapObj,
                                                    XYZ        = [X,Y,Z])

        refPoint = self.findVoxRefPoint(symOps,[X,Y,Z],density)

        props = self.decideWhichVoxToKeep(refPoint = refPoint,
                                          symOps   = symOps,
                                          X        = X,
                                          dens     = density)
        keptPoints    = props[0]
        symOrOriginal = props[1]
        samexyzs      = props[2]

        if self.findClusters:
            (negClustInfo,posClustInfo) = self.findPosNegClusters(keptPoints)

            if self.foundNegClust and self.foundPosClust:
                vecDic = self.findNegToPosDensShift(negClustInfo,posClustInfo)

                if self.writeVecsToFile:
                    self.writeVecStartStopToPDBFile(mainNegCentroid = negClustInfo['main centroid'],
                                                    weightedVectors = vecDic['weighted vectors'],
                                                    densMapObj      = self.densMapObj)

        if self.makePlot:
            self.make3dScatterPlot(keptPoints    = keptPoints,
                                   negClustInfo  = negClustInfo,
                                   posClustInfo  = posClustInfo,
                                   symOrOriginal = symOrOriginal,
                                   vectorDic     = vecDic)

    def addToOutput(self,
                    prop = []):

        # append to the list of desired outputs

        self.output.append(prop)

    def getxyzPerAtom(self,
                      xyzsPerAtomDic = {}):

        # get list of xyz positions of voxels for atom.
        # 'xyzsPerAtomDic' is taken from readAtomMap.py

        if self.atmNum not in xyzsPerAtomDic.keys():
            return False

        X = np.array([xyz[0] for xyz in xyzsPerAtomDic[self.atmNum]])
        Y = np.array([xyz[1] for xyz in xyzsPerAtomDic[self.atmNum]])
        Z = np.array([xyz[2] for xyz in xyzsPerAtomDic[self.atmNum]])

        return X,Y,Z

    def getVxlDensitiesPerAtom(self,
                               vxlsPerAtomDic = {}):

        # get list of density values of voxels for atom.
        # 'vxlsPerAtomDic' is taken from readAtomMap.py.
        # Keep densities here to 3dp since sym copies may have 
        # differing densities to higher accuracy than this

        return np.round(np.array(vxlsPerAtomDic[self.atmNum]),3) 

    def getbasicSymRelatedPoints(self,
                                 densMapObj = [],
                                 XYZ        = []):

        # retrieve the symmetry operations for the map file.
        # 'densMapObj' is taken from readAtomMap.py (self.densmap).
        # For any sym op that takes a point outside of unit cell,
        # move translate this back into the unit cell of interest

        symOps = densMapObj.getSymOps(XYZ[0],XYZ[1],XYZ[2])
        symOps = np.mod(symOps,1).tolist()

        return symOps

    def getExtendedSymRelatedPoints(self,
                                    densMapObj = [],
                                    XYZ        = [],
                                    threshold  = 0.01):

        # retrieve the symmetry operations for the map file.
        # 'densMapObj' is taken from readAtomMap.py (self.densmap).
        # for case where original xyz points are close to a unit
        # cell boundary, include additional adjacent unit cells

        symOps = self.getbasicSymRelatedPoints(densMapObj,XYZ)

        translateUnitCell = []

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
            if self.report:
                print 'Necessary to consider {} adjacent unit cell(s)'.format(l)

        return symOps

    def findVoxRefPoint(self,
                        symOps = [],
                        XYZ    = [],
                        dens   = []):

        # find suitable reference voxel point (out of those assigned to the atom)
        # that is nearest to centre of asym unit. This will be the reference point
        # from to which the distance to over voxels will be determined

        centreOfInterest = np.mean(symOps[0], axis = 1)
        originalXYZs = np.array([[symOps[0][j][ii] for j in range(3)] for ii in range(len(XYZ[0]))])
        dists = list(np.linalg.norm(originalXYZs - np.array(centreOfInterest),axis = 1))
        minInd = dists.index(min(dists))
        refPoint = [XYZ[0][minInd],XYZ[1][minInd],XYZ[2][minInd]] + [dens[minInd]]

        return refPoint

    def decideWhichVoxToKeep(self,
                             refPoint  = [],
                             symOps    = [],
                             X         = [],
                             dens      = [],
                             threshold = 1e-10):

        # loop through all the located voxels assigned to the atom
        # and decide which ones form a single well defined cluster
        # in space. This will then be what is treated as the 
        # voxel set of this particular atom

        keptPoints    = [refPoint]
        symOrOriginal = [0]
        samexyzs      = 0
        dupcount      = 0

        for ii in range(len(X)):

            symxyzs   = [[symOp[j][ii] for j in range(3)] for symOp in symOps]
            symPoints = [sym + [dens[ii]] for sym in symxyzs]
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

            if min(minDist) < threshold:
                samexyzs += 1 

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

        self.decideWhichVoxToKeepReport(keptPoints,symOrOriginal,dupcount,samexyzs)

        return keptPoints,symOrOriginal,samexyzs

    def decideWhichVoxToKeepReport(self,
                                   keptPoints    = [],
                                   symOrOriginal = [],
                                   dupcount      = [], 
                                   samexyzs      = []):

        # print to command line how the previous method did

        if self.report:
            print 'Atom Identity: {}'.format(self.atmId)
            print 'Number points found: {}'.format(len(keptPoints))
            print '\tNumber that are original: {}'.format(len(symOrOriginal) - np.sum(symOrOriginal))
            print '\tNumber that are sym copies: {}'.format(np.sum(symOrOriginal))
            print '\tNumber of duplicate positions: {}'.format(dupcount + samexyzs)
            print '\t\tidentical: {}'.format(dupcount)
            print '\t\tsame xyz: {}'.format(samexyzs)

    def findPosNegClusters(self,
                           keptPoints = []):

        # now that a complete set of voxels have been found around an atom
        # determine the location of clusters of positive and negative density

        posVals = [pt for pt in keptPoints if pt[3] >= 0]
        negVals = [pt for pt in keptPoints if pt[3] < 0]

        negClustInfo = self.findNegClusters(negVals = negVals)
        posClustInfo = self.findPosClusters(posVals = posVals)

        return negClustInfo,posClustInfo

    def findNegClusters(self,
                        n_clusters = 2,
                        negVals    = []):

        # find negative point clusters

        negClustInfo = {}

        if len(negVals) > 1:
            self.foundNegClust = True

            negKmeans    = KMeans(n_clusters = n_clusters)
            negClusters  = negKmeans.fit_predict(negVals)
            negCluster1  = [val for i,val in enumerate(negVals) if negClusters[i] == 0]
            negCluster2  = [val for i,val in enumerate(negVals) if negClusters[i] == 1]
            mean1        = np.mean([c[3] for c in negCluster1])
            mean2        = np.mean([c[3] for c in negCluster2])
            centroids    = negKmeans.cluster_centers_

            if mean2 > mean1:
                mainClust   = negCluster1
                mainCentrID = 0
                clusters    = np.mod(np.array(negClusters)-1,2)
                self.addToOutput(mean1)

            else:
                mainClust   = negCluster2
                mainCentrID = 0
                clusters    = negClusters
                self.addToOutput(mean2)   

            negClustInfo['clusters']      = clusters
            negClustInfo['main cluster']  = mainClust
            negClustInfo['centroids']     = centroids
            negClustInfo['main centroid'] = centroids[mainCentrID, :]
            negClustInfo['mean 1']        = mean1
            negClustInfo['mean 2']        = mean2

        else:
            self.foundNegClust       = False
            negClustInfo['clusters'] = [1]

        return negClustInfo

    def findPosClusters(self,
                        numPosClusts = 10,
                        posVals      = []):

        # find positive point clusters

        posClustInfo = {}

        if len(posVals) > 1:

            # reduce number of clusters if too many specified
            if numPosClusts > len(posVals):
                numPosClusts = len(posVals) 

            self.foundPosClust = True
            posKmeans     = KMeans(n_clusters = numPosClusts)
            clusters   = posKmeans.fit_predict([pt[:3] for pt in posVals])
            posClustList  = []

            for j in range(numPosClusts):
                posClustList.append([val for i,val in enumerate(posVals) if clusters[i] == j])   

            posClustMeans  = [np.mean([c[3] for c in clust]) for clust in posClustList]
            maxPosInd      = posClustMeans.index(max(posClustMeans))

            posClustInfo['clusters']     = clusters
            posClustInfo['centroids']    = posKmeans.cluster_centers_
            posClustInfo['main cluster'] = maxPosInd
            posClustInfo['cluster list'] = posClustList
            posClustInfo['means']        = posClustMeans

        else:
            self.foundPosClust = False
            posClustInfo['clusters'] = [1]

        return posClustInfo

    def findNegToPosDensShift(self,
                              negClustInfo = {},
                              posClustInfo = {},
                              printText    = False):

        # now that clusters have been located, find the shift vectors 
        # which determine the magnitude and direction in which density
        # appears to have shifted. A choice can be made between only 
        # considering one main positive cluster of interest or all of them

        if self.mainVecOnly:
            ind = posClustInfo['main cluster']
            clusterLoop = [posClustInfo['cluster list'][ind]]
        else:
            clusterLoop = posClustInfo['cluster list']

        XYZ = ['x','y','z']
        UVW = ['u','v','w']

        resultantVectors = []
        vectorDic = {}

        for clust in clusterLoop:

            quiverPts  = { i : [] for i in XYZ + UVW }

            for negPt in negClustInfo['main cluster']:

                densDiffs  = []
                directions = []

                for posPt in clust:
                    densDiffs.append(posPt[3]-negPt[3])
                    vectorChange = np.array(posPt[:3])-np.array(negPt[:3])
                    directions.append(vectorChange/np.linalg.norm(vectorChange))

                totalComponent = np.dot(densDiffs,directions)

                for n,xyz in enumerate(XYZ):
                    quiverPts[xyz].append(negPt[n])

                for n,uvw in enumerate(UVW):
                    quiverPts[uvw].append(totalComponent[n])

            resultantVectors.append(np.mean([quiverPts[uvw] for uvw in UVW], axis = 1))

        vMagnitudes = [np.linalg.norm(vec) for vec in resultantVectors]
        
        vectorDic['normed vectors']   = [vec/length for vec,length in zip(resultantVectors,vMagnitudes)]
        vectorDic['weighted vectors'] = [0.03*vec/20 for vec in resultantVectors]
        vectorDic['magnitudes'] = vMagnitudes

        netDensShift = np.sum(resultantVectors,0)

        self.addToOutput(vectorDic['normed vectors'])
        self.addToOutput(vectorDic['magnitudes'])
        self.addToOutput(netDensShift)

        if printText:
            for k in vectorDic.keys():
                print '{}: {}'.format(k,np.around(vectorDic[k],decimals = 2))

        return vectorDic

    def writeVecStartStopToPDBFile(self,
                                   mainNegCentroid = [],
                                   outputPDBname   = 'Dshift-vectors.pdb',
                                   densMapObj      = [],
                                   weightedVectors = []):

        # write located density shift vector start and stop xyz positions 
        # to a pdb file so that they can then be visualised in pymol (or otherwise).
        # 'densMapObj' is taken from readAtomMap.py (self.densmap).

        from classHolder import singlePDB
        from PDBFileManipulation import writePDBline

        vecStart = densMapObj.abc2xyz(asymIndices = mainNegCentroid[:3],
                                      fracInput   = True)

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

        vectorPDBfile = open(outputPDBname,'a')
        
        atmIDparts = self.atmId.split('-')

        vecPosAsAtom = singlePDB(atomnum    = lineCount,
                                 residuenum = atmIDparts[1],
                                 atomtype   = atmIDparts[3],
                                 basetype   = 'V0',
                                 chaintype  = atmIDparts[0],
                                 X_coord    = vecStart[0],
                                 Y_coord    = vecStart[1],
                                 Z_coord    = vecStart[2])

        pdbLine = writePDBline(vecPosAsAtom,1)
        vectorPDBfile.write(pdbLine + '\n')

        for i,vec in enumerate(weightedVectors):
            vecEnd   = densMapObj.abc2xyz(asymIndices = np.array(mainNegCentroid[:3]) + vec,
                                          fracInput   = True)

            # NOTE THAT PYMOL ARROWS SHOULD BE OF THE FORM BELOW
            # print 'cgo_arrow [{}],[{}],hlength=0.25,radius=0.05,hradius=0.075'.
            # format(','.join([str(round(v,2)) for v in vecStart]),','.join([str(round(v,2)) for v in vecEnd]))

            vecPosAsAtom = singlePDB(atomnum    = lineCount + i + 1,
                                     residuenum = atmIDparts[1],
                                     atomtype   = atmIDparts[3],
                                     basetype   = 'V' + str(i+1),
                                     chaintype  = atmIDparts[0],
                                     X_coord    = vecEnd[0],
                                     Y_coord    = vecEnd[1],
                                     Z_coord    = vecEnd[2])

            pdbLine = writePDBline(vecPosAsAtom,1)
            vectorPDBfile.write(pdbLine + '\n')

        vectorPDBfile.close()

    def make3dScatterPlot(self,
                          keptPoints    = [],
                          negClustInfo  = {},
                          posClustInfo  = {},
                          symOrOriginal = [],
                          vectorDic     = {},
                          showSymm     = False,
                          showDens     = False,
                          transparency = 0.6):

        # make a 3d scatter plot of the located voxels around an atom
        # in xyz cartesian space. The colour scheme can take multiple 
        # options (see below for details)

        import matplotlib.pyplot as plt 
        from mpl_toolkits.mplot3d import Axes3D

        xs = [xyz[0] for xyz in keptPoints]
        ys = [xyz[1] for xyz in keptPoints]
        zs = [xyz[2] for xyz in keptPoints]
        ds = [xyz[3] for xyz in keptPoints]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')

        colors = []
        cn,cp = 0,0
        for d in ds:
            if d < 0:
                cmap = plt.cm.Reds
                colors.append(cmap((negClustInfo['clusters'][cn]+1)/float(4)))
                cn += 1
            else:
                cmap = plt.cm.Greens
                colors.append(cmap((posClustInfo['clusters'][cp]+1)/float(len(posClustInfo['cluster list'])+2)))
                cp += 1

        if self.showSymm:
            ax.scatter(xs, ys, zs, c = symOrOriginal, s = 10, alpha = transparency)

        elif self.showDens:
            ax.scatter(xs, ys, zs, c = colors, s = 1000*np.abs(ds), alpha = transparency)

            if self.findClusters:
                if self.plotCentroids:
                    if self.foundNegClust:
                        ax.scatter(negClustInfo['centroids'][:3, 0],
                                   negClustInfo['centroids'][:3, 1],
                                   negClustInfo['centroids'][:3, 2],
                                   marker     = 'x',
                                   s          = 100,
                                   linewidths = 3,
                                   color      = 'k')

                    if self.foundPosClust:
                        ax.scatter(posClustInfo['centroids'][:, 0],
                                   posClustInfo['centroids'][:, 1],
                                   posClustInfo['centroids'][:, 2],
                                   marker     = 'd',
                                   s          = 100,
                                   linewidths = 3,
                                   color      = 'b')

                if self.foundNegClust and self.foundPosClust:
                    for vec,length in zip(vectorDic['normed vectors'],vectorDic['magnitudes']):
                        ax.quiver(negClustInfo['main centroid'][0],
                                  negClustInfo['main centroid'][1],
                                  negClustInfo['main centroid'][2],
                                  vec[0],
                                  vec[1],
                                  vec[2],
                                  length = 0.001*length,
                                  pivot  = 'tail',
                                  color  = 'k',
                                  arrow_length_ratio = 0.1)
    
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.title(self.atmId)
        plt.show()


