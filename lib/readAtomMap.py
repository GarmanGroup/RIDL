# -*- coding: utf-8 -*-
from residueFormatter import densper_resatom_NOresidueclass,densper_res
from vxlsPerAtmAnalysisPlots import plotVxlsPerAtm
from densityAnalysisPlots import edens_scatter
from deleteListIndices import multi_delete 
from savevariables import save_objectlist
from itertools import izip as zip, count
from PDBFileManipulation import PDBtoList
from map2VoxelClassList import readMap
from progbar import progress
from logFile import logFile
import numpy as np
import sys  
import time
import math

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
    # assign values within a density map (mapFileName2) to specific atoms, using
    # the an atom-tagged map (mapFileName1) to determine which regions of space
    # are to be assigned to each atom

    def __init__(self,
                 filesIn    = '',
                 filesOut   = '',
                 pdbName    = '',
                 atomTagMap = '',
                 densityMap = '',
                 FCmap      = '',
                 plot       = False):

        self.filesIn  = filesIn
        self.filesOut = filesOut # output directory
        self.pdbName  = pdbName
        self.map1     = atomTagMap # atom-tagged map
        self.map2     = densityMap # density map (typically Fo-Fo)
        self.map3     = FCmap # FC map
        self.plot     = plot

    def maps2atmdensity(self):
        self.printTitle()

        # write a log file for this eTrack run
        logName = '{}{}_log-mapProcessing.txt'.format(self.filesIn,self.pdbName)
        self.log = logFile(fileName      = logName,
                           fileDir       = self.filesIn,
                           printToScreen = True)

        self.lgwrite(ln='eTrack run - map processing\n')

        self.lgwrite(ln='input directory: {}'.format(self.filesIn), strip = False)        
        self.lgwrite(ln='output directory: {}\n'.format(self.filesOut), strip = False)        

        self.readPDBfile()
        self.readAtomMap()
        self.readDensityMap()
        self.reportDensMapInfo()
        self.checkMapCompatibility()

        self.readFCMap()


        self.createVoxelList()
        self.plotDensHistPlots()
        self.calculateDensMetrics()

        if self.plot == True:
            self.plotDensScatterPlots()
            self.plotPerResidueBoxPlots()

        self.pickleAtomList()

    def readPDBfile(self):
        # read in pdb file info here
        logfile = open('{}{}_log.txt'.format(self.filesOut,self.pdbName),'a')

        self.startTimer()
        self.lgwrite(ln='Reading in pdb file...')
        self.lgwrite(ln='pdb name: {}{}.pdb'.format(self.filesIn,self.pdbName))

        # read in the pdb file to fill list of atom objects
        self.PDBarray = PDBtoList('{}{}.pdb'.format(self.filesIn,self.pdbName))
        self.success()
        self.stopTimer()

        # make sure array of atoms ordered by atom number
        self.PDBarray.sort(key=lambda x: x.atomnum)
           
        # need to get VDW radius for each atom:
        for atom in self.PDBarray:
            atom.VDW_get()  
        
    def readAtomMap(self):
        # read in the atom map
        self.startTimer()
        self.fillerLine()
        self.lgwrite(ln='Reading in atom-tagged map file...')
        self.lgwrite(ln='Atom map name: {}{}'.format(self.filesIn,self.map1))

        self.atmmap,self.atomIndices = readMap(dirIn   = self.filesIn,
                                               dirOut  = self.filesOut,
                                               mapName = self.map1,
                                               mapType = 'atom_map',
                                               log     = self.log)  
        self.success()
        self.stopTimer()

        # find number of atoms in structure
        num_atoms = len(self.PDBarray)

        # find atom numbers present in list (repeated atom numbers removed)
        seen = set()
        seen_add = seen.add
        uniq_atms = [x for x in self.atmmap.vxls_val if not (x in seen or seen_add(x))] 
        
        # find set of atoms numbers not present (i.e atoms not assigned to voxels)
        Atms_notpres = set(range(1,num_atoms+1)) - set(uniq_atms)
        self.lgwrite(ln='Number of atoms not assigned to voxels: {}'.format(len(Atms_notpres)))

    def readDensityMap(self):
        # read in the density map
        self.startTimer()
        self.fillerLine()
        self.lgwrite(ln='Reading in Density map file...')
        self.lgwrite(ln='Density map name: {}{}'.format(self.filesIn, self.map2))
        
        self.densmap = readMap(dirIn    = self.filesIn,
                               dirOut   = self.filesOut,
                               mapName  = self.map2,
                               mapType  = 'density_map',
                               atomInds = self.atomIndices,
                               log      = self.log)  
        self.success()
        self.stopTimer()

    def readFCMap(self):
        # read in the FC density map
        self.startTimer()
        self.fillerLine()
        self.lgwrite(ln='Reading in FC density map file...')
        self.lgwrite(ln='Density map name: {}{}'.format(self.filesIn, self.map2))

        self.FCmap = readMap(dirIn    = self.filesIn,
                             dirOut   = self.filesOut,
                             mapName  = self.map3,
                             mapType  = 'density_map',
                             atomInds = self.atomIndices,
                             log      = self.log)  

        self.success()
        self.stopTimer()

    def reportDensMapInfo(self):
        # print density map summary information to command line
        totalNumVxls     = np.product(self.atmmap.nxyz.values())
        structureNumVxls = len(self.densmap.vxls_val)
        totalMean        = self.densmap.density['mean']
        structureMean    = np.mean(self.densmap.vxls_val)
        solvNumVxls      = totalNumVxls - structureNumVxls
        solvMean         = (totalNumVxls*totalMean - structureNumVxls*structureMean)/solvNumVxls

        self.lgwrite(ln='For voxels assigned to structure:')
        self.lgwrite(ln='\tmean structure density : {}'.format(round(structureMean,4)))
        self.lgwrite(ln='\tmax structure density : {}'.format(round(max(self.densmap.vxls_val)),4))
        self.lgwrite(ln='\tmin structure density : {}'.format(round(min(self.densmap.vxls_val)),4))
        self.lgwrite(ln='\tstd structure density : {}'.format(round(np.std(self.densmap.vxls_val)),4))
        self.lgwrite(ln='\t# voxels included : {}'.format(structureNumVxls))
        self.lgwrite(ln='\nFor voxels assigned to solvent:')
        self.lgwrite(ln='\tmean solvent-region density : {}'.format(round(solvMean),4))
        self.lgwrite(ln='\t# voxels included : {}'.format(solvNumVxls))

    def checkMapCompatibility(self):
        # check that atom-tagged and density map can be combined successfully
        self.fillerLine()
        # print 'Checking that maps have same dimensions and sampling properties...' 
        self.lgwrite(ln='Checking that maps have same dimensions and sampling properties...' )

        self.startTimer()
        # Check that the maps have the same dimensions, grid sampling,..
        if (self.atmmap.axis != self.densmap.axis or
            self.atmmap.gridsamp != self.densmap.gridsamp or 
            self.atmmap.start != self.densmap.start or
            self.atmmap.nxyz != self.densmap.nxyz or
            self.atmmap.type != self.densmap.type):

            self.lgwrite(ln='Incompatible map properties --> terminating script')
            sys.exit()

        elif self.atmmap.celldims != self.densmap.celldims:
            self.lgwrite(ln='Not exact same map grid dimensions..')

            # now check if grid dims same to a specific dp and consider continuing
            stop = True
            for i in list(reversed(range(7))):
                count = 0
                for key in self.atmmap.celldims.keys():
                    if np.round(self.atmmap.celldims[key],i) == np.round(self.densmap.celldims[key],i):
                        count += 1
                if count == 6:
                    self.lgwrite(ln='Map grid dimensions same to {}dp --> continuing with processing anyway'.format(i))
                    stop = False
                    break
            if stop == True:
                self.lgwrite(ln='Map grid dimensions still not same to 0dp --> terminating script')
                sys.exit()

        else:
            self.success()
            self.lgwrite(ln='The atom and density map are of compatible format!')
        self.stopTimer()

        self.fillerLine()
        self.lgwrite(ln='Total number of voxels assigned to atoms: {}'.format(len(self.atmmap.vxls_val)))

    def createVoxelList(self):
        # create dictionary of voxels with atom numbers as keys 
        self.startTimer()
        self.fillerLine()
        self.lgwrite(ln='Combining voxel density and atom values...')
        vxl_list = {atm:[] for atm in self.atmmap.vxls_val}
        for atm,dens in zip(self.atmmap.vxls_val,self.densmap.vxls_val):
            vxl_list[atm].append(dens)
        self.vxlsPerAtom = vxl_list

        # delete atmmap and densmap now to save memory
        self.densmap,self.atmmap =[],[]
        self.stopTimer()
 
    def plotDensHistPlots(self):
        # histogram & kde plots of number of voxels per atom
        self.startTimer()
        self.lgwrite(ln='Plotting histogram plots of voxels per atom...')
        self.lgwrite(ln='Plots written to "{}plots"'.format(self.filesOut))
        for p in ('histogram','kde'):
            plotVxlsPerAtm(pdbName=self.pdbName,
                           where=self.filesOut,
                           vxlsPerAtom=self.vxlsPerAtom,
                           plotType=p)
        self.stopTimer()

    def calculateDensMetrics(self):
        # determine density summary metrics per atom
        self.fillerLine()
        self.startTimer()
        self.lgwrite(ln='Calculating electron density statistics per atom...')
        for atom in self.PDBarray:
            atomVxls = self.vxlsPerAtom[atom.atomnum]
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
        self.success()
        self.stopTimer()

        # delete the vxlsPerAtom list now to save memory
        del self.vxlsPerAtom

        # get additional metrics per atom
        for atom in self.PDBarray:
            atom.getAdditionalMetrics()

    def plotDensScatterPlots(self):
        # plot scatter plots for density metrics 
        self.startTimer()
        self.fillerLine(style='line')
        # print 'Plotting scatter plots for electron density statistics...'
        self.lgwrite(ln='Plotting scatter plots for electron density statistics...')
        plotVars = (['mean','max'],
                    ['mean','median'],
                    ['mean','min'],
                    ['min','max'],
                    ['mean','std'],
                    ['std','rsd'],
                    ['min','min90tile'],
                    ['max','max90tile'],
                    ['min90tile','min95tile'],[
                    'max90tile','max95tile'],
                    ['std','range'],
                    ['mean','range'])
        for pVars in plotVars:
            edens_scatter(self.filesOut,pVars,self.PDBarray,self.pdbName)

    def plotPerResidueBoxPlots(self):
        # perform residue analysis for datatset, outputting boxplots for each atom specific
        # to each residue, and also a combined boxplot across all residues in structures.
        self.fillerLine(style='line')
        for densMet in ('mean','min','max'):
            resArray = densper_resatom_NOresidueclass(where    = self.filesOut,
                                                      PDBarray = self.PDBarray,
                                                      plot     = True,
                                                      densMet  = densMet,
                                                      pdbName  = self.pdbName)
        densper_res(where    = self.filesOut,
                    resArray = resArray,
                    pdbName  = self.pdbName)

        # remove residueArray now to save memory 
        resArray = []
        self.stopTimer()

    def pickleAtomList(self):
        self.pklFileName = save_objectlist(self.PDBarray,self.pdbName)

    def startTimer(self):
        self.timeStart = time.time()

    def stopTimer(self):
        elapsedTime = time.time() - self.timeStart
        self.lgwrite(ln='section time: {}s\n'.format(round(elapsedTime,3)))
        sys.stdout.flush()

    def success(self):
        self.lgwrite(ln='---> success')

    def fillerLine(self,style='stars'):
        if style == 'stars':
            print '\n***'
        elif style == 'line':
            print '\n'+'-'*30
        elif style == 'blank':
            print '\n'

    def printTitle(self,title='eTrack run'):
        print '\n'+'='*50
        print '-'*50
        print '|||'+' '*15+title+' '*15+'|||'
        print '-'*50
        print '='*50+'\n'

    def lgwrite(self,ln='',strip=True):
        # write line to log file
        self.log.writeToLog(str=ln, strip=strip)




