# -*- coding: utf-8 -*-

from deleteListIndices import multi_delete 
from map2VoxelClassList import readMap
import sys   
from PDBFileManipulation import PDBtoList
from densityAnalysisPlots import edens_scatter
from residueFormatter import densper_resatom_NOresidueclass,densper_res
import numpy as np
from itertools import izip as zip, count
from vxlsPerAtmAnalysisPlots import plotVxlsPerAtm
import math
from progbar import progress
import time
from savevariables import save_objectlist
from logFile import logFile

class voxel_density():
    # A class for .map file voxel
    def __init__(self,density=0,atmnum=0):
        self.density = density
        self.atmnum = atmnum

class vxls_of_atm():
    # A class for collecting voxels per atom
    def __init__(self,vxls=[],atmnum = 0):
        self.vxls = vxls
        self.atmnum = atmnum

class maps2DensMetrics():
    def __init__(self,filesIn='',filesOut='',pdbName='',mapFileName1='',
                 mapType1='atom_map',mapFileName2='',mapType2='density_map',plot=False):
        self.filesIn = filesIn
        self.filesOut = filesOut # output directory
        self.pdbName = pdbName
        self.map1 = {'filename':mapFileName1,'type':mapType1}
        self.map2 = {'filename':mapFileName2,'type':mapType2}
        self.plot = plot

    def maps2atmdensity(self):
        self.printTitle()
        # write a log file for this eTrack run
        self.log = logFile('{}{}_log-mapProcessing.txt'.format(self.filesOut,self.pdbName))
        self.log.writeToLog('eTrack run - map processing\n')
        self.log.writeToLog('input directory: {}'.format(self.filesIn))        
        self.log.writeToLog('output directory: {}\n'.format(self.filesOut))        

        self.readPDBfile()
        self.readAtomMap()
        self.readDensityMap()
        self.reportDensMapInfo()
        self.checkMapCompatibility()
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
        self.log.writeToLog('Reading in pdb file...')
        self.log.writeToLog('pdb name: {}{}.pdb'.format(self.filesIn,self.pdbName))

        # next read in the pdb structure file:
        # run function to fill PDBarray list with atom objects from structure
        self.PDBarray = PDBtoList('{}{}.pdb'.format(self.filesIn,self.pdbName),[])
        self.success()
        self.stopTimer()

        # want to make sure array of structure atoms ordered by atom number
        # before reading through them
        self.PDBarray.sort(key=lambda x: x.atomnum)
           
        # need to get VDW radius for each atom:
        for atom in self.PDBarray:
            atom.VDW_get()  
        
    def readAtomMap(self):
        # read in the atom map
        self.startTimer()
        self.fillerLine()
        self.log.writeToLog('Reading in Atom map file...')
        self.log.writeToLog('Atom map name: {}{}'.format(self.filesIn,self.map1['filename']))
        self.atmmap,self.atom_indices = readMap(self.filesIn,self.filesOut,self.pdbName,
                                                self.map1['filename'],self.map1['type'],[],self.log)  

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
        self.log.writeToLog('Number of atoms not assigned to voxels: {}'.format(len(Atms_notpres)))

        
    def readDensityMap(self):
        # read in the density map
        self.startTimer()
        self.fillerLine()
        self.log.writeToLog('Reading in Density map file...')
        self.log.writeToLog('Density map name: {}{}'.format(self.filesIn,self.map2['filename']))
        
        self.densmap = readMap(self.filesIn,self.filesOut,self.pdbName,
                               self.map2['filename'],self.map2['type'],
                               self.atom_indices,self.log)  
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

        self.log.writeToLog('For voxels assigned to structure:')
        self.log.writeToLog('mean structure density : {}'.format(structureMean))
        self.log.writeToLog('max structure density : {}'.format(max(self.densmap.vxls_val)))
        self.log.writeToLog('min structure density : {}'.format(min(self.densmap.vxls_val)))
        self.log.writeToLog('std structure density : {}'.format(np.std(self.densmap.vxls_val)))
        self.log.writeToLog('# voxels included : {}'.format(structureNumVxls))
        self.log.writeToLog('For voxels assigned to solvent:')
        self.log.writeToLog('mean solvent-region density : {}'.format(solvMean))
        self.log.writeToLog('# voxels included : {}'.format(solvNumVxls))

    def checkMapCompatibility(self):
        # check that atom-tagged and density map can be combined successfully
        self.fillerLine()
        # print 'Checking that maps have same dimensions and sampling properties...' 
        self.log.writeToLog('Checking that maps have same dimensions and sampling properties...' )

        self.startTimer()
        # Check that the maps have the same dimensions, grid sampling,..
        if (self.atmmap.axis != self.densmap.axis or
            self.atmmap.gridsamp != self.densmap.gridsamp or 
            self.atmmap.start != self.densmap.start or
            self.atmmap.nxyz != self.densmap.nxyz or
            self.atmmap.type != self.densmap.type):

            self.log.writeToLog('Incompatible map properties --> terminating script')
            sys.exit()

        elif self.atmmap.celldims != self.densmap.celldims:
            self.log.writeToLog('Not exact same map grid dimensions..')

            # now check if grid dims same to a specific dp and consider continuing
            stop = True
            for i in list(reversed(range(7))):
                count = 0
                for key in self.atmmap.celldims.keys():
                    if np.round(self.atmmap.celldims[key],i) == np.round(self.densmap.celldims[key],i):
                        count += 1
                if count == 6:
                    self.log.writeToLog('Map grid dimensions same to {}dp --> continuing with processing anyway'.format(i))
                    stop = False
                    break
            if stop == True:
                self.log.writeToLog('Map grid dimensions still not same to 0dp --> terminating script')
                sys.exit()

        else:
            self.success()
            self.log.writeToLog('The atom and density map are of compatible format!')
        self.stopTimer()

        self.fillerLine()
        self.log.writeToLog('Total number of voxels assigned to atoms: {}'.format(len(self.atmmap.vxls_val)))

    def createVoxelList(self):
        # create dictionary of voxels with atom numbers as keys 
        self.startTimer()
        self.fillerLine()
        self.log.writeToLog('Combining voxel density and atom values...')
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
        self.log.writeToLog('Plotting histogram plots of voxels per atom...')
        self.log.writeToLog('Plots written to "{}plots"'.format(self.filesOut))
        for plotType in ('histogram','kde'):
            plotVxlsPerAtm(self.pdbName,self.filesOut,self.vxlsPerAtom,plotType)
        self.stopTimer()

    def calculateDensMetrics(self):
        # determine density summary metrics per atom
        self.fillerLine()
        self.startTimer()
        self.log.writeToLog('Calculating electron density statistics per atom...')
        for atom in self.PDBarray:
            atomVxls = self.vxlsPerAtom[atom.atomnum]
            if len(atomVxls) != 0:
                atom.meandensity    = np.mean(atomVxls)
                atom.mediandensity  = np.median(atomVxls)
                atom.mindensity     = min(atomVxls)
                atom.maxdensity     = max(atomVxls)
                atom.stddensity     = np.std(atomVxls)
                atom.min90tile      = np.percentile(atomVxls,10)
                atom.max90tile      = np.percentile(atomVxls,90)
                atom.min95tile      = np.percentile(atomVxls,5)
                atom.max95tile      = np.percentile(atomVxls,95)
                atom.numvoxels      = len(atomVxls)
        
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
        self.fillerLine()
        # print 'Plotting scatter plots for electron density statistics...'
        self.log.writeToLog('Plotting scatter plots for electron density statistics...')
        plotVars = (['mean','max'],['mean','median'],['mean','min'],['min','max'],
                    ['mean','std'],['std','rsd'],['min','min90tile'],['max','max90tile'],
                    ['min90tile','min95tile'],['max90tile','max95tile'],
                    ['std','range'],['mean','range'])
        for pVars in plotVars:
            edens_scatter(self.filesOut,pVars,self.PDBarray,self.pdbName)

    def plotPerResidueBoxPlots(self):
        # perform residue analysis for datatset, outputting boxplots for each atom specific
        # to each residue, and also a combined boxplot across all residues in structures.
        for densMet in ('mean','min','max'):
            residueArray = densper_resatom_NOresidueclass(self.filesOut,self.PDBarray,'y',densMet,self.pdbName)

        minresnum = 0
        sideormain = ['sidechain','mainchain']
        densper_res(self.filesOut,residueArray,minresnum,sideormain,'min',self.pdbName)

        # remove residueArray now to save memory 
        residueArray = []
        self.stopTimer()

    def pickleAtomList(self):
        self.pklFileName = save_objectlist(self.PDBarray,self.pdbName)

    def startTimer(self):
        self.timeStart = time.time()

    def stopTimer(self):
        elapsedTime = time.time() - self.timeStart
        self.log.writeToLog('section time: {}s\n'.format(round(elapsedTime,3)))

        sys.stdout.flush()

    def success(self):
        self.log.writeToLog('---> success')

    def fillerLine(self):
        print '\n------------------------------------------------'

    def printTitle(self):
        print '\n================================================'
        print '------------------------------------------------'
        print '|||              eTrack run                  |||'
        print '------------------------------------------------'
        print '================================================\n'






