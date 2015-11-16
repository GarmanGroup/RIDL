# -*- coding: utf-8 -*-

from deleteListIndices import multi_delete 
from map2VoxelClassList import readMap
import sys   
from PDBFileManipulation import PDBtoCLASSARRAY_v2 as pdb2list
from densityAnalysisPlots import edens_scatter
from residueFormatter import densper_resatom_NOresidueclass,densper_res
import numpy as np
from itertools import izip as zip, count
from densAndAtomMaps2VxlList import combinevxl_atmanddensvals
from vxlsPerAtmAnalysisPlots import plotVxlsPerAtm
import math
from progbar import progress
import time
from savevariables import save_objectlist

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
    def __init__(self,where,pdbname,mapfilname1,maptype1,mapfilname2,maptype2,plot):
        self.where = where # output directory
        self.pdbname = pdbname
        self.map1 = {'filename':mapfilname1,'type':maptype1}
        self.map2 = {'filename':mapfilname2,'type':maptype2}
        self.plot = plot

    def maps2atmdensity(self):
        print '\n================================================'
        print '------------------------------------------------'
        print '|||              eTrack run                  |||'
        print '------------------------------------------------'
        print '================================================\n'

        # write a log file for this eTrack run
        logfile = open('{}output/{}_log.txt'.format(self.where,self.pdbname),'w')
        logfile.write('eTrack run log file\n')
        logfile.write('Date: '+ str(time.strftime("%d/%m/%Y"))+'\n')
        logfile.write('Time: '+ str(time.strftime("%H:%M:%S"))+'\n')

        self.readPDBfile()
        self.readAtomMap()
        self.readDensityMap()
        self.checkMapCompatibility()
        self.createVoxelList()
        self.assignVoxels2Atoms()
        self.plotDensHistPlots()
        self.calculateDensMetrics()
        if self.plot == True:
            self.plotDensScatterPlots()
            self.plotPerResidueBoxPlots()
        self.pickleAtomList()

    def readPDBfile(self):
        # read in pdb file info here
        logfile = open('{}output/{}_log.txt'.format(self.where,self.pdbname),'a')
        self.startTimer()
        print 'Reading in pdb file...'
        print 'pdb name: {}{}.pdb'.format(self.where,self.pdbname)
        logfile.write('pdb name: {}{}.pdb\n'.format(self.where,self.pdbname))

        # next read in the pdb structure file:
        # run function to fill PDBarray list with atom objects from structure
        self.PDBarray = pdb2list('{}{}.pdb'.format(self.where,self.pdbname),[])
        self.success()
        self.stopTimer()

        # want to make sure array of structure atoms ordered by atomnumber
        # before reading through them
        self.PDBarray.sort(key=lambda x: x.atomnum)
           
        # need to get VDW radius for each atom:
        for atom in self.PDBarray:
            atom.VDW_get()  
        
    def readAtomMap(self):
        # read in the atom map
        logfile = open('{}output/{}_log.txt'.format(self.where,self.pdbname),'a')

        self.startTimer()
        self.fillerLine()
        print 'Reading in Atom map file...'
        print 'Atom map name: {}{}'.format(self.where,self.map1['filename'])
        logfile.write('atom map name: {}{}\n'.format(self.where,self.map1['filename']))

        self.atmmap,self.atom_indices = readMap(self.where,self.pdbname,self.map1['filename'],self.map1['type'],[])  

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
        print 'Number of atoms not assigned to voxels: %s' %str(len(Atms_notpres))

        # append to log file for this eTrack run
        logfile.write('Number of atoms not assigned to voxels: %s\n' %str(len(Atms_notpres)))
        
    def readDensityMap(self):
        # read in the density map
        logfile = open('{}output/{}_log.txt'.format(self.where,self.pdbname),'a')

        self.startTimer()
        self.fillerLine()
        print 'Reading in Density map file...'
        print 'Density map name: {}{}'.format(self.where,self.map2['filename'])
        logfile.write('density map name: {}{}\n'.format(self.where,self.map2['filename']))

        self.densmap = readMap(self.where,self.pdbname,self.map2['filename'],self.map2['type'],self.atom_indices)  

        self.success()
        self.stopTimer()

    def checkMapCompatibility(self):
        # check that atom-tagged and density map can be combined successfully
        logfile = open('{}output/{}_log.txt'.format(self.where,self.pdbname),'a')

        self.fillerLine()
        print 'Checking that maps have same dimensions and sampling properties...' 
        self.startTimer()
        # Check that the maps have the same dimensions, grid sampling,..
        if (self.atmmap.axis != self.densmap.axis or
            self.atmmap.celldims != self.densmap.celldims or
            self.atmmap.gridsamp != self.densmap.gridsamp or 
            self.atmmap.start != self.densmap.start or
            self.atmmap.nxyz != self.densmap.nxyz or
            self.atmmap.type != self.densmap.type):
        
            print 'Incompatible map properties --> terminating script'
            logfile.write('Incompatible map properties --> terminating script\n')
            sys.exit()
        else:
            self.success()
            print 'The atom and density map are of compatible format!'
            logfile.write('The atom and density map are of compatible format!\n')
        self.stopTimer()

        self.fillerLine()
        print 'Total number of voxels assigned to atoms: %s' %str(len(self.atmmap.vxls_val))
        logfile.write('Total number of voxels assigned to atoms: %s\n' %str(len(self.atmmap.vxls_val)))
        logfile.close()

    def createVoxelList(self):
        # create list of voxel objects in class voxel_density 
        self.startTimer()
        self.fillerLine()
        print 'Combining voxel density and atom values...'
        self.vxlList = combinevxl_atmanddensvals(self.atmmap.vxls_val,self.densmap.vxls_val)
        
        # delete atmmap and densmap now to save memory
        self.densmap,self.atmmap =[],[]
        self.stopTimer()

    def assignVoxels2Atoms(self):
        # next assign voxels to atom objects in PDBarray list
        self.startTimer()
        self.fillerLine()
        print 'Assigning voxels to corresponding atoms...'
        # first order voxels by atom number
        self.vxlList.sort(key=lambda x: x.atmnum)
        
        # initialise list of voxels for each atom before can append voxels. Note 
        # that voxels are saved in list attribute of local vxls_of_atom class objects
        self.PDBarray.sort(key=lambda x: x.atomnum)   
        self.stopTimer()

        self.startTimer()
        self.vxlsPerAtom = {}
        for atm in self.PDBarray:
          self.vxlsPerAtom[atm.atomnum] = []
        self.stopTimer()

        # assign voxels to each atom within a dictionary
        self.startTimer()
        for vxl in self.vxlList:
          self.vxlsPerAtom[vxl.atmnum].append(vxl)

        # can remove vxlList variable now to save memory
        self.vxlList = []   
        self.stopTimer()
 
    def plotDensHistPlots(self):
        # histogram & kde plots of number of voxels per atom
        for plotType in ('histogram','kde'):
            plotVxlsPerAtm(self.pdbname,self.where,self.vxlsPerAtom,plotType)

    def calculateDensMetrics(self):
        # determine density summary metrics per atom, including:
        # max, min, mean, median, standard deviation, 90-tile min,
        # 90-tile max, 95-tile min, 95-tile max, mode (why not!), 
        # relative standard deviation (rsd = std/mean)
        self.startTimer()
        self.fillerLine()
        print 'Calculating electron density statistics per atom...'
        for atom in self.PDBarray:
            if len(self.vxlsPerAtom[atom.atomnum]) != 0:
                atomVxls = [vxl.density for vxl in self.vxlsPerAtom[atom.atomnum]]
                atom.meandensity = np.mean(atomVxls)
                atom.mediandensity = np.median(atomVxls)
                atom.mindensity = min(atomVxls)
                atom.maxdensity = max(atomVxls)
                atom.stddensity = np.std(atomVxls)
                atom.min90tile = np.percentile(atomVxls,10)
                atom.max90tile = np.percentile(atomVxls,90)
                atom.min95tile = np.percentile(atomVxls,5)
                atom.max95tile = np.percentile(atomVxls,95)
                atom.modedensity = max(set(atomVxls), key=atomVxls.count)
                atom.rsddensity = float(atom.stddensity)/atom.meandensity
                atom.rangedensity = np.linalg.norm(atom.maxdensity - atom.mindensity)
                atom.numvoxels = len(atomVxls)
        
        # check that the last step has worked
        # (may not be necessary to have but good to check when testing!)
        for atom in self.PDBarray:
            for vxl in self.vxlsPerAtom[atom.atomnum]:
                if atom.atomnum != vxl.atmnum:
                    print 'error!'
                    return
         
        self.success()
        self.stopTimer()

        # delete the vxlsPerAtom list now to save memory
        del self.vxlsPerAtom

    def plotDensScatterPlots(self):
        # plot scatter plots for density metrics 
        self.startTimer()
        self.fillerLine()
        print 'Plotting scatter plots for electron density statistics...'
        plotVars = (['mean','max'],['mean','median'],['mean','min'],['min','max'],
                    ['mean','std'],['std','rsd'],['min','min90tile'],['max','max90tile'],
                    ['min90tile','min95tile'],['max90tile','max95tile'],['mean','mode'],
                    ['std','range'],['mean','range'])
        for pVars in plotVars:
            edens_scatter(self.where,pVars,self.PDBarray,self.pdbname)

    def plotPerResidueBoxPlots(self):
        # perform residue analysis for datatset, outputting boxplots for each atom specific
        # to each residue, and also a combined boxplot across all residues in structures.
        for densMet in ('mean','min','max'):
            residueArray = densper_resatom_NOresidueclass(self.where,self.PDBarray,'y',densMet,self.pdbname)

        minresnum = 0
        sideormain = ['sidechain','mainchain']
        densper_res(self.where,residueArray,minresnum,sideormain,'min',self.pdbname)

        # remove residueArray now to save memory 
        residueArray = []
        self.stopTimer()

    def pickleAtomList(self):
        self.pklFileName = save_objectlist(self.PDBarray,self.pdbname)

    def startTimer(self):
        self.timeStart = time.time()

    def stopTimer(self):
        elapsedTime = time.time() - self.timeStart
        print elapsedTime
        sys.stdout.flush()

    def success(self):
        print '---> success'

    def fillerLine(self):
        print '\n------------------------------------------------'






