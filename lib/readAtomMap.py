from residueFormatter import densper_resatom_NOresidueclass,densper_res
from vxlsPerAtmAnalysisPlots import plotVxlsPerAtm
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
                 calcFCmap   = False):

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
        vxl_list = {atm:[] for atm in self.atmmap.vxls_val}

        for atm,dens in zip(self.atmmap.vxls_val,self.densmap.vxls_val):
            vxl_list[atm].append(dens)
        self.vxlsPerAtom = vxl_list

        if useFCmap is True:
            vxl_list2 = {atm:[] for atm in self.atmmap.vxls_val}
            for atm,dens in zip(self.atmmap.vxls_val,self.FCmap.vxls_val):
                vxl_list2[atm].append(dens)
            self.FCperAtom = vxl_list2

        # delete atmmap and densmap now to save memory
        self.densmap,self.atmmap,self.FCmap =[],[],[]
        self.stopTimer()
 
    def plotDensHistPlots(self,
                          getVoxelStats = False):

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

        self.stopTimer()

    def calculateDensMetrics(self,
                             FCweighted = True,
                             plotDistn  = False):

        # determine density summary metrics per atom

        self.startTimer()
        self.printStepNumber()
        self.lgwrite(ln = 'Calculating electron density statistics per atom...')

        for atom in self.PDBarray:

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
                    if plotDistn is True:
                        self.plotFCdistnPlot(atomOfInterest    = atom,
                                             atomFCvals        = atomFCvals,
                                             atomFCvalsMaxNorm = atomFCvalsMaxNormalised,
                                             FCatMin           = FCatMin,
                                             reliability       = reliability)

        self.success()
        self.stopTimer()

        # delete the vxlsPerAtom list now to save memory
        del self.vxlsPerAtom

        # get additional metrics per atom
        for atom in self.PDBarray:
            atom.getAdditionalMetrics()

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
                fig.savefig('testDistnPlot-{}{}'.format(atomOfInterest.getAtomID(),plotType))

    def plotDensScatterPlots(self,
                             printText = False):

        # plot scatter plots for density metrics 

        self.startTimer()
        self.fillerLine(style = 'line')
        str = 'Plotting scatter plots for electron density statistics...'
        self.lgwrite(ln = str,forcePrint = printText)

        plotVars = (['meandensity','maxdensity'],
                    ['meandensity','mediandensity'],
                    ['meandensity','mindensity'],
                    ['mindensity','maxdensity'],
                    ['meandensity','stddensity'],
                    ['mindensity','min90tile'],
                    ['maxdensity','max90tile'],
                    ['min90tile','min95tile'],
                    ['max90tile','max95tile'])

        for pVars in plotVars:
            logStr = edens_scatter(outputDir = self.filesOut,
                                   metrics   = pVars,
                                   PDBarray  = self.PDBarray,
                                   pdbName   = self.pdbName)
            self.lgwrite(ln=logStr)

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




