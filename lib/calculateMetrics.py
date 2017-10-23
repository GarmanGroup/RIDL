from combinedAtomList import combinedAtomList
from savevariables import retrieve_objectlist
from savevariables import saveGenericObject
from PDBFileManipulation import PDBtoList
from mapsToDensityMetrics import maps2DensMetrics
from shutil import move
from os import path, makedirs, listdir, remove, rmdir


class calculateMetrics(object):

    # a class for retrieving the RIDL input text file information and running
    # the RIDL pipeline to calculate per-atom damage metrics for a specified
    # damage series. This code requires previously generated atom-tagged and
    # density maps (typically Fobs(n) - Fobs(1) Fourier difference maps) to
    # have been generated for the current damage series (as specified within
    # the input .txt file parsed below). If run as part of the full RIDL
    # pipeline (by running 'python runRIDL.py -i [inputfile.txt] -pc') then
    # this will automatically run directly after the suitable map files have
    # been generated, with no need to explicitly write a new input file for
    # this class to work.

    def __init__(self,
                 mapDir='./', outDir='./', pklFiles=[],
                 initialPDB="", seriesName="untitled-series", pklDataFile="",
                 doses=[], plot='no', output='simple', logFile='',
                 inclFCmets=True, laterDatasets=[], autoRun=True,
                 RIDLinputFile='untitled.txt'):

        # the input map file directory
        self.mapDir = mapDir

        # the output file directory
        self.outDir = outDir

        # list of pkl files from map_processing
        self.pklFiles = pklFiles

        # the first dataset pdb code
        self.initialPDB = initialPDB

        # the list of later dataset pdb files
        self.laterDatasets = laterDatasets.split(',')

        # the general series name
        self.seriesName = seriesName

        # the combined series pkl file from post_processing
        self.pklDataFile = pklDataFile

        # list of increasing doses
        self.doses = [float(d) for d in doses.split(',')]

        # input file for main RIDL job (this can be appended to here)
        self.RIDLinputFile = RIDLinputFile

        # (bool) decide to plot per-residue summary plots per dataset
        self.plot = plot

        # decide whether to plot v large heatmap for metrics (time consuming)
        self.plotHeatMaps = False

        # the amount of output to provide (either 'simple' for
        # just Dloss or 'full' larger selection of output files)
        self.output = output

        # log file for the current RIDL job
        self.logFile = logFile

        # generate metrics which require FC maps to be made
        self.inclFCmets = inclFCmets

        # if number of initial datasets given doesn't match
        # number of later datasets, assume same initial dataset
        # used for every later dataset (fix as first one given)
        initialPDBs = self.initialPDB.split(',')
        numLaterDsets = len(self.laterDatasets)
        if len(initialPDBs) != numLaterDsets:
            initialPDBs = [initialPDBs[0]]*numLaterDsets
        l = []
        for pdb in initialPDBs:
            if not pdb.endswith('.pdb'):
                l.append(pdb+'.pdb')
            else:
                l.append(pdb)
        self.initialPDB = l

        if autoRun:
            self.runPipeline()

    def runPipeline(self,
                    map_process=True, post_process=True):

        # the following function reads in the above functions
        # one by one in a scripted pipeline. Takes inputs to
        # specify which parts of above pipeline are to be
        # included. For each function input specify True if
        # this part is to be performed and False otherwise.

        if map_process or post_process:
            success = self.checkInOutDirExist()
        if not success:
            return

        self.setOutputDirs()

        if map_process:
            self.map_processing()
        else:
            self.logFile.writeToLog(str='Map processing task not chosen...')
        self.fillerLine()

        if post_process:
            self.post_processing()

        else:
            self.logFile.writeToLog(str='Post processing job not chosen...')

        self.fillerLine(blank=True)

    def checkInOutDirExist(self,
                           checkMapdir=True, checkOutdir=True):

        # check that an input/output directories have been
        # found and make subdirectories if present

        if checkMapdir:
            if not path.isdir(self.mapDir):
                self.logFile.writeToLog(
                    str='Input file location: {} does '.format(self.mapDir) +
                        'not exist. Please select an appropriate directory.')
                return False

        if checkOutdir:
            if not path.isdir(self.outDir):
                self.logFile.writeToLog(
                    str='Output file location: {} does '.format(self.outDir) +
                        'not exist. Please select an appropriate directory.')
                return False

        return True

    def makeOutputDir(self,
                      dirName='./'):

        # if the above sub directory does not exist, make it

        if not path.exists(dirName):
            makedirs(dirName)
            self.logFile.writeToLog(
                str='New sub directory ' +
                    '"{}" '.format(dirName.replace(self.outputDir, '')) +
                    'created to contain output files')

    def setOutputDirs(self):

        # set the locations of the output directories

        self.outputDir = '{}RIDL-metrics/'.format(self.outDir)

        # add pkl file names as attribute if specified in input file
        if len(self.pklFiles) != 0:
            self.pklFiles = [self.outputDir+f for f in self.pklFiles]

    def map_processing(self):

        # combine the density map and atom-tagged map for a given dataset,
        # to calculate per-atom density metrics for each refined atom

        self.logFile.writeToLog(
            str='Combining density maps and atom-tagged maps to calculate ' +
                'per-atom density metrics for each refined atom.\n')

        self.logFile.writeToLog(
            str='input directory:  {}\n'.format(self.mapDir) +
                'output directory: {}'.format(self.outputDir),
            strip=False)

        # create additional subdirectories
        self.makeOutputDir(dirName=self.outputDir)

        # make pklFiles and dir to move all generated
        # per-dataset pkl files to this
        pklFileDir = 'pklFiles-perDataset/'
        self.makeOutputDir(dirName='{}{}'.format(self.outputDir, pklFileDir))

        pklFileNames = []
        i = 0
        for d, initPDB in zip(self.laterDatasets, self.initialPDB):
            i += 1
            # derive per-atom density metrics from maps
            mapName1 = '{}_atoms.map'.format(d)
            mapName2 = '{}_density.map'.format(d)
            # mapName3 = initPDB.replace('.pdb','_FC.map')
            mapName3 = '{}_FC.map'.format(d)

            self.logFile.writeToLog(
                str='\n---------------------------------\n' +
                    'Higher dose dataset {} starts here'.format(i))

            maps2DensMets = maps2DensMetrics(filesIn=self.mapDir,
                                             filesOut=self.outputDir,
                                             pdbName=d, atomTagMap=mapName1,
                                             densityMap=mapName2,
                                             FCmap=mapName3,
                                             plotHist=self.plot,
                                             logFile=self.logFile,
                                             calcFCmap=self.inclFCmets)

            maps2DensMets.maps2atmdensity()

            # move pkl file to working output directory
            pklFileName = maps2DensMets.pklFileName

            move(pklFileName,
                 '{}{}{}'.format(self.outputDir, pklFileDir, pklFileName))

            pklFileNames.append(
                '{}{}{}'.format(self.outputDir, pklFileDir, pklFileName))

        self.pklFiles = pklFileNames

    def post_processing(self):

        # group the per-atom density metrics for each dataset together

        self.logFile.writeToLog(
            str='Combining density metric information for each dataset ' +
                'together within the damage series')

        txt = 'Input pkl files for post processing chosen from input file:'
        for file in self.pklFiles:
            txt += '\n\t{}'.format(file.replace(self.outDir, ""))
        self.logFile.writeToLog(str=txt)

        # next read in the pdb structure file as list of atom objects
        initialPDBlist = PDBtoList(pdbFileName=self.get1stDsetPDB())

        # retrieve object lists of atoms for each damage set
        ln = '\nReading in pkl files for higher dataset structures...'
        self.logFile.writeToLog(str=ln)

        dList = []
        for pkl_filename in self.pklFiles:
            ln = 'Damage file number: {}'.format(len(dList)+1)
            self.logFile.writeToLog(str=ln)
            PDB_ret = retrieve_objectlist(fileName=pkl_filename,
                                          logFile=self.logFile)
            # remove pkl file since no longer needed
            remove(pkl_filename)

            # add new retrieved damage set list to dList
            dList.append(PDB_ret)
        pklDir = '/'.join(pkl_filename.split('/')[:-1])+'/'

        # remove directory if it is now empty
        if listdir(pklDir) == []:
            rmdir(pklDir)

        # create a list of atom objects with attributes as lists varying over
        # dose range, only including atoms present in ALL damage datasets
        self.logFile.writeToLog(
            str='New list of atoms over full dose range calculated...')
        combinedAtoms = combinedAtomList(datasetList=dList,
                                         numLigRegDsets=len(dList),
                                         doseList=self.doses,
                                         initialPDBList=initialPDBlist,
                                         outputDir=self.outputDir,
                                         seriesName=self.seriesName,
                                         inclFCmetrics=self.inclFCmets)

        combinedAtoms.getMultiDoseAtomList()

        # calculate 'average' variant Dloss metrics
        combinedAtoms.calcAdditionalMetrics(newMetric='average')

        # calculate Calpha normalised metrics, if Calpha atoms exist
        if self.checkCalphasPresent(atomObjList=combinedAtoms):

            metricsOfInterest = ['loss', 'mean', 'gain', 'Bfactor']

            if self.inclFCmets:
                metricsOfInterest += ['density_weighted_mean_negOnly',
                                      'density_weighted_loss']

            for m in metricsOfInterest:
                combinedAtoms.calcAdditionalMetrics(metric=m)

        self.combinedAtoms = combinedAtoms

        # save metric data to pkl file
        pklDataFile = saveGenericObject(obj=self.combinedAtoms,
                                        fileName=self.seriesName)

        move(pklDataFile, '{}{}'.format(self.outputDir, pklDataFile))
        self.pklDataFile = pklDataFile

        inputfile = open(self.RIDLinputFile, 'a')
        inputfile.write('\npklDataFile ' + pklDataFile)
        inputfile.close()

    def checkCalphasPresent(self,
                            atomObjList=[]):

        # check whether structure contains any Calpha
        # protein backbone atoms within it

        return atomObjList.checkCalphaAtomsExist()

    def get1stDsetPDB(self):

        # retrieve name of first dataset pdb coordinate file.
        # If multiple initial datasets input, take the first only

        pdbFile = self.mapDir + self.initialPDB[0]

        return pdbFile

    def fillerLine(self,
                   blank=False):

        # print a filler line to command line

        if not blank:
            ln = '\n----------------------------------------------------------'
        else:
            ln = '\n'
        self.logFile.writeToLog(str=ln)
