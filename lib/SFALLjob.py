from ccp4Job import ccp4Job, checkInputsExist
from mapTools import mapTools
from errors import error


class SFALLjob():

    def __init__(self,
                 inputPDBfile='', outputDir='./', VDWR=1, symmetrygroup='P1',
                 inputMapfile='', gridDimensions=[], task='atom map',
                 mapoutType='ATMMOD', outputMapFile='', outputMtzFile='',
                 runLog=''):

        self.inputPDBfile = inputPDBfile
        self.inputMapfile = inputMapfile
        self.outputDir = outputDir
        self.symGroup = symmetrygroup
        self.VDWR = VDWR
        self.gridDims = gridDimensions
        self.task = task
        # atom-map 'ATMMOD' or solvent-map 'SOLVMAP' if task is 'atom map'
        self.mapoutType = mapoutType
        self.runLog = runLog

        if task == 'atom map':
            if outputMapFile == '':
                self.outputMapFile = '{}_SFALL.map'.format(
                    inputPDBfile.split('.pdb')[0])
            else:
                self.outputMapFile = outputMapFile
            self.outputFile = self.outputMapFile
        elif task in ('mtz from pdb', 'mtz from map'):
            # specify whether to generate mtz from pdb or map
            if outputMtzFile == '':
                self.outputMtzFile = '{}_SFALL.mtz'.format(
                    inputPDBfile.split('.pdb')[0])
            else:
                self.outputMtzFile = outputMtzFile
            self.outputFile = self.outputMtzFile

    def run(self):
        inputFiles = [self.inputPDBfile]
        if not checkInputsExist(inputFiles, self.runLog):
            return False
        self.runSFALL()
        if self.jobSuccess:
            self.provideFeedback()
            return True
        else:
            error(
                text='Job did not run successfully, see job log file ' +
                     '"{}"'.format(self.outputLogfile),
                log=self.runLog, type='error')
            return False

    def runSFALL(self):

        # run SFALL job using the external ccp4Job class

        self.printPurpose()
        title = 'run of sfall'

        if self.task == 'atom map':
            # use to make an atom-tagged map from a coordinate model
            self.commandInput1 = 'sfall ' +\
                     'XYZIN {} '.format(self.inputPDBfile) +\
                     'ATOMSF atomsf.lib ' +\
                     'MAPOUT {} '.format(self.outputMapFile) +\
                     'SYMINFO syminfo.lib '

            if len(self.gridDims) == 0:
                gridDims = ''
            else:
                gridDims = 'GRID {} {} {}\n'.format(
                    self.gridDims[0], self.gridDims[1], self.gridDims[2])

            self.commandInput2 = 'MODE ATMMAP {}\n'.format(self.mapoutType) +\
                                 'SYMMETRY {}\n'.format(self.symGroup) +\
                                 'VDWR {}\n'.format(self.VDWR) +\
                                 'title {}\n'.format(title) +\
                                 '{}'.format(gridDims) +\
                                 'END'

        elif self.task == 'mtz from pdb':
            # make a new Fcalc column in a new mtz file
            self.commandInput1 = 'sfall ' +\
                     'XYZIN {} '.format(self.inputPDBfile) +\
                     'ATOMSF atomsf.lib ' +\
                     'HKLOUT {} '.format(self.outputMtzFile) +\
                     'SYMINFO syminfo.lib '

            self.commandInput2 = 'title {}\n'.format(title) +\
                                 'labout  FC=FCalc PHIC=PHICalc\n' +\
                                 'NAME -\nPROJECT RIDL -\nCRYSTAL RIDL -\n' +\
                                 'DATASET RIDL\nMODE SFCALC -\nXYZIN\n' +\
                                 'SYMMETRY {}\n'.format(self.symGroup) +\
                                 'BADD 0.0\nVDWR 2.5\nEND'

        elif self.task == 'mtz from map':
            # otherwise make a new Fcalc column in an input map file
            self.commandInput1 = 'sfall ' +\
                     'MAPIN {} '.format(self.inputMapfile) +\
                     'ATOMSF atomsf.lib ' +\
                     'HKLOUT {} '.format(self.outputMtzFile) +\
                     'SYMINFO syminfo.lib '

            self.commandInput2 = 'title {}\n'.format(title) +\
                                 'labout  FC=FCalc PHIC=PHICalc\n' +\
                                 'NAME -\nPROJECT RIDL -\nCRYSTAL RIDL -\n' +\
                                 'DATASET RIDL\nMODE SFCALC -\MAPIN\n' +\
                                 'SYMMETRY {}\n'.format(self.symGroup) +\
                                 'BADD 0.0\nVDWR 2.5\nEND'

        self.outputLogfile = 'SFALLlogfile.txt'

        # run SFALL job
        job = ccp4Job(jobName='SFALL',
                      commandInput1=self.commandInput1,
                      commandInput2=self.commandInput2,
                      outputDir=self.outputDir,
                      outputLog=self.outputLogfile,
                      outputFile=self.outputFile)

        self.jobSuccess = job.checkJobSuccess(self.runLog)

    def provideFeedback(self,
                        includeDir=False):

        # provide some feedback

        if not includeDir:
            fileIn = self.inputPDBfile.split('/')[-1]
            fileOut = self.outputFile.split('/')[-1]
        else:
            fileIn = self.inputPDBfile
            fileOut = self.outputFile

        self.runLog.writeToLog(
            'SFALL Summary:\nInput pdb file: {}\nOutput file: {}'.format(
                fileIn, fileOut))

        if self.task == 'atom map':
            Map = mapTools(mapName=self.outputMapFile, logFile=self.runLog)
            Map.printMapInfo()

    def printPurpose(self,
                     include=True):

        # provide a summary of what this does
        # (within RIDL) to the command line

        if self.task == 'atom map':
            self.runLog.writeToLog(
                'Creating atom-tagged .map file over unit cell for model ' +
                '"{}"'.format(self.inputPDBfile.split('/')[-1]))
        elif self.task == 'mtz from pdb':
            self.runLog.writeToLog(
                'generating new mtz with Fcalc column derived from model ' +
                '"{}"'.format(self.inputPDBfile.split('/')[-1]))
        elif self.task == 'mtz from map':
            self.runLog.writeToLog(
                'generating new mtz from an input map ')