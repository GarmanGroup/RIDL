from ccp4Job import ccp4Job, checkInputsExist
from mapTools import mapTools


class MAPMASKjob():

    def __init__(self,
                 mapFile1='', mapFile2='', outputDir='./', runLog=''):

        self.inputMapFile = mapFile1
        self.inputMapFile2 = mapFile2
        self.outputDir = outputDir
        self.runLog = runLog

    def defineCorrectOutputMap(self,
                               switch=False, factorMultiply=False):

        # give correct naming scheme to output map file

        if switch:
            self.outputMapFile = self.inputMapFile.split('.map')[0]
            +'_switchedAxes.map'
        elif not factorMultiply:
            self.outputMapFile = self.inputMapFile.split('.map')[0]
            +'_cropped.map'
        else:
            self.outputMapFile = self.inputMapFile

    def defineCommandInput(self):

        # write the first part of command line input for mapmask run

        ci1 = 'mapmask MAPIN {} '.format(self.inputMapFile)
        if self.inputMapFile2 != '':
            ci1 += 'MAPLIM {} '.format(self.inputMapFile2)
        ci1 += 'MAPOUT {} '.format(self.outputMapFile) +\
               'SYMINFO syminfo.lib '
        self.commandInput1 = ci1

    def switchAxisOrder(self,
                        order=[], symGroup="", includeDir=False):

        # switch the axis order of an input .map file.
        # order = [1,2,3] for example

        self.defineCorrectOutputMap(switch=True)
        xyz = {'1': 'X', '2': 'Y', '3': 'Z'}
        axisOrder = [xyz[str(i)] for i in order]

        self.printPurpose(mode='switch axes', axisOrder=axisOrder)

        inputFiles = [self.inputMapFile]
        if not checkInputsExist(inputFiles, self.runLog):
            return False

        self.defineCommandInput()
        self.commandInput2 = 'SYMMETRY {}\nAXIS {}\nEND'.format(
            symGroup, ' '.join(axisOrder))
        self.outputLogfile = 'MAPMASKlogfile.txt'

        # run MAPMASK job
        job = ccp4Job(jobName='MAPMASK_switchAxisOrder',
                      commandInput1=self.commandInput1,
                      commandInput2=self.commandInput2,
                      outputDir=self.outputDir,
                      outputLog=self.outputLogfile,
                      outputFile=self.outputMapFile)

        self.jobSuccess = job.checkJobSuccess()
        success = self.provideFeedback()
        return success

    def crop2AsymUnit(self,
                      includeDir=False):

        # Crop map 1 to asymmetric unit

        # fillerLine()
        self.printPurpose(mode='crop to asym')
        self.defineCorrectOutputMap()
        inputFiles = [self.inputMapFile]
        if not checkInputsExist(inputFiles, self.runLog):
            return False

        self.defineCommandInput()
        self.commandInput2 = 'EXTEND\nXYZLIM ASU\nEND'
        self.outputLogfile = 'MAPMASKlogfile.txt'

        # run MAPMASK job
        job = ccp4Job(jobName='MAPMASK_crop2AsymUnit',
                      commandInput1=self.commandInput1,
                      commandInput2=self.commandInput2,
                      outputDir=self.outputDir,
                      outputLog=self.outputLogfile,
                      outputFile=self.outputMapFile)

        self.jobSuccess = job.checkJobSuccess()
        success = self.provideFeedback()
        return success

    def cropMap2Map(self,
                    includeDir=False):

        # Crop map 1 to map 2

        # fillerLine()
        self.printPurpose(mode='crop to map')
        self.defineCorrectOutputMap()
        inputFiles = [self.inputMapFile, self.inputMapFile2]
        if not checkInputsExist(inputFiles, self.runLog):
            return False

        maps = [self.inputMapFile, self.inputMapFile2]
        if not includeDir:
            maps = [m.split('/')[-1] for m in maps]

        self.defineCommandInput()
        self.commandInput2 = 'EXTEND\nXYZLIM MATCH\nEND'
        self.outputLogfile = 'MAPMASKlogfile.txt'

        # run MAPMASK job
        job = ccp4Job(jobName='MAPMASK_cropMap2Map',
                      commandInput1=self.commandInput1,
                      commandInput2=self.commandInput2,
                      outputDir=self.outputDir,
                      outputLog=self.outputLogfile,
                      outputFile=self.outputMapFile)

        self.jobSuccess = job.checkJobSuccess()
        success = self.provideFeedback()
        return success

    def multipleByFactor(self,
                         includeDir=False, factor=-1.0, symGroup=""):

        # multiple all points in a density
        # map file by a constant factor

        self.printPurpose(mode='Multiply by factor')
        self.defineCorrectOutputMap(factorMultiply=True)
        inputFiles = [self.inputMapFile]
        if not checkInputsExist(inputFiles, self.runLog):
            return False

        self.defineCommandInput()

        self.commandInput2 = 'SYMMETRY {}\nSCALE -\nFACTOR {} 0.0\nEND'.format(
            symGroup, factor)
        self.outputLogfile = 'MAPMASKlogfile.txt'

        # run MAPMASK job
        job = ccp4Job(jobName='MAPMASK_multipyByFactor',
                      commandInput1=self.commandInput1,
                      commandInput2=self.commandInput2,
                      outputDir=self.outputDir,
                      outputLog=self.outputLogfile,
                      outputFile=self.outputMapFile)

        self.jobSuccess = job.checkJobSuccess()
        success = self.provideFeedback()
        return success

    def provideFeedback(self,
                        includeMapInfo=False):

        # get feedback if successful job

        if self.jobSuccess:
            if includeMapInfo:
                self.feedback()
            return True
        else:
            err = 'Job did not run successfully, see job log file "{}"'.format(
                self.outputLogfile)
            self.runLog.writeToLog(err)
            return False

    def feedback(self,
                 includeDir=False):

        # provide some feedback

        if not includeDir:
            fileIn1 = self.inputMapFile.split('/')[-1]
            if self.inputMapFile2 != '':
                fileIn2 = self.inputMapFile2.split('/')[-1]
            fileOut = self.outputMapFile.split('/')[-1]
        else:
            fileIn1 = self.inputMapFile
            if self.inputMapFile2 != '':
                fileIn2 = self.inputMapFile2
            fileOut = self.outputMapFile

        txt = 'MAPMASK Summary:\n' +\
              'Input map file: {}\n'.format(fileIn1)
        if self.inputMapFile2 != '':
            txt += 'Input map 2 file: {}\n'.format(fileIn2)
        txt += 'Output map file: {}'.format(fileOut)
        self.runLog.writeToLog(txt)

        Map = mapTools(mapName=self.outputMapFile,
                       logFile=self.runLog)
        Map.printMapInfo()

    def printPurpose(self,
                     include=True, mode='switch axes',
                     includeDir=False, axisOrder=[]):

        # provide a summary of what this does
        # (within RIDL) to the command line

        maps = [self.inputMapFile, self.inputMapFile2]

        if not includeDir:
            maps = [m.split('/')[-1] for m in maps]

        if mode == 'switch axes':
            ln = 'Switching map "{}" axis ordering to {}'.format(
                maps[0], axisOrder)
        if mode == 'crop to asym':
            ln = 'Cropping map "{}" to asym unit'.format(maps[0])
        if mode == 'crop to map':
            ln = 'Cropping map "{}" to map "{}"'.format(*maps)
        if mode == 'Multiply by factor':
            ln = 'Multiple map by a constant factor'
        self.runLog.writeToLog(ln)
