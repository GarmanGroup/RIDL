import os.path
from CADjob import CADjob
from SCALEITjob import SCALEITjob
from logFile import logFile
from SIGMAAjob import SIGMAAjob
import shutil


class pipeline():

    # class to run CAD job to combine F and SIGF columns from
    # two merged mtz files, then to scale the 2nd datasets F
    # structure factors against the 1st datasets

    def __init__(self,
                 outputDir='', jobName='untitled-job', log='',
                 mtzIn1='./untitled.mtz', Mtz1LabelName='FP',
                 Mtz1SIGFPlabel='SIGFP', RfreeFlag1='FREE',
                 Mtz1LabelRename='D1', mtzIn2='./untitled2.mtz',
                 Mtz2LabelName='FP',  Mtz2SIGFPlabel='SIGFP',
                 Mtz2LabelRename='D2', mtzIn3='./untitled3.mtz',
                 Mtz3phaseLabel='PHIC', Mtz3FcalcLabel='FC',
                 Mtz3LabelRename='DP', inputPDBfile='./untitled.pdb',
                 densMapType='DIFF', scaleType='ANISOTROPIC',
                 deleteMtzs=True, FOMweight='NONE'):

        # specify where output files should be written
        self.outputDir = outputDir
        self.makeOutputDir(dirName=self.outputDir)
        self.findFilesInDir()

        self.jobName = jobName
        self.mtzIn1 = mtzIn1
        self.Mtz1LabelName = Mtz1LabelName
        self.Mtz1SIGFPlabel = Mtz1SIGFPlabel
        self.RfreeFlag1 = RfreeFlag1
        self.Mtz1LabelRename = Mtz1LabelRename
        self.mtzIn2 = mtzIn2
        self.Mtz2LabelName = Mtz2LabelName
        self.Mtz2SIGFPlabel = Mtz2SIGFPlabel
        self.Mtz2LabelRename = Mtz2LabelRename
        self.mtzIn3 = mtzIn3
        self.Mtz3phaseLabel = Mtz3phaseLabel
        self.Mtz3FcalcLabel = Mtz3FcalcLabel
        self.Mtz3LabelRename = Mtz3LabelRename
        self.inputPDBfile = inputPDBfile
        self.densMapType = densMapType
        self.scaleType = scaleType
        self.deleteMtzs = deleteMtzs
        self.FOMweight = FOMweight

        if log == '':
            f = '{}{}_runLog1.log'.format(self.outputDir+'RIDL-log', jobName)
            self.runLog = logFile(fileName=f,
                                  fileDir=self.outputDir)
        else:
            self.runLog = log

        # specify output files for parts of pipeline
        self.CADoutputMtz = '{}{}_CADcombined.mtz'.format(
            self.outputDir, self.jobName)
        self.SCALEITinputMtz = self.CADoutputMtz
        self.SCALEIToutputMtz = '{}{}_SCALEITcombined.mtz'.format(
            self.outputDir, self.jobName)

    def makeOutputDir(self,
                      dirName='./'):

        # if the above sub directory does not exist, make it

        if not os.path.exists(dirName):
            os.makedirs(dirName)
            self.runLog.writeToLog(
                str='New sub directory "{}" '.format(dirName) +
                    'created to contain output files')

    def runPipeline(self):

        # run the current subroutine to parse an input
        # file and run CAD and SCALEIT to combine the
        # mtz information for a low and high dose
        # dataset within a damage series

        # copy input mtz files to working directory and rename
        self.moveInputMtzs()

        # run SIGMAA job if required to generate a new FOM weight column
        if self.FOMweight == 'recalculate':

            if self.densMapType == '2FOFC':
                mtzLbls_in = self.Mtz2LabelName
                mtzLbls_out = self.Mtz2LabelRename

            else:
                mtzLbls_in = self.Mtz1LabelName
                mtzLbls_out = self.Mtz1LabelName

            self.printStepNumber()
            sigmaa = SIGMAAjob(inputMtz=self.SIGMAAinputMtz,
                               MtzLabelNameIn=mtzLbls_in,
                               MtzLabelNameOut=mtzLbls_out,
                               RfreeFlag=self.RfreeFlag1,
                               inputPDB=self.inputPDBfile,
                               outputDir=self.outputDir,
                               runLog=self.runLog)
            success = sigmaa.run()

            if not success:
                return False

            # if 2FO-FC map required, use FWT column from
            # sigmaa-output mtz (we are done here)
            if self.densMapType == '2FOFC':
                self.cleanUpDir()
                return True

            self.CADinputMtz1 = sigmaa.outputMtz
        else:
            self.CADinputMtz1 = self.SIGMAAinputMtz

        # run CAD job
        self.printStepNumber()
        cad = CADjob(inputMtz1=self.CADinputMtz1,
                     inputMtz2=self.CADinputMtz2,
                     inputMtz3=self.CADinputMtz3,
                     Mtz1LabelName=self.Mtz1LabelName,
                     Mtz2LabelName=self.Mtz2LabelName,
                     Mtz3phaseLabel=self.Mtz3phaseLabel,
                     Mtz3FcalcLabel=self.Mtz3FcalcLabel,
                     Mtz1LabelRename=self.Mtz1LabelRename,
                     Mtz2LabelRename=self.Mtz2LabelRename,
                     Mtz3LabelRename=self.Mtz3LabelRename,
                     outputMtz=self.CADoutputMtz,
                     outputDir=self.outputDir,
                     runLog=self.runLog,
                     FOMWeight=self.FOMweight)
        success = cad.run()

        if not success:
            return False

        # run SCALEIT job
        if self.scaleType != 'NONE':
            self.printStepNumber()
            scaleit = SCALEITjob(inputMtz=self.SCALEITinputMtz,
                                 outputMtz=self.SCALEIToutputMtz,
                                 Mtz1Label=self.Mtz1LabelRename,
                                 Mtz2Label=self.Mtz2LabelRename,
                                 outputDir=self.outputDir,
                                 scaling=self.scaleType,
                                 runLog=self.runLog)
            success = scaleit.run()

            if not success:
                return False

        # end of pipeline reached
        self.cleanUpDir()
        return True

    def moveInputMtzs(self):

        # move input mtz files to working directory and rename as suitable

        if self.densMapType == '2FOFC':
            self.SIGMAAinputMtz = '{}{}.mtz'.format(
                self.outputDir, self.Mtz2LabelRename.strip())
            shutil.copy2(self.mtzIn2, self.SIGMAAinputMtz)

        else:
            self.SIGMAAinputMtz = '{}{}.mtz'.format(
                self.outputDir, self.Mtz1LabelRename.strip())
            self.CADinputMtz2 = '{}{}.mtz'.format(
                self.outputDir, self.Mtz2LabelRename.strip())
            self.CADinputMtz3 = '{}{}.mtz'.format(
                self.outputDir, self.Mtz3LabelRename.strip())
            shutil.copy2(self.mtzIn1, self.SIGMAAinputMtz)
            shutil.copy2(self.mtzIn2, self.CADinputMtz2)
            shutil.copy2(self.mtzIn3, self.CADinputMtz3)

    def deleteNonFinalMtzs(self):

        # delete all non-final mtz files within run
        # (such as those output by CAD before SCALEIT etc).
        # Not current used at runtime!
        return

        # give option to delete all mtz files within
        # output directory except the final resulting
        # mtz for job - used to save room if necessary
        if self.deleteMtzs.lower() != 'true':
            return
        if self.densMapType == '2FOFC':
            fileEnd = 'sigmaa.mtz'
        else:
            fileEnd = 'SCALEITcombined.mtz'
        for f in os.listdir(self.outputDir):
            if ((f.endswith('.mtz') and not f.endswith(fileEnd)) or
               f.endswith('.tmp')):
                os.remove(self.outputDir+f)

    def cleanUpDir(self):

        # give option to clean up working directory

        # delete non-final mtz files
        self.runLog.writeToLog(str='\nCleaning up working directory...')

        self.deleteNonFinalMtzs()

        # move txt files to subdir
        self.makeOutputDir(dirName='{}txtFiles/'.format(self.outputDir))

        for file in os.listdir(self.outputDir):
            if file.endswith('.txt') and file not in self.filesInDir:
                args = [self.outputDir, file]
                shutil.move('{}{}'.format(*args),
                            '{}txtFiles/{}'.format(*args))

    def findFilesInDir(self):

        # find files initially in working directory

        self.filesInDir = os.listdir(self.outputDir)

    def printStepNumber(self):

        # print a string indicating the current pipeline
        # step number directory to the command line

        try:
            self.stepNumber
        except AttributeError:
            self.stepNumber = 1
        self.runLog.writeToLog(
            str='\n_______\nSTEP {})'.format(self.stepNumber))
        self.stepNumber += 1
