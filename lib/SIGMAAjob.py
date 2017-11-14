from ccp4Job import ccp4Job, checkInputsExist, fillerLine
from errors import error


class SIGMAAjob():
    # run SIGMAA job to combine generate FOM weights
    # if not present in input mtz file
    def __init__(self,
                 inputMtz='', MtzLabelNameIn='', MtzLabelNameOut='',
                 RfreeFlag='', inputPDB='', outputDir='./', runLog='',
                 noScale=True):

        self.inputMtz = inputMtz
        self.LabelName = MtzLabelNameIn
        self.LabelRename = MtzLabelNameOut
        self.RfreeFlag = RfreeFlag
        self.tmpMtz = outputDir+'/'+(
            inputMtz.split('/')[-1]).split('.mtz')[0]+'.tmp'
        self.outputMtz = self.tmpMtz.split('.tmp')[0]+'_sigmaa.mtz'
        self.inputPDB = inputPDB
        self.outputDir = outputDir
        self.runLog = runLog
        self.noScale = noScale

        self.runLog.writeToLog('Running SIGMAA job')

    def run(self):
        inputFiles = [self.inputMtz, self.inputPDB]
        if not checkInputsExist(inputFiles, self.runLog):
            return False

        success = self.getSpaceGroup()
        if not success:
            return False

        self.runSIGMAA()
        if self.jobSuccess:
            self.provideFeedback()
            self.runLog.writeToLog('Output files:')
            self.runLog.writeToLog('{}'.format(self.outputMtz))
            return True
        else:
            error(text='Job did not run successfully, see job log file ' +
                       '"{}"'.format(self.outputLogfile),
                  log=self.runLog, type='error')
            return False

    def runSIGMAA(self):
        title = 'run of sigma-a'
        fillerLine()

        cmd1 = 'sfall ' +\
               'HKLIN {} '.format(self.inputMtz) +\
               'HKLOUT {} '.format(self.tmpMtz) +\
               'XYZIN {}'.format(self.inputPDB)

        cmd2 = 'title {}\n'.format(title) +\
               'LABIN  FP=F{} SIGFP=SIGF{} FREE={}\n'.format(
                    self.LabelName, self.LabelName, self.RfreeFlag) +\
               'labout -\nFC=FC PHIC=PHIC\nMODE SFCALC -\n' +\
               'XYZIN -\nHKLIN\nsymmetry {}\n'.format(self.spaceGroup)

        if self.noScale is True:
            self.commandInput2 += 'NOSCALE\nend'
        else:
            self.commandInput2 += 'end'

        self.outputLogfile = 'SIGMAAlogfile1.txt'

        # run SFALL job within SIGMAA routine
        job = ccp4Job(jobName='SFALL', commandInput1=cmd1, commandInput2=cmd2,
                      outputDir=self.outputDir, outputLog=self.outputLogfile,
                      outputFile=self.tmpMtz)

        self.jobSuccess = job.checkJobSuccess(self.runLog)
        if not self.jobSuccess:
            return

        cmd1 = 'sigmaa ' +\
               'HKLIN {} '.format(self.tmpMtz) +\
               'HKLOUT {} '.format(self.outputMtz)

        cmd2 = 'title {}\n'.format(title) +\
               'LABIN  FP=F{} SIGFP=SIGF{} FC=FC PHIC=PHIC\n'.format(
                *([self.LabelName]*2)) +\
               'labout -\nDELFWT=DELFWT{} FWT=FWT{} WCMB=FOM{}\n'.format(
                *([self.LabelRename]*3)) +\
               'ranges 20 -\n1000\nsymmetry "{}"\n'.format(self.spaceGroup) +\
               'partial -\ndamp 1.0\nend'

        self.outputLogfile = 'SIGMAAlogfile2.txt'

        # run SIGMAA next
        job = ccp4Job(jobName='SIGMAA', commandInput1=cmd1, commandInput2=cmd2,
                      outputDir=self.outputDir, outputLog=self.outputLogfile,
                      outputFile=self.outputMtz)

        self.jobSuccess = job.checkJobSuccess()

    def provideFeedback(self, includeDir=False):

        # provide some feedback

        if not includeDir:
            fileIn = self.inputMtz.split('/')[-1]
            fileOut = self.outputMtz.split('/')[-1]
        else:
            fileIn = self.inputMtz
            fileOut = self.outputMtz

        self.runLog.writeToLog(
            'SIGMAA Summary:\nInput mtz file: {}\nOutput mtz file: {}'.format(
                fileIn, fileOut))

    def getSpaceGroup(self):
        pdbin = open(self.inputPDB, 'r')
        for line in pdbin.readlines():
            if line.split()[0] == 'CRYST1':
                self.spaceGroup = line[55:66].replace(' ', '')
                self.runLog.writeToLog(
                    'Retrieving space group from file: "{}"'.format(
                        self.inputPDB))
                self.runLog.writeToLog(
                    'Space group determined to be {}'.format(self.spaceGroup))
        try:
            self.spaceGroup
        except AttributeError:
            self.runLog.writeToLog(
                'Cannot find space group from file: {}'.format(self.inputPDB))
            return False
        return True
