from ccp4Job import ccp4Job, checkInputsExist


class CADjob():

    # run CAD job to combine F and SIGF columns
    # from two merged mtz files

    def __init__(self,
                 inputMtz1='', inputMtz2='', inputMtz3='',
                 Mtz1FPlabel='FP', Mtz1SIGFPlabel='SIGFP',
                 Mtz2FPlabel='FP', Mtz2SIGFPlabel='SIGFP',
                 Mtz1LabelName='', Mtz2LabelName='',
                 Mtz3phaseLabel='', Mtz3FcalcLabel='',
                 Mtz1LabelRename='', Mtz2LabelRename='',
                 Mtz3LabelRename='', outputMtz='',
                 outputDir='./', runLog='', FOMWeight='False',
                 ignoreSIGFs=False, ignoreDset1=False):

        # mtz containing initial dataset Fs
        self.inputMtz1 = inputMtz1

        # mtz containing later dataset Fs
        self.inputMtz2 = inputMtz2

        # mtz containing phases
        self.inputMtz3 = inputMtz3

        # labels in input mtz files
        self.Mtz1FPlabel = Mtz1FPlabel
        self.Mtz1SIGFPlabel = Mtz1SIGFPlabel
        self.Mtz2FPlabel = Mtz2FPlabel
        self.Mtz2SIGFPlabel = Mtz2SIGFPlabel
        self.Mtz3phaseLabel = Mtz3phaseLabel
        self.Mtz3FcalcLabel = Mtz3FcalcLabel
        self.Mtz1LabelRename = Mtz1LabelRename
        self.Mtz2LabelRename = Mtz2LabelRename
        self.Mtz3LabelRename = Mtz3LabelRename

        # provide option to ignore initial dataset in merge
        # to allow case where only later dataset is specified
        self.ignoreDset1 = ignoreDset1

        # the output mtz name
        self.outputMtz = outputMtz

        # the output location where files should be written
        self.outputDir = outputDir

        # the run log file (this is the same log as the overall RIDL one)
        self.runLog = runLog

        # provide option to not include SIGF labels
        self.ignoreSIGFs = ignoreSIGFs

        # the figure of merit column naming if specified
        self.FOMWeight = FOMWeight
        self.FOMtag = {'in': '', 'out': '', 'type': ''}

    def run(self):

        # a new labelling convention for output files

        self.renameLabels = [self.Mtz1LabelRename, self.Mtz2LabelRename,
                             self.Mtz3LabelRename]

        inputFiles = [self.inputMtz2, self.inputMtz3]
        if not self.ignoreDset1:
            inputFiles += [self.inputMtz1]

        if self.FOMWeight != 'False':
            self.FOMtag['out'] = '- \nE3 = FOM_{}'.format(self.renameLabels[0])
            self.FOMtag['type'] = '- \nE3 = W'

        if self.FOMWeight == 'recalculate':
            self.FOMtag['in'] = '- \nE3 = FOM{}'.format(self.labels[0])

        elif 'preset' in self.FOMWeight:
            lbl = self.FOMWeight.split(',')[-1]
            self.FOMtag['in'] = '- \nE3 = FOM{}'.format(lbl)

        if not checkInputsExist(inputFiles, self.runLog):
            return False

        self.runCAD()
        if self.jobSuccess:
            self.provideFeedback()
            return True

        else:
            self.runLog.writeToLog('Job did not run successfully, see job' +
                                   'log file "{}"'.format(self.outputLogfile))
            return False

    def runCAD(self):

        self.printPurpose()

        if not self.ignoreDset1:

            cmd1 = 'cad ' +\
                   'HKLIN1 {} '.format(self.inputMtz1) +\
                   'HKLIN2 {} '.format(self.inputMtz2) +\
                   'HKLIN3 {} '.format(self.inputMtz3) +\
                   'HKLOUT {}'.format(self.outputMtz)

            if not self.ignoreSIGFs:
                cmd2 = 'title CAD JOB\n' +\
                       'monitor BRIEF\n' +\
                       'labin file 1 - \n' +\
                       'E1 = {} - \n'.format(self.Mtz1FPlabel) +\
                       'E2 = {} '.format(self.Mtz1SIGFPlabel) +\
                       '{} \n'.format(self.FOMtag['in']) +\
                       'labout file 1 - \n' +\
                       'E1 = FP_{} - \n'.format(self.renameLabels[0]) +\
                       'E2 = SIGFP_{} '.format(self.renameLabels[0]) +\
                       '{} \n'.format(self.FOMtag['out']) +\
                       'ctypin file 1 - \n' +\
                       'E1 = F - \n' +\
                       'E2 = Q ' +\
                       '{} \n'.format(self.FOMtag['type']) +\
                       'labin file 2 - \n' +\
                       'E1 = {} - \n'.format(self.Mtz2FPlabel) +\
                       'E2 = {} \n'.format(self.Mtz2SIGFPlabel) +\
                       'labout file 2 - \n' +\
                       'E1 = FP_{} - \n'.format(self.renameLabels[1]) +\
                       'E2 = SIGFP_{} \n'.format(self.renameLabels[1]) +\
                       'ctypin file 2 - \n' +\
                       'E1 = F - \n' +\
                       'E2 = Q \n' +\
                       'labin file 3 - \n' +\
                       'E1 = {} -\n'.format(self.Mtz3phaseLabel) +\
                       'E2 = {} \n'.format(self.Mtz3FcalcLabel) +\
                       'labout file 3 - \n' +\
                       'E1 = PHIC_{} - \n'.format(self.renameLabels[2]) +\
                       'E2 = FC_{} \n'.format(self.renameLabels[2]) +\
                       'ctypin file 3 - \n' +\
                       'E1 = P -\n' +\
                       'E2 = F \n'
            else:
                cmd2 = 'title CAD JOB\n' +\
                       'monitor BRIEF\n' +\
                       'labin file 1 - \n' +\
                       'E1 = {} \n'.format(self.Mtz1FPlabel) +\
                       '{} \n'.format(self.FOMtag['in']) +\
                       'labout file 1 - \n' +\
                       'E1 = FP_{} \n'.format(self.renameLabels[0]) +\
                       '{} \n'.format(self.FOMtag['out']) +\
                       'ctypin file 1 - \n' +\
                       'E1 = F \n' +\
                       '{} \n'.format(self.FOMtag['type']) +\
                       'labin file 2 - \n' +\
                       'E1 = {} \n'.format(self.Mtz2FPlabel) +\
                       'labout file 2 - \n' +\
                       'E1 = FP_{} \n'.format(self.renameLabels[1]) +\
                       'ctypin file 2 - \n' +\
                       'E1 = F \n' +\
                       'labin file 3 - \n' +\
                       'E1 = {} -\n'.format(self.Mtz3phaseLabel) +\
                       'E2 = {} \n'.format(self.Mtz3FcalcLabel) +\
                       'labout file 3 - \n' +\
                       'E1 = PHIC_{} - \n'.format(self.renameLabels[2]) +\
                       'E2 = FC_{} \n'.format(self.renameLabels[2]) +\
                       'ctypin file 3 - \n' +\
                       'E1 = P -\n' +\
                       'E2 = F \n'

        else:
            cmd1 = 'cad ' +\
                 'HKLIN1 {} '.format(self.inputMtz2) +\
                 'HKLIN2 {} '.format(self.inputMtz3) +\
                 'HKLOUT {}'.format(self.outputMtz)

            if not self.ignoreSIGFs:
                cmd2 = 'title CAD JOB\n' +\
                       'monitor BRIEF\n' +\
                       'labin file 1 - \n' +\
                       'E1 = {} - \n'.format(self.Mtz2FPlabel) +\
                       'E2 = {} \n'.format(self.Mtz2SIGFPlabel) +\
                       'labout file 1 - \n' +\
                       'E1 = FP_{} - \n'.format(self.renameLabels[1]) +\
                       'E2 = SIGFP_{} \n'.format(self.renameLabels[1]) +\
                       'ctypin file 1 - \n' +\
                       'E1 = F - \n' +\
                       'E2 = Q \n' +\
                       'labin file 2 - \n' +\
                       'E1 = {} -\n'.format(self.Mtz3phaseLabel) +\
                       'E2 = {} \n'.format(self.Mtz3FcalcLabel) +\
                       'labout file 2 - \n' +\
                       'E1 = PHIC_{} - \n'.format(self.renameLabels[2]) +\
                       'E2 = FC_{} \n'.format(self.renameLabels[2]) +\
                       'ctypin file 2 - \n' +\
                       'E1 = P -\n' +\
                       'E2 = F \n'
            else:
                cmd2 = 'title CAD JOB\n' +\
                       'monitor BRIEF\n' +\
                       'labin file 1 - \n' +\
                       'E1 = {} \n'.format(self.Mtz2FPlabel) +\
                       'labout file 1 - \n' +\
                       'E1 = FP_{} \n'.format(self.renameLabels[1]) +\
                       'ctypin file 1 - \n' +\
                       'E1 = F \n' +\
                       'labin file 2 - \n' +\
                       'E1 = {} -\n'.format(self.Mtz3phaseLabel) +\
                       'E2 = {} \n'.format(self.Mtz3FcalcLabel) +\
                       'labout file 2 - \n' +\
                       'E1 = PHIC_{} - \n'.format(self.renameLabels[2]) +\
                       'E2 = FC_{} \n'.format(self.renameLabels[2]) +\
                       'ctypin file 2 - \n' +\
                       'E1 = P -\n' +\
                       'E2 = F \n'

        self.outputLogfile = 'CADlogfile.txt'

        # run CAD job
        job = ccp4Job(jobName='CAD', commandInput1=cmd1,
                      commandInput2=cmd2, outputDir=self.outputDir,
                      outputLog=self.outputLogfile,
                      outputFile=self.outputMtz)

        self.jobSuccess = job.checkJobSuccess(self.runLog)

    def provideFeedback(self,
                        includeDir=False):

        # provide some feedback

        if not includeDir:
            fileIn1 = self.inputMtz1.split('/')[-1]
            fileIn2 = self.inputMtz2.split('/')[-1]
            fileIn3 = self.inputMtz3.split('/')[-1]
            fileOut = self.outputMtz.split('/')[-1]
        else:
            fileIn1 = self.inputMtz1
            fileIn2 = self.inputMtz2
            fileIn3 = self.inputMtz3
            fileOut = self.outputMtz

        txt = 'CAD Summary:\n' +\
              'Input mtz files: {}\n{}{}\n{}{}\n'.format(
                fileIn1, ' '*17, fileIn2, ' '*17, fileIn3) +\
              'Output mtz file: {}'.format(fileOut)
        self.runLog.writeToLog(txt)

    def printPurpose(self,
                     include=True):

        # provide a summary of what this does
        # (within RIDL) to the command line

        txt = 'Combining relevant columns from multiple ' +\
              'input mtz files into 1 single mtz file'
        self.runLog.writeToLog(txt)
