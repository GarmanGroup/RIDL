from ccp4Job import ccp4Job, checkInputsExist
from mapTools import mapTools


class FFTjob():

    def __init__(self,
                 mapType='DIFF-NO-SIGMA',  mapTag='', FOMweight='false',
                 pdbFile='', mtzFile='', outputDir='./', axes=[1, 2, 3],
                 gridSamps=[0, 0, 0], labels1=['', '', '', ''],
                 labels2=['', '', '', ''], F1Scale=1.0, F2Scale=1.0,
                 highResCutoff='', lowResCutoff='', runLog='',
                 outputMapFile='', useSigLabs=True):

        self.mapType = mapType
        self.FOMweight = FOMweight
        self.inputMtzFile = mtzFile
        self.F1Scale = F1Scale
        self.F2Scale = F2Scale
        self.useSigLabs = useSigLabs
        self.labels1 = labels1
        self.labels2 = labels2
        self.mapTag = mapTag
        self.pdbFile = pdbFile
        self.outputMapFile = outputMapFile
        self.outputDir = outputDir
        XYZ = ['X', 'Y', 'Z']
        self.fastAxis = XYZ[int(axes[0])-1]
        self.medAxis = XYZ[int(axes[1])-1]
        self.slowAxis = XYZ[int(axes[2])-1]
        self.gridSamp1 = gridSamps[0]
        self.gridSamp2 = gridSamps[1]
        self.gridSamp3 = gridSamps[2]
        self.runLog = runLog
        self.highResCutoff = highResCutoff
        self.lowResCutoff = lowResCutoff

    def run(self):

        inputFiles = [self.inputMtzFile]
        if not checkInputsExist(inputFiles, self.runLog):
            return False

        # can change unique identifier for map type in file name
        if self.mapTag == '':
            self.mapTag = self.mapType

        if self.outputMapFile == '':
            self.outputMapFile = '{}-{}_FFT.map'.format(
                self.pdbFile.split('.pdb')[0], self.mapTag)

        self.runFFT()

        if self.jobSuccess:
            self.provideFeedback()
            self.runLog.writeToLog('Output files:')
            self.runLog.writeToLog('{}'.format(self.outputMapFile))
            return True

        else:
            self.runLog.writeToLog('Job did not run successfully, see job ' +
                                   'log file "{}"'.format(self.outputLogfile))
            return False

    def runFFT(self):

        # run FFT job using the external ccp4Job class

        F1, F2 = self.labels1[0], self.labels2[0]
        SIG1, SIG2 = self.labels1[1], self.labels2[1]
        FOM1, FOM2 = self.labels1[2], self.labels2[2]
        PHI1, PHI2 = self.labels1[3], self.labels2[3]

        self.printPurpose()
        cmd1 = 'fft ' +\
               'HKLIN {} '.format(self.inputMtzFile) +\
               'MAPOUT {} '.format(self.outputMapFile) +\
               'SYMINFO syminfo.lib '

        # if FOM is specified the weighting is applied to the map
        if self.FOMweight.lower() == 'recalculate' or 'preset,' in self.FOMweight.lower():
            FOMstring = 'W={}'.format(FOM2)
            FOMstring2 = 'W2={}'.format(FOM2)
        else:
            FOMstring, FOMstring2 = '', ''

        # DIFF is the default map option. It calculates a F-F
        # difference map between a low and high dataset
        if self.mapType == 'DIFF':
            if self.useSigLabs:
                labinStr = 'LABIN F1={} SIG1={} F2={} SIG2={} '.format(
                    F1, SIG1, F2, SIG2)
            else:
                labinStr = 'LABIN F1={} F2={} '.format(F1, F2)
            phiStr = 'PHI={} {} PHI2={} {}\n'.format(
                PHI2, FOMstring, PHI2, FOMstring2)

            scaleStr = 'SCALE F1 {} 0.0 F2 {} 0.0\n'.format(
                self.F1Scale, self.F2Scale)

        # SIMPLE sets the low dataset weight to 0. This is essentially
        # map generated only from the high dataset.
        if self.mapType == 'SIMPLE':
            if self.useSigLabs:
                labinStr = 'LABIN F1={} SIG1={} F2={} SIG2={} '.format(
                    F1, SIG1, F2, SIG2)
            else:
                labinStr = 'LABIN F1={} F2={} '.format(F1, F2)
            phiStr = 'PHI={} {} PHI2={} {}\n'.format(
                PHI2, FOMstring, PHI2, FOMstring2)
            scaleStr = 'SCALE F1 {} 0.0 F2 0.0 0.0\n'.format(
                self.F1Scale, self.F2Scale)

        # this option is very similar to SIMPLE above, however the
        # low dataset is neglected entirely.
        # It is designed to superseed SIMPLE now and doesn't require
        # a valid low dataset F column to be present
        if self.mapType == 'HIGHONLY':
            if self.useSigLabs:
                labinStr = 'LABIN F1={} SIG1={} '.format(F1, SIG1)
            else:
                labinStr = 'LABIN F1={} '.format(F1)
            phiStr = 'PHI={} {}\n'.format(PHI2, FOMstring)
            scaleStr = 'SCALE F1 {} 0.0 '.format(self.F1Scale)

        # the option below is distinct from DIFF above, since it allows
        # a difference map to be generated after SIGF labels have been
        # used here to scale datasets, however does NOT use SIGFs in the
        # map calculation
        if self.mapType == 'DIFF-NO-SIGMA':
            labinStr = 'LABIN F1={} F2={} '.format(F1, F2)
            phiStr = 'PHI={} {} PHI2={} {}\n'.format(
                PHI1, FOMstring, PHI2, FOMstring2)
            scaleStr = 'SCALE F1 {} 0.0 F2 {} 0.0\n'.format(
                self.F1Scale, self.F2Scale)

        # choose sigmaa-derived FWT and PHIC for 2FOFC map setting
        if self.mapType == '2FOFC':
            labinStr = 'LABIN F1={}'.format(F1)
            phiStr = 'PHI={}'.format(PHI1)
            scaleStr = 'SCALE F1 1.0 0.0'

        # for case of FC map generation
        if self.mapType == 'FC':
            labinStr = 'LABIN F1={}'.format(F1)
            phiStr = 'PHI={}'.format(PHI1)
            scaleStr = 'SCALE F1 {} 0.0'.format(self.F1Scale)

        # include a resolution cutoff if specified (not used for Fcalc maps)
        resStr = ''
        if self.mapType != 'FC':
            if self.highResCutoff != '':
                resStr = 'RESOLUTION {}'.format(self.highResCutoff)
                if self.lowResCutoff != '':
                    resStr += ' {}'.format(self.lowResCutoff)
                resStr += '\n'

        cmd2 = 'AXIS {} {} {}\n'.format(
            self.fastAxis, self.medAxis, self.slowAxis) +\
            'title FTT DENSITY MAP RUN\n' +\
            'grid {} {} {}\n'.format(
            self.gridSamp1, self.gridSamp2, self.gridSamp3) +\
            'xyzlim 0 1 0 1 0 1\n' +\
            '{}'.format(resStr) +\
            '{} {}\n{}\nEND'.format(labinStr, phiStr, scaleStr)

        self.outputLogfile = 'FFTlogfile.txt'

        # run FFT job
        job = ccp4Job(jobName='FFT', commandInput1=cmd1,
                      commandInput2=cmd2,
                      outputDir=self.outputDir, outputLog=self.outputLogfile,
                      outputFile=self.outputMapFile)

        self.jobSuccess = job.checkJobSuccess(self.runLog)

    def provideFeedback(self,
                        includeDir=False):

        # provide some feedback

        if not includeDir:
            fileIn = self.inputMtzFile.split('/')[-1]
            fileOut = self.outputMapFile.split('/')[-1]
        else:
            fileIn = self.inputMtzFile
            fileOut = self.outputMapFile

        self.runLog.writeToLog('FFT Summary:\n' +
                               'Input mtz file: {}\n'.format(fileIn) +
                               'Output map file: {}'.format(fileOut))

        Map = mapTools(mapName=self.outputMapFile, logFile=self.runLog)
        Map.printMapInfo()

    def printPurpose(self,
                     include=True):

        # provide a summary of what this does
        # (within RIDL) to the command line

        if self.mapTag == 'DIFF':
            ln = 'Generating Fourier difference map over crystal unit cell\n'
        if self.mapTag == 'HIGHONLY':
            ln = 'Generating simple map over crystal unit cell using ' +\
                 'high dataset info only\n'
        if self.mapTag == 'SIMPLE':
            ln = 'Generating Fobs density map over crystal unit cell\n'
        if self.mapTag == '2FOFC':
            ln = 'Generating 2Fo-Fc density map over crystal unit cell\n'
        if self.mapTag == 'FC':
            ln = 'Generating Fcalc density map over crystal unit cell\n'
        ln += 'with same grid dimensions as SFALL-output atom-tagged .map file'
        self.runLog.writeToLog(ln)
