import sys
sys.path.insert(0, './lib')
from calculateMetrics import calculateMetrics
import os


class run():

    # run the metric calculation step of the RIDL pipeline
    #  to calculate per-atom density metrics for a series
    # of increasing doses

    def __init__(self,
                 calculate=True, inputFileLoc='',
                 logFile='', skipToSumFiles=False,
                 writeSumFiles=False, inclFCmets=False):

        self.inputFileName = inputFileLoc + 'metricCalc_inputfile.txt'
        self.logFile = logFile

        if calculate:
            self.runCalculateMetrics(sumFiles=writeSumFiles,
                                     inclFCmets=inclFCmets)

        if skipToSumFiles:
            self.runCalculateMetrics(mapProcess=False, postProcess=False,
                                     retrieve=True, sumFiles=True,
                                     inclFCmets=inclFCmets)

    def runCalculateMetrics(self,
                            mapProcess=True, postProcess=True,
                            retrieve=False, sumFiles=False,
                            inclFCmets=False):

        # run the metric calculation scripts for
        # the currently defined input file
        self.logFile.writeToLog(str='\n\n**** METRIC CALCULATIONS ****\n')

        exists = self.checkInputFileExists()
        if not exists:
            return

        eT = calculateMetrics(logFile=self.logFile, sumFiles=sumFiles,
                              inclFCmets=inclFCmets)

        eT.runPipeline(map_process=mapProcess, post_process=postProcess,
                       retrieve=retrieve, inputFileName=self.inputFileName)
        self.et = eT

    def defineDoseList(self, doses, names, version):

        # do not include first dataset dose if difference maps chosen

        if version == 'DIFF':
            if len(doses.split(',')[1:]) > 1:
                dosesOut = ','.join(doses.split(',')[1:])
                namesOut = ','.join(names.split(',')[1:])
            else:
                dosesOut = doses.split(',')[1]
                namesOut = names.split(',')[1]
            return dosesOut, namesOut
        else:
            return doses, names

    def writeBlankInputFile(self):

        # Need to create input file for a generic damage
        # series to be completed by the user

        args = ['inDir <input file location>',
                'outDir <output file location>',
                'damageset_name <consistent name of series, ' +
                'e.g. for TRAP1.pdb, TRAP2.pdb, this is TRAP>',
                'damageset_num <dataset numbers, e.g. 1,2 for ' +
                'TRAP. Can be letters corresponding to pdb series>',
                'initialPDB <initial dataset pdb file e.g. TRAP1.pdb>',
                'doses <list of doses, length must match length ' +
                'of damageset_num above>',
                'PKLMULTIFILE <if already processed in RIDL, can ' +
                'specify single output .pkl for series>']
        self.writeInputFile(*args)

    def writeInputFile(self,
                       inDir='', outDir='', damSetName='',
                       laterDatasets='', initialPDB='', doses='',
                       PKLMULTIFILE='', outputGraphs='yes'):

        # write a generic input file for a damage series here

        inputString = 'inDir {}\n'.format(inDir) +\
                      'outDir {}\n'.format(outDir) +\
                      'damageset_name {}\n'.format(damSetName) +\
                      'laterDatasets {}\n'.format(laterDatasets) +\
                      'initialPDB {}\n'.format(initialPDB) +\
                      'doses {}\n'.format(doses)

        if PKLMULTIFILE != '':
            inputString += 'PKLMULTIFILE {}'.format(PKLMULTIFILE)

        if outputGraphs == 'yes':
            inputString += 'plot'
        elif outputGraphs == 'slim':
            inputString += 'slim-plot'

        inputFile = open(self.inputFileName, 'w')
        inputFile.write(inputString)
        inputFile.close()

    def checkInputFileExists(self):

        # check whether input file present for the metric
        # calculation step of the RIDL pipeline

        if (os.path.isfile(self.inputFileName) and
                os.access(self.inputFileName, os.R_OK)):
            return True
        else:
            print "Either input file is missing or is not readable"
            print 'Check whether "{}" is a '.format(self.inputFileName) +\
                  'suitable input file.'
            return False
