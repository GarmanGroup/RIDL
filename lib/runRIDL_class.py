import os
import time
from errors import error
from processFiles import processFiles
from logFile import logFile


class process():

    # a class to overview the generation of atom-tagged and
    # density maps from pdb & mtz files for a specified damage
    # series. This class contains some methods to write input
    # files, and can call the class processFiles, in which the
    # actual map-generation (and subsequent calculation of
    # per-atom damage metrics) can be performed

    def __init__(self,
                 inputFile='fullInput.txt', makeMaps=True,
                 makeMetrics=False, cleanUpFinalFiles=False,
                 printverboseOutput=False, printOutput=True,
                 makeSummaryFile=False, keepMapDir=True):

        self.inputFile = inputFile
        self.makeMaps = makeMaps
        self.makeMetrics = makeMetrics
        self.makeSummaryFile = makeSummaryFile
        self.cleanUpFinalFiles = cleanUpFinalFiles
        self.printverboseOutput = printverboseOutput
        self.printOutput = printOutput
        self.keepMapDir = keepMapDir

    def setInputFile(self, name):

        # specify input file name as required for the run

        self.inputFile = name
        print('Input file name set as "{}"'.format(self.inputFile))

    def checkInputFileExists(self):

        # check if input file exists and is readable

        if (os.path.isfile(self.inputFile) and
                os.access(self.inputFile, os.R_OK)):
            return True
        else:
            error(text="Either input file is missing or is not readable")
            return False

    def run(self):

        # run to process files with the specified input file

        success = self.checkInputFileExists()
        if not success:
            return success

        self.startLogFile()
        self.titleCaption()
        self.info()

        processFiles(
            inputFile=self.inputFile, makeMaps=self.makeMaps,
            makeMetrics=self.makeMetrics, makeSummaryFile=self.makeSummaryFile,
            cleanFinalFiles=self.cleanUpFinalFiles, logFileObj=self.logFile,
            keepMapDir=self.keepMapDir)

    def printInputFile(self):

        # print the contents of the specified input file

        fileIn = open(self.inputFile, 'r')
        for line in fileIn.readlines():
            if len(line.strip()) != 0:
                print(line.strip())
        fileIn.close()

    def writeTemplateInputFile(self,
                               numHigherDoseDatasets=1):

        # write a template input file to current directory to be completed

        f = open(self.inputFile, 'w')
        v = [','.join(['???']*numHigherDoseDatasets)]*2
        string = 'dir ???\n\nINITIALDATASET\nname1 ???\nmtz1 ???' +\
                 '\nmtzlabels1 ???\npdb1 ???\nRfreeFlag1 ???\n' +\
                 '\nLATERDATASET\nname2 {}\nmtz2 {}\nmtzlabels2 '.format(*v) +\
                 '\nPHASEDATASET\nname3 ???\nmtz3 ???\nphaseLabel ' +\
                 '???\nFcalcLabel ???'
        f.write(string)
        f.close()
        print('Template input file "{}" written'.format(self.inputFile))

    def howToWriteInputFile(self,
                            printStr=True):

        # some information on how write the input file

        infoString = """
*** Information on how to successfully write an input file for job ***

* Note:
* If multiple higher-dose datasets are to be processed successively in a single batch,
* then please input a comma-separated list for each argument within the LATERDATASET
* section below:

FILELOCATION
    dir         : full path to output working directory

INITIALDATASET
    name1       : assign a name to your low dose damage set
    mtz1        : full path to the low dose mtz file
    mtzlabels1  : F & SIGF column labels (look in low dose mtz file), if "FP_X" then type "FP_X"
    pdb1        : full path to the low dose pdb file
    RfreeFlag1  : the Rfree flag label within the low dose mtz file (e.g. "FreeR_flag")

LATERDATASET
    name2       : assign a name to your high dose damage set (e.g. a pdb code "1qid")
    mtz2        : full path to the high dose mtz file
    mtzlabels2  : F & SIGF column labels (look in high dose mtz file), if "FP_X" then type "FP_X"

PHASEDATASET
    name3       : assign a name to your low dose damage set (same as INITIALDATASET above)
    mtz3        : full path to the low dose mtz file (same as INITIALDATASET above)
    phaseLabel  : PHI phase column labels (look in low dose mtz file), if "PHIC_X" then type "PHIC_X"
    FcalcLabel  : calculated structure amplitudes from same file as above, if "FC_X" then type "FC_X"

MAPINFO
    For difference map analysis, do not change these parameters, contact the author for further
    information before varying these
    """
        if printStr:
            print(infoString)

        return infoString

    def fillerLine(self):

        # a simple filler line to print to command line

        ln = '---------------------------------------------------------------'
        self.logFile.writeToLog(str=ln)

    def titleCaption(self):

        # a title line to print to command line

        self.logFile.writeToLog(
          str='\n----------------------------------------------------\n' +
              ' '*24+'RIDL\n' +
              ' '*4+'(Radiation-Induced Density Loss Analysis)\n' +
              '----------------------------------------------------')

    def info(self):

        # print a small ammount of information to the
        # command line (date and email)

        txt = 'date: {}\n'.format(time.strftime("%c")) +\
              'email queries to csbury@me.com\n'
        self.logFile.writeToLog(str=txt)

    def quickParseInputFile(self):

        # parse input file to retrieve directory name for output

        f = open(self.inputFile, 'r')
        for l in f.readlines():
            if len(l) == 0:
                continue
            if 'dir' == l.split()[0]:
                self.outputDir = l.split()[1]
                if self.outputDir[-1] not in ('/'):
                    self.outputDir += '/'
                return

    def startLogFile(self):

        # create log file for the current job

        self.quickParseInputFile()
        self.checkOutputDirExists(printToScreen=self.printOutput)

        if self.outputDir[-1] == '/':
            logDir = self.outputDir[:-1]+'/RIDL-log/'
        else:
            logDir = self.outputDir+'/RIDL-log/'

        if not os.path.exists(logDir):
            os.makedirs(logDir)

        name = '{}RIDLjob.log'.format(logDir)

        logName = lambda x: name.replace('.log', '_{}{}'.format(x, '.log'))
        i = 1
        while os.path.isfile(logName(i)):
            i += 1
        uniqLogName = logName(i)

        log = logFile(fileName=uniqLogName, fileDir=self.outputDir,
                      printToScreenMajor=self.printOutput,
                      printToScreenMinor=self.printverboseOutput)
        self.logFile = log

    def checkOutputDirExists(self,
                             printToScreen=True):

        # check whether output directory exists and make if not

        if not self.outputDir.endswith('/'):
            if printToScreen:
                print('Working directory specified in input ' +
                      'file must end in "/" - appending.')
            self.outputDir += '/'

        if not os.path.exists(self.outputDir):
            if printToScreen:
                print('Output directory "{}" '.format(self.outputDir) +
                      'not found, making directory')
            os.makedirs(self.outputDir)
