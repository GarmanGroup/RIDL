import os
import shutil
from errors import error


class ccp4Job():

    def __init__(self,
                 jobName='untitled', commandInput1='',
                 commandInput2='', outputDir='./',
                 outputLog='', outputFile=''):

        self.jobName = jobName
        self.commandInput1 = commandInput1
        self.commandInput2 = commandInput2
        self.outputDir = outputDir
        self.outputLogfile = outputLog
        self.outputFile = outputFile

        # automatically run ccp4 program
        self.runCCP4program()

    def runCCP4program(self):
        # generic method to run a ccp4 program on command line

        # write commandInput2 to txt file
        textinput = open('{}inputfile.txt'.format(self.jobName), 'w')
        textinput.write(self.commandInput2)
        textinput.close()

        # run ccp4 program job
        os.system('{} < {}inputfile.txt > {}'.format(self.commandInput1,
                                                     self.jobName,
                                                     self.outputLogfile))

        # move ccp4 job input and log files to working sub directory
        shutil.move('{}inputfile.txt'.format(self.jobName),
                    '{}{}inputfile.txt'.format(self.outputDir,
                                               self.jobName))

        shutil.move(self.outputLogfile,
                    '{}{}'.format(self.outputDir,
                                  self.outputLogfile))

    def checkJobSuccess(self, runLog):

        # job success checked, based on:
        # (a) the program log file has completed
        # (b) whether output files exist

        fIn = open(self.outputDir+self.outputLogfile, 'r')
        logFinished = False
        for l in fIn.readlines():
            if 'normal termination' in l.lower():
                logFinished = True
                break
        fIn.close()

        if not logFinished:
            error(
                text='{} did not proceed to completion'.format(self.jobName),
                log=runLog)
            return False

        if not os.path.isfile(self.outputFile):
            error(
                text='{} did not proceed to completion'.format(self.jobName),
                log=runLog)
            return False
        else:
            return True


def checkInputsExist(inputFiles, runLog):
    # check if input files exist for a ccp4 job
    for fileName in inputFiles:
        if not os.path.isfile(fileName):
            runLog.writeToLog(
                str='Failed to find required input file "{}"'.format(fileName))
            return False
    else:
        return True


def fillerLine(longLine=True, linebreak=True):
    if longLine:
        ln = '--------------------------'*2
    else:
        ln = '--------------------------'
    if linebreak:
        ln = '\n'+ln
    print ln
