
import imp
import os


class checkDependencies():

    def __init__(self,
                 checkAll=False,
                 optionals=False):

        pythonPackageList = ['sys', 'argparse', 'os',  'time',
                             'shutil', 'numpy', 'imp', 'matplotlib',
                             'math', 'scipy', 'string',
                             'pandas', 'operator', 'warnings',
                             'struct', 'mmap', 'random']

        if optionals:
            pythonPackageList += ['numexpr']

        allFound = True
        if checkAll:
            found = self.checkSeaborn()
            allFound *= found
            found = self.checkCCP4()
            allFound *= found
            for pkg in pythonPackageList:
                found = self.checkPythonPackage(packageName=pkg)
                allFound *= found

            if allFound == 1:
                print 'All dependencies have been successfully located'
            else:
                print 'Warning:\nOne or more required dependencies' +\
                      ' not successfully located\nSee above for details.'

    def checkCCP4(self,
                  printText=True,
                  logFile=''):

        # check whether ccp4 program suite
        # is present and flag if not.

        if printText:
            self.printOrWriteToLog(logFile=logFile,
                                   txt='Checking whether CCP4 ' +
                                       'program suite accessible...')

        if 'ccp4' not in os.environ["PATH"].lower():
            self.printOrWriteToLog(logFile=logFile,
                                   txt='CCP4 program suite not ' +
                                   'found in PATH..\n')
            return False
        else:
            if printText:
                self.printOrWriteToLog(logFile=logFile,
                                       txt='---> success!\n')

        return True

    def checkPythonPackage(self,
                           printText=True, logFile='',
                           packageName='seaborn'):

        # check whether a python package is present and
        # flag if not. logFile is '' (do not write to log
        # file) or a logFile.py class object

        if printText:
            self.printOrWriteToLog(logFile=logFile,
                                   txt='Checking if "{}"'.format(packageName) +
                                       ' python package present...')

        try:
            imp.find_module(packageName)
        except ImportError:
            self.printOrWriteToLog(logFile=logFile,
                                   txt='"{}" package'.format(packageName) +
                                       ' not found..\nTry \n' +
                                       '"pip install {}"'.format(packageName) +
                                       ' to install?\n')
            return False

        if printText:
            self.printOrWriteToLog(logFile=logFile,
                                   txt='---> success!\n')

        return True

    def checkSeaborn(self,
                     printText=True, logFile=''):

        # check whether seaborn is present and flag if not.
        # logFile is '' (do not write to log file file) or
        # a logFile.py class object

        if printText:

            self.printOrWriteToLog(logFile=logFile,
                                   txt='Checking whether seaborn ' +
                                       'plotting library present...')

        try:
            imp.find_module('seaborn')
        except ImportError:
            self.printOrWriteToLog(logFile=logFile,
                                   txt='Plotting library seaborn not\n' +
                                       'found.. Will not create any plots' +
                                       ' for current run.\nUse "pip ' +
                                       'install seaborn" to install for ' +
                                       'future use.\n')
            return False

        if printText:
            self.printOrWriteToLog(logFile=logFile,
                                   txt='---> success!\n')

        return True

    def printOrWriteToLog(self,
                          logFile='', txt=''):

        # print to command line or write to log file

        if logFile == '':
            print txt
        else:
            logFile.writeToLog(str=txt)
