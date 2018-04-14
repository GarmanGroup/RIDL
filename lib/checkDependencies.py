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
                             'struct', 'mmap', 'random', 'seaborn']

        if optionals:
            pythonPackageList += ['numexpr']

        allFound = True
        if checkAll:
            found = self.checkCCP4()
            allFound *= found
            self.printOrWriteToLog(txt='Checking required python packages...')
            pyPkgsFound = True
            for pkg in pythonPackageList:
                found = self.checkPyPackage(packageName=pkg)
                allFound *= found
                pyPkgsFound *= found
            if pyPkgsFound == 1:
                self.printOrWriteToLog(txt='---> success!\n')

            if allFound == 1:
                self.printOrWriteToLog(
                  txt='All dependencies have been successfully located')
            else:
                self.printOrWriteToLog(
                  txt='Warning:\nOne or more required dependencies' +
                      ' not successfully located.\nSee above for details.')

    def checkCCP4(self,
                  logFile=''):

        # check whether ccp4 program suite
        # is present and flag if not.

        self.printOrWriteToLog(logFile=logFile,
                               txt='Checking whether CCP4 ' +
                                   'program suite accessible...')

        if 'ccp4' not in os.environ["PATH"].lower():
            self.printOrWriteToLog(logFile=logFile,
                                   txt='CCP4 program suite not ' +
                                   'found in PATH..\n')
            return False
        else:
            self.printOrWriteToLog(logFile=logFile,
                                   txt='---> success!\n')

        return True

    def checkPyPackage(self,
                       logFile='', packageName='seaborn'):

        # check whether a python package is present and
        # flag if not. logFile is '' (do not write to log
        # file) or a logFile.py class object

        try:
            imp.find_module(packageName)
        except ImportError:
            self.printOrWriteToLog(logFile=logFile,
                                   txt='"{}" package'.format(packageName) +
                                       ' not found..\nTry ' +
                                       '"pip install {}"'.format(packageName) +
                                       ' to install?\n')

            if packageName == 'seaborn':
                self.printOrWriteToLog(
                  logFile=logFile,
                  txt='Warning:\nPlotting library seaborn not found.. Will ' +
                      'not create any plots for current run.\nNote that the ' +
                      'HTML-format summary file will not be available if ' +
                      'seaborn is not present.\nNote that the metrics can ' +
                      'still be generated using "python runRIDL.py -i ' +
                      'input.txt -pc"\n')

            return False

        return True

    def checkSeaborn(self,
                     logFile=''):

        # check whether seaborn is present and flag if not.
        # logFile is '' (do not write to log file file) or
        # a logFile.py class object

        self.printOrWriteToLog(logFile=logFile,
                               txt='\nChecking whether seaborn ' +
                                   'plotting library present...')

        try:
            imp.find_module('seaborna')
        except ImportError:
            self.printOrWriteToLog(logFile=logFile,
                                   txt='Warning:\nPlotting library seaborn ' +
                                       'not found.. Will not create any ' +
                                       'plots for current run.\nUse "pip ' +
                                       'install seaborn" to install for ' +
                                       'future use.\nNote that the HTML-' +
                                       'format summary file will not be ' +
                                       'available if seaborn is not present\n')
            return False

        self.printOrWriteToLog(logFile=logFile, txt='---> success!\n')

        return True

    def printOrWriteToLog(self,
                          logFile='', txt=''):

        # print to command line or write to log file

        if logFile == '':
            print(txt)
        else:
            logFile.writeToLog(str=txt)
