
import imp
import os

class checkDependencies():

	def __init__(self,
				 checkAll = False):

		pythonPackageList = ['sys',
							 'argparse',
							 'os',
							 'time',
							 'shutil',
							 'numpy',
							 'imp',
							 'matplotlib',
							 'math',
							 'scipy',
							 'numexpr',
							 'string',
							 'pandas',
							 'operator',
							 'warnings',
							 'struct',
							 'mmap',
							 'random']

		# pythonPackageList += ['sklearn']

		allFound = True
		if checkAll:
			found = self.checkSeaborn()
			allFound *= found
			found = self.checkCCP4()
			allFound *= found
			for pkg in pythonPackageList:
				found = self.checkPythonPackage(packageName = pkg)
				allFound *= found

			if allFound == 1:
				print 'All required dependencies have been successfully located'
			else:
				print 'Warning:\nOne or more required dependencies'+\
					  ' not successfully located\nSee above for details.'

	def checkCCP4(self,
				  printText = False,
				  logFile   = ''):

		# check whether ccp4 program suite 
		# is present and flag if not.

		if printText is True:

			ln = 'Checking whether CCP4 program suite accessible...'
			self.printOrWriteToLog(logFile = logFile,
							  	   txt     = ln)

		t = 'testccp4works.txt'
		os.system('aimless &> {}'.format(t))
		l = open(t,'r').readlines()[0]
		if l == 'sh: aimless: command not found\n':
			err = 'CCP4 program suite not found..\n'
			self.printOrWriteToLog(logFile = logFile,
						  	       txt     = err)
			return False

		if printText is True:	
			ln = '---> success!\n'
			self.printOrWriteToLog(logFile = logFile,
							  	   txt     = ln)

		return True

	def checkPythonPackage(self,
					       printText   = False,
					       logFile     = '',
					       packageName = 'seaborn'):

		# check whether a python package is present and 
		# flag if not. logFile is '' (do not write to log 
		# file) or a logFile.py class object

		if printText is True:

			ln = 'Checking whether "{}" python package present...'.format(packageName)
			self.printOrWriteToLog(logFile = logFile,
							  	   txt     = ln)

		try:
			imp.find_module(packageName)
		except ImportError:
			err = 'Python package {} not found..\n'.format(packageName)+\
			 	  'Recommendation:\n'+\
			      'Try "pip install {}" to install?\n'.format(packageName)
			self.printOrWriteToLog(logFile = logFile,
						  	  	   txt     = err)
			return False

		if printText is True:	
			ln = '---> success!\n'
			self.printOrWriteToLog(logFile = logFile,
							       txt     = ln)

		return True

	def checkSeaborn(self,
					 printText = False,
					 logFile   = ''):

		# check whether seaborn is present and flag if not.
		# logFile is '' (do not write to log file file) or
		# a logFile.py class object

		if printText is True:

			ln = 'Checking whether seaborn plotting library present...'
			self.printOrWriteToLog(logFile = logFile,
							  	   txt     = ln)

		try:
			imp.find_module('seaborn')
		except ImportError:
			err = 'Plotting library seaborn not found..\n'+\
			 	  'Will not create any plots for current run.\n'+\
			      'Use "pip install seaborn" to install for future use.\n'
			self.printOrWriteToLog(logFile = logFile,
						  	  	   txt     = err)
			return False

		if printText is True:	
			ln = '---> success!\n'
			self.printOrWriteToLog(logFile = logFile,
							       txt     = ln)

		return True

	def printOrWriteToLog(self,
						  logFile = '',
						  txt     = ''):

		# print to command line or write to log file

		if logFile == '':
			print txt
		else:
			logFile.writeToLog(str = txt)