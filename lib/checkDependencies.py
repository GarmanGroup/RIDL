
import imp
import os

class checkDependencies():

	def __init__(self,
				 checkAll = False):

		if checkAll:
			self.checkSeaborn(printText = True)
			self.checkCCP4(printText = True)

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
			if printText is True:

				err = 'CCP4 program suite not found..\n'

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
			if printText is True:

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