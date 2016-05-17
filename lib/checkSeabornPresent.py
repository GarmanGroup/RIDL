import imp

def checkSeabornPresent(printText = False,
						logFile   = ''):

	# check whether seaborn is present and flag if not.
	# logFile is '' (do not write to log file file) or
	# a logFile.py class object

	if printText is True:

		ln = 'Checking whether seaborn plotting library present...'
		printOrWriteToLog(logFile = logFile,
						  txt     = ln)

	try:
		imp.find_module('seaborn')
	except ImportError:
		if printText is True:

			err = 'Plotting library seaborn not found..\n'+\
			 	  'Will not create any plots for current run.\n'+\
			      'Use "pip install seaborn" to install for future use.\n'

			printOrWriteToLog(logFile = logFile,
						  	  txt     = err)
		return False

	if printText is True:	
		ln = '---> success!\n'
		printOrWriteToLog(logFile = logFile,
						  txt     = ln)

	return True

def printOrWriteToLog(logFile = '',
					  txt     = ''):

	# print to command line or write to log file

	if logFile == '':
		print txt
	else:
		logFile.writeToLog(str = txt)
