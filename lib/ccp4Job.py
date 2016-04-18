import os

class ccp4Job():

	def __init__(self,
				 jobName       = 'untitled',
				 commandInput1 = '',
				 commandInput2 = '',
				 outputDir     = './',
				 outputLog     = '',
				 outputFile    = ''):

		self.jobName 		= jobName
		self.commandInput1 	= commandInput1
		self.commandInput2 	= commandInput2
		self.outputDir 		= outputDir
		self.outputLogfile 	= outputLog
		self.outputFile 	= outputFile # used to check success of job

		# automatically run ccp4 program
		self.runCCP4program()

	def runCCP4program(self):
		# generic method to run a ccp4 program on command line

		# write commandInput2 to txt file
		textinput = open('{}inputfile.txt'.format(self.jobName),'w')
		textinput.write(self.commandInput2)
		textinput.close()

		# run ccp4 program job
		os.system('{} < {}inputfile.txt > {}'.format(self.commandInput1,
												     self.jobName,
												     self.outputLogfile))

		# move ccp4 job input and log files to working sub directory
		os.system('mv {}inputfile.txt {}{}inputfile.txt'.format(self.jobName,
			    												self.outputDir,
			    												self.jobName))
		os.system('mv {} {}{}'.format(self.outputLogfile,
									  self.outputDir,
									  self.outputLogfile))

	def checkJobSuccess(self):
		# job success checked, based on whether output files exist
		ErrorString = '{} did not proceed to completion'.format(self.jobName)
		if os.path.isfile(self.outputFile) is False:
			print ErrorString
			return False
		else:
			return True

def checkInputsExist(inputFiles,runLog):
	# check if input files exist for a ccp4 job
	for fileName in inputFiles:
		if os.path.isfile(fileName) is False:
			str = 'Failed to find required input file "{}"'.format(fileName)
			runLog.writeToLog(str)
			return False
	else:
		runLog.writeToLog('Input files:')
		for fileName in inputFiles:
			runLog.writeToLog('{}'.format(fileName))	
		return True

def fillerLine(long=True,linebreak=True):
	if long is True:
		ln = '--------------------------'*2
	else:
		ln = '--------------------------'
	if linebreak is True:
		ln = '\n'+ln
	print ln


