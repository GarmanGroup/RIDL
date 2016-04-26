from time import gmtime, strftime

class logFile():
	def __init__(self,fileName='untitled-log',fileDir='',printToScreen=False):
		self.logFile = fileName
		self.createLogFile()
		self.fileDir = fileDir # where majority of files come from
		self.allocateDir()
		self.printToScreen = printToScreen

	def createLogFile(self):
		log = open(self.logFile,'w')
		log.write('Created at: {}\n'.format(self.getTime()))
		log.close()

	def allocateDir(self):
		# this is the specified directory where most files have come from.
		# this is written at the top of the log file to avoid repetition.
		log = open(self.logFile,'a')
		log.write('All files come from the following directory unless otherwise given:\n"{}"\n'.format(self.fileDir))
		log.close()

	def writeToLog(self,str='',strip=True,forcePrint=False):
		# write string to current log file

		# strip away the common directory name 
		if strip is True:
			logstring = str.replace(self.fileDir,'')
		else:
			logstring = str

		currentTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
		with open(self.logFile, "a") as logfile:
			logfile.write('{}\t{}\n'.format(self.getTime(),logstring))

			if self.printToScreen is True or forcePrint is True:
				if logstring.startswith('Running'):
					print '\n'+logstring
				else:
					print logstring

	def getTime(self):
		return strftime("%Y-%m-%d %H:%M:%S", gmtime())

