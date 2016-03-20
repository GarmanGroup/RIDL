from time import gmtime, strftime

class logFile():
	def __init__(self,fileName):
		self.logFile = fileName
		self.createLogFile()

	def createLogFile(self):
		log = open(self.logFile,'w')
		log.write('Created at: {}\n'.format(self.getTime()))
		log.close()

	def writeToLog(self,logstring):
		# write string to current log file
		currentTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
		with open(self.logFile, "a") as logfile:
			logfile.write('{}\t{}\n'.format(self.getTime(),logstring))
			if logstring.startswith('Running'):
				print '\n'+logstring
			else:
				print logstring

	def getTime(self):
		return strftime("%Y-%m-%d %H:%M:%S", gmtime())
