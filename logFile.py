from time import gmtime, strftime

class logFile():
	def __init__(self,fileName):
		self.logFile = fileName

		fileName = open(self.logFile,'w')
		fileName.write('Created at: {}'.format(self.getTime()))

	def writeToLog(self,logstring):
		# write string to current log file
		currentTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
		with open(self.logFile, "a") as logfile:
			logfile.write('{}\t{}\n'.format(self.getTime(),logstring))

	def getTime(self):
		return strftime("%Y-%m-%d %H:%M:%S", gmtime())
