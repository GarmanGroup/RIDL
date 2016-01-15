from ccp4Job import ccp4Job,checkInputsExist
from mapTools import mapTools

class MAPMASKjob():

	def __init__(self,inputMapFile,inputMapFile2,outputDir,runLog):
		self.inputMapFile 	= inputMapFile 
		self.inputMapFile2 	= inputMapFile2
		self.outputMapFile 	= inputMapFile.split('.map')[0]+'_cropped.map'
		self.outputDir 		= outputDir
		self.runLog			= runLog
		self.runLog.writeToLog('Running MAPMASK job')

	def crop2AsymUnit(self):

		inputFiles = [self.inputMapFile]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False

		self.commandInput1 = 'mapmask '+\
				 			 'MAPIN {} '.format(self.inputMapFile)+\
						 	 'MAPOUT {} '.format(self.outputMapFile)+\
						 	 'SYMINFO /Applications/ccp4-6.5/lib/data/syminfo.lib '
		self.commandInput2 = 'EXTEND\nXYZLIM ASU\nEND'
		self.outputLogfile = 'MAPMASKlogfile.txt'

		# run MAPMASK job
		job = ccp4Job('MAPMASK_crop2AsymUnit',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.outputMapFile)
		self.jobSuccess = job.checkJobSuccess()

		# get feedback if successful job
		if self.jobSuccess is True:
			self.provideFeedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMapFile))
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

	def cropMap2Map(self):

		inputFiles = [self.inputMapFile,self.inputMapFile2]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False

		self.commandInput1 = 'mapmask '+\
				 			 'MAPIN {} '.format(self.inputMapFile)+\
				 			 'MAPLIM {} '.format(self.inputMapFile2)+\
						 	 'MAPOUT {} '.format(self.outputMapFile)+\
						 	 'SYMINFO /Applications/ccp4-6.5/lib/data/syminfo.lib '
		self.commandInput2 = 'EXTEND\nXYZLIM MATCH\nEND'
		self.outputLogfile = 'MAPMASKlogfile.txt'

		# run MAPMASK job
		job = ccp4Job('MAPMASK_cropMap2Map',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.outputMapFile)
		self.jobSuccess = job.checkJobSuccess()

		# get feedback if successful job
		if self.jobSuccess is True:
			self.provideFeedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMapFile))
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

	def provideFeedback(self):
		# provide some feedback
		print '--------------------------'
		print 'MAPMASK Summary:'
		print 'Input map file: {}'.format(self.inputMapFile)
		if self.inputMapFile2 != '':
			print 'Input map 2 file: {}'.format(self.inputMapFile2)
		print 'Output map file: {}'.format(self.outputMapFile)
		Map = mapTools(self.outputMapFile)
		Map.printMapInfo()
		print '--------------------------'
