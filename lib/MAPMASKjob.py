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

		ci1 = 'mapmask MAPIN {} '.format(self.inputMapFile)
		if self.inputMapFile2 != '': 
			ci1 += 'MAPLIM {} '.format(self.inputMapFile2)
		ci1 += 'MAPOUT {} '.format(self.outputMapFile)+\
			   'SYMINFO /Applications/ccp4-6.5/lib/data/syminfo.lib '
		self.commandInput1 = ci1

	def switchAxisOrder(self,order,symGroup):
		# switch the axis order of an input .map file. order = [1,2,3] for example
		xyz = {'1':'X','2':'Y','3':'Z'}
		axisOrder = [xyz[str(i)] for i in order]
		self.outputMapFile 	= inputMapFile.split('.map')[0]+'_switchedAxes.map'

		inputFiles = [self.inputMapFile]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False

		self.runLog.writeToLog('Switching map "{}" axis ordering to {}'.format(self.inputMapFile,axisOrder))
		self.commandInput2 = 'SYMMETRY {}\nAXIS {}\nEND'.format(symGroup,' '.join(axisOrder))
		self.outputLogfile = 'MAPMASKlogfile.txt'

		# run MAPMASK job
		job = ccp4Job('MAPMASK_switchAxisOrder',self.commandInput1,self.commandInput2,
					  self.outputDir,self.outputLogfile,self.outputMapFile)
		self.jobSuccess = job.checkJobSuccess()
		self.provideFeedback()

	def crop2AsymUnit(self):
		# Crop map 1 to asymmetric unit
		inputFiles = [self.inputMapFile]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False

		self.runLog.writeToLog('Cropping map "{}" to asym unit'.format(self.inputMapFile))
		self.commandInput2 = 'EXTEND\nXYZLIM ASU\nEND'
		self.outputLogfile = 'MAPMASKlogfile.txt'

		# run MAPMASK job
		job = ccp4Job('MAPMASK_crop2AsymUnit',self.commandInput1,self.commandInput2,
					  self.outputDir,self.outputLogfile,self.outputMapFile)
		self.jobSuccess = job.checkJobSuccess()
		self.provideFeedback()

	def cropMap2Map(self):
		# Crop map 1 to map 2
		inputFiles = [self.inputMapFile,self.inputMapFile2]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False

		self.runLog.writeToLog('Cropping map "{}" to map "{}"'.format(self.inputMapFile,self.inputMapFile2))
		self.commandInput2 = 'EXTEND\nXYZLIM MATCH\nEND'
		self.outputLogfile = 'MAPMASKlogfile.txt'

		# run MAPMASK job
		job = ccp4Job('MAPMASK_cropMap2Map',self.commandInput1,self.commandInput2,
					  self.outputDir,self.outputLogfile,self.outputMapFile)
		self.jobSuccess = job.checkJobSuccess()
		self.provideFeedback()

	def provideFeedback(self):
		# get feedback if successful job
		if self.jobSuccess is True:
			self.feedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMapFile))
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

	def feedback(self):
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
