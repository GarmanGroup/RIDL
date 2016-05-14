from ccp4Job import ccp4Job,checkInputsExist,fillerLine
from mapTools import mapTools

class MAPMASKjob():

	def __init__(self,
				 mapFile1  = '',
				 mapFile2  = '',
				 outputDir = './',
				 runLog    = ''):

		self.inputMapFile 	= mapFile1 
		self.inputMapFile2 	= mapFile2
		self.outputDir 		= outputDir
		self.runLog			= runLog
		
	def defineCorrectOutputMap(self,switch):

		# give correct naming scheme to output map file

		if switch is True:
			self.outputMapFile = self.inputMapFile.split('.map')[0]+'_switchedAxes.map'
		else:
			self.outputMapFile 	= self.inputMapFile.split('.map')[0]+'_cropped.map'

	def defineCommandInput(self):

		# write the first part of command line input for mapmask run

		ci1 = 'mapmask MAPIN {} '.format(self.inputMapFile)
		if self.inputMapFile2 != '': 
			ci1 += 'MAPLIM {} '.format(self.inputMapFile2)
		ci1 += 'MAPOUT {} '.format(self.outputMapFile)+\
			   'SYMINFO syminfo.lib '
		self.commandInput1 = ci1

	def switchAxisOrder(self,
					    order      = [],
					    symGroup   = "",
					    includeDir = False):

		# switch the axis order of an input .map file. 
		# order = [1,2,3] for example

		self.printPurpose(mode = 'switch axes')
		self.defineCorrectOutputMap(True)
		xyz = {'1':'X','2':'Y','3':'Z'}
		axisOrder = [xyz[str(i)] for i in order]

		inputFiles = [self.inputMapFile]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False

		if includeDir is False:
			map1 = self.inputMapFile.split('/')[-1]
		else:
			map1 = self.inputMapFile

		self.defineCommandInput()
		self.commandInput2 = 'SYMMETRY {}\nAXIS {}\nEND'.format(symGroup,' '.join(axisOrder))
		self.outputLogfile = 'MAPMASKlogfile.txt'

		# run MAPMASK job
		job = ccp4Job(jobName       = 'MAPMASK_switchAxisOrder',
					  commandInput1 = self.commandInput1,
					  commandInput2 = self.commandInput2,
					  outputDir     = self.outputDir,
					  outputLog     = self.outputLogfile,
					  outputFile    = self.outputMapFile)

		self.jobSuccess = job.checkJobSuccess()
		success = self.provideFeedback()
		return success

	def crop2AsymUnit(self,
					  includeDir = False):

		# Crop map 1 to asymmetric unit

		# fillerLine()
		self.printPurpose(mode = 'crop to asym')
		self.defineCorrectOutputMap(False)
		inputFiles = [self.inputMapFile]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False

		if includeDir is False:
			map1 = self.inputMapFile.split('/')[-1]
		else:
			map1 = self.inputMapFile

		self.defineCommandInput()
		self.commandInput2 = 'EXTEND\nXYZLIM ASU\nEND'
		self.outputLogfile = 'MAPMASKlogfile.txt'

		# run MAPMASK job
		job = ccp4Job(jobName       = 'MAPMASK_crop2AsymUnit',
					  commandInput1 = self.commandInput1,
					  commandInput2 = self.commandInput2,
					  outputDir     = self.outputDir,
					  outputLog     = self.outputLogfile,
					  outputFile    = self.outputMapFile)

		self.jobSuccess = job.checkJobSuccess()
		success = self.provideFeedback()
		return success

	def cropMap2Map(self,
					includeDir = False):

		# Crop map 1 to map 2

		# fillerLine()
		self.printPurpose(mode = 'crop to map')
		self.defineCorrectOutputMap(False)
		inputFiles = [self.inputMapFile,self.inputMapFile2]
		if checkInputsExist(inputFiles,self.runLog) is False:
			 return False

		maps = [self.inputMapFile,self.inputMapFile2]
		if includeDir is False:
			maps = [m.split('/')[-1] for m in maps]

		self.defineCommandInput()
		self.commandInput2 = 'EXTEND\nXYZLIM MATCH\nEND'
		self.outputLogfile = 'MAPMASKlogfile.txt'

		# run MAPMASK job
		job = ccp4Job(jobName       = 'MAPMASK_cropMap2Map',
					  commandInput1 = self.commandInput1,
					  commandInput2 = self.commandInput2,
					  outputDir     = self.outputDir,
					  outputLog     = self.outputLogfile,
					  outputFile    = self.outputMapFile)

		self.jobSuccess = job.checkJobSuccess()
		success = self.provideFeedback()
		return success

	def provideFeedback(self,
						includeMapInfo = False):

		# get feedback if successful job

		if self.jobSuccess is True:
			if includeMapInfo is True:
				self.feedback()
			return True
		else:
			err = 'Job did not run successfully, see job log file "{}"'.format(self.outputLogfile)
			self.runLog.writeToLog(err)
			return False

	def feedback(self,
				 includeDir = False):

		# provide some feedback

		if includeDir is False:
			fileIn1  = self.inputMapFile.split('/')[-1]
			if self.inputMapFile2 != '':
				fileIn2  = self.inputMapFile2.split('/')[-1]
			fileOut = self.outputMapFile.split('/')[-1]
		else:
			fileIn1  = self.inputMapFile
			if self.inputMapFile2 != '':
				fileIn2  = self.inputMapFile2
			fileOut = self.outputMapFile

		txt = 'MAPMASK Summary:\n'+\
		      'Input map file: {}\n'.format(fileIn1)
		if self.inputMapFile2 != '':
			txt += 'Input map 2 file: {}\n'.format(fileIn2)
		txt += 'Output map file: {}'.format(fileOut)
		self.runLog.writeToLog(txt)

		Map = mapTools(mapName = self.outputMapFile,
					   logFile = self.runLog)
		Map.printMapInfo()

	def printPurpose(self,
					 include    = True,
					 mode       = 'switch axes',
					 includeDir = False):

		# provide a summary of what this does 
		# (within ETRACK) to the command line
		
		maps = [self.inputMapFile,
				self.inputMapFile2]

		if includeDir is False:
			maps = [m.split('/')[-1] for m in maps]

		if mode == 'switch axes':
			ln = 'Switching map "{}" axis ordering to {}'.format(maps[0],axisOrder)
		if mode == 'crop to asym':
			ln = 'Cropping map "{}" to asym unit'.format(maps[0])
		if mode == 'crop to map':
			ln = 'Cropping map "{}" to map "{}"'.format(*maps)
		self.runLog.writeToLog(ln)



