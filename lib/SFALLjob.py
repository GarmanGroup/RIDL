from ccp4Job import ccp4Job,checkInputsExist,fillerLine
from mapTools import mapTools

class SFALLjob():

	def __init__(self,
				 inputPDBfile   = '',
				 outputDir      = './',
				 VDWR           = 1,
				 symmetrygroup  = 'P1',
				 gridDimensions = [],
				 mapoutType     = 'ATMMOD',
				 runLog         = ''):

		self.inputPDBfile  = inputPDBfile
		self.outputMapFile = '{}_sfall.map'.format(inputPDBfile.split('_reordered.pdb')[0])
		self.outputDir 	   = outputDir
		self.symGroup 	   = symmetrygroup
		self.VDWR 		   = VDWR
		self.gridDims 	   = gridDimensions
		self.mapoutType    = mapoutType # atom-map ATMMOD or solvent-map SOLVMAP
		self.runLog 	   = runLog

	def run(self):
		inputFiles = [self.inputPDBfile]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False
		self.runSFALL()
		if self.jobSuccess is True:
			self.provideFeedback()
			return True
		else:
			ln = 'Job did not run successfully, see job log file "{}"'.format(self.outputLogfile)
			self.runLog.writeToLog(ln)
			return False

	def runSFALL(self):

		# run SFALL job using the external ccp4Job class

		self.printPurpose()
		title = 'run of sfall'

		self.commandInput1 = 'sfall '+\
				 'XYZIN {} '.format(self.inputPDBfile)+\
				 'ATOMSF atomsf.lib '+\
			 	 'MAPOUT {} '.format(self.outputMapFile)+\
				 'SYMINFO syminfo.lib '

		if len(self.gridDims) == 0:
			self.commandInput2 = 'MODE ATMMAP {}\n'.format(self.mapoutType)+\
							 	 'SYMMETRY {}\n'.format(self.symGroup)+\
							 	 'VDWR {}\n'.format(self.VDWR)+\
							 	 'title {}\n'.format(title)+\
				   			 	 'END'
		else:
			self.commandInput2 = 'MODE ATMMAP {}\n'.format(self.mapoutType)+\
							 	 'SYMMETRY {}\n'.format(self.symGroup)+\
							 	 'VDWR {}\n'.format(self.VDWR)+\
							 	 'title {}\n'.format(title)+\
							 	 'GRID {} {} {}\n'.format(self.gridDims[0],self.gridDims[1],self.gridDims[2])+\
							 	 'END'

		self.outputLogfile = 'SFALLlogfile.txt'

		# run SFALL job
		job = ccp4Job(jobName       = 'SFALL',
					  commandInput1 = self.commandInput1,
					  commandInput2 = self.commandInput2,
					  outputDir     = self.outputDir,
					  outputLog     = self.outputLogfile,
					  outputFile    = self.outputMapFile)

		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self,
						includeDir = False):

		# provide some feedback

		if includeDir is False:
			fileIn  = self.inputPDBfile.split('/')[-1]
			fileOut = self.outputMapFile.split('/')[-1]
		else:
			fileIn  = self.inputPDBfile
			fileOut = self.outputMapFile

		txt = 'SFALL Summary:\n'+\
			  'Input pdb file: {}\n'.format(fileIn)+\
			  'Output map file: {}'.format(fileOut)
		self.runLog.writeToLog(txt)

		Map = mapTools(mapName = self.outputMapFile,
					   logFile = self.runLog)
		Map.printMapInfo()

	def printPurpose(self,
					 include = True):

		# provide a summary of what this does 
		# (within ETRACK) to the command line

		ln = 'Creating atom-tagged .map file over unit cell for model "{}"'.format(self.inputPDBfile.split('/')[-1])
		self.runLog.writeToLog(ln)


