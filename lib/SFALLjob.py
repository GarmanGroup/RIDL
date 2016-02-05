from ccp4Job import ccp4Job,checkInputsExist
from mapTools import mapTools

class SFALLjob():

	def __init__(self,inputPDBfile,outputDir,VDWR,symmetrygroup,gridDimensions,mapoutType,runLog):
		self.inputPDBfile 	= inputPDBfile
		self.outputMapFile 	= '{}_sfall.map'.format(inputPDBfile.split('_reordered.pdb')[0])
		self.outputDir 		= outputDir
		self.symGroup 		= symmetrygroup
		self.VDWR 			= VDWR
		self.gridDims 		= gridDimensions
		self.mapoutType 	= mapoutType # atom-map ATMMOD or solvent-map SOLVMAP
		self.runLog 		= runLog
		self.runLog.writeToLog('Running SFALL job')

	def run(self):
		inputFiles = [self.inputPDBfile]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False
		self.runSFALL()
		if self.jobSuccess is True:
			self.provideFeedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMapFile))
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

	def runSFALL(self):
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
		job = ccp4Job('SFALL',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.outputMapFile)
		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self):
		# provide some feedback
		print '--------------------------'
		print 'SFALL Summary:'
		print 'Input pdb file: {}'.format(self.inputPDBfile)
		print 'Output map file: {}'.format(self.outputMapFile)
		Map = mapTools(self.outputMapFile)
		Map.printMapInfo()
		print '--------------------------'
