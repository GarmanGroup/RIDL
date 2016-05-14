from ccp4Job import ccp4Job,checkInputsExist,fillerLine
import os

class SCALEITjob():
	# run SCALEIT job to scale 2nd dataset Fs against 1st datasets	
	def __init__(self,
				 inputMtz  = '',
				 outputMtz = 'scaleit.mtz',
				 Mtz1Label = '',
				 Mtz2Label = '',
				 outputDir = '',
				 runLog    = ''):

		self.inputMtz 	= inputMtz
		self.outputMtz 	= outputMtz
		self.Mtz1Label 	= Mtz1Label
		self.Mtz2Label 	= Mtz2Label
		self.outputDir 	= outputDir
		self.runLog		= runLog
		self.runLog.writeToLog('Running SCALEIT job')

	def run(self):
		inputFiles = [self.inputMtz]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False
		self.runSCALEIT()
		if self.jobSuccess is True:
			self.provideFeedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMtz))
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

	def runSCALEIT(self):
		# run SCALEIT job to scale 2nd dataset Fs against 1st datasets
		# fillerLine()
		self.printPurpose()
		title = 'run of scaleit'

		# run SCALEIT from command line
		self.commandInput1 = 'scaleit '+\
							 'HKLIN {} '.format(self.inputMtz)+\
							 'HKLOUT {}'.format(self.outputMtz)
		self.commandInput2 = 'title SCALEIT JOB\n'+\
							 'NOWT\n'+\
							 'converge NCYC 4\n'+\
							 'converge ABS 0.001\n'+\
							 'converge TOLR -7\n'+\
							 'REFINE ANISOTROPIC\n'+\
							 'auto\n'+\
							 'LABIN  FP=FP_{} SIGFP=SIGFP_{} -\n'.format(self.Mtz1Label,
							 										 	 self.Mtz1Label)+\
							 ' FPH1=FP_{} SIGFPH1=SIGFP_{}\n'.format(self.Mtz2Label,
							  									 	 self.Mtz2Label)+\
							 'END'

		self.outputLogfile = 'SCALEITlogfile.txt'

		# run SCALEIT job
		job = ccp4Job(jobName       = 'SCALEIT',
					  commandInput1 = self.commandInput1,
					  commandInput2 = self.commandInput2,
					  outputDir     = self.outputDir,
					  outputLog     = self.outputLogfile,
					  outputFile    = self.outputMtz)

		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self,includeDir=False):
		# provide some feedback
		if includeDir is False:
			fileIn  = self.inputMtz.split('/')[-1]
			fileOut = self.outputMtz.split('/')[-1]
		else:
			fileIn  = self.inputMtz
			fileOut = self.outputMtz

		print 'SCALEIT Summary:'
		print 'Input mtz file: {}'.format(fileIn)
		print 'Output mtz file: {}'.format(fileOut)

	def printPurpose(self,include=True):
		# provide a summary of what this does (within ETRACK) to the command line
		str = 'Scaling 2nd dataset Fobs column against 1st dataset (low dose) Fobs'
		print str


