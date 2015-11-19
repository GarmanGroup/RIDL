from ccp4Job import ccp4Job
import os

class SCALEITjob():
	# run SCALEIT job to scale 2nd dataset Fs against 1st datasets	
	def __init__(self,inputMtz,outputMtz,Mtz1Label,Mtz2Label,outputDir,runLog):
		self.inputMtz 	= inputMtz
		self.outputMtz 	= outputMtz
		self.Mtz1Label 	= Mtz1Label
		self.Mtz2Label 	= Mtz2Label
		self.outputDir 	= outputDir
		self.runLog		= runLog
		self.runLog.writeToLog('Running SCALEIT job')

	def run(self):
		if self.checkInputsExist() is False:
			return False
		self.runSCALEIT()
		if self.jobSuccess is True:
			self.provideFeedback()
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

	def checkInputsExist(self):
		# check if input mtz files exist
		if os.path.isfile(self.inputMtz) is False:
			self.runLog.writeToLog('Failed to find required input file')
			return False
		else:
			self.runLog.writeToLog('Input files:')
			self.runLog.writeToLog('{}'.format(self.inputMtz))	
			self.runLog.writeToLog('Output files:')
			self.runLog.writeToLog('{}'.format(self.outputMtz))	
			return True

	def runSCALEIT(self):
		# run SCALEIT job to scale 2nd dataset Fs against 1st datasets
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
		job = ccp4Job('SCALEIT',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.outputMtz)
		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self):
		# provide some feedback
		print '--------------------------'
		print 'SCALEIT Summary:'
		print 'Input mtz files: {}'.format(self.inputMtz)
		print 'Output mtz file: {}'.format(self.outputMtz)
		print '--------------------------'
