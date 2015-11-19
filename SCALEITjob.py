from ccp4Job import ccp4Job

class SCALEITjob():
	# run SCALEIT job to scale 2nd dataset Fs against 1st datasets	
	def __init__(self,SCALEITinputMtz,SCALEIToutputMtz,Mtz1Label,Mtz2Label,outputDir,runLog):
		self.SCALEITinputMtz 	= SCALEITinputMtz
		self.SCALEIToutputMtz 	= SCALEIToutputMtz
		self.Mtz1Label 			= Mtz1Label
		self.Mtz2Label 			= Mtz2Label
		self.outputDir 			= outputDir
		self.runLog				= runLog

		runLog.writeToLog('Running SCALEIT job')

	def checkInputsExist(self):
		# check if input mtz files exist
		if self.checkFileExists(self.SCALEITinputMtz) is False:
			string = 'Failed to find required input file'
			print string
			runLog.writeToLog(string)
			return False
		else:
			runLog.writeToLog('Input files:')
			runLog.writeToLog('{}'.format(self.SCALEITinputMtz))	
			return True

	def run(self):
		if self.checkInputsExist() is False:
			return
		self.runSCALEIT()
		if self.jobSuccess is True:
			self.provideFeedback()
		else:
			string = 'Job did not run successfully, see job log file "{}"',format(self.outputLogfile)
			runLog.writeToLog(string)
			return

	def runSCALEIT(self):
		# run SCALEIT job to scale 2nd dataset Fs against 1st datasets
		title = 'run of scaleit'

		# run SCALEIT from command line
		self.commandInput1 = 'scaleit '+\
							 'HKLIN {} '.format(self.SCALEITinputMtz)+\
							 'HKLOUT {}'.format(self.SCALEIToutputMtz)
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
		job = ccp4Job('SCALEIT',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.SCALEIToutputMtz)
		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self):
		# provide some feedback
		print '--------------------------'
		print 'CAD Summary:'
		print 'Input mtz files: {}'.format(self.SCALEITinputMtz)
		print 'Output mtz file: {}'.format(self.SCALEIToutputMtz)
		print '--------------------------'
