from ccp4Job import ccp4Job,checkInputsExist
from mapTools import mapTools

class FFTjob():

	def __init__(self,mapType,FOMweight,laterPDB,inputMtzFile,outputDir,axes,gridSamps,labels1,labels2,runLog):
		self.mapType 		= mapType
		self.FOMweight		= FOMweight
		self.inputMtzFile 	= inputMtzFile
		self.outputMapFile 	= '{}/{}_fft.map'.format(outputDir,laterPDB)
		self.outputDir		= outputDir
		XYZ 				= ['X','Y','Z']
		self.fastAxis 		= XYZ[int(axes[0])-1]
		self.medAxis 		= XYZ[int(axes[1])-1]
		self.slowAxis 		= XYZ[int(axes[2])-1]
		self.gridSamp1 		= gridSamps[0]
		self.gridSamp2 		= gridSamps[1]
		self.gridSamp3 		= gridSamps[2]
		self.F1 			= labels1[0]
		self.F2 			= labels2[0]
		self.SIG1 			= labels1[1]
		self.SIG2 			= labels2[1]
		self.FOM1 			= labels1[2]
		self.FOM2 			= labels2[2]
		self.PHI1 			= labels1[3]
		self.PHI2 			= labels2[3]
		self.runLog 		= runLog
		self.runLog.writeToLog('Running FFT job')

	def run(self):
		inputFiles = [self.inputMtzFile]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False
		self.runFFT()
		if self.jobSuccess is True:
			self.provideFeedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMapFile))
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False
			
	def runFFT(self):

		self.commandInput1 = 'fft '+\
				 			 'HKLIN {} '.format(self.inputMtzFile)+\
			 	 			 'MAPOUT {} '.format(self.outputMapFile)+\
			 	 			 'SYMINFO /Applications/ccp4-6.5/lib/data/syminfo.lib '

		# can distinguish here between SIMPLE and DIFFERENCE map types
		if self.mapType == 'SIMPLE':
			F2Scale = 0.0
		else:
			F2Scale = 1.0

		#if a FOM is specified exists the weighting is applied to the map
		if self.FOMweight in ('True','TRUE','true'):
			FOMstring  = 'W={}'.format(self.FOM2)
			FOMstring2 = 'W2={}'.format(self.FOM2)
		else:
			FOMstring = ''

		self.commandInput2 = 	'AXIS {} {} {}\n'.format(self.fastAxis,self.medAxis,self.slowAxis)+\
				 			 	'title FTT DENSITY MAP RUN\n'+\
				 				'grid {} {} {}\n'.format(self.gridSamp1,self.gridSamp2,self.gridSamp3)+\
				 				'xyzlim 0 1 0 1 0 1\n'+\
								'LABIN F1={} SIG1={} F2={} SIG2={} '.format(self.F1,self.SIG1,self.F2,self.SIG2)+\
								'PHI={} {} PHI2={} {}\n'.format(self.PHI2,FOMstring,self.PHI2,FOMstring2)+\
								'SCALE F1 1.0 0.0 F2 {} 0.0\n'.format(F2Scale)+\
								'END'

		self.outputLogfile = 'FFTlogfile.txt'

		# run FFT job
		job = ccp4Job('FFT',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.outputMapFile)
		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self):
		# provide some feedback
		print '--------------------------'
		print 'FFT Summary:'
		print 'Input mtz file: {}'.format(self.inputMtzFile)
		print 'Output map file: {}'.format(self.outputMapFile)
		Map = mapTools(self.outputMapFile)
		Map.printMapInfo()
		print '--------------------------'
