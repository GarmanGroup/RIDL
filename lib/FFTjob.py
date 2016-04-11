from ccp4Job import ccp4Job,checkInputsExist,fillerLine
from mapTools import mapTools

class FFTjob():
	def __init__(self,mapType,FOMweight,inputPDBname,inputMtzFile,
				 outputDir,axes,gridSamps,labels1,labels2,runLog):

		self.mapType 		= mapType
		self.FOMweight		= FOMweight
		self.inputMtzFile 	= inputMtzFile
		self.outputMapFile 	= '{}_fft.map'.format(inputPDBname.split('_reordered.pdb')[0])
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
		# run FFT job using the external ccp4Job class
		fillerLine()
		self.printPurpose()
		self.commandInput1 = 'fft '+\
				 			 'HKLIN {} '.format(self.inputMtzFile)+\
			 	 			 'MAPOUT {} '.format(self.outputMapFile)+\
			 	 			 'SYMINFO syminfo.lib '

		# if FOM is specified the weighting is applied to the map
		if self.FOMweight == 'recalculate' or 'preset,' in self.FOMweight:
			FOMstring  = 'W={}'.format(self.FOM2)
			FOMstring2 = 'W2={}'.format(self.FOM2)
		else:
			FOMstring,FOMstring2 = '',''

		# if DIFFERENCE or SIMPLE map types map type chosen
		if self.mapType in ('DIFF','SIMPLE'):
			labinStr 	= 'LABIN F1={} SIG1={} F2={} SIG2={} '.format(self.F1,self.SIG1,self.F2,self.SIG2)
			phiStr  	= 'PHI={} {} PHI2={} {}\n'.format(self.PHI2,FOMstring,self.PHI2,FOMstring2)
			# can distinguish here between SIMPLE and DIFFERENCE map types
			if self.mapType == 'SIMPLE': 
				F2Scale = 0.0
			else: 
				F2Scale = 1.0
			scaleStr 	= 'SCALE F1 1.0 0.0 F2 {} 0.0\n'.format(F2Scale)

		# choose sigmaa-derived FWT and PHIC for 2FOFC map setting
		if self.mapType == '2FOFC':
			labinStr 	= 'LABIN F1={}'.format(self.F1)
			phiStr 		= 'PHI={}'.format(self.PHI1)
			scaleStr 	= 'SCALE F1 1.0 0.0'

		self.commandInput2 = 	'AXIS {} {} {}\n'.format(self.fastAxis,self.medAxis,self.slowAxis)+\
				 			 	'title FTT DENSITY MAP RUN\n'+\
				 				'grid {} {} {}\n'.format(self.gridSamp1,self.gridSamp2,self.gridSamp3)+\
				 				'xyzlim 0 1 0 1 0 1\n'+\
								'{} {}\n{}\n'.format(labinStr,phiStr,scaleStr)+\
								'END'

		self.outputLogfile = 'FFTlogfile.txt'

		# run FFT job
		job = ccp4Job('FFT',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.outputMapFile)
		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self,includeDir=False):
		# provide some feedback
		if includeDir is False:
			fileIn  = self.inputMtzFile.split('/')[-1]
			fileOut = self.outputMapFile.split('/')[-1]
		else:
			fileIn  = self.inputMtzFile
			fileOut = self.outputMapFile

		print 'FFT Summary:'
		print 'Input mtz file: {}'.format(fileIn)
		print 'Output map file: {}'.format(fileOut)
		Map = mapTools(self.outputMapFile)
		Map.printMapInfo()

	def printPurpose(self,include=True):
		# provide a summary of what this does (within ETRACK) to the command line
		if self.mapType == 'DIFF':
			str = 'Generating Fourier difference map over crystal unit cell,\n'
		if self.mapType == 'SIMPLE':
			str = 'Generating Fobs electron density map over crystal unit cell,\n'
		if self.mapType == '2FOFC':
			str = 'Generating 2Fo-Fc electron density map over crystal unit cell,\n'
		str += 'with same grid sampling dimensions as SFALL-output atom-tagged .map file'
		print str



