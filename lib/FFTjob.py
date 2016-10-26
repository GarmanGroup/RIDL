from ccp4Job import ccp4Job,checkInputsExist,fillerLine
from mapTools import mapTools

class FFTjob():
	
	def __init__(self,
				 mapType   = 'DIFF',
				 mapTag    = '',
				 FOMweight = False,
				 pdbFile   = '',
				 mtzFile   = '',
				 outputDir = './',
				 axes      = [1,2,3],
				 gridSamps = [0,0,0],
				 labels1   = ['','','',''],
				 labels2   = ['','','',''],
				 F1Scale   = 1.0,
				 F2Scale   = 1.0,
				 runLog    = ''):

		self.mapType 	  = mapType
		self.FOMweight	  = FOMweight
		self.inputMtzFile = mtzFile
		self.F1Scale      = F1Scale
		self.F2Scale      = F2Scale

		# can change unique identifier for map type in file name 
		if mapTag == '':
			self.mapTag = mapType
		else:
			self.mapTag = mapTag

		self.outputMapFile 	= '{}-{}-fft.map'.format(pdbFile.split('_reordered.pdb')[0],self.mapTag)
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

		self.highResCutoff  = ''
		self.lowResCutoff   = ''

	def run(self):

		inputFiles = [self.inputMtzFile]
		if not checkInputsExist(inputFiles,self.runLog):
			return False
		self.runFFT()

		if self.jobSuccess:
			self.provideFeedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMapFile))
			return True

		else:
			err = 'Job did not run successfully, see job log file "{}"'.format(self.outputLogfile)
			self.runLog.writeToLog(err)
			return False
			
	def runFFT(self):

		# run FFT job using the external ccp4Job class

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
			FOMstring, FOMstring2 = '', ''

		# if DIFFERENCE or SIMPLE map types map type chosen
		if self.mapType in ('DIFF','SIMPLE'):
			labinStr 	= 'LABIN F1={} SIG1={} F2={} SIG2={} '.format(self.F1,self.SIG1,self.F2,self.SIG2)
			phiStr  	= 'PHI={} {} PHI2={} {}\n'.format(self.PHI2,FOMstring,self.PHI2,FOMstring2)
			# can distinguish here between SIMPLE and DIFFERENCE map types
			if self.mapType == 'SIMPLE': 
				self.F2Scale = 0.0

			scaleStr 	= 'SCALE F1 {} 0.0 F2 {} 0.0\n'.format(self.F1Scale,self.F2Scale)

		# choose sigmaa-derived FWT and PHIC for 2FOFC map setting
		if self.mapType == '2FOFC':
			labinStr 	= 'LABIN F1={}'.format(self.F1)
			phiStr 		= 'PHI={}'.format(self.PHI1)
			scaleStr 	= 'SCALE F1 1.0 0.0'

		# for case of FC map generation
		if self.mapType == 'FC':
			labinStr 	= 'LABIN F1={}'.format(self.F1)
			phiStr 		= 'PHI={}'.format(self.PHI1)
			scaleStr 	= 'SCALE F1 {} 0.0'.format(self.F1Scale)

		# include a resolution cutoff if specified (not used for Fcalc maps)
		resStr = ''
		if self.mapType != 'FC':
			if self.highResCutoff != '':
				resStr = 'RESOLUTION {}'.format(self.highResCutoff)
				if self.lowResCutoff != '':
					resStr += ' {}'.format(self.lowResCutoff)
				resStr += '\n'

		self.commandInput2 = 'AXIS {} {} {}\n'.format(self.fastAxis,self.medAxis,self.slowAxis)+\
				 			 'title FTT DENSITY MAP RUN\n'+\
			 				 'grid {} {} {}\n'.format(self.gridSamp1,self.gridSamp2,self.gridSamp3)+\
			 				 'xyzlim 0 1 0 1 0 1\n'+\
			 				 '{}'.format(resStr)+\
							 '{} {}\n{}\n'.format(labinStr,phiStr,scaleStr)+\
							 'END'

		self.outputLogfile = 'FFTlogfile.txt'

		# run FFT job
		job = ccp4Job(jobName       = 'FFT',
					  commandInput1 = self.commandInput1,
					  commandInput2 = self.commandInput2,
					  outputDir     = self.outputDir,
					  outputLog     = self.outputLogfile,
					  outputFile    = self.outputMapFile)

		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self,
						includeDir = False):

		# provide some feedback

		if not includeDir:
			fileIn  = self.inputMtzFile.split('/')[-1]
			fileOut = self.outputMapFile.split('/')[-1]
		else:
			fileIn  = self.inputMtzFile
			fileOut = self.outputMapFile

		txt = 'FFT Summary:\n'+\
			  'Input mtz file: {}\n'.format(fileIn)+\
			  'Output map file: {}'.format(fileOut)
		self.runLog.writeToLog(txt)

		Map = mapTools(mapName = self.outputMapFile,
					   logFile = self.runLog)
		Map.printMapInfo()

	def printPurpose(self,
					 include = True):

		# provide a summary of what this does 
		# (within RIDL) to the command line

		if self.mapTag == 'DIFF':
			ln = 'Generating Fourier difference map over crystal unit cell,\n'
		if self.mapTag == 'SIMPLE':
			ln = 'Generating Fobs electron density map over crystal unit cell,\n'
		if self.mapTag == '2FOFC':
			ln = 'Generating 2Fo-Fc electron density map over crystal unit cell,\n'
		if self.mapTag == 'FC':
			ln = 'Generating Fcalc electron density map over crystal unit cell,\n'
		ln += 'with same grid sampling dimensions as SFALL-output atom-tagged .map file'
		self.runLog.writeToLog(ln)

