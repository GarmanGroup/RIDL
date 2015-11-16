from ccp4Job import ccp4Job
from mapTools import mapTools

class FFTjob():

	def __init__(self,mapType,laterPDB,inputMtzFile,outputDir,axes,gridSamps,labels1,labels2):
		self.mapType 		= mapType
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

	def run(self):
		self.runFFT()
		if self.jobSuccess is True:
			self.provideFeedback()
		else:
			return
			
	def runFFT(self):

		self.commandInput1 = '/Applications/ccp4-6.4.0/bin/fft '+\
				 				  'HKLIN {} '.format(self.inputMtzFile)+\
			 	 				  'MAPOUT {} '.format(self.outputMapFile)+\
			 	 				  'SYMINFO /Applications/ccp4-6.4.0/lib/data/syminfo.lib '

		# if FOM damaged exists, sub this line in below:
		# "LABIN F1=%s SIG1=%s F2=%s SIG2=%s " %(str(Fobs_dam),str(SIGobs_dam),str(Fobs_init),str(SIGobs_init)) +\
		#       "PHI=%s W=%s W2=%s\n" %(str(PHIC_dam),str(FOM_dam),str(FOM_init)) +\

		# MODIFIED HERE TO DISTINGUISH BETWEEN SIMPLE AND DIFFERENCE MAP TYPES
		if self.mapType == 'SIMPLE':
			F2Scale = 0.0
		else:
			F2Scale = 1.0

		self.commandInput2 = 	'AXIS {} {} {}\n'.format(self.fastAxis,self.medAxis,self.slowAxis)+\
				 			 	'title FTT DENSITY MAP RUN\n'+\
				 				'grid {} {} {}\n'.format(self.gridSamp1,self.gridSamp2,self.gridSamp3)+\
				 				'xyzlim 0 1 0 1 0 1\n'+\
								'LABIN F1={} SIG1={} F2={} SIG2={} '.format(self.F1,self.SIG1,self.F2,self.SIG2)+\
								'PHI={}\n'.format(self.PHI1)+\
								'SCALE F1 1.0 0.0 F2 {} 0.0\n'.format(F2Scale)+\
								'END'
		# 'PHI={} W={} W2 ={}\n'.format(self.PHI1,self.FOM1,self.FOM2)+\
		# LINE ABOVE ONLY FOR CASE WHERE FOM LINES IN FILE (replace with other line above)

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
