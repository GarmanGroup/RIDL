import os.path
from time import gmtime, strftime

class pipeline():
	# class to run CAD job to combine F and SIGF columns from
	# two merged mtz files, then to scale the 2nd datasets F 
	# structure factors against the 1st datasets
	def __init__(self,where):

		# specify where output files should be written
		self.outputDir 			= where
		# if the above sub directory does not exist, make it
		if not os.path.exists(self.outputDir):
			os.makedirs(self.outputDir)
			print 'New sub directory "{}" made to contain output files'.format(self.outputDir)

		# specify output files for parts of pipeline
		self.CADoutputMtz 		= self.outputDir+'/CADcombined.mtz'
		self.SCALEITinputMtz 	= self.outputDir+'/CADcombined.mtz'
		self.SCALEIToutputMtz 	= self.outputDir+'/SCALEITcombined.mtz'
		self.txtInputFile 		= 'Inputs.txt'
		self.runLogfile 		= self.outputDir+'/runLog.txt'

		# these are indicators that specific jobs have been performed
		self.inputFileRead 		= False
		self.CADbeenRun 		= False

	def writeToLog(self,logstring):
		# write any steps performed by pipeline to log file
		currentTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
		with open(self.runLogfile, "a") as logfile:
			logfile.write('{}\t{}\n'.format(currentTime,logstring))

	def runPipeline(self):
		# run through pipeline in sequential order

		self.readInputs()		# read input file
		self.runCAD()			# run CAD job 
		self.runSCALEIT()		# run SCALEIT job

	def readInputs(self):
		# open input file and parse inputs for CAD job
		self.writeToLog('Reading Inputs from {}'.format(self.txtInputFile))

		# if Input.txt not found, flag error
		if self.checkFileExists(self.txtInputFile) is False:
			return

		# parse input file
		inputFile = open(self.txtInputFile,'r')
		for line in inputFile.readlines():
			if line.split()[0] == 'END':
				break
			elif line[0] == '#':
				continue
			elif line.split()[0] == 'filename1':
				self.CADinputMtz1 	= line.split()[1]
			elif line.split()[0] == 'labels1':
				self.Mtz1LabelName 	= line.split()[1]
			elif line.split()[0] == 'filename2':
				self.CADinputMtz2 	= line.split()[1]
			elif line.split()[0] == 'labels2':
				self.Mtz2LabelName 	= line.split()[1]
			elif line.split()[0] == 'filename3':
				self.CADinputMtz3 	= line.split()[1]
			elif line.split()[0] == 'labels3':
				self.Mtz3LabelName 	= line.split()[1]
			elif line.split()[0] == 'label1rename':
				self.Mtz1LabelRename 	= line.split()[1]
			elif line.split()[0] == 'label2rename':
				self.Mtz2LabelRename 	= line.split()[1]
			elif line.split()[0] == 'label3rename':
				self.Mtz3LabelRename 	= line.split()[1]
		inputFile.close()
		self.inputFileRead = True

	def runCAD(self):
		# run CAD job to combine F and SIGF columns from two merged 
		# mtz files
		self.jobName = 'CAD'
		self.writeToLog('Running CAD job')

		# check if input txt file read in 
		if self.inputFileRead is False:
			ErrorString = 'Need to read input file {} first'.format(self.txtInputFile)
			print ErrorString
			self.writeToLog(ErrorString)
			return 

		# check if input mtz files exist
		if (self.checkFileExists(self.CADinputMtz1) is False or
		    self.checkFileExists(self.CADinputMtz2) is False):
			return
		else:
			self.writeToLog('Input files:')
			self.writeToLog('{}'.format(self.CADinputMtz1))	
			self.writeToLog('{}'.format(self.CADinputMtz2))												  

		# run CAD from command line
		self.commandInput1 	= 	'/Applications/ccp4-6.5/bin/cad '+\
								'HKLIN1 {} '.format(self.CADinputMtz1)+\
								'HKLIN2 {} '.format(self.CADinputMtz2)+\
								'HKLIN3 {} '.format(self.CADinputMtz3)+\
								'HKLOUT {}'.format(self.CADoutputMtz)

		labels = [self.Mtz1LabelName,self.Mtz2LabelName,self.Mtz3LabelName]
		renameLabels = [self.Mtz1LabelRename,self.Mtz2LabelRename,self.Mtz3LabelRename]
		self.commandInput2 	=	'title CAD JOB\n'+\
								'monitor BRIEF\n'+\
								'labin file 1 - \n'+\
								'E1 = F_{} - \n'.format(labels[0])+\
								'E2 = SIGF_{} \n'.format(labels[0])+\
								'labout file 1 - \n'+\
								'E1 = FP_{} - \n'.format(renameLabels[0])+\
								'E2 = SIGFP_{} \n'.format(renameLabels[0])+\
								'ctypin file 1 - \n'+\
								'E1 = F - \n'+\
								'E2 = Q \n'+\
								'labin file 2 - \n'+\
								'E1 = F_{} - \n'.format(labels[1])+\
								'E2 = SIGF_{} \n'.format(labels[1])+\
								'labout file 2 - \n'+\
								'E1 = FP_{} - \n'.format(renameLabels[1])+\
								'E2 = SIGFP_{} \n'.format(renameLabels[1])+\
								'ctypin file 2 - \n'+\
								'E1 = F - \n'+\
								'E2 = Q \n'+\
								'labin file 3 - \n'+\
								'E1 = PHI{} \n'.format(labels[2])+\
								'labout file 3 - \n'+\
								'E1 = PHIC_{} \n'.format(renameLabels[2])+\
								'ctypin file 3 - \n'+\
								'E1 = P \n'

		self.outputLogfile = 'CADlogfile.txt'

		# run CAD job
		self.runCCP4program()
		self.CADbeenRun = True

	def runSCALEIT(self):
		# run SCALEIT job to scale 2nd dataset Fs against 1st datasets
		self.jobName = 'SCALEIT'
		self.writeToLog('Running SCALEIT job')

		# check that CAD job has been run prior to this
		if self.CADbeenRun is False:
			ErrorString = 'Need to run CAD job first'
			print 'Need to run CAD job first'
			self.writeToLog(ErrorString)
			return 

		# check if input mtz file exist
		if self.checkFileExists(self.SCALEITinputMtz) is False:
			return
		else:
			self.writeToLog('Input files:')
			self.writeToLog('{}'.format(self.SCALEITinputMtz))

		# run SCALEIT from command line
		self.commandInput1 = '/Applications/ccp4-6.5/bin/scaleit '+\
						'HKLIN {} '.format(self.SCALEITinputMtz)+\
						'HKLOUT {}'.format(self.SCALEIToutputMtz)
		self.commandInput2 = 'title SCALEIT JOB\n'+\
							 'NOWT\n'+\
							 'converge NCYC 4\n'+\
							 'converge ABS 0.001\n'+\
							 'converge TOLR -7\n'+\
							 'REFINE ANISOTROPIC\n'+\
							 'auto\n'+\
							 'LABIN  FP=FP_{} SIGFP=SIGFP_{} -\n'.format(self.Mtz1LabelRename,
							 										 	 self.Mtz1LabelRename)+\
							 ' FPH1=FP_{} SIGFPH1=SIGFP_{}\n'.format(self.Mtz2LabelRename,
							  									 	 self.Mtz2LabelRename)+\
							 'END'

		self.outputLogfile = 'SCALEITlogfile.txt'

		# run SCALEIT job
		self.runCCP4program()

	def runCCP4program(self):
		# generic method to run a ccp4 program on command line

		# write commandInput2 to txt file
		textinput = open('{}inputfile.txt'.format(self.jobName),'w')
		textinput.write(self.commandInput2)
		textinput.close()

		# run ccp4 program job
		os.system('{} < {}inputfile.txt > {}'.format(self.commandInput1,
												   self.jobName,
												   self.outputLogfile))
		self.writeToLog('Program input file written: {}inputfile.txt'.format(self.jobName))
		self.writeToLog('Program log file written: {}'.format(self.outputLogfile))

		# move ccp4 job input and log files to working sub directory
		os.system('mv {}inputfile.txt {}/{}inputfile.txt'.format(self.jobName,self.outputDir,self.jobName))
		os.system('mv {} {}/{}'.format(self.outputLogfile,self.outputDir,self.outputLogfile))

	def checkFileExists(self,filename):
		# method to check if file exists
		if os.path.isfile(filename) is False:
			ErrorString = 'File {} not found'.format(filename)
			print ErrorString
			self.writeToLog(ErrorString)
			return False
		else:
			return True


