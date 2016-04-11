import os.path
from time import gmtime, strftime
from CADjob import CADjob
from SCALEITjob import SCALEITjob
from logFile import logFile
from SIGMAAjob import SIGMAAjob

class pipeline():
	# class to run CAD job to combine F and SIGF columns from
	# two merged mtz files, then to scale the 2nd datasets F 
	# structure factors against the 1st datasets
	def __init__(self,outputDir='',inputFile='',jobName='untitled-job'):

		# specify where output files should be written
		self.outputDir 			= outputDir
		self.makeOutputDir()
		self.findFilesInDir() 	
		self.txtInputFile 		= inputFile
		self.jobName 			= jobName
		self.runLog 			= logFile(fileName='{}{}_runLog1.log'.format(self.outputDir,jobName),
										  fileDir=self.outputDir)

		# specify output files for parts of pipeline
		self.CADoutputMtz 		= '{}{}_CADcombined.mtz'.format(self.outputDir,self.jobName)
		self.SCALEITinputMtz 	= self.CADoutputMtz
		self.SCALEIToutputMtz 	= '{}{}_SCALEITcombined.mtz'.format(self.outputDir,self.jobName)

	def makeOutputDir(self):
		# if the above sub directory does not exist, make it
		if not os.path.exists(self.outputDir):
			os.makedirs(self.outputDir)
			print 'New sub directory "{}" made to contain output files'.format(self.outputDir)

	def runPipeline(self):

		# read input file
		success = self.readInputs()	
		if success is False:
			return 1

		# copy input mtz files to working directory and rename
		self.moveInputMtzs()

		# run SIGMAA job if required to generate a new FOM weight column
		if self.FFTmapWeight == 'recalculate':

			if self.densMapType == '2FOFC':
				mtzLbls_in  = self.Mtz2LabelName
				mtzLbls_out = self.Mtz2LabelRename
			else: 
				mtzLbls_in  = self.Mtz1LabelName
				mtzLbls_out = self.Mtz1LabelName

			sigmaa = SIGMAAjob(self.SIGMAAinputMtz,mtzLbls_in,mtzLbls_out,self.RfreeFlag1,
							   self.inputPDBfile,self.outputDir,self.runLog)	
			success = sigmaa.run()
			if success is False: 
				return 2

			# if 2FO-FC map required, use FWT column from sigmaa-output mtz (we are done here)
			if self.densMapType == '2FOFC':
				self.cleanUpDir()
				return 0

			self.CADinputMtz1 = sigmaa.outputMtz
		else:
			self.CADinputMtz1 = self.SIGMAAinputMtz

		# run CAD job 
		cad = CADjob(self.CADinputMtz1,self.CADinputMtz2,self.CADinputMtz3,
					 self.Mtz1LabelName,self.Mtz2LabelName,self.Mtz3LabelName,
					 self.Mtz1LabelRename,self.Mtz2LabelRename,self.Mtz3LabelRename,
					 self.CADoutputMtz,self.outputDir,self.runLog,self.FFTmapWeight)
		success = cad.run()
		if success is False:
			return 3

 		# run SCALEIT job 
		scaleit = SCALEITjob(self.SCALEITinputMtz,self.SCALEIToutputMtz,
							 self.Mtz1LabelRename,self.Mtz2LabelRename,
							 self.outputDir,self.runLog)
		success = scaleit.run()
		if success is False:
			return 4

		# end of pipeline reached
		self.cleanUpDir()	
		return 0

	def readInputs(self):
		# open input file and parse inputs for CAD job.
		# if Input.txt not found, flag error
		if self.checkFileExists(self.txtInputFile) is False:
			self.runLog.writeToLog(str='Required input file {} not found..'.format(self.txtInputFile))
			return False

		self.runLog.writeToLog(str='Reading inputs from {}'.format(self.txtInputFile))

		# parse input file
		inputFile = open(self.txtInputFile,'r')
		for l in inputFile.readlines():
			try:
				l.split()[1]
			except IndexError:
				continue # ignore blank lines
			if l.strip()[0] == '#':
				continue # ignore commented lines
			if l.split()[0] == 'END':
				break
			setattr(self,l.split()[0],l.split()[1])
		inputFile.close()

		# check that all required properties have been found
		requiredProps = ['mtzIn1','Mtz1LabelName','RfreeFlag1','Mtz1LabelRename',
						 'mtzIn2','Mtz2LabelName','Mtz2LabelRename',
						 'mtzIn3','Mtz3LabelName','Mtz3LabelRename',
						 'inputPDBfile','densMapType','deleteMtzs']

		for prop in requiredProps:
			try:
				getattr(self,prop)
			except AttributeError:
				print 'Necessary input not found: {}'.format(prop)
				return False
		return True

	def moveInputMtzs(self):
		# move input mtz files to working directory and rename as suitable
		if self.densMapType == '2FOFC':
			self.SIGMAAinputMtz  = '{}{}.mtz'.format(self.outputDir,self.Mtz2LabelRename.strip())
			os.system('cp {} {}'.format(self.mtzIn2,self.SIGMAAinputMtz))
		else:
			self.SIGMAAinputMtz  = '{}{}.mtz'.format(self.outputDir,self.Mtz1LabelRename.strip())
			self.CADinputMtz2 	 = '{}{}.mtz'.format(self.outputDir,self.Mtz2LabelRename.strip())
			self.CADinputMtz3    = '{}{}.mtz'.format(self.outputDir,self.Mtz3LabelRename.strip())
			os.system('cp {} {}'.format(self.mtzIn1,self.SIGMAAinputMtz))
			os.system('cp {} {}'.format(self.mtzIn2,self.CADinputMtz2))
			os.system('cp {} {}'.format(self.mtzIn3,self.CADinputMtz3))

	def deleteNonFinalMtzs(self):
		return
		# give option to delete all mtz files within output directory except the final 
		# resulting mtz for job - used to save room if necessary
		if self.deleteMtzs.lower() != 'true': 
			return
		if self.densMapType == '2FOFC':
			fileEnd = 'sigmaa.mtz'
		else:
			fileEnd = 'SCALEITcombined.mtz'
		for f in os.listdir(self.outputDir): 
			if (f.endswith('.mtz') and not f.endswith(fileEnd)) or f.endswith('.tmp'):
				os.system('rm {}{}'.format(self.outputDir,f))

	def cleanUpDir(self):
		# give option to clean up working directory 
		# delete non-final mtz files
		print 'Cleaning up working directory...\n'
		self.deleteNonFinalMtzs()
		# move txt files to subdir
		os.system('mkdir {}txtFiles/'.format(self.outputDir))
		for file in os.listdir(self.outputDir): 
			if file.endswith('.txt') and file not in self.filesInDir:
				os.system('mv {}{} {}txtFiles/{}'.format(self.outputDir,file,self.outputDir,file))

	def findFilesInDir(self):
		# find files initially in working directory
		self.filesInDir = os.listdir(self.outputDir)

	def checkFileExists(self,filename):
		# method to check if file exists
		if os.path.isfile(filename) is False:
			ErrorString = 'File {} not found'.format(filename)
			print ErrorString
			self.runLog.writeToLog(str=ErrorString)
			return False
		else:
			return True
