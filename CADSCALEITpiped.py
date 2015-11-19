import os.path
from time import gmtime, strftime
from CADjob import CADjob
from SCALEITjob import SCALEITjob
from logFile import logFile

class pipeline():
	# class to run CAD job to combine F and SIGF columns from
	# two merged mtz files, then to scale the 2nd datasets F 
	# structure factors against the 1st datasets
	def __init__(self,where,inputFile):

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
		self.txtInputFile 		= self.outputDir+'/'+inputFile
		self.pipelineLog 		= logFile(self.outputDir+'/runLog.txt')

	def runPipeline(self):
		# read input file
		success = self.readInputs()	
		if success is False:
			return 0	

		# run CAD job 
		cad = CADjob(self.CADinputMtz1,self.CADinputMtz2,self.CADinputMtz3,
		self.Mtz1LabelName,self.Mtz2LabelName,self.Mtz3LabelName,
		self.Mtz1LabelRename,self.Mtz2LabelRename,self.Mtz3LabelRename,
		self.CADoutputMtz,self.outputDir,self.pipelineLog)
		success = cad.run()
		if success is False:
			return 1

 		# run SCALEIT job 
		scaleit = SCALEITjob(self.SCALEITinputMtz,self.SCALEIToutputMtz,
							 self.Mtz1LabelRename,self.Mtz2LabelRename,
							 self.outputDir,self.pipelineLog)
		success = scaleit.run()
		if success is False:
			return 2

	def readInputs(self):
		# open input file and parse inputs for CAD job
		self.pipelineLog.writeToLog('Reading Inputs from {}'.format(self.txtInputFile))

		# if Input.txt not found, flag error
		if self.checkFileExists(self.txtInputFile) is False:
			print 'Required input file {} not found..'.format(self.txtInputFile)
			return False

		# parse input file
		inputFile = open(self.txtInputFile,'r')
		for line in inputFile.readlines():
			if line.split()[0] == 'END':
				break
			elif line[0] == '#':
				continue
			elif line.split()[0] == 'filename1':
				self.CADinputMtz1 		= line.split()[1]
			elif line.split()[0] == 'labels1':
				self.Mtz1LabelName 		= line.split()[1]
			elif line.split()[0] == 'filename2':
				self.CADinputMtz2 		= line.split()[1]
			elif line.split()[0] == 'labels2':
				self.Mtz2LabelName 		= line.split()[1]
			elif line.split()[0] == 'filename3':
				self.CADinputMtz3 		= line.split()[1]
			elif line.split()[0] == 'labels3':
				self.Mtz3LabelName 		= line.split()[1]
			elif line.split()[0] == 'label1rename':
				self.Mtz1LabelRename 	= line.split()[1]
			elif line.split()[0] == 'label2rename':
				self.Mtz2LabelRename 	= line.split()[1]
			elif line.split()[0] == 'label3rename':
				self.Mtz3LabelRename 	= line.split()[1]
		inputFile.close()
		return True

	def checkFileExists(self,filename):
		# method to check if file exists
		if os.path.isfile(filename) is False:
			ErrorString = 'File {} not found'.format(filename)
			print ErrorString
			self.pipelineLog.writeToLog(ErrorString)
			return False
		else:
			return True
