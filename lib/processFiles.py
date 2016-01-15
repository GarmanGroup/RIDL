from CADSCALEITpiped import pipeline as pipe1
from SFALLFFTpiped import pipeline as pipe2
import os

class processFiles():
	def __init__(self,inputFile):
		self.inputFile = inputFile

	def runProcessing(self):
		success = self.readMainInputFile()
		if success == False:
			return
		self.checkOutputDirExists()
		self.findFilesInDir()
		self.writePipeline1Inputs()
		self.writePipeline2Inputs()
		success = self.runPipeline1()
		if success == True:
			success = self.runPipeline2()
			if success == True:
				success = self.cleanUpDirectory()
		return success

	def readMainInputFile(self):
		# split the input file into separate input files for 
		# each section of the pipeline

		# if Input.txt not found, flag error
		if os.path.isfile(self.inputFile) is False:
			print 'Required input file {} not found..'.format(self.inputFile)
			return

		fileIn = open(self.inputFile,'r')

		for line in fileIn.readlines():
			try:
				line.split()[1]
			except IndexError:
				continue # ignore blank lines
			if line.strip()[0] == '#':
				continue # ignore commented lines
			setattr(self,line.split()[0],line.split()[1])
		fileIn.close()

		# check that all required properties have been found
		requiredProps = ['name1','mtz1','mtzlabels1','pdb1','RfreeFlag1',
						 'name2','mtz2','mtzlabels2','pdb2',
						 'name3','mtz3','mtzlabels3',''
						 'sfall_VDWR','densMapType','dir','FFTmapWeight']

		for prop in requiredProps:
			try:
				getattr(self,prop)
			except AttributeError:
				print 'Necessary input not found: {}'.format(prop)
				return False
		print 'All necessary inputs found in input file'
		return True

	def checkOutputDirExists(self):
		# check whether output directory exists and make if not
		if not os.path.exists(self.dir):
			print 'Output directory "{}" not found, making directory'.format(self.dir)
			os.makedirs(self.dir)

	def writePipeline1Inputs(self):

		self.pipe1FileName 	= self.dir+'/'+self.inputFile.split('.')[0]+'_cadscaleit.txt'

		fileOut1 = open(self.pipe1FileName,'w')
		inputString = 'filename1 {}\n'.format(self.mtz1)+\
					  'labels1 {}\n'.format(self.mtzlabels1)+\
					  'label1rename {}\n'.format(self.name1)+\
					  'RfreeFlag1 {}\n'.format(self.RfreeFlag1)+\
					  'filename2 {}\n'.format(self.mtz2)+\
					  'labels2 {}\n'.format(self.mtzlabels2)+\
					  'label2rename {}\n'.format(self.name2)+\
					  'filename3 {}\n'.format(self.mtz3)+\
					  'labels3 {}\n'.format(self.mtzlabels3)+\
					  'label3rename {}\n'.format(self.name3)+\
					  'inputPDBfile {}\n'.format(self.pdb1)+\
					  'END'
		fileOut1.write(inputString)
		fileOut1.close()
		self.jobName = '{}-{}'.format(self.name2,self.name1)

	def writePipeline2Inputs(self):

		self.pipe2FileName 	= self.dir+'/'+self.inputFile.split('.')[0]+'_sfallfft.txt'

		fileOut2 = open(self.pipe2FileName,'w')
		inputString = 'pdbIN {}\n'.format(self.pdb2)+\
					  'initialPDB {}\n'.format(self.name1)+\
					  'laterPDB {}\n'.format(self.name2)+\
					  'sfall_VDWR {}\n'.format(self.sfall_VDWR)+\
					  'mtzIN {}/{}_SCALEITcombined.mtz\n'.format(self.dir,self.jobName)+\
					  'densMapType {}\n'.format(self.densMapType)+\
					  'FFTmapWeight {}\n'.format(self.FFTmapWeight)+\
					  'END'
		fileOut2.write(inputString)
		fileOut2.close()

	def runPipeline1(self):
		self.p1 = pipe1(self.dir,self.pipe1FileName,self.jobName)
		outcome = self.p1.runPipeline()
		if outcome == 0:
			print 'Pipeline ran to completion'
			return True
		else:
			print 'Pipeline failed to run to completion'
			return False

	def runPipeline2(self):
		self.p2 = pipe2(self.dir,self.pipe2FileName,self.jobName)
		outcome = self.p2.runPipeline()
		if outcome == 0:
			print 'Pipeline ran to completion'
			return True
		else:
			print 'Pipeline failed to run to completion'
			return False

	def findFilesInDir(self):
		# find files initially in working directory
		self.filesInDir = os.listdir(self.dir)

	def cleanUpDirectory(self):
		# after successful completion of pipeline clean up working directory
		print 'Cleaning up directory'

		# distinguish between FFT and END map output formats
		if self.densMapType in ('SIMPLE','DIFF'): densMapProg = 'fft'
		elif self.densMapType in ('END'): densMapProg = 'END'

		params 			= [self.dir,self.name2]
		keyLogFiles 	= [self.p1.runLog.logFile,self.p2.runLog.logFile]
		mapFiles 		= ['{}/{}_{}_cropped_cropped.map'.format(self.dir,self.name2,densMapProg),
						   '{}/{}_sfall_cropped.map'.format(*params)]
		pdbFiles 		= ['{}/{}_reordered.pdb'.format(*params)]
		outputFiles 	= mapFiles + pdbFiles
		renameFiles 	= ['{}/{}_density.map'.format(*params),
					  	   '{}/{}_atoms.map'.format(*params),
					   	   '{}/{}.pdb'.format(*params)]

		subdir = '{}/{}_additionalFiles'.format(self.dir,self.jobName)
		os.system('mkdir {}'.format(subdir))

		for file in os.listdir(self.dir): 
			if file not in self.filesInDir:
				fileName = '{}/{}'.format(self.dir,file)
				if fileName not in keyLogFiles+mapFiles+pdbFiles+[subdir]:
					os.system('mv {} {}/{}'.format(fileName,subdir,file))

		os.system('tar -zcvf {}.tar.gz {}'.format(subdir,subdir))
		os.system('rm -rf {}'.format(subdir))

		# rename final map files
		# params1 = [self.dir,self.name2,densMapProg,self.dir,self.name2]
		# params2 = [self.dir,self.name2,self.dir,self.name2]
		# os.system('mv {}/{}_{}_cropped_cropped.map {}/{}_density.map'.format(*params1))
		# os.system('mv {}/{}_sfall_cropped.map {}/{}_atoms.map'.format(*params2))
		# os.system('mv {}/{}_reordered.pdb {}/{}.pdb'.format(*params2))

		# rename final map files
		for i in range(3): os.system('mv {} {}'.format(outputFiles[i],renameFiles[i]))

		# check that resulting files are found
		self.findFilesInDir()
		for f in renameFiles:
			if f.split('/')[-1] not in self.filesInDir:
				print 'Not all key output files found'
				return False
		return True






