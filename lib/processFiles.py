from CADSCALEITpiped import pipeline as pipe1
from SFALLFFTpiped import pipeline as pipe2
import os

class processFiles():
	def __init__(self,inputFile):
		self.inputFile = inputFile

	def runProcessing(self):
		self.titleCaption('ETRACK file preparation')
		success = self.readMainInputFile()
		if success == False:
			return
		self.checkOutputDirExists()
		self.findFilesInDir()
		self.writePipeline1Inputs()
		success = self.runPipeline1()
		if success == True:
			self.writePipeline2Inputs()
			success = self.runPipeline2()
			if success == True:
				success = self.cleanUpDir()
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
						 'sfall_VDWR','densMapType','dir',
						 'FFTmapWeight','deleteMtzs']

		for prop in requiredProps:
			try:
				getattr(self,prop)
			except AttributeError:
				print 'Necessary input not found: {}'.format(prop)
				return False
		print 'All necessary inputs found in input file'

		# create dataset name here, used to distinguish input file naming scheme 
		# (if name2 not the same as pdb2 then this is required)
		self.dsetName = str((self.pdb2).split('/')[-1]).strip('.pdb')

		return True

	def checkOutputDirExists(self):
		# check whether output directory exists and make if not
		if not os.path.exists(self.dir):
			print 'Output directory "{}" not found, making directory'.format(self.dir)
			os.makedirs(self.dir)

	def writePipeline1Inputs(self):

		self.pipe1FileName 	= self.dir+'/'+self.inputFile.split('.')[0]+'_cadscaleit.txt'

		fileOut1 = open(self.pipe1FileName,'w')
		inputString = 'mtzIn1 {}\n'.format(self.mtz1)+\
					  'Mtz1LabelName {}\n'.format(self.mtzlabels1)+\
					  'Mtz1LabelRename {}\n'.format(self.name1)+\
					  'RfreeFlag1 {}\n'.format(self.RfreeFlag1)+\
					  'mtzIn2 {}\n'.format(self.mtz2)+\
					  'Mtz2LabelName {}\n'.format(self.mtzlabels2)+\
					  'Mtz2LabelRename {}\n'.format(self.name2)+\
					  'mtzIn3 {}\n'.format(self.mtz3)+\
					  'Mtz3LabelName {}\n'.format(self.mtzlabels3)+\
					  'Mtz3LabelRename {}\n'.format(self.name3)+\
					  'inputPDBfile {}\n'.format(self.pdb2)+\
					  'densMapType {}\n'.format(self.densMapType)+\
					  'deleteMtzs {}\n'.format(self.deleteMtzs)+\
					  'END'

		fileOut1.write(inputString)
		fileOut1.close()
		self.jobName = '{}-{}'.format(self.name2,self.name1)

	def writePipeline2Inputs(self):

		self.pipe2FileName 	= self.dir+'/'+self.inputFile.split('.')[0]+'_sfallfft.txt'

		if self.densMapType != '2FOFC':
			mtzIn = '{}_SCALEITcombined.mtz'.format(self.jobName)
		else:
			mtzIn = '{}_sigmaa.mtz'.format(self.name2)

		fileOut2 = open(self.pipe2FileName,'w')
		inputString = 'pdbIN {}\n'.format(self.pdb2)+\
					  'initialPDB {}\n'.format(self.name1)+\
					  'laterPDB {}\n'.format(self.name2)+\
					  'sfall_VDWR {}\n'.format(self.sfall_VDWR)+\
					  'mtzIN {}/{}\n'.format(self.dir,mtzIn)+\
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

	def cleanUpDir(self):
		# after successful completion of pipeline clean up working directory
		print '\nCleaning up working directory'

		# distinguish between FFT and END map output formats depending on program used (FFT/END shellscript)
		if self.densMapType == 'END':
			densMapProg = 'END_switchedAxes'
		else:
			densMapProg = 'fft'

		params 			= [self.dir,self.dsetName]
		renameParams	= [self.dir,self.name2]
		keyLogFiles 	= [self.p1.runLog.logFile,self.p2.runLog.logFile]
		mapFiles 		= ['{}/{}_{}_cropped_cropped.map'.format(self.dir,self.dsetName,densMapProg),
						   '{}/{}_sfall_cropped.map'.format(*params)]
		pdbFiles 		= ['{}/{}_reordered.pdb'.format(*params)]
		outputFiles 	= mapFiles + pdbFiles
		renameFiles 	= ['{}/{}_density.map'.format(*renameParams),
					  	   '{}/{}_atoms.map'.format(*renameParams),
					   	   '{}/{}.pdb'.format(*renameParams)]

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
		for i in range(3): os.system('mv {} {}'.format(outputFiles[i],renameFiles[i]))

		# check that resulting files are found
		self.findFilesInDir()
		for f in renameFiles:
			if f.split('/')[-1] not in self.filesInDir:
				print 'Not all key output files found'
				return False
		return True

	def titleCaption(self,title):
		print '\n\n||========================== {} ==========================||'.format(title)


