from CADSCALEITpiped import pipeline as pipe1
from SFALLFFTpiped import pipeline as pipe2
import os

class processFiles():
	def __init__(self,inputFile):
		self.inputFile = inputFile

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
		requiredProps = ['name1','mtz1','mtzlabels1','pdb1',
						 'name2','mtz2','mtzlabels2','pdb2',
						 'name3','mtz3','mtzlabels3',''
						 'sfall_VDWR','FFTmapType',
						 'dir','FFTmapWeight']

		for prop in requiredProps:
			try:
				getattr(self,prop)
			except attributeError:
				print 'Necessary input not found: {}'.format(prop)
				return False
		print 'All necessary inputs found in input file'
		return True

	def writePipeline1Inputs(self):

		self.pipe1FileName 	= self.dir+'/'+self.inputFile.split('.')[0]+'_cadscaleit.txt'

		fileOut1 = open(self.pipe1FileName,'w')
		inputString = 'filename1 {}\n'.format(self.mtz1)+\
					  'labels1 {}\n'.format(self.mtzlabels1)+\
					  'label1rename {}\n'.format(self.name1)+\
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
					  'runname {}\n'.format(self.jobName)+\
					  'initialPDB {}\n'.format(self.name1)+\
					  'laterPDB {}\n'.format(self.name2)+\
					  'foldername {}\n'.format(self.dir)+\
					  'sfall_VDWR {}\n'.format(self.sfall_VDWR)+\
					  'mtzIN {}/{}_SCALEITcombined.mtz\n'.format(self.dir,self.jobName)+\
					  'FFTmapType {}\n'.format(self.FFTmapType)+\
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
		self.p2 = pipe2(self.dir,self.pipe2FileName)
		outcome = self.p2.runPipeline()
		if outcome == 0:
			print 'Pipeline ran to completion'
			return True
		else:
			print 'Pipeline failed to run to completion'
			return False

	def cleanUpDirectory(self):
		# after successful completion of pipeline clean up working directory
		print 'Cleaning up directory'
		originalFiles 	= [self.mtz1,self.mtz2,self.mtz3,self.pdb1,self.pdb2]
		keyLogFiles 	= [self.p1.runLog.logFile,self.p2.runLog.logFile]
		mapFiles 		= ['{}/{}_fft_cropped_cropped.map'.format(self.dir,self.name2),
						   '{}/{}_sfall_cropped.map'.format(self.dir,self.name2)]

		subdir = '{}/additionalFiles'.format(self.dir)
		os.system('mkdir {}'.format(subdir))

		for file in os.listdir(self.dir):
			fileName = '{}/{}'.format(self.dir,file)
			if fileName not in originalFiles+keyLogFiles+mapFiles+[subdir]:
				os.system('mv {} {}/{}'.format(fileName,subdir,file))

		os.system('tar -zcvf {}.tar.gz {}'.format(subdir,subdir))
		os.system('rm -rf {}'.format(subdir))

		# rename final map files
		params = [self.dir,self.name2,self.dir,self.name2]
		os.system('mv {}/{}_fft_cropped_cropped.map {}/{}_dens.map'.format(*params))
		os.system('mv {}/{}_sfall_cropped.map {}/{}_atoms.map'.format(*params))




