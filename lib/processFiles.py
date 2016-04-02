from CADSCALEITpiped import pipeline as pipe1
from SFALLFFTpiped import pipeline as pipe2
import os
import shutil

class processFiles():
	def __init__(self,inputFile):
		self.inputFile = inputFile

	def runProcessing(self):
		success = self.readMainInputFile()
		if success == False:
			return
		self.checkOutputDirExists()
		self.findFilesInDir()

		self.checkForMultipleDatasets()

		try:
			self.multiDatasets
		except AttributeError:
			return #don't proceed if error in input file

		success = self.checkCorrectInputFormats()
		if success is False:
			return success

		if self.multiDatasets is False:
			self.getCurrentInputParams()
			success = self.runPipelines()
		else:
			for i in range(self.numDsets):
				self.getCurrentInputParams(jobNumber=i)
				success = self.runPipelines()
		return success

	def runPipelines(self):
		# run two processing pipelines sequentially
		self.writePipeline1Inputs()
		success = self.runPipeline1()
		if success == True:
			self.writePipeline2Inputs()
			success = self.runPipeline2()
			if success == True:
				if self.deleteIntermediateFiles.lower() == 'true':
					success = self.cleanUpDir()
				else:
					success = self.cleanUpDir(removeMtzs=False,
											  removeMaps=False,
											  removePdbs=False)
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

		success = self.checkAllRequiredInputsFound()
		return success

	def checkAllRequiredInputsFound(self):
		# check that all required properties have been found
		requiredProps = ['name1','mtz1','mtzlabels1','pdb1','RfreeFlag1',
						 'name2','mtz2','mtzlabels2','pdb2',
						 'name3','mtz3','mtzlabels3',''
						 'sfall_VDWR','densMapType','dir',
						 'FFTmapWeight','deleteIntermediateFiles']

		for prop in requiredProps:
			try:
				getattr(self,prop)
			except AttributeError:
				print 'Necessary input not found: {}'.format(prop)
				return False
		print 'All necessary inputs found in input file'
		return True

	def checkCorrectInputFormats(self):
		# check that the required properties have correct formats
		props = ['mtz1','pdb1','mtz2','pdb2','mtz3']
		for p in props:
			if self.multiDatasets is False:
				success = self.checkFileExists(fileName=getattr(self,p))
				if success is False:
					return False
			else:
				for f in getattr(self,p).split(','):
					success = self.checkFileExists(fileName=f)
					if success is False:
						return False

		if self.multiDatasets is False:
			if self.name1 == self.name2:
				print '"name1" and "name2" inputs must be different - CAD will fail otherwise'
				return False
		else:
			error = False
			if self.repeatedFile1InputsUsed is True:
				for n in self.name2.split(','):
					if self.name1 == n:
						error = True
			else:
				if self.name1 == self.name2:
					error = True
			if error is True:
				print '"name1" and "name2" inputs must be different for each job in batch'
				return False

		if self.densMapType not in ('DIFF','SIMPLE','2FOFC','END'):
			print '"densMapType" input of incompatible format, default is DIFF'
			return False
		if self.deleteIntermediateFiles.lower() not in ('true','false'):
			print '"deleteIntermediateFiles" input of incompatible format (true,false) - case insensitive'
			return False
		try:
			float(self.sfall_VDWR)
		except ValueError:
			print '"sfall_VDWR" input must be a float'
			return False
		if float(self.sfall_VDWR) <= 0:
			print '"sfall_VDWR" input must be a positive float'
			return False

		print 'All input file parameters appear to be of suitable format'
		return True

	def checkFileExists(self,fileName=''):
		if os.path.isfile(fileName) and os.access(fileName,os.R_OK):
			return True
		else:
			print 'File "{}" could not be located..'.format(fileName)
			return False

	def checkForMultipleDatasets(self):
		# multiple processing jobs can be run sequentially if input file contains comma
		# separated lists for particular files. Determine whether this is the
		# case and check that correctly formatted.
		props = ['name2','mtz2','mtzlabels2','pdb2']
		lengths = []
		for p in props:
			val = getattr(self,p)
			lengths.append(len(val.split(',')))
		self.numDsets  = lengths[0]
		if lengths == [1]*4:
			print 'Only single dataset located and will be processed'
			self.multiDatasets = False
		elif lengths[1:] != lengths[:-1]:
			errorStr = 'Error! Input file properties ({}) '.format(','.join(props))+\
				  	   'must have same number of comma-separated inputs'
			print errorStr

		else:
			print 'Multiple datasets located in input file to be processed'
			props = ['name1','mtz1','mtzlabels1','pdb1','RfreeFlag1']
			lengths = []
			for p in props:
				val = getattr(self,p)
				lengths.append(len(val.split(',')))

			if lengths not in ([1]*5,[self.numDsets]*5):
				errorStr ='Error! Input file properties ({})'.format(','.join(props))+\
					  	  'must all each have either 1 or 5 comma-separated inputs'
				print errorStr
			else:
				s = False
				if lengths == [1]*5:
					s = True
				self.repeatedFile1InputsUsed = s
				self.multiDatasets = True

	def checkOutputDirExists(self):
		# check whether output directory exists and make if not
		if not os.path.exists(self.dir):
			print 'Output directory "{}" not found, making directory'.format(self.dir)
			os.makedirs(self.dir)

	def getCurrentInputParams(self,jobNumber=0):
		# select correct input parameters if multiple jobs within main input file.
		# used when writing separate input files for each part of pipeline.
		props1 = ['name2','mtz2','mtzlabels2','pdb2']
		props2 = ['name1','mtz1','mtzlabels1','pdb1','RfreeFlag1']

		if self.multiDatasets == False:
			for prop in props1+props2:
				setattr(self,prop+'_current',getattr(self,prop))
			self.createDatasetName()
			return
		for prop in props1:
			setattr(self,prop+'_current',getattr(self,prop).split(',')[jobNumber])
		for prop in props2:
			if self.repeatedFile1InputsUsed == True:
				setattr(self,prop+'_current',getattr(self,prop))
			else:
				setattr(self,prop+'_current',getattr(self,prop).split(',')[jobNumber])
		self.createDatasetName()

	def createDatasetName(self):
		# create dataset name here, used to distinguish input file naming scheme 
		# (if name2 not the same as pdb2 then this is required)
		self.dsetName = str((self.pdb2_current).split('/')[-1]).strip('.pdb')

	def writePipeline1Inputs(self):

		self.pipe1FileName 	= self.dir+'/'+self.inputFile.split('.')[0]+'_cadscaleit.txt'

		props = {'mtzIn1':'mtz1',
				 'Mtz1LabelName':'mtzlabels1_current',
				 'Mtz1LabelRename':'name1_current',
				 'RfreeFlag1':'RfreeFlag1_current',
				 'mtzIn2':'mtz2_current',
				 'Mtz2LabelName':'mtzlabels2_current',
				 'Mtz2LabelRename':'name2_current',
				 'mtzIn3':'mtz3',
				 'Mtz3LabelName':'mtzlabels3',
				 'Mtz3LabelRename':'name3',
				 'inputPDBfile':'pdb2_current',
				 'densMapType':'densMapType',
				 'deleteMtzs':'deleteIntermediateFiles',
				 'FFTmapWeight':'FFTmapWeight'}

		fileOut1 = open(self.pipe1FileName,'w')
		inputString = ''
		for p in props.keys():
			inputString += '{} {}\n'.format(p,getattr(self,props[p]))
		inputString += 'END'
		fileOut1.write(inputString)
		fileOut1.close()

		self.jobName = '{}-{}'.format(self.name2_current,self.name1_current)

	def writePipeline2Inputs(self):

		self.pipe2FileName 	= self.dir+'/'+self.inputFile.split('.')[0]+'_sfallfft.txt'

		if self.densMapType != '2FOFC':
			self.mtzIn = '{}/{}_SCALEITcombined.mtz'.format(self.dir,self.jobName)
		else:
			self.mtzIn = '{}/{}_sigmaa.mtz'.format(self.dir,self.name2_current)

		props = {'pdbIN':'pdb2_current',
				 'initialPDB':'name1_current',
				 'laterPDB':'name2_current',
				 'sfall_VDWR':'sfall_VDWR',
				 'mtzIN':'mtzIn',
				 'densMapType':'densMapType',
				 'FFTmapWeight':'FFTmapWeight'}

		fileOut2 = open(self.pipe2FileName,'w')
		inputString = ''
		for p in props.keys():
			inputString += '{} {}\n'.format(p,getattr(self,props[p]))
		inputString += 'END'
		fileOut2.write(inputString)
		fileOut2.close()

	def runPipeline1(self):
		self.p1 = pipe1(outputDir=self.dir,
						inputFile=self.pipe1FileName,
						jobName=self.jobName)

		outcome = self.p1.runPipeline()
		if outcome == 0:
			print 'Pipeline ran to completion'
			return True
		else:
			print 'Pipeline failed to run to completion'
			return False

	def runPipeline2(self):
		self.p2 = pipe2(outputDir=self.dir,
						inputFile=self.pipe2FileName,
						jobName=self.jobName)


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

	def cleanUpDir(self,removeMtzs=True,removeMaps=True,removePdbs=True):
		# after successful completion of pipeline clean up working directory
		print '\nCleaning up working directory'

		# distinguish between FFT and END map output formats depending on program used (FFT/END shellscript)
		if self.densMapType == 'END':
			densMapProg = 'END_switchedAxes'
		else:
			densMapProg = 'fft'

		params 			= [self.dir,self.dsetName]
		renameParams	= [self.dir,self.name2_current]
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

		for file in os.listdir(subdir):
			remove = False
			if removeMtzs is True and file.endswith('.mtz'):
					remove = True
			if removeMaps is True and file.endswith('.map'):
					remove = True
			if removePdbs is True and file.endswith('.pdb'):
					remove = True
			if remove is True:
				os.system('rm {}/{}'.format(subdir,file))

		shutil.make_archive(subdir, 'zip', subdir)
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