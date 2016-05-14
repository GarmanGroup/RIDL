from CADSCALEITpiped import pipeline as pipe1
from SFALLFFTpiped import pipeline as pipe2
from cleanUpFiles import cleanUpFinalFiles
from logFile import logFile
from errors import error
import shutil
import os

class processFiles():

	def __init__(self,
				 inputFile       = '',
				 proceedToETRACK = False,
				 skipToETRACK	 = False,
				 outputGraphs    = True,
				 eTrackInputName = 'e_Track_inputfile.txt',
				 cleanFinalFiles = False):

		# class to read an input file and generate a set of density 
		# and atom-tagged maps for a damage series. Can handle a single
		# low & high dose dataset pair, or a low dose dataset and a series
		# of high dose datasets. Can optionally proceed directly on to
		# the metric calculation stage of the pipeline with 'proceedToETRACK'
		# or skip directly to it (if the maps have already been created) 
		# with 'skipToETRACK'.

		self.inputFile 		 = inputFile
		self.proceedToETRACK = proceedToETRACK
		self.outputGraphs    = outputGraphs
		self.eTrackInputName = eTrackInputName

		if skipToETRACK is False:
			success = self.runProcessing()
			self.jobSuccess = success
		else:
			success = self.skipToETRACK(inputName = eTrackInputName)
			self.jobSuccess = success

		if self.jobSuccess is True and cleanFinalFiles is True:
				if skipToETRACK is True or proceedToETRACK is True:
					cleanUpFinalFiles(outputDir = self.dir)

	def skipToETRACK(self,inputName = 'input.txt'):

		# do not generate maps for job, but proceed 
		# directly to metric calculations, if correct 
		# input file exists in working directory

			success = self.readMainInputFile()
			if success is False:
				return False

			self.checkOutputDirExists(makeProcessDir = False)
			self.findFilesInDir(mapProcessDir = False)

			if self.eTrackInputName not in self.filesInDir:
				err = 'Unable to find input file '+\
					  '"{}" in "{}"\n'.format(self.eTrackInputName,self.dir)+\
					  'Ensure "python ETRACK -i <input.txt> '+\
					  '-p" has been run prior to -c'
				error(text = err)
				return False

			success = self.writeETRACKInputFile(run          = True,
							 		            imports      = True,
							 		            skipToETRACK = True)
			return True

	def runProcessing(self):

		success = self.readMainInputFile()
		if success is False:
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

		success = self.checkMtzLabelsExist()
		if success is False:
			return success

		if self.multiDatasets is False:
			self.getCurrentInputParams()
			success = self.runMapGenerationPipelines()
		else:
			for i in range(self.numDsets):
				self.getCurrentInputParams(jobNumber = i)
				success = self.runMapGenerationPipelines()
				if success is False:
					return success

		if success is True:
			if self.proceedToETRACK is True:
				self.writeETRACKInputFile()

		return success

	def runMapGenerationPipelines(self):

		# run two map generation pipelines sequentially 
		# for a single dataset

		self.setJobName()
		self.startLogFile()
		self.writePipeline1Inputs()
		success = self.runPipeline1(logFile = self.logFile)
		if success is True:
			self.writePipeline2Inputs()
			success = self.runPipeline2(logFile = self.logFile)
			if success is True:
				if self.deleteIntermediateFiles.lower() == 'true':
					success = self.cleanUpIntermediateFiles()
				else:
					success = self.cleanUpIntermediateFiles(removeMtzs = False,
											  				removeMaps = False,
											  				removePdbs = False)
		return success

	def readMainInputFile(self):

		# split the input file into separate input 
		# files for each section of the pipeline

		# if Input.txt not found, flag error
		if os.path.isfile(self.inputFile) is False:
			err = 'Required input file {} not found..'.format(self.inputFile)
			error(text = err)
			return False

		fileIn = open(self.inputFile,'r')

		for line in fileIn.readlines():
			try:
				line.split()[1]
			except IndexError:
				continue # ignore blank lines
			if line.strip()[0] == '#':
				continue # ignore commented lines
			inputPart = ''.join(line.split()[1:])
			setattr(self,line.split()[0],inputPart)
		fileIn.close()

		success = self.checkAllRequiredInputsFound()
		return success

	def checkAllRequiredInputsFound(self):

		# check that all required properties have been found

		requiredProps = ['name1',
						 'mtz1',
						 'mtzlabels1',
						 'pdb1',
						 'RfreeFlag1',
						 'name2',
						 'mtz2',
						 'mtzlabels2',
						 'pdb2',
						 'name3',
						 'mtz3',
						 'phaseLabel',
						 'FcalcLabel',
						 'sfall_VDWR',
						 'densMapType',
						 'dir',
						 'FFTmapWeight',
						 'deleteIntermediateFiles',
						 'dose1',
						 'dose2']

		for prop in requiredProps:
			try:
				getattr(self,prop)
			except AttributeError:
				err = 'Necessary input not found: {}'.format(prop)
				error(text = err)
				return False

		print 'All necessary inputs found in input file'
		return True

	def checkCorrectInputFormats(self):

		# check that the required properties have correct formats

		props = ['mtz1',
				 'pdb1',
				 'mtz2',
				 'pdb2',
				 'mtz3']

		for p in props:
			if self.multiDatasets is False:
				success = self.checkFileExists(fileName = getattr(self,p))
				if success is False:
					return False
			else:
				for f in getattr(self,p).split(','):
					success = self.checkFileExists(fileName=f)
					if success is False:
						return False

		if self.multiDatasets is False:
			if self.name1 == self.name2:
				err = '"name1" and "name2" inputs must be different, '+\
					  'otherwise CAD will fail. Both currently set '+\
					  'as "{}".'.format(self.name1)
				error(text = err)
				return False

		else:
			isError = False
			if self.repeatedFile1InputsUsed is True:
				for n in self.name2.split(','):
					if self.name1 == n:
						isError = True
						sameName = n
						break
			else:
				for i,(n1, n2) in enumerate(zip(self.name1.split(','), self.name2.split(','))):
					if n1 == n2:
						isError = True
						sameName = n1
						break

			if isError is True:
				err = '"name1" and "name2" inputs must be different '+\
					  'for each job in batch, otherwise CAD will fail. '+\
					  'Currently for one batch, the "initial" and "later" '+\
					  'datasets are both currently called "{}".'.format(sameName)
				error(text = err)
				return False

		# add to name3 to identify it from other datasets
		self.name3 += '-phases'

		if self.densMapType not in ('DIFF','SIMPLE','2FOFC','END'):
			err = '"densMapType" input of incompatible format, (default is "DIFF"). '+\
				  'Currently set as "{}" in input file.'.format(self.densMapType)
			error(text = err)
			return False

		if self.deleteIntermediateFiles.lower() not in ('true','false'):
			err = '"deleteIntermediateFiles" input of incompatible '+\
				  'format ("true","false"), case insensitive. '+\
				  'Currently set as "{}" in input file'.format(self.deleteIntermediateFiles)
			error(text = err)
			return False
		try:
			float(self.sfall_VDWR)
		except ValueError:
			err = '"sfall_VDWR" input must be a float. '+\
			      'Currently set as "{}" in input file.'.format(self.sfall_VDWR)
			error(text = err)
			return False

		if float(self.sfall_VDWR) <= 0:
			err = '"sfall_VDWR" input must be a positive float. '+\
				  'Currently set as "{}" in input file.'.format(self.sfall_VDWR)

			error(text = err)
			return False

		if 'preset' not in self.FFTmapWeight and self.FFTmapWeight not in ('recalculate','False'):
			err = '"FFTmapWeight" input must take value in ("recalculate", '+\
				  '"False","preset,x") where "x" is the FOM weight label '+\
				  '(of the form FOMx) of a preset FOM column within the '+\
				  'input .mtz file. Currently set as "{}" in input file'.format(self.FFTmapWeight)
			error(text = err)
			return False

		if self.FFTmapWeight.startswith('preset') is True:
			if self.FFTmapWeight.replace('preset','')[0] != ',':
				err = 'If "FFTmapWeight" input is specified by "preset", '+\
					  'then this must be followed a comma and then the FOM '+\
					  'label name within the input .mtz file (i.e. "preset,x" '+\
					  'for column label "FOMx"). If column label is simply '+\
					  '"FOM" then use "preset,". Currently set as "{}" in input file'.format(self.FFTmapWeight)
				error(text = err)
				return False

		doseStr = 'Doses can be calculated using RADDOSE-3D (visit '+\
				  'www.raddo.se for more details). If no doses have '+\
				  'been calculated but you still wish to run this program, '+\
				  'please set dose inputs to NOTCALCULATED within input file.\n'

		if self.dose1 != 'NOTCALCULATED':
			err = '"dose1" input must be a positive float. '+\
				  'Currently set as "{}" in input file.\n{}'.format(self.dose1,doseStr)
			try:
				float(self.dose1)
			except ValueError:
				error(text = err)
				return False
			if float(self.dose1) < 0:
				error(text = err)
				return False

		if self.dose2 != 'NOTCALCULATED':
			if self.multiDatasets is False:
				err = '"dose2" input must be a positive float. '+\
				      'Currently set as "{}" in input file.\n{}'.format(self.dose2,doseStr)				
				try:
					float(self.dose2)
				except ValueError:
					error(text = err)
					return False
				if float(self.dose2) < 0:
					error(text = err)
					return False

			else:
				for dose in self.dose2.split(','):
					err = 'All doses in "dose2" input must be positive floats. '+\
				          'Currently set as "{}" in input file.\n{}'.format(self.dose2,doseStr)
					try:
						float(dose)
					except ValueError:
						error(text = err)
						return False
					if float(dose) < 0:
						error(text = err)
						return False

		print 'All input file parameters appear to be of suitable format.'
		return True

	def checkMtzLabelsExist(self):

		# for each mtz file input, check that the correct labels as 
		# specified within the txt input file are successfully found

		# initial dataset mtz files checked here
		if self.repeatedFile1InputsUsed is True:
			mtzFiles  = [self.mtz1]
			mtzLabels = [self.mtzlabels1]
		else:
			mtzFiles  = self.mtz1.split(',')
			mtzLabels = self.mtzlabels1.split(',')

		for i,(f, lab) in enumerate(zip(mtzFiles, mtzLabels)):
			foundLabels = self.getLabelsFromMtz(fileName = f)
			labels = ['F'+lab,
					  'SIGF'+lab]
			for lab in labels:
				if lab not in foundLabels:
					self.mtzLabelNotFound(mtzFile = f,
										  label   = lab)
					return False

		# later dataset mtz files checked here
		if self.multiDatasets is False:
			mtzFiles  = [self.mtz2]
			mtzLabels = [self.mtzlabels2]
		else:
			mtzFiles  = self.mtz2.split(',')
			mtzLabels = self.mtzlabels2.split(',')

		for i,(f, lab) in enumerate(zip(mtzFiles, mtzLabels)):
			foundLabels = self.getLabelsFromMtz(fileName = f)
			labels = ['F'+lab,
					  'SIGF'+lab]
			for lab in labels:
				if lab not in foundLabels:
					self.mtzLabelNotFound(mtzFile = f,
										  label   = lab)
					return False

		# phase information dataset mtz file checked here
		foundLabels = self.getLabelsFromMtz(fileName = self.mtz3)
		labels = [self.phaseLabel,
				  self.FcalcLabel]
		for lab in labels:
			if lab not in foundLabels:
				self.mtzLabelNotFound(mtzFile = f,
									  label   = lab)
				return False

		return True

	def mtzLabelNotFound(self,
						 mtzFile = 'untitled.mtz',
						 label   = 'unspecified'):

		# if an mtz label has not been found, print an error message
		err = 'Column "{}" not found in file "{}"'.format(label,mtzFile.split('/')[-1])+\
			  'Please check the .txt input file used for the current job'+\
		      'Also check the mtz file to ensure the column name "{}" exists'.format(label)
		error(text = err)

	def getLabelsFromMtz(self,
						 fileName = 'untitled.mtz'):

		# get list of column names from an mtz file

		# write commands for mtzdump to run
		inFile = open('mtzdumpInput.txt','w')
		inFile.write('header\nend')
		inFile.close()

		# run mtzdump and parse output to find column labels
		os.system('mtzdump hklin {} < mtzdumpInput.txt > mtzdump.tmp'.format(fileName))
		output = open('mtzdump.tmp','r')
		labelsFound = False
		for l in output.readlines():
			if l.replace('\n','') == '':
				continue
			if '* Column Types :' in l:
				break
			if labelsFound is True:
				labels = l.split()
			if '* Column Labels :' in l:
				labelsFound = True
				continue
		output.close()
		os.remove('mtzdump.tmp')

		return labels

	def checkFileExists(self, 
					    fileName = ''):

		# check that file exists and return bool 

		if os.path.isfile(fileName) and os.access(fileName,os.R_OK):
			return True
		else:
			err = 'File "{}" could not be located..'.format(fileName)
			error(text = err)
			return False

	def checkForMultipleDatasets(self):

		# multiple processing jobs can be run sequentially if input 
		# file contains comma separated lists for particular files. 
		# Determine whether this is the case and check that correctly 
		# formatted.

		props = ['name2',
				 'mtz2',
				 'mtzlabels2',
				 'pdb2']

		if self.dose2 != 'NOTCALCULATED':
			props.append('dose2')

		lengths = []
		for p in props:
			val = getattr(self,p)
			lengths.append(len(val.split(',')))

		self.numDsets  = lengths[0]
		if lengths == [1]*4:
			print 'Only single dataset located and will be processed'
			self.multiDatasets = False
		elif lengths[1:] != lengths[:-1]:
			err = 'Error! Input file properties ({}) '.format(','.join(props))+\
				  'must have same number of comma-separated inputs'
			error(text = err)

		else:
			print 'Multiple datasets located in input file to be processed'
			props = ['name1',
					 'mtz1',
					 'mtzlabels1',
					 'pdb1',
					 'RfreeFlag1']
			lengths = []
			for p in props:
				val = getattr(self,p)
				lengths.append(len(val.split(',')))

			if lengths not in ([1]*5,[self.numDsets]*5):
				err = 'Error! Input file properties ({})'.format(','.join(props))+\
					  'must all each have either 1 or 5 comma-separated inputs'
				error(text = err)
			else:
				s = False
				if lengths == [1]*5:
					s = True
				self.repeatedFile1InputsUsed = s
				self.multiDatasets = True

	def checkOutputDirExists(self,
							 makeProcessDir = True):

		# check whether output directory exists and make if not

		if self.dir.endswith('/') is False:
			print 'Working directory specified in input '+\
				  'file must end in "/" - appending'
			self.dir += '/'

		if not os.path.exists(self.dir):
			print 'Output directory "{}" not found, making directory'.format(self.dir)
			os.makedirs(self.dir)
		print 'Working directory set to "{}"'.format(self.dir)

		self.mapProcessDir = self.dir + 'ETRACK-mapProcessing/'
		if makeProcessDir is True:
			if 'ETRACK-mapProcessing' not in os.listdir(self.dir):
				os.makedirs(self.mapProcessDir)

	def getCurrentInputParams(self, 
							  jobNumber = 0):

		# select correct input parameters if multiple jobs within 
		# main input file. used when writing separate input files 
		# for each part of pipeline.

		props1 = ['name2',
				  'mtz2',
				  'mtzlabels2',
				  'pdb2']

		props2 = ['name1',
				  'mtz1',
				  'mtzlabels1',
				  'pdb1',
				  'RfreeFlag1']

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

		# create dataset name here, used to distinguish 
		# input file naming scheme (if name2 not the same 
		# as pdb2 then this is required)

		self.dsetName = str((self.pdb2_current).split('/')[-1]).strip('.pdb')

	def writePipeline1Inputs(self):

		# write the input file for the first subroutine
		# (run of CAD and SCALEIT)

		self.pipe1FileName 	= '{}{}_cadscaleit.txt'.format(self.mapProcessDir,
														   self.inputFile.split('.')[0])

		props = {'mtzIn1'          : 'mtz1',
				 'Mtz1LabelName'   : 'mtzlabels1_current',
				 'Mtz1LabelRename' : 'name1_current',
				 'RfreeFlag1'      : 'RfreeFlag1_current',
				 'mtzIn2'          : 'mtz2_current',
				 'Mtz2LabelName'   : 'mtzlabels2_current',
				 'Mtz2LabelRename' : 'name2_current',
				 'mtzIn3'          : 'mtz3',
				 'Mtz3phaseLabel'  : 'phaseLabel',
				 'Mtz3FcalcLabel'  : 'FcalcLabel',
				 'Mtz3LabelRename' : 'name3',
				 'inputPDBfile'    : 'pdb2_current',
				 'densMapType'     : 'densMapType',
				 'deleteMtzs'      : 'deleteIntermediateFiles',
				 'FFTmapWeight'    : 'FFTmapWeight'}

		fileOut1 = open(self.pipe1FileName,'w')
		inputString = ''
		for p in props.keys():
			inputString += '{} {}\n'.format(p,getattr(self,props[p]))
		inputString += 'END'
		fileOut1.write(inputString)
		fileOut1.close()

	def setJobName(self):

		# set a name for the current job

		self.jobName = '{}-{}'.format(self.name2_current,self.name1_current)

	def writePipeline2Inputs(self):

		# write input file for the second subroutine 
		# (run of SFALL and FFT)

		self.pipe2FileName 	= '{}{}_sfallfft.txt'.format(self.mapProcessDir,
														 self.inputFile.split('.')[0])

		if self.densMapType != '2FOFC':
			self.mtzIn = '{}{}_SCALEITcombined.mtz'.format(self.mapProcessDir,self.jobName)
		else:
			self.mtzIn = '{}{}_sigmaa.mtz'.format(self.mapProcessDir,self.name2_current)

		props = {'pdbIN'        : 'pdb2_current',
				 'initialPDB'   : 'name1_current',
				 'laterPDB'     : 'name2_current',
				 'phaseDataset' : 'name3',
				 'sfall_VDWR'   : 'sfall_VDWR',
				 'mtzIN'        : 'mtzIn',
				 'densMapType'  : 'densMapType',
				 'FFTmapWeight' : 'FFTmapWeight'}

		fileOut2 = open(self.pipe2FileName,'w')
		inputString = ''
		for p in props.keys():
			inputString += '{} {}\n'.format(p,getattr(self,props[p]))
		inputString += 'END'
		fileOut2.write(inputString)
		fileOut2.close()

	def startLogFile(self):

		# create log file for full map processing pipeline

		fName = '{}{}_run.log'.format(self.mapProcessDir,self.jobName)
		log = logFile(fileName = fName,
					  fileDir  = self.mapProcessDir)
		self.logFile = log

	def runPipeline1(self,
				     logFile = ''):

		# run the first subroutine (CAD and SCALEIT run)

		self.p1 = pipe1(outputDir = self.mapProcessDir,
						inputFile = self.pipe1FileName,
						jobName   = self.jobName,
						log       = logFile)

		outcome = self.p1.runPipeline()
		if outcome == 0:
			print 'Pipeline ran to completion'
			return True
		else:
			err = 'Pipeline failed to run to completion'
			error(text = err)
			return False

	def runPipeline2(self,
					 logFile = ''):

		# run the second subroutine (SFALL and FFT etc)

		self.p2 = pipe2(outputDir = self.mapProcessDir,
						inputFile = self.pipe2FileName,
						jobName   = self.jobName,
						log       = logFile)

		outcome = self.p2.runPipeline()
		if outcome == 0:
			print 'Pipeline ran to completion'
			return True
		else:
			print 'Pipeline failed to run to completion'
			return False

	def findFilesInDir(self,
					   mapProcessDir = True):

		# find files initially in working directory

		if mapProcessDir is True:
			self.filesInDir = os.listdir(self.mapProcessDir)
		else:
			self.filesInDir = os.listdir(self.dir)

	def cleanUpIntermediateFiles(self,
				   				 removeMtzs   = True,
				   				 removeMaps   = True,
				   				 removePdbs   = True,
				   				 includeFCmap = True):

		# after successful completion of the map processing 
		# part of the pipeline for each SINGLE dataset
		# clean up working directory

		print '\nCleaning up working directory'

		# distinguish between FFT and END map output formats 
		# depending on program used (FFT/END shellscript)
		if self.densMapType == 'END':
			densMapProg = 'END_switchedAxes'
		else:
			densMapProg = 'fft'

		params 		  = [self.mapProcessDir,self.dsetName]
		renameParams  = [self.mapProcessDir,self.name2_current]
		renameParams2 = [self.mapProcessDir,self.name1_current]

		keyLogFiles  = [self.p1.runLog.logFile,self.p2.runLog.logFile]

		if self.densMapType != 'END':
			mapFiles = ['{}{}-{}-{}_cropped_cropped.map'.format(self.mapProcessDir,
																self.dsetName,
																self.densMapType,
																densMapProg)]
		else:
			mapFiles = ['{}{}-{}_cropped_cropped.map'.format(self.mapProcessDir,
															 self.dsetName,
															 densMapProg)]
		mapFiles += ['{}{}_sfall_cropped.map'.format(*params)]

		renameMaps = ['{}{}_density.map'.format(*renameParams),
					   '{}{}_atoms.map'.format(*renameParams)]

		if includeFCmap is True:
			mapFiles += ['{}{}-FC-{}_cropped_cropped.map'.format(self.mapProcessDir,
																 self.dsetName,
																 densMapProg)]
			renameMaps += ['{}{}_FC.map'.format(*renameParams2)]

		pdbFiles  = ['{}{}_reordered.pdb'.format(*params)]
		renamePDB = ['{}{}.pdb'.format(*renameParams)]

		outputFiles = mapFiles + pdbFiles
		renameFiles = renameMaps + renamePDB

		subdir = '{}{}_additionalFiles/'.format(self.mapProcessDir,
												self.jobName)
		self.makeOutputDir(dirName = subdir)

		keyFiles = keyLogFiles + mapFiles + pdbFiles
		for file in os.listdir(self.mapProcessDir): 
			if file.endswith('_additionalFiles'):
				continue
			if file not in self.filesInDir:
				fileName = '{}{}'.format(self.mapProcessDir,file)
				if fileName not in keyFiles:
					shutil.move(fileName,'{}{}'.format(subdir,file))
		for file in os.listdir(subdir):
			remove = False
			if removeMtzs is True and file.endswith('.mtz'):
				remove = True
			if removeMaps is True and file.endswith('.map'):
				remove = True
			if removePdbs is True and file.endswith('.pdb'):
				remove = True
			if remove is True:
				os.remove(subdir+file)

		shutil.make_archive(subdir, 'zip', subdir)
		shutil.rmtree(subdir)

		# rename final map & pdb files
		for i in range(len(outputFiles)): 
			shutil.move(outputFiles[i],renameFiles[i])

		# move initial dataset pdb files to working directory
		self.moveInitialPDBfile()

		# check that resulting files are found
		self.findFilesInDir()
		for f in renameFiles:
			if f.split('/')[-1] not in self.filesInDir:
				err = 'Not all key output files found'
				error(text = err)
				return False
		return True

	def moveInitialPDBfile(self):

		# copy the initial dataset pdb file to the working 
		# directory useful if further ETRACK jobs to be run)

		if self.multiDatasets == False or self.repeatedFile1InputsUsed == True:
			shutil.copy2(self.pdb1,
						 '{}{}.pdb'.format(self.mapProcessDir,
						 			   self.name1))
		else:
			for pdb in self.pdb1.split(','):
				shutil.copy2(pdb,'{}{}.pdb'.format(self.mapProcessDir,self.name1))

	def writeETRACKInputFile(self,
							 write        = True,
							 run          = True,
							 imports      = True,
							 skipToETRACK = False):

		# writes an input file for the run of ETRACK after this processing 
		# pipeline has completed. Allows the user to immediately run the 
		# rest of the program, IF the user did specify all required datasets 
		# within a damage series. If the user only processed a subset of 
		# required datasets within a damage series, another ETRACK input 
		# file will need to be manually created to run the rest of the 
		# program with ALL datasets within the damage series included.

		if imports is True:
			from runETRACK import run as runETRACK

		if skipToETRACK is False:
			if write is False:
				return True

			r = runETRACK(runAll = False)
			seriesName = self.dir.split('/')[-2]

			if self.multiDatasets is True and self.repeatedFile1InputsUsed is False:
				err =  'Must have single INITIALDATASET inputs in input '+\
					   'file to run ETRACK immediately here'
				error(text = err)
				return False

			else:
				if self.dose2 == 'NOTCALCULATED':
					doses = ','.join(map(str,range(len(self.name2.split(',')))))
				else:
					doses = self.dose2

			r.writeInputFile(inDir         = self.mapProcessDir,
							 outDir        = self.dir,
							 damSetName    = seriesName,
							 laterDatasets = self.name2,
							 initialPDB    = self.name1,
							 doses         = doses,
							 outputGraphs  = self.outputGraphs)

			shutil.move(r.inputFileName,
						'{}/{}'.format(self.dir,r.inputFileName))

		if run is True:
			r = runETRACK(inputFileLoc = self.dir)

		return True

	def makeOutputDir(self,
					  dirName = './'):

		# if the above sub directory does not exist, make it

		if not os.path.exists(dirName):
			os.makedirs(dirName)
			print 'New sub directory "{}"'.format(dirName)+\
				  'created to contain output files'




