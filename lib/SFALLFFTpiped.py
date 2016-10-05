from MAPMASKjob import MAPMASKjob
from PDBCURjob import PDBCURjob
from SFALLjob import SFALLjob
from mapTools import mapTools
from logFile import logFile
from FFTjob import FFTjob
from ENDjob import ENDjob
import shutil
import struct
import sys
import os

class pipeline():

	# class to run SFALL job to generate atom-tagged map, and FFT to 
	# generate a corresponding density map (typically a Fobs(n)-Fobs(1)
	# Fourier difference map). The maps are then both restricted to the 
	# crystal asymmetric unit, with identical grid sampling properties 
	# (i.e. they are 'compatible'). This subroutine runs after the 
	# CAD-SCALEIT subroutine (see processFiles.py for how these are called).


	def __init__(self,
				 outputDir = '',
				 inputFile = '',
				 jobName   = 'untitled-job',
				 log       = ''):

		self.outputDir	= outputDir
		self.inputFile 	= inputFile
		self.jobName 	= jobName
		self.findFilesInDir() # find files initially in working dir

		# create log file
		if log == '':
			self.runLog = logFile(fileName = '{}{}_runLog1.log'.format(self.outputDir,jobName),
								  fileDir  = self.outputDir)
		else:
			self.runLog = log

	def runPipeline(self):

		# run the current subroutine within this class

		# read input file first
		success = self.readInputFile()
		if success is False:
			return 1

		# run pdbcur job 
		self.printStepNumber()
		pdbcur = PDBCURjob(inputPDBfile = self.pdbcurPDBinputFile,
						   outputDir    = self.outputDir,
						   runLog       = self.runLog)
		success = pdbcur.run()

		if success is False:
			return 2

		self.PDBCURoutputFile = pdbcur.outputPDBfile

		# reorder atoms in PDB file
		self.renumberPDBFile()

		# get space group from PDB file
		success = self.getSpaceGroup()
		if success is False:
			return 3

		# run SFALL job
		self.printStepNumber()
		sfall = SFALLjob(inputPDBfile   = self.reorderedPDBFile,
						 outputDir      = self.outputDir,
						 VDWR           = self.sfall_VDWR,
						 symmetrygroup  = self.spaceGroup,
						 gridDimensions = self.sfall_GRID,
						 runLog         = self.runLog)
		success = sfall.run()
		if success is False:
			return 4

		# run FFT job
		sfallMap  = mapTools(mapName = sfall.outputMapFile)
		axes 	  = [sfallMap.fastaxis,
					 sfallMap.medaxis,
					 sfallMap.slowaxis]

		gridSamps = [sfallMap.gridsamp1,
					 sfallMap.gridsamp2,
					 sfallMap.gridsamp3]

		if self.densMapType in ('DIFF','SIMPLE'):
			tags = ['FP_',
					'SIGFP_',
					'FOM_']
			labelsInit 	= [i+self.initPDB for i in tags] + ['PHIC_'+self.phaseDataset]
			labelsLater = [i+self.laterPDB for i in tags] + ['PHIC_'+self.phaseDataset]

		if self.densMapType == '2FOFC':
			labelsInit 	= ['']*4
			labelsLater = ['FWT_{}'.format(self.laterPDB),'','','PHIC']
		
		if self.densMapType != 'END':
			self.printStepNumber()
			fft = FFTjob(mapType   = self.densMapType,
						 FOMweight = self.FOMweight,
						 pdbFile   = self.reorderedPDBFile,
						 mtzFile   = self.inputMtzFile,
					     outputDir = self.outputDir,
					     axes      = axes,
					     gridSamps = gridSamps,
					     labels1   = labelsLater,
					     labels2   = labelsInit,
					     runLog    = self.runLog)
			success = fft.run()
		else:
			# run END job if required (may take time to run!!)
			self.printStepNumber()
			endInputPDB = self.pdbcurPDBinputFile
			endInputMTZ = ''.join(endInputPDB.split('.')[:-1]+['.mtz'])
			endInputEFF = ''.join(endInputPDB.split('.')[:-1]+['.eff'])
			
			end = ENDjob(pdbFile   = endInputPDB,
						 mtzFile   = endInputMTZ,
						 effFile   = endInputEFF,
						 outputDir = self.outputDir,
						 gridSamps = gridSamps,
						 runLog    = self.runLog)
			success = end.run()

		if success is False:
			return 5

		# generate FC map using FFT
		if self.includeFCmaps() is True:
			self.printStepNumber()

			fft_FC = FFTjob(mapType   = 'FC',
						    FOMweight = self.FOMweight,
						    pdbFile   = self.reorderedPDBFile,
						    mtzFile   = self.inputMtzFile,
					        outputDir = self.outputDir,
					        axes      = axes,
					        gridSamps = gridSamps,
					        labels1   = ['FC_{}'.format(self.phaseDataset),'','','PHIC_'+self.phaseDataset],
					        runLog    = self.runLog)
			success = fft_FC.run()

			if success is False: 
				return 6

		# crop atom-tagged map to asymmetric unit:
		self.printStepNumber()
		mapmask1 = MAPMASKjob(mapFile1  = sfall.outputMapFile,
							  outputDir = self.outputDir,
							  runLog    = self.runLog)
		success = mapmask1.crop2AsymUnit()

		if success is False:
			return 7

		# choose correct density map to include in MAPMASK cropping below
		if self.densMapType != 'END':
			inputDensMap = fft.outputMapFile
		else: 
			inputDensMap = end.outputMapFile

		# switch END map axes to match SFALL atom-tagged map if required
		if self.densMapType == 'END':
			mapmask_END = MAPMASKjob(mapFile1  = inputDensMap,
									 outputDir = self.outputDir,
									 runLog    = self.runLog)
			success = mapmask_END.switchAxisOrder(order = axes,
												  symGroup = self.spaceGroup)
			
			if success is False:
				return 8.0
			else: 
				inputDensMap = mapmask_END.outputMapFile

		# crop density map to atom-tagged map over asym unit:
		newMap = self.cropDensmapToSFALLmap(mapType = self.densMapType,
								            densMap = inputDensMap,
								            atmMap  = mapmask1.outputMapFile)

		# crop FC map to atom-tagged map over asym unit:
		if self.includeFCmaps() is True:
			newMap = self.cropDensmapToSFALLmap(mapType = 'FC',
									   			densMap = fft_FC.outputMapFile,
									   			atmMap  = mapmask1.outputMapFile)

		# print the map info (grid size, num voxels)
		ln =  'All generated maps should now be restricted to asym unit '+\
			  'with properties: '
		self.runLog.writeToLog(str = ln)
		Map = mapTools(mapName = mapmask1.outputMapFile,
					   logFile = self.runLog)
		Map.printMapInfo()

		# perform map consistency check between cropped fft and sfall maps
		fftMap = mapTools(newMap)
		fftMap.readHeader()
		sfallMap = mapTools(mapmask1.outputMapFile)
		sfallMap.readHeader()
		success = self.mapConsistencyCheck(sfallMap,fftMap)
		if success is False:
			return 10
		else:
			self.cleanUpDir()
			return 0

	def readInputFile(self):

		# read in input file for this subroutine

		# if Input.txt not found, flag error
		if os.path.isfile(self.inputFile) is False:
			err = 'Required input file {} not found..'.format(self.inputFile)
			self.runLog.writeToLog(str = err)
			return False

		inputFile = open(self.inputFile,'r')
		props = {'pdbIN'           : 'pdbcurPDBinputFile',
				 'runname'         : 'runName',
				 'sfall_VDWR'      : 'sfall_VDWR',
				 'mtzIN'           : 'inputMtzFile',
				 'foldername'      : 'outputDir',
				 'initialPDB'      : 'initPDB',
				 'laterPDB'        : 'laterPDB',
				 'phaseDataset'    : 'phaseDataset',
				 'densMapType'     : 'densMapType',
				 'FFTmapWeight'    : 'FOMweight',
				 'calculateFCmaps' : 'FCmaps'}

		self.sfall_GRID = []
		for l in inputFile.readlines():
			if 'END' == l.split()[0]:
				break
			elif l.split()[0] in props.keys():
				setattr(self,props[l.split()[0]],l.split()[1])
			elif l.split()[0] == 'sfall_GRID':
				self.sfall_GRID = l.split()[1:4]
		return True

	def includeFCmaps(self):

		# determine whether FCalc maps should be generated

		if self.FCmaps.upper() == 'TRUE':
			return True
		else:
			return False

	def renumberPDBFile(self):

		# reorder atoms in pdb file since some may be 
		# missing now after pdbcur has been run

		ln = 'Renumbering input pdb file: {}'.format(self.PDBCURoutputFile)
		self.runLog.writeToLog(str = ln)

		self.reorderedPDBFile = (self.PDBCURoutputFile).split('_pdbcur.pdb')[0]+'_reordered.pdb'

		pdbin = open(self.PDBCURoutputFile,'r')
		pdbout = open(self.reorderedPDBFile,'w')

		counter = 0
		for line in pdbin.readlines():
			if ('ATOM' not in line[0:4] and 'HETATM' not in line[0:6]):
				pdbout.write(line)
			else:
				counter += 1
				pdbout.write(line[0:6])
				new_atomnum = " "*(5-len(str(counter))) + str(counter) #has length 5
				pdbout.write(new_atomnum)
				pdbout.write(line[11:80]+'\n')
		pdbin.close()
		pdbout.close()
		self.runLog.writeToLog(str = 'Output pdb file: {}'.format(self.reorderedPDBFile))

	def getSpaceGroup(self):

		# parse the space group from the input pdb file

		pdbin = open(self.reorderedPDBFile,'r')
		for line in pdbin.readlines():
			if line.split()[0] == 'CRYST1':
				self.spaceGroup = line[55:66].replace(' ','')
				ln = 'Retrieving space group from file: {}.\n'.format(self.PDBCURoutputFile)+\
					 'Space group determined to be {}'.format(self.spaceGroup)
				self.runLog.writeToLog(str = ln)

		try:
			self.spaceGroup
		except attributeError:
			err = 'Unable to find space group from file: {}'.format(self.PDBCURoutputFile)
			self.runLog.writeToLog(str = err)

			return False
		return True

	def cropDensmapToSFALLmap(self,
							  mapType = 'DIFF',
							  densMap = '',
							  atmMap  = ''):

		# crop the density map to exact same dimensions
		# as SFALL atom-tagged map
		
		# run MAPMASK job to crop fft density map to asym unit
		mapmask2 = MAPMASKjob(mapFile1  = densMap,
							  outputDir = self.outputDir,
							  runLog    = self.runLog)
		success = mapmask2.crop2AsymUnit()

		if success is False:
			return 8.1

		# run MAPMASK job to crop fft density map to same grid 
		# sampling dimensions as SFALL atom map
		mapmask3 = MAPMASKjob(mapFile1  = mapmask2.outputMapFile,
							  mapFile2  = atmMap,
							  outputDir = self.outputDir,
							  runLog    = self.runLog)
		success = mapmask3.cropMap2Map()

		if success is False:
			return 9
		else:
			return mapmask3.outputMapFile

	def mapConsistencyCheck(self,sfallMap,fftMap):

		# this function determines whether the atom map and density 
		# map calculated using SFALL and FFT are compatible - meaning 
		# the grid dimensions/filtering are the same and the ordering
		# of the fast, medium, and slow axes are identical.

		err = 'Checking that atom map (SFALL) and density map (FFT) are compatible...'
		self.runLog.writeToLog(str = err)

		if (sfallMap.gridsamp1 != fftMap.gridsamp1 or
			sfallMap.gridsamp2 != fftMap.gridsamp2 or
			sfallMap.gridsamp3 != fftMap.gridsamp3):
			err = 'Incompatible grid sampling found...'
			self.runLog.writeToLog(str = err)

			return False

		if (sfallMap.fastaxis != fftMap.fastaxis or
			sfallMap.medaxis != fftMap.medaxis or
			sfallMap.slowaxis != fftMap.slowaxis):
			err = 'Incompatible fast,med,slow axes ordering found...'
			self.runLog.writeToLog(str = err)
			return False

		if (sfallMap.numCols != fftMap.numCols or
			sfallMap.numRows != fftMap.numRows or
			sfallMap.numSecs != fftMap.numSecs):
			err = 'Incompatible number of rows, columns and sections...'
			self.runLog.writeToLog(str = err)
			return False

		if sfallMap.getMapSize() != fftMap.getMapSize():
			err = 'Incompatible map file sizes'
			self.runLog.writeToLog(str = err)
			return False

		self.runLog.writeToLog(str = '---> success!')
		return True

	def cleanUpDir(self):

		# give option to clean up working directory 

		ln = '\nCleaning up working directory...'
		self.runLog.writeToLog(str = ln)

		# move txt files to subdir
		self.makeOutputDir(dirName = '{}txtFiles/'.format(self.outputDir))
		
		for file in os.listdir(self.outputDir): 
			if file.endswith('.txt') and file not in self.filesInDir:
				args = [self.outputDir,file]
				shutil.move('{}{}'.format(*args),'{}txtFiles/{}'.format(*args))

	def findFilesInDir(self):

		# find files initially in working directory

		self.filesInDir = os.listdir(self.outputDir)

	def makeOutputDir(self,
					  dirName = './'):

		# if the above sub directory does not exist, make it

		if not os.path.exists(dirName):
			os.makedirs(dirName)
			ln = 'New sub directory "{}" created to contain output files.'.format(dirName)
			self.runLog.writeToLog(str = ln)

	def printStepNumber(self):

		# print a string indicating the current pipeline 
		# step number directory to the command line

		try:
			self.stepNumber
		except AttributeError:
			self.stepNumber = 3
		txt = '\n_______'+\
			  '\nSTEP {})'.format(self.stepNumber)
		self.runLog.writeToLog(str = txt)

		self.stepNumber += 1



