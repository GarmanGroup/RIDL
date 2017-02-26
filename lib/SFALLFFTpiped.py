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
				 outputDir      = '',
				 inputFile      = '',
				 jobName        = 'untitled-job',
				 log            = '',
				 useUnscaledMtz = False):

		self.outputDir	= outputDir
		self.inputFile 	= inputFile
		self.jobName 	= jobName
		self.findFilesInDir() # find files initially in working dir

		# create log file
		if log == '':
			f = '{}{}_runLog1.log'.format(self.outputDir+'RIDL-log',jobName)
			self.runLog = logFile(fileName = f,
								  fileDir  = self.outputDir)
		else:
			self.runLog = log

		# provide option to use non-scaleit scaled Fobs(n) for
		# density map calculations (retrieved directly from CAD instead)
		self.useUnscaledMtz = useUnscaledMtz

	def runPipeline(self):

		# run the current subroutine within this class

		success = self.readInputFile()
		if not success:
			return False

		if self.useUnscaledMtz:
			self.inputMtzFile = self.inputMtzFile.replace('_SCALEIT','_CAD')

		success = self.curatePdbFile()
		if not success:
			return False

		self.renumberPDBFile()

		success = self.getSpaceGroup()
		if not success:
			return False

		success = self.getAtomTaggedMap()
		if not success:
			return 4

		########################################################
		# TEST BIT  - DO NOT RUN FOR NORMAL USE. 
		# THIS WILL USE PHENIX TO GENERATE A FO-FO DIFF MTZ 
		# AND THEN CONVERT TO MAP USING FFT. MTZ FROM CAD IS 
		# USED AS INPUT TO PHENIX, SO NO SCALEIT SCALING DONE
		makeFoFoMapWithPhenix = False
		if makeFoFoMapWithPhenix:
			mtzInput = self.inputMtzFile.replace('_SCALEIT','_CAD')
			cmd = 	'phenix.fobs_minus_fobs_map '+\
					'f_obs_1_file_name={} '.format(mtzInput)+\
					'f_obs_2_file_name={} '.format(mtzInput)+\
					'phase_source={} '.format(self.reorderedPDBFile)+\
					'f_obs_1_label=FP_{} '.format(self.laterPDB)+\
					'f_obs_2_label=FP_{}'.format(self.initPDB)
			os.system(cmd+'>phenix.log')
			self.fcalcMtz   = self.inputMtzFile
			self.densMapMtz = self.inputMtzFile.replace('_SCALEITcombined','phenixFoFo')
			shutil.move('FoFoPHFc.mtz',self.densMapMtz)
			shutil.move('phenix.log','{}phenix.log'.format(self.outputDir))
			self.densMapType = 'FC'
			labelsLater  = ['FoFo','','','PHFc']
			labelsInit   = ['','','','']
			self.mapTag  = 'DIFF'
			success = self.generateDensMap(labelsInit,labelsLater)
		else:
			self.mapTag = ''
			self.fcalcMtz   = self.inputMtzFile
			self.densMapMtz = self.inputMtzFile
			success = self.generateDensMap()
        ########################################################

		if not success:
			return False

		if self.includeFCmaps():
			success = self.generateFcalcMap()
			if not success:
				return False

		success = self.cropAtmTaggedMapToAsymUnit()
		if not success:
			return False

		if self.densMapType == 'END':
			success = self.ensureSameMapAxesOrder()
			if not success:
				return False

		success = self.cropMapToAtomTaggedMap(densMap = self.densityMap)
		if not success:
			return False

		if self.includeFCmaps():
			success = self.cropMapToAtomTaggedMap(densMap = self.FcalcMap)
			if not success:
				return False

		self.reportCroppedMapInfo()
		success = self.mapConsistencyCheck()
		if not success:
			return False

		self.cleanUpDir()
		return True

	def curatePdbFile(self):

		# run pdbcur job

		self.printStepNumber()
		pdbcur = PDBCURjob(inputPDBfile = self.pdbcurPDBinputFile,
						   outputDir    = self.outputDir,
						   runLog       = self.runLog)
		success = pdbcur.run()
		self.PDBCURoutputFile = pdbcur.outputPDBfile
		return success

	def getAtomTaggedMap(self):

		# run SFALL job to generate atom tagged map

		self.printStepNumber()
		sfall = SFALLjob(inputPDBfile   = self.reorderedPDBFile,
						 outputDir      = self.outputDir,
						 VDWR           = self.sfall_VDWR,
						 symmetrygroup  = self.spaceGroup,
						 gridDimensions = self.sfall_GRID,
						 runLog         = self.runLog)
		success = sfall.run()

		sfallMap  = mapTools(mapName = sfall.outputMapFile)
		self.axes 	   = [sfallMap.fastaxis,
					 	  sfallMap.medaxis,
					      sfallMap.slowaxis]
		self.gridSamps = [sfallMap.gridsamp1,
					 	  sfallMap.gridsamp2,
					 	  sfallMap.gridsamp3]

		self.atomTaggedMap = sfall.outputMapFile

		return success

	def generateDensMap(self,
						labelsInit  = [],
						labelsLater = []):

		# run FFT job to generate density map

		self.printStepNumber()
		if self.densMapType in ('DIFF','SIMPLE'):
			tags = ['FP_',
					'SIGFP_',
					'FOM_']
			labelsInit 	= [i+self.initPDB for i in tags] + ['PHIC_'+self.phaseDataset]
			labelsLater = [i+self.laterPDB for i in tags] + ['PHIC_'+self.phaseDataset]
			# labelsLater = [i+self.laterPDB for i in tags] + ['PHIC_2'] # FOR specific test case only

		if self.densMapType == '2FOFC':
			labelsInit 	= ['']*4
			labelsLater = ['FWT_{}'.format(self.laterPDB),'','','PHIC']
		
		if self.densMapType != 'END':
			fft = FFTjob(mapType   = self.densMapType,
						 mapTag    = self.mapTag,
						 FOMweight = self.FOMweight,
						 pdbFile   = self.reorderedPDBFile,
						 mtzFile   = self.densMapMtz,
					     outputDir = self.outputDir,
					     axes      = self.axes,
					     gridSamps = self.gridSamps,
					     labels1   = labelsLater,
					     labels2   = labelsInit,
					     runLog    = self.runLog)
			success = fft.run()
			self.densityMap = fft.outputMapFile

		else:
			# run END job if required (may take time to run!!)
			endInputPDB = self.pdbcurPDBinputFile
			endInputMTZ = ''.join(endInputPDB.split('.')[:-1]+['.mtz'])
			endInputEFF = ''.join(endInputPDB.split('.')[:-1]+['.eff'])
			
			end = ENDjob(pdbFile   = endInputPDB,
						 mtzFile   = endInputMTZ,
						 effFile   = endInputEFF,
						 outputDir = self.outputDir,
						 gridSamps = self.gridSamps,
						 runLog    = self.runLog)
			success = end.run()
			self.densityMap = end.outputMapFile
		return success

	def generateFcalcMap(self):

		# generate FC map using FFT

		self.printStepNumber()
		fft_FC = FFTjob(mapType   = 'FC',
					    FOMweight = self.FOMweight,
					    pdbFile   = self.reorderedPDBFile,
					    mtzFile   = self.fcalcMtz,
				        outputDir = self.outputDir,
				        axes      = self.axes,
				        gridSamps = self.gridSamps,
				        labels1   = ['FC_{}'.format(self.phaseDataset),'','','PHIC_'+self.phaseDataset],
				        runLog    = self.runLog)
		success = fft_FC.run()
		self.FcalcMap = fft_FC.outputMapFile

		return success

	def cropAtmTaggedMapToAsymUnit(self):

		# crop atom-tagged map to asymmetric unit:

		self.printStepNumber()
		mapmask1 = MAPMASKjob(mapFile1  = self.atomTaggedMap,
							  outputDir = self.outputDir,
							  runLog    = self.runLog)
		success = mapmask1.crop2AsymUnit()
		self.atomTaggedMap = mapmask1.outputMapFile

		return success

	def ensureSameMapAxesOrder(self):

		# switch map axes to match SFALL atom-tagged map if 
		# required(only typically required for END maps)

			mapmask = MAPMASKjob(mapFile1  = self.densityMap,
								 outputDir = self.outputDir,
								 runLog    = self.runLog)
			success = mapmask.switchAxisOrder(order    = self.axes,
											  symGroup = self.spaceGroup)
			
			self.densityMap = mapmask.outputMapFile

			return success

	def readInputFile(self):

		# read in input file for this subroutine

		# if Input.txt not found, flag error
		if not os.path.isfile(self.inputFile):
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

	def cropMapToAtomTaggedMap(self,
							   densMap = 'untitled.map'):

		# crop the density map to exact same dimensions
		# as SFALL atom-tagged map
		
		# run MAPMASK job to crop fft density map to asym unit
		mapmask2 = MAPMASKjob(mapFile1  = densMap,
							  outputDir = self.outputDir,
							  runLog    = self.runLog)
		success = mapmask2.crop2AsymUnit()
		if not success:
			return False

		# run MAPMASK job to crop fft density map to same
		# grid sampling dimensions as SFALL atom map
		mapmask3 = MAPMASKjob(mapFile1  = mapmask2.outputMapFile,
							  mapFile2  = self.atomTaggedMap,
							  outputDir = self.outputDir,
							  runLog    = self.runLog)
		success = mapmask3.cropMap2Map()
		self.croppedDensityMap = mapmask3.outputMapFile

		return success

	def reportCroppedMapInfo(self):

		# print the map info (grid size, num voxels)

		ln =  'All generated maps should now be restricted to asym unit '+\
			  'with properties: '
		self.runLog.writeToLog(str = ln)
		Map = mapTools(mapName = self.atomTaggedMap,
					   logFile = self.runLog)
		Map.printMapInfo()

	def mapConsistencyCheck(self):

		# this function determines whether the atom map and density 
		# map calculated using SFALL and FFT are compatible - meaning 
		# the grid dimensions/filtering are the same and the ordering
		# of the fast, medium, and slow axes are identical.

		fftMap = mapTools(self.croppedDensityMap)
		fftMap.readHeader()
		sfallMap = mapTools(self.atomTaggedMap)
		sfallMap.readHeader()

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



