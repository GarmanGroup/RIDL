print 'Importing packages...'
import time
t0 = time.time()
from combinedAtomList import combinedAtomList
from savevariables import retrieve_objectlist,save_objectlist,saveGenericObject,retrieveGenericObject
from PDBFileManipulation import PDBtoList,writePDBline
from readAtomMap import maps2DensMetrics
from time import gmtime, strftime
from numpy import isnan
from shutil import move
from os import path,makedirs,listdir,system,remove,rmdir
from ridlFeedback import provideFeedback
from furtherOutput import furtherAnalysis
print 'Total import time: {}s'.format(round(time.time()-t0,3))

class calculateMetrics(object):

	# a class for retrieving the RIDL input text file information and running 
	# the RIDL pipeline to calculate per-atom damage metrics for a specified 
	# damage series. This code requires previously generated atom-tagged and
	# density maps (typically Fobs(n) - Fobs(1) Fourier difference maps) to 
	# have been generated for the current damage series (as specified within 
	# the input .txt file parsed below). If run as part of the full RIDL 
	# pipeline (by running 'python runRIDL.py -i [inputfile.txt] -pc') then 
	# this will automatically run directly after the suitable map files have
	# been generated, with no need to explicitly write a new input file for
	# this class to work.

	def __init__(self,
				 inDir      = './',
				 outDir     = './',
				 pdbNames   = [],
				 pklFiles   = [],
				 initialPDB = "",
				 seriesName = "untitled-series",
				 pklSeries  = "",
				 doses      = [],
				 plot       = 'no',
				 output     = 'simple',
				 logFile    = '',
				 sumFiles   = True,
				 inclFCmets = True):

		self.inDir 		= inDir 		# the input file directory
		self.outDir  	= outDir 		# the output file directory
		self.pdbNames 	= pdbNames 		# the list of pdb codes for series
		self.pklFiles 	= pklFiles 		# list of pkl files from map_processing
		self.initialPDB = initialPDB 	# the first dataset pdb code
		self.seriesName = seriesName 	# the general series name 
		self.pklSeries	= pklSeries 	# the combined series pkl file from post_processing
		self.doses 		= doses 		# list of increasing doses
		self.plot       = plot 			# (bool) decide to plot per-residue summary plots per dataset
		self.output     = output        # the amount of output to provide (either 'simple' for just Dloss
										# info or 'full' larger selection of output files)
		self.logFile    = logFile       # log file for the current RIDL job
		self.sumFiles   = sumFiles      # write summary html files and graphs (bool)
		self.inclFCmets = inclFCmets    # generate metrics which require FC maps to be made
		
	def runPipeline(self,
					map_process   = True,
					post_process  = True,
					retrieve      = True,
					inputFileName = ''):

		# the following function reads in the above functions one by one in a scripted 
		# pipeline. Takes inputs to specify which parts of above pipeline are to be
		# included. For each function input specify True if this part is to be performed
		# and False otherwise.

		self.inputFileName = inputFileName

		# check whether valid inputs to function
		valid = self.checkValidInputs(map_process,post_process,retrieve)
		if not valid: 
			return

		# first need to run function above to read in input file containing info
		# on where input files are and where output files should be written
		ln = 'Reading input file: {}'.format(self.inputFileName)
		self.logFile.writeToLog(str = ln)

		success = self.readInputFile()
		if not success: 
			return

		success = self.checkInOutDirExist()
		if not success: 
			return

		self.setOutputDirs()

		if map_process:
			self.map_processing()
		else:
			ln = 'Map processing task not chosen...'
			self.logFile.writeToLog(str = ln)
		self.fillerLine()

		if post_process:
			self.post_processing()

			# save metric data to pkl file
			pklSeries = saveGenericObject(obj      = self.combinedAtoms,
										  fileName = self.seriesName)

			move(pklSeries,
				 '{}{}'.format(self.outputDir,pklSeries))
			self.pklSeries = pklSeries

			inputfile = open(self.inputFileName,'a')
			inputfile.write('\npklDataFile ' + pklSeries)
			inputfile.close()

			self.feedback(csvOnly = not self.sumFiles)

		else: 
			ln = 'Post processing job not chosen...'
			self.logFile.writeToLog(str = ln)

		if retrieve:
			self.PDBmulti_retrieve()
			self.feedback()

		self.fillerLine(blank = True)

	def checkValidInputs(self,map_process,post_process,retrieve):

		# check the runPipeline inputs to make sure that they are valid

		for var in [map_process,post_process,retrieve]:
			if var not in (True,False):
				print 'Only Boolian parameters expected here'
				return False
		try:
			f = open(self.inputFileName,'r')
			f.close()
		except IOError:
			err = 'metric calculation input file "{}" not found'.format(self.inputFileName)
			self.logFile.writeToLog(str = err)
			return False
		return True

	def readInputFile(self,
					  printText = True):

		# read input file e_Track_input.txt to specify location 
		# of input files and where to write output files

		props = {'inDir'          : 'inDir',
				 'outDir'         : 'outDir',
				 'damageset_name' : 'seriesName',
				 'initialPDB'     : 'initialPDB',
				 'pklDataFile'    : 'pklSeries',
				 'laterDatasets'  : 'laterDatasets'}

		inputfile = open(self.inputFileName,'r')
		pklFiles = []
		for line in inputfile.readlines():
			l = line.split()
			if '#' == line[0]: 
				continue
			elif l[0] in props.keys():
				setattr(self,props[l[0]],l[1])
			elif 'damageset_num' in l[0]:
				datasetNums = l[1]
			elif 'PKLFILE' == l[0]:
				self.pklFiles.append(l[1])
			elif 'doses' == l[0]:
				self.doses 	= [float(d) for d in l[1].split(',')]
			elif 'plot' == l[0]:
				self.plot         = True
				self.plotHeatMaps = True
			elif 'slim-plot' == l[0]:
				self.plot         = True
				self.plotHeatMaps = False
		inputfile.close()

		# if number of initial datasets given doesn't match
		# number of later datasets, assume same initial dataset
		# used for every later dataset (fix as first one given)
		initialPDBs   = self.initialPDB.split(',')
		numLaterDsets = len(self.laterDatasets.split(','))
		if len(initialPDBs) != numLaterDsets:
			initialPDBs = [initialPDBs[0]]*numLaterDsets
		l = []
		for pdb in initialPDBs:
			if not pdb.endswith('.pdb'):
				l.append(pdb+'.pdb')
			else:
				l.append(pdb)
		self.initialPDB = l

		# locate the correct format for the list of datasets within damage 
		# series. Currently two formats acceptable: (a) series-name + dataset-id
		# (per dataset), (b) input list of full dataset names (recommended).
		found = True
		try:
			datasetNums
		except UnboundLocalError:
			found = False
		if found:
			self.pdbNames = [self.seriesName+num for num in datasetNums.split(',')]
			return True
		found = True
		try:
			self.laterDatasets
		except AttributeError:
			found = False
		if found:
			self.pdbNames = self.laterDatasets.split(',')
			return True
		else:
			if printText:
				err = 'Error! Unable to extract list of dataset names from input file'
				self.logFile.writeToLog(str = err)
			return False

	def checkInOutDirExist(self):

		# check that an input/output directories have been 
		# found and make subdirectories if present

		for dir in ([[self.inDir,'Input'],[self.outDir,'Output']]):
			if not path.isdir(dir[0]):
				err = '{} file location: {} does not exist. '.format(dir[1],dir[0])+\
					  'Please select an appropriate directory'
				self.logFile.writeToLog(str = err)
				return False
		return True

	def makeOutputDir(self,
					  dirName   = './'):

		# if the above sub directory does not exist, make it

		if not path.exists(dirName):
			makedirs(dirName)
			ln = 'New sub directory "{}" created to contain output files'.format(dirName.replace(self.outputDir,''))
			self.logFile.writeToLog(str = ln)

	def setOutputDirs(self):

		# set the locations of the output directories

		self.outputDir 		= '{}RIDL-metrics/'.format(self.outDir)
		self.outputPlotDir 	= '{}plots/'.format(self.outputDir)

		# add pkl file names as attribute if specified in input file
		if len(self.pklFiles) != 0: 
			self.pklFiles = [self.outputDir+f for f in self.pklFiles]

	def map_processing(self):

		# combine the density map and atom-tagged map for a given dataset,
		# to calculate per-atom density metrics for each refined atom

		txt = 'Combining density maps and atom-tagged maps to calculate '+\
			  'per-atom density metrics for each refined atom in structure.\n'
		self.logFile.writeToLog(str = txt)

		txt = 'input directory:  {}\n'.format(self.inDir)+\
			  'output directory: {}'.format(self.outputDir)
		self.logFile.writeToLog(str   = txt,
								strip = False)

		# create additional subdirectories
		for oDir in (self.outputDir,self.outputPlotDir):
			self.makeOutputDir(dirName = oDir)

		# make pklFiles and dir to move all generated per-dataset pkl files to this
		pklFileDir = 'pklFiles-perDataset/'
		self.makeOutputDir(dirName = '{}{}'.format(self.outputDir,pklFileDir))

		txt = '{} higher dose datasets located within input file.'.format(len(self.pdbNames))+\
			  'Calculating per-atom density metrics for each dataset individually.\n'

		pklFileNames = []
		i = 0
		for d,initPDB in zip(self.pdbNames,self.initialPDB):
			i += 1
			# derive per-atom density metrics from maps
			mapName1 = '{}_atoms.map'.format(d)
			mapName2 = '{}_density.map'.format(d)
			mapName3 = initPDB.replace('.pdb','_FC.map')

			txt = '\n---------------------------------\n'+\
				  'Higher dose dataset {} starts here'.format(i)
			self.logFile.writeToLog(str = txt)

			maps2DensMets = maps2DensMetrics(filesIn     = self.inDir,
										     filesOut    = self.outputDir,
										     pdbName     = d,
										     atomTagMap  = mapName1,
										     densityMap  = mapName2,
										     FCmap 	     = mapName3,
										     plotHist    = self.plot,
										     logFile     = self.logFile,
										     calcFCmap   = self.inclFCmets)

   			maps2DensMets.maps2atmdensity()

			# move pkl file to working output directory
			pklFileName = maps2DensMets.pklFileName

			move(pklFileName,
				 '{}{}{}'.format(self.outputDir,
								 pklFileDir,
								 pklFileName))

			pklFileNames.append('{}{}{}'.format(self.outputDir,
												pklFileDir,
												pklFileName))

		self.pklFiles = pklFileNames

	def post_processing(self):

		# group the per-atom density metrics for each dataset together

		txt = 'Combining density metric information for each dataset '+\
			  'together within the damage series'
		self.logFile.writeToLog(str = txt)

		txt = 'Input pkl files for post processing chosen from input file:'
		for file in self.pklFiles: 
			txt += '\n\t{}'.format(file.replace(self.outDir,""))
		self.logFile.writeToLog(str = txt)

		# next read in the pdb structure file as list of atom objects
		initialPDBlist = PDBtoList(pdbFileName = self.get1stDsetPDB())

		# retrieve object lists of atoms for each damage set
		ln = '\nReading in pkl files for higher dataset structures...'
		self.logFile.writeToLog(str = ln)

		dList = [] # dataset list
		for pkl_filename in self.pklFiles:
			ln = 'Damage file number: {}'.format(len(dList)+1)
			self.logFile.writeToLog(str = ln)
			PDB_ret = retrieve_objectlist(fileName = pkl_filename,
										  logFile  = self.logFile)
			# remove pkl file since no longer needed
			remove(pkl_filename)

			# add new retrieved damage set list to dList
			dList.append(PDB_ret)
		pklDir = '/'.join(pkl_filename.split('/')[:-1])+'/'

		# remove directory if it is now empty
		if listdir(pklDir) == []:
			rmdir(pklDir)

		# create a list of atom objects with attributes as lists varying over 
		# dose range, only including atoms present in ALL damage datasets
		ln = 'New list of atoms over full dose range calculated...'
		self.logFile.writeToLog(str = ln)
		combinedAtoms = combinedAtomList(datasetList    = dList,
										 numLigRegDsets = len(dList),
										 doseList       = self.doses,
										 initialPDBList = initialPDBlist,
										 outputDir      = self.outputDir,
										 seriesName     = self.seriesName,
										 inclFCmetrics  = self.inclFCmets)

		combinedAtoms.getMultiDoseAtomList()

		# calculate 'average' variant Dloss metrics
		combinedAtoms.calcAdditionalMetrics(newMetric = 'average')

		# calculate Calpha normalised metrics, if Calpha atoms exist
		if self.checkCalphasPresent(atomObjList = combinedAtoms):
			for m in ('loss','mean','gain','Bfactor'):
				combinedAtoms.calcAdditionalMetrics(metric = m)
		
		self.combinedAtoms = combinedAtoms

	def PDBmulti_retrieve(self):

		# retrieve list of atom objects from .pkl file

		self.fillerLine(blank = True)

		txt = 'Retrieving per-atom damage metric information from '+\
			  ' .pkl file for each dataset within damage series.'
		self.logFile.writeToLog(str = txt)

		txt = 'Input pkl file for data retrieval chosen from input file:\n'+\
			  '\t{}'.format(self.pklSeries)
		self.logFile.writeToLog(str = txt)

		# retrieve the combinedAtoms object from the pkl file
		self.combinedAtoms = retrieveGenericObject(fileName = self.outputDir+self.pklSeries)

	def feedback(self,
				 csvOnly      = False,
				 includeTests = False):

		# write feedback files for current RIDL job.
		# if csvOnly is True then ONLY csv files will
		# be output from the run (i.e. no html summary
		# file and no plots)

		if not includeTests:
			provideFeedback(csvOnly       = csvOnly,
							atmsObjs      = self.combinedAtoms,
							logFile       = self.logFile,
							outputDir     = self.outputDir,
							outputPlotDir = self.outputPlotDir,
							pklSeries     = self.pklSeries,
							plotHeatMaps  = self.plotHeatMaps,
							doses         = self.doses,
							pdbNames      = self.pdbNames,
							inputDir      = self.inDir,
							initialPDB    = self.initialPDB)
		else:
			furtherAnalysis(csvOnly       = csvOnly,
							atmsObjs      = self.combinedAtoms,
							logFile       = self.logFile,
							outputDir     = self.outputDir,
							outputPlotDir = self.outputPlotDir,
							pklSeries     = self.pklSeries,
							plotHeatMaps  = self.plotHeatMaps,
							doses         = self.doses,
							pdbNames      = self.pdbNames,
							inputDir      = self.inDir,
							initialPDB    = self.initialPDB)

	def checkCalphasPresent(self,
							atomObjList = []):

		# check whether structure contains any Calpha 
		# protein backbone atoms within it

		return atomObjList.checkCalphaAtomsExist()

	def get1stDsetPDB(self):

		# retrieve name of first dataset pdb coordinate file.
		# If multiple initial datasets input, take the first only

		pdbFile = self.inDir + str(self.initialPDB[0])

		return pdbFile

	def fillerLine(self,
				   blank = False):

		# print a filler line to command line

		if not blank:
			ln = '\n---------------------------------------------------------------'	
		else:
			ln = '\n'
		self.logFile.writeToLog(str = ln)
