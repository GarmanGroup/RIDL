# -*- coding: utf-8 -*-
from savevariables import retrieve_objectlist,save_objectlist,saveGenericObject,retrieveGenericObject
from PDBFileManipulation import PDBtoList,writePDBline
from combinedAtomList import combinedAtomList
from readAtomMap import maps2DensMetrics
from time import gmtime, strftime
import numpy as np
import shutil
import os

print 'Checking whether seaborn plotting library present...'
try:
	import seaborn as sns
	seabornFound = True
except ImportError:
	print 'Plotting library seaborn not found..'
	print 'Will not create any plots for current run.'
	print 'Use "pip install seaborn" to install for future use'
	seabornFound = False

class eTrack(object):
	# A class for retrieving the eTrack input text file information and running 
	# the eTrack pipeline functions separately or together in full pipeline

	def __init__(self,
				 inDir      = './',
				 outDir     = './',
				 pdbNames   = [],
				 pklFiles   = [],
				 initialPDB = "",
				 seriesName = "untitled-series",
				 pklSeries  = "",
				 doses      = [],
				 plot       = False,
				 output     = 'simple'):

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
		
		if seabornFound is False:
			self.plot = False 

	def runPipeline(self,
					map_process   = True,
					post_process  = True,
					retrieve      = True,
					inputFileName = ''):

		# the following function reads in the above functions one by one in a scripted 
		# pipeline. Takes inputs to specify which parts of above pipeline are to be
		# included. For each function input specify True if this part is to be performed
		# and False otherwise.

		self.titleCaption(title='ETRACK pipeline')
		self.inputFileName = inputFileName

		# check whether valid inputs to function
		valid = self.checkValidInputs(map_process,post_process,retrieve)
		if valid is False: 
			return

		# first need to run function above to read in input file containing info
		# on where input files are and where output files should be written
		print 'Reading input file: {}'.format(self.inputFileName)
		success = self.readInputFile()
		if success is False: 
			return
		success = self.checkInOutDirExist()
		if success is False: 
			return
		self.setOutputDirs()

		if map_process is True:
			self.map_processing()
		else:
			print 'Map processing task not chosen...'
		self.fillerLine()

		if post_process is True:
			self.post_processing()

			# save PDBmulti as pkl file
			pklSeries = saveGenericObject(obj      = self.combinedAtoms,
										  fileName = self.seriesName)

			shutil.move(pklSeries,
						'{}{}'.format(self.outputDir,pklSeries))
			self.pklSeries = pklSeries

			# provide summary txt file on Dloss metric metric per-dataset
			self.summaryFile(normType = 'Standard') 

			if self.CalphaPresent is True:
				self.summaryFile(normType='Calpha normalised') 
			
			self.writeDamSitesToFile()

		else: 
			print 'Post processing job not chosen...'

		if retrieve is True:
			self.PDBmulti_retrieve()
		else:
			print 'PDBmulti retrieval from pkl file not chosen...'
		self.fillerLine(blank=True)

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
			print 'ETRACK input file "{}" not found'.format(self.inputFileName)
			return False
		return True

	def readInputFile(self, printText = True):

		# read input file e_Track_input.txt to specify location 
		# of input files and where to write output files

		props = {'inDir'          : 'inDir',
				 'outDir'         : 'outDir',
				 'damageset_name' : 'seriesName',
				 'initialPDB'     : 'initialPDB',
				 'PKLMULTIFILE'   : 'pklSeries',
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
				self.plot = True
		inputfile.close()

		if not self.initialPDB.endswith('.pdb'):
			self.initialPDB += '.pdb'

		# locate the correct format for the list of datasets within damage 
		# series. Currently two formats acceptable: (a) series-name + dataset-id
		# (per dataset), (b) input list of full dataset names (recommended).
		found = True
		try:
			datasetNums
		except UnboundLocalError:
			found = False
		if found is True:
			self.pdbNames = [self.seriesName+num for num in datasetNums.split(',')]
			return True
		found = True
		try:
			self.laterDatasets
		except AttributeError:
			found = False
		if found is True:
			self.pdbNames = self.laterDatasets.split(',')
			return True
		else:
			if printText is True:
				print 'Error! Unable to extract list of dataset names from input file'
			return False

	def checkInOutDirExist(self,
						   printText = True):

		# check that an input/output directories have been 
		# found and make subdirectories if present

		for dir in ([[self.inDir,'Input'],[self.outDir,'Output']]):
			if os.path.isdir(dir[0]) == False:
				str = '{} file location: {} does not exist. '.format(dir[1],dir[0])+\
					  'Please select an appropriate directory'
				if printText is True:
					print str
				return False
		return True

	def makeOutputDir(self, dirName = './'):

		# if the above sub directory does not exist, make it

		if not os.path.exists(dirName):
			os.makedirs(dirName)
			print 'New sub directory "{}" created to contain output files'.format(dirName)

	def setOutputDirs(self):

		# set the locations of the output directories

		self.outputDir 		= '{}ETRACK-output/'.format(self.outDir)
		self.outputPlotDir 	= '{}plots/'.format(self.outputDir)

		# add pkl file names as attribute if specified in input file
		if len(self.pklFiles) != 0: 
			self.pklFiles = [self.outputDir+f for f in self.pklFiles]

	def map_processing(self):
		self.titleCaption(title='Map Processing')

		# create additional subdirectories
		for oDir in (self.outputDir,self.outputPlotDir):
			self.makeOutputDir(dirName = oDir)

		# make pklFiles and dir to move all generated per-dataset pkl files to this
		pklFileDir = 'pklFiles-perDataset/'
		self.makeOutputDir(dirName = '{}{}'.format(self.outputDir,pklFileDir))

		pklFileNames = []
		for dataset in self.pdbNames:
			# derive per-atom density metrics from maps
			mapName1 = '{}_atoms.map'.format(dataset)
			mapName2 = '{}_density.map'.format(dataset)
			mapName3 = '{}_FC.map'.format(self.initialPDB.replace('.pdb',''))

			maps2DensMets 	= maps2DensMetrics(filesIn    = self.inDir,
											   filesOut   = self.outputDir,
											   pdbName    = dataset,
											   atomTagMap = mapName1,
											   densityMap = mapName2,
											   FCmap 	  = mapName3,
											   plotScatter = False,
											   plotHist   = self.plot,
											   plotBar    = False)

   			maps2DensMets.maps2atmdensity()

			# move pkl file to working output directory
			pklFileName = maps2DensMets.pklFileName

			shutil.move(pklFileName,
						'{}{}{}'.format(self.outputDir,
										pklFileDir,
										pklFileName))

			pklFileNames.append('{}{}{}'.format(self.outputDir,
												pklFileDir,
												pklFileName))

		self.pklFiles = pklFileNames

	def post_processing(self):
		self.titleCaption(title = 'Post Processing')
		print 'Input pkl files for post processing chosen from input file:'
		for file in self.pklFiles: 
			print '\t{}'.format(file.replace(self.outDir,""))

		# next read in the pdb structure file as list of atom objects
		print 'Reading in initial pdb file...'
		initialPDBlist = PDBtoList(pdbFileName = self.inDir+self.initialPDB)

		# retrieve object lists of atoms for each damage set
		self.fillerLine()
		print 'Reading in damaged pkl files...'
		dList = [] # dataset list
		for pkl_filename in self.pklFiles:
			print 'Damage file number: {}'.format(len(dList)+1)
			PDB_ret = retrieve_objectlist(pkl_filename)

			# add new retrieved damage set list to dList
			dList.append(PDB_ret)

		# create a list of atom objects with attributes as lists varying over 
		# dose range, only including atoms present in ALL damage datasets
		print 'New list of atoms over full dose range calculated...'
		combinedAtoms = combinedAtomList(datasetList    = dList,
										 numLigRegDsets = len(dList),
										 doseList       = self.doses,
										 initialPDBList = initialPDBlist,
										 outputDir      = self.outputDir,
										 seriesName     = self.seriesName)

		combinedAtoms.getMultiDoseAtomList()

		# calculate 'average' variant Dloss metrics
		combinedAtoms.calcAdditionalMetrics(newMetric='average')

		# calculate Calpha normalised metrics, if Calpha atoms exist
		self.CalphaPresent = combinedAtoms.checkCalphaAtomsExist()
		if self.CalphaPresent is True:
			for m in ('loss','mean','gain','Bfactor'):
				combinedAtoms.calcAdditionalMetrics(metric=m)
		
		self.combinedAtoms = combinedAtoms
		self.writeCsvFiles(output = self.output)

	def writeCsvFiles(self,
					  move   = True,
					  output = 'simple'):

		# write atom numbers and density metrics to simple 
		# csv files,one for each density metric separately

		self.fillerLine()
		print 'Writing .csv file for per-atom density metric:'

		if output == 'simple':
			metrics = [['loss','Standard'],
					   ['loss','Calpha normalised'],
					   ['loss','reliability']]
		else:
			metrics = self.combinedAtoms.getDensMetrics()

		for densMet in metrics:
			print '\tmetric: {}, normalisation: {}'.format(*densMet)
			self.combinedAtoms.writeMetric2File(where    = self.outputDir,
										        metric   = densMet[0],
										        normType = densMet[1])

		for groupBy in ('residue','atomtype'):
			self.combinedAtoms.writeMetric2File(where   = self.outputDir,
										   		groupBy = groupBy)
			if self.CalphaPresent is True:
				self.combinedAtoms.writeMetric2File(where    = self.outputDir,
										  	        groupBy  = groupBy,
										            normType = 'Calpha normalised')

		# make csvFiles dir and move all generated csv files to this
		if move is True:
			self.makeOutputDir(dirName = '{}csvFiles/'.format(self.outputDir))
			self.makeOutputDir(dirName = '{}csvFiles/Calpha-normalised/'.format(self.outputDir))
			self.makeOutputDir(dirName = '{}csvFiles/Standard/'.format(self.outputDir))

			for file in os.listdir(self.outputDir):
				if file.endswith(".csv"):
					if '-Calphanormalised' in file:
						loc = 'Calpha-normalised/'
					else:
						loc = 'Standard/'
					shutil.move('{}{}'.format(self.outputDir,file),
								'{}csvFiles/{}{}'.format(self.outputDir,loc,file))

	def PDBmulti_retrieve(self):
		self.fillerLine(blank = True)
		self.titleCaption(title='Atom List Retrieval')
		print 'Input pkl file for data retrieval chosen from input file:'
		print '\t{}'.format(self.pklSeries)

		# retrieve the combinedAtoms object from the pkl file
		self.combinedAtoms = retrieveGenericObject(fileName = self.outputDir+self.pklSeries)

	def summaryFile(self,
					fileType = 'html',
					metric   = 'loss',
					normType = 'Standard'):

		if fileType == 'html':
			self.summaryHTML(metric        = metric,
							 normType      = normType,
							 includeGraphs = self.plot)
		else:
			self.summaryTxt(metric   = metric,
				            normType = normType)

	def summaryTxt(self, 
				   metric   = 'loss', 
				   normType = 'Standard'):
		# produce a selection of per-dataset summary statistics
		# by default set metric as 'loss' for Dloss metric

		numDsets = self.getNumDatasets()
		summaryFile = open('{}summaryFile-D{}-{}.txt'.format(self.outputDir,metric,normType.replace(' ','-')),'w')
		summaryFile.write('D{} ({}) eTrack summary file\n'.format(metric,normType))
		summaryFile.write('Created: {}\n'.format(strftime("%Y-%m-%d %H:%M:%S", gmtime())))
		summaryFile.write('Summary information derived from {}\n'.format(self.pklSeries))
		summaryFile.write('Email charles.bury@dtc.ox.ac.uk for queries\n')

		summaryFile.write('\n-----------------------\n')
		summaryFile.write('Number of datasets reported in file: {}\n'.format(numDsets))
		summaryFile.write('Providing analysis on each dataset individually below\n')

		summaryFile.write('\n-----------------------\n')
		summaryFile.write('Key of reported quantities:\n')
		summaryFile.write('mean: average of metric calculated over all atoms of a specified type\n')
		summaryFile.write('std: standard deviation of metric calculated over all atoms of a specified type\n')
		summaryFile.write('#atoms: total number of atoms of a specified type\n')
		summaryFile.write('outliers: assuming a symmetric distn around the mode, number of atoms that fall outside this domain\n')
		summaryFile.write('skew: skewness of metric distribution for atoms of specified type\n')
		summaryFile.write('ratio: net ratio of distn values either side of metric distn mode\n\n')

		av,std = self.combinedAtoms.getAverageMetricVals(metric,normType) # get structure-wide metric average & std dev

		for i in range(numDsets):
			summaryFile.write('-----------------------'*3+'\n')
			summaryFile.write('Dataset info:\nNumber in series: {}\n'.format(i+1))
			summaryFile.write('Name: {}\n'.format(self.pdbNames[i]))
			summaryFile.write('Dose (MGy): {}\n'.format(self.doses[i]))

			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Structure-wide D{} summary:\n'.format(metric))
			summaryFile.write('Average D{}: {}\n'.format(metric,round(av[i],3)))
			summaryFile.write('Std dev in D{}: {}\n'.format(metric,round(std[i],3)))

			summaryFile.write('# atoms with {} Dloss metric above N std dev of structure-wide mean:\n'.format(normType))
			summaryFile.write('N\t\t#atoms\n')
			n = self.combinedAtoms.numAtmsWithMetricAboveLevel(dataset   = i,
															   metric    = metric,
															   normType  = normType,
															   threshold = 0.5)
			summaryFile.write('{}\t\t{}\n'.format(0.5,n))
			t = 0
			while n != 0:
				t += 1
				n = self.combinedAtoms.numAtmsWithMetricAboveLevel(dataset   = i,
															   	   metric    = metric,
															       normType  = normType,
															   	   threshold = t)
				summaryFile.write('{}\t\t{}\n'.format(t,n))

			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-atom Statistics:\n')
			summaryFile.write('Top hits ranked by D{} metric:\n'.format(metric))
			statsOut = self.combinedAtoms.getTopNAtomsString(metric,normType,i,25)
			summaryFile.write(statsOut)

			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-atom-type Statistics:\n')
			summaryFile.write('D{} metric ranked by mean value:\n'.format(metric))
			statsOut = self.combinedAtoms.getPerAtmtypeStats(metric,normType,i,'mean',25)
			summaryFile.write(statsOut[0]+'\n')

			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-atom Statistics:\n')
			summaryFile.write('Top D{} hits grouped by residue type:\n'.format(metric))
			n = self.combinedAtoms.getNumAtoms()*0.05 # take top 5% of atoms here
			statsOut = self.combinedAtoms.breakdownTopNatomsBy(metric   = metric,
															   normType = normType,
															   dataset  = i,
															   n        = n)

			summaryFile.write(statsOut[0]+statsOut[1]+'\n')

			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-residue Statistics:\n')
			summaryFile.write('D{} metric ranked by mean value:\n'.format(metric))
			statsOut = self.combinedAtoms.getPerResidueStats(metric,normType,i,'mean','all')
			summaryFile.write(statsOut[0]+'\n')

			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-chain Statistics:\n')
			summaryFile.write('D{} metric ranked by mean value:\n'.format(metric))
			statsOut = self.combinedAtoms.getPerChainStats(metric,normType,i,'mean','all')
			summaryFile.write(statsOut[0]+'\n')
		summaryFile.close()

	def metricBarplots(self):
		# plot barplots of damage metric for each susceptible residue type
		numDsets = self.getNumDatasets()
		for i in range(numDsets):
			for set in [1,2]:
				for b in ('Box','Bar'):
					self.combinedAtoms.susceptAtmComparisonBarplot(metric,normType,i,set,b)

	def summaryHTML(self, 
					metric        = 'loss', 
					normType      = 'Standard',
					includeGraphs = True):

		# produce a selection of per-dataset summary statistics

		numDsets = self.getNumDatasets()
		summaryFile = open('{}summaryFile-D{}-{}.html'.format(self.outputDir,metric,normType.replace(' ','-')),'w')
		summaryFile.write('<!DOCTYPE html>\n<html>\n<head>\n')

		headString = '<html>\n<head>\n<title>eTrack summary file</title>\n'+\
					 '<style>\ntable, th, td {\nborder: 1px solid black;\nborder-collapse: collapse;\n}\n'+\
					 'th, td {\npadding: 5px;\ntext-align: center;\n}\n</style>\n</head>\n'
		summaryFile.write(headString)

		bodyString = '<body>\n<h1>D{} ({}) eTrack summary file</h1>\n'.format(metric,normType)+\
					 'Created: {}<br>\n'.format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))+\
					 'Summary information derived from {}<br>\n'.format(self.pklSeries)+\
					 'Email charles.bury@dtc.ox.ac.uk for queries<br>\n'+\
					 'Number of datasets reported in file: {}<br>\n'.format(numDsets)+\
					 'Providing analysis on each dataset individually below<br>\n'+\
					 '<h3>Key of reported quantities</h3>\n'+\
					 '<ul><li>mean: average of metric calculated over all atoms of a specified type</li>\n'+\
					 '<li>std: standard deviation of metric calculated over all atoms of a specified type</li>\n'+\
					 '<li>#atoms: total number of atoms of a specified type</li>\n'+\
					 '<li>outliers: assuming a symmetric distn around the mode, number of atoms that fall outside this domain</li>\n'+\
					 '<li>skew: skewness of metric distribution for atoms of specified type</li>\n'+\
					 '<li>normality: p-value for null hypothesis that distribution of metric values is normally distributed. If not enough atoms are present to perform this test, "n/a" is given</li>\n</ul>'
		summaryFile.write(bodyString)

		# get structure-wide metric average & std dev
		av,std = self.combinedAtoms.getAverageMetricVals(densMet  = metric,
														 normType = normType)
		for i in range(numDsets):
			summaryFile.write('<hr>')
			bodyString = '<h3>Dataset info</h3>\n'+\
						 'Number in series: {}<br>\n'.format(i+1)+\
						 'Dose (MGy): {}<br>\n'.format(self.doses[i])+\
						 '<h3>Structure-wide D{} summary</h3>\n'.format(metric)+\
						 'Average D{}: {}<br>\n'.format(metric,round(av[i],3))+\
						 'Std dev in D{}: {}<br>\n'.format(metric,round(std[i],3))

			if normType == 'Calpha normalised':
				CAweights = self.combinedAtoms.retrieveCalphaWeight(metric = metric)
				bodyString += 'Calpha weight for current dataset: {}<br>\n'.format(round(CAweights.weight[metric][i],3))
			summaryFile.write(bodyString)

			if includeGraphs is True:
				failedPlots = self.makeDistnPlots(densMet  = metric,
											      normType = normType,
											      plotSet  = 4)
				plotString = '<h3>Distribution of Dloss for all refined atoms within structure</h3>\n'+\
							 '<img src="plots/DistnPlot_Residues-all_Metric-D{}_Normalisation-{}_Dataset-{}.png"><br>'.format(metric,normType.replace(' ',''),i)
				summaryFile.write(plotString)

			bodyString = '# atoms with {} Dloss metric above N std dev of structure-wide mean:<br><br>\n'.format(normType)
			summaryFile.write(bodyString)

			n = self.combinedAtoms.numAtmsWithMetricAboveLevel(dataset   = i,
															   metric    = metric,
															   normType  = normType,
															   threshold = 0.5)
			tableString = '<table border="1" style="width:10%">\n'+\
						  '<tr><th>N</th><th>#atoms</th></tr>\n'+\
						  '<tr><td>{}</td><td>{}</td></tr>\n'.format(0.5,n)
			t = 0
			while n != 0:
				t += 1
				n = self.combinedAtoms.numAtmsWithMetricAboveLevel(dataset   = i,
															       metric    = metric,
															       normType  = normType,
															       threshold = t)
				tableString += '<tr><td>{}</td><td>{}</td></tr>\n'.format(t,n)
			tableString += '</table>\n'
			summaryFile.write(tableString)

			statsOut = self.combinedAtoms.getTopNAtomsString(metric,normType,i,25)
			infoString = '<h3>Per-atom Statistics</h3>\n'+\
						 'Top hits ranked by D{} metric.<br>Dmean and Dgain are mean '.format(metric)+\
						 'and maximum voxel difference density values, respectively, '+\
						 'assigned within a local region around each atom.<br>Proximity (from 0 '+\
						 'to 1) is a measure of the closeness of the voxel exhibiting the '+\
						 'maximum density loss Dloss value from the specified atom (higher '+\
						 'values indicate smaller distances):<br><br>\n'
			infoString += self.convertPlainTxtTable2html(statsOut,width='50%')
			summaryFile.write(infoString)

			statsOut = self.combinedAtoms.getPerAtmtypeStats(metric   = metric,
															 normType = normType,
															 dataset  = i)
			infoString = '<h3>Per-atom-type Statistics</h3>\n'+\
						 'D{} metric ranked by mean value:<br><br>\n'.format(metric)
			infoString += self.convertPlainTxtTable2html(statsOut[0],width='60%')
			summaryFile.write(infoString)

			n = self.combinedAtoms.getNumAtoms()*0.05 # take top 5% of atoms here
			statsOut = self.combinedAtoms.breakdownTopNatomsBy(metric   = metric,
															   normType = normType,
															   dataset  = i,
															   n        = n)

			infoString = '<h3>Per-atom Statistics</h3>\n'+\
						 '{}<br><br>\n'.format(statsOut[0])
			infoString += self.convertPlainTxtTable2html(statsOut[1])
			summaryFile.write(infoString)

			statsOut = self.combinedAtoms.getPerResidueStats(metric   = metric,
															 normType = normType,
															 dataset  = i)
			infoString = '<h3>Per-residue Statistics</h3>\n'+\
						 'D{} metric ranked by mean value:<br><br>\n'.format(metric)
			infoString += self.convertPlainTxtTable2html(statsOut[0],width='60%')
			summaryFile.write(infoString)

			statsOut = self.combinedAtoms.getPerChainStats(metric   = metric,
														   normType = normType,
														   dataset  = i,
														   n        = 'all')
			infoString = '<h3>Per-chain Statistics</h3>\n'+\
						 'D{} metric ranked by mean value:<br><br>\n'.format(metric)
			infoString += self.convertPlainTxtTable2html(statsOut[0],width='60%')
			summaryFile.write(infoString)

			if includeGraphs is True:
				failedPlots = self.makeDistnPlots(densMet  = metric,
											      normType = normType,
											      plotSet  = 1)
				if i not in failedPlots['GLUASPCYSMETTYR']:
					plotString = '<h3>Distribution of Dloss for known susceptible residue types</h3>\n'+\
							 '<img src="plots/DistnPlot_Residues-GLUASPCYSMETTYR_Metric-D{}_Normalisation-{}_Dataset-{}.png"><br>'.format(metric,normType.replace(' ',''),i)
					summaryFile.write(plotString)

				failedPlots = self.makeDistnPlots(densMet  = metric,
												  normType = normType,
												  plotSet  = 3)
				if i not in failedPlots['DADCDGDT']:
					plotString = '<h3>Distribution of Dloss for known susceptible nucleotide types</h3>\n'+\
								 '<img src="plots/DistnPlot_Residues-DADCDGDT_Metric-D{}_Normalisation-{}_Dataset-{}.png"><br>'.format(metric,normType.replace(' ',''),i)
					summaryFile.write(plotString)

			infoString = '<h3>Suspicious Atoms</h3>\n'
			suspAtomLens = []
			for t in [6,5,4,3,2]:
				suspAtoms = self.combinedAtoms.detectSuspiciousAtoms(dataset   = i,
							  										 metric    = metric,
							  										 normType  = normType,
							  										 threshold = t)

				infoString+= '{} atoms found with unusually high/low D{} '.format(len(suspAtoms),metric)+\
							 'values compared to average of that atom type '+\
							 '(over {} standard deviations from average)'.format(t)

				if len(suspAtoms) > 0:
					infoString += ':<br>'
					suspAtomLens.append(len(suspAtoms))
					for s in suspAtoms:
						infoString += '{}<br>'.format(s)
				infoString += '<br>'

				if len(suspAtomLens) > 1:
					if suspAtomLens[-1] > 0 and suspAtomLens[-2] > 0:
						break

			summaryFile.write(infoString)

		summaryFile.write('</body>\n</html>')
		summaryFile.close()

		# TEST CASES - NOT NORMALLY PLOTTED
		# self.plotLinePlot(restype   = 'GLU',
		# 				  errorBars = 'NONE',
		# 				  atomType  = '')
		# self.plotLinePlot(restype   = 'GLU',
		# 				  errorBars = 'ATOMTYPE',
		# 				  atomType  = '')
		# self.plotLinePlot(restype   = 'CYS',
		# 	              errorBars = 'NONE',
		# 	              atomType  = 'SG')

	def plotLinePlot(self,
					 restype   = 'GLU',
					 errorBars = 'NONE',
					 atomType  = ''):

		self.combinedAtoms.graphMetric(atomType  = atomType,
									   restype   = restype,
									   errorBars = errorBars,
									   outputDir = self.outputPlotDir,
									   saveFig   = True)

	def convertPlainTxtTable2html(self,
								  plainText      = '',
								  numHeaderLines = 1,
								  width          = '40%',
								  border         = '1'):

		htmlTable = '<table border="{}" style="width:{}">\n'.format(border,width)
		for i, l in enumerate(plainText.split('\n')):
			if i < numHeaderLines:
				newStr = '<tr><th>'+' '.join(l.split()).replace(' ','</th><th>')+'</th></tr>'
			else:
				newStr = '<tr><td>'+' '.join(l.split()).replace(' ','</td><td>')+'</td></tr>'
			htmlTable += newStr+'\n'
		htmlTable += '</table>\n'
		return htmlTable

	def getNumDatasets(self):
		# get number of datasets in damage series
		numDsets = self.combinedAtoms.atomList[0].getNumDatasets()
		return numDsets

	def writeDamSitesToFile(self, 
							metric      = 'loss', 
							normType    = 'Standard', 
							numDamSites = 25):
		# write top damage sites to .pdb file for each dataset

		self.damSitesPDB = []
		for i in range(self.getNumDatasets()): 
			damPDB = self.combinedAtoms.getTopNAtomsPDBfile(metric   = metric,
															normType = normType,
															dataset  = i,
															n        = numDamSites,
															pdbFile  = self.inDir+self.initialPDB)

			self.damSitesPDB.append(damPDB)

	def colorByMetric(self, 
					  metric    = 'loss', 
					  normType  = 'Standard', 
					  dataset   = 0, 
					  singleRes = ''):

		# for the initial pdb file in damage series, convert the Bfactor 
		# column to values for the specified metric. If 'singleRes' is 
		# specified (use 3-letter residue code) then the average density 
		# metric for each atom within this residue is calculated within 
		# the structure and only this is output to resulting PDB file, 
		# otherwise use ''

		if normType == 'Calpha normalised': 
			self.combinedAtoms.calcAdditionalMetrics(metric=metric,normType=normType)

		pdbIn = open(self.inDir + self.initialPDB,'r')
		fileOut = self.outputDir + self.initialPDB.strip('.pdb')+'_{}D{}_{}.pdb'.format(normType.replace(" ",""),metric,dataset)
		if singleRes != '':
			fileOut = fileOut.strip('.pdb')+'-{}.pdb'.format(singleRes)

		pdbOut = open(fileOut,'w')
		pdbOut.write('REMARK\tBfactor column replaced by {} D{} metric values\n'.format(normType,metric))
		for l in pdbIn.readlines():
			if l.split()[0] in ('CRYST1','SCALE1','SCALE2','SCALE3'):
				pdbOut.write(l)
			elif 'ATOM' in l.split()[0]:
				break
		pdbIn.close()

		if singleRes == '':
			for atm in self.combinedAtoms.atomList:
				dens = atm.densMetric[metric][normType]['values'][dataset]
				if not np.isnan(dens): # don't include atoms for which not calculated properly
					l = writePDBline(atm,dens)
					pdbOut.write(l+'\n')
		else:
			atmDic = self.combinedAtoms.getAvMetricPerAtmInRes(singleRes,metric,normType,dataset)
			for atm in self.combinedAtoms.atomList:
				if atm.basetype == singleRes:
					resNum = atm.residuenum # save the residue number now
					break
			for atm in self.combinedAtoms.atomList:
				if atm.residuenum == resNum:
					dens = atmDic[atm.atomtype]
					l = writePDBline(atm,dens)
					pdbOut.write(l+'\n')
		pdbOut.write('END')
		pdbOut.close()

	def visualiseDamSites(self, 
						  dataset  = 0, 
						  metric   = 'loss', 
						  software = 'pymol', 
						  size     = 1,
						  savePic  = True):

		# open coot/pymol to view top damage sites. size is the 
		# density metric scale factor used when visualising damage 
		# sites as spheres of vdw = size*{density metric} within pymol

		if software not in ('coot','pymol'):
			print 'Damage sites can be visualised in either "coot" or "pymol"'
			print 'Please make sure the paths to these programs are '+\
				  'correctly set up before proceeding'
			return
		try: self.damSitesPDB
		except AttributeError:
			return "Must run .writeDamSitesToFile() before damage sites can be read in {}".format(software)
		if software == 'coot':
			os.system('coot -pdb {} -pdb2 {}'.format(self.inDir+self.initialPDB,self.damSitesPDB[dataset]))
		else:
			# need to write script for pymol to run with
			damSitesTag = (self.damSitesPDB[dataset].split('/')[-1]).strip('.pdb')
			structureTag = self.initialPDB.strip('.pdb')
			scriptName = self.outputDir+'runPymolScript.pml'
			pymolScript = open(scriptName,'w')

			pymolStr = 'load {}\n'.format(self.inDir+self.initialPDB)+\
					   'load {}\n'.format(self.damSitesPDB[dataset])+\
					   'hide lines\nshow cartoon\n'+\
					   'set antialias, 1\n'+\
					   'set ray_trace_mode, 0\n'+\
					   'set cartoon_fancy_helices, 1\n'+\
					   'set cartoon_side_chain_helper, on\n'+\
					   'set ray_opaque_background, 0\n'+\
					   'show sphere, {}\n'.format(damSitesTag)+\
					   'alter {}, vdw=b*{}\n'.format(damSitesTag,size)+\
					   'rebuild\n'+\
					   'select nearDamage, {} w. 4 of {}\n'.format(structureTag,damSitesTag)+\
					   'show sticks, nearDamage\n'+\
					   'set stick_transparency, 0.6\n'+\
					   'color red, {}\n'.format(damSitesTag)+\
					   'color white, {}\n'.format(structureTag)

			if savePic is True:
				# DOESN'T WORK CURRENTLY
				pymolStr += 'ray\n png {}\nquit'.format(self.damSitesPDB[dataset].replace('.pdb','.png'))

			pymolScript.write(pymolStr)
			pymolScript.close()
			os.system('pymol -q {}'.format(scriptName))

	def makeDistnPlots(self,
					   densMet  = 'loss',
					   normType = 'Standard',
					   plotType = 'both',
					   plotSet  = 1):

		# retrieve the distribution for damage metric values for all 
		# atoms of specified residues (see 'resList' below), and then 
		# plot a kde plot for each

		if plotSet == 1:
			resList = [['GLU','ASP','CYS','MET','TYR']]
		elif plotSet == 2:
			resList = [['GLU','GLN'],
					   ['ASP','ASN'],
					   ['ILE','LEU'],
					   ['TYR','PHE'],
				   	   ['GLU','GLY'],
				   	   ['ASP','GLY'],
				   	   ['TYR','GLY'],
				   	   ['PHE','GLY'],
				   	   ['CYS','GLY'],
				   	   ['MET','GLY'],
				   	   ['TYR','ASP','GLU'],
				   	   ['TYR','PHE','GLY'],
				   	   ['GLU','ASP','CYS','MET','TYR']]
		elif plotSet == 3:
			resList = [['DA','DC','DG','DT']]

		elif plotSet == 4:
			resList = ['all']

		for resGroup in resList:
			k = ''.join(resGroup)
			failedPlots = {k:[]}
			if len(resGroup) > 6:
				plotType = 'kde'
			else:
				plotType = 'both'
			for i in range(self.getNumDatasets()):
				data = self.combinedAtoms.graphMetricDistn(metric    = densMet,
														   normType  = normType,
														   valType   = i,
														   plotType  = plotType,
														   resiType  = resGroup,
														   outputDir = self.outputPlotDir,
														   printText = False)	
				if data == {}:
					failedPlots[k].append(i)
		return failedPlots

	def sensAtomPlots(self,
					  metric   = 'loss',
					  normType = 'Standard'):

		# set of plots investigating damage progression 
		# for sensitive atoms within structure.

		pairs = {'susceptRes': [['GLU','CD','CG'],
							    ['GLU','CD','OE1'],
							    ['ASP','CG','CB'],
							    ['ASP','CG','OD1'],
							    ['TYR','OH','CZ'],
							    ['TYR','CZ','CE2'],
							    ['MET','SD','CG'],
							    ['MET','SD','CE'],
							    ['CYS','SG','CB'],
							    ['CYS','CA','CB']],
				'PHEvTYR'    : [['TYR','CZ','OH'],
					            ['TYR','CZ','CE2'],
					            ['TYR','CZ','CE1'],
							    ['TYR','CZ','CG'],
							    ['PHE','CZ','CE2'],
							    ['PHE','CZ','CE1'],
							    ['PHE','CZ','CG']],
				'TYR'        : [['TYR','CB','CG'],
					            ['TYR','CB','CD1'],
							    ['TYR','CB','CD2'],
							    ['TYR','CB','CE1'],
							    ['TYR','CB','CE2'],
							    ['TYR','CB','CZ'],
							    ['TYR','CB','OH']],
				'PHE'        : [['PHE','CB','CG'],
							    ['PHE','CB','CD1'],
							    ['PHE','CB','CD2'],
							    ['PHE','CB','CE1'],
							    ['PHE','CB','CE2'],
							    ['PHE','CB','CZ']]}

		# determine ratios of atoms within susceptible residue 
		# types and write output txt file and plot
		for met in ('distance','ratio'):

			for key in pairs.keys():
				self.combinedAtoms.findMetricRatioKeyResidues(metric   = 'loss',
															  normType = 'Standard',
															  rType    = met,
															  pairs    = pairs[key],
															  title    = key)
		atomtypes = [[['GLU','CD'],
					  	['ASP','CG'],
					  		['TYR','OH'],
					  			['CYS','SG'],
					  				['MET','SD']],
					 [['GLU','CD'],
					  	['ASP','CG'],
					  		['TYR','OH'],
					  			['CYS','SG'],
					  				['MET','SD'],
					 					['PHE','CZ'],
					  						['TYR','CZ']],
					 [['GLU','CG'],
					  	['ASP','CB'],
					  		['TYR','OH'],
					  			['CYS','CB'],
					  				['MET','CG'],
					  					['MET','CE'],
					  						['PHE','CZ'],
					  							['TYR','CZ']],
					 [['TYR','OH'],
					  	['PHE','CZ'],
					  		['TYR','CZ'],
					  			['PHE','CE1'],
					  				['TYR','CE1'],
					  					['PHE','CE2'],
					  						['TYR','CE2']]]

		for aTypes in atomtypes:
			self.combinedAtoms.plotSusceptibleAtoms(densMet = metric,
													sus = aTypes,
													normType  = normType)

			self.combinedAtoms.plotSusceptibleAtoms(densMet   = metric,
													normType  = normType,
													errorbars = False,
													susAtms   = aTypes)

	def getSpaceGroup(self):
		pdbFile = self.inDir+self.initialPDB
		pdbin = open(self.inDir+self.initialPDB,'r')
		for line in pdbin.readlines():
			if line.split()[0] == 'CRYST1':
				self.spaceGroup = line[55:66].replace(' ','')
		pdbin.close()
		try: 
			self.spaceGroup
		except attributeError:
			print 'Unable to find space group from file: {}'.format(pdbFile)
			return False
		return self.spaceGroup

	def fillerLine(self,blank=False):
		if blank is False:
			print '---------------------------------------------------------------'	
		else:
			print '\n'

	def titleCaption(self,title='unspecified title'):
		print '||========================== {} ==========================||'.format(title)
