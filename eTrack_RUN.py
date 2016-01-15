# -*- coding: utf-8 -*-

from readAtomMap import maps2DensMetrics
from savevariables import retrieve_objectlist,save_objectlist,saveGenericObject,retrieveGenericObject
from PDBFileManipulation import PDBtoList
from combinedAtomList import combinedAtomList
import os
 
class eTrack(object):
	# A class for retrieving the eTrack input text file information and running 
	# the eTrack pipeline functions separately or together in full pipeline
	def __init__(self, where = "", pdbNames = [], pklFiles = [], initialPDB = "",
				 seriesName = "", pklSeries = "",doses=[],plot = False):

		self.where 			= where 		# the input and output file directory
		self.pdbNames 		= pdbNames 		# the list of pdb codes for series
		self.pklFiles 		= pklFiles 		# list of pkl files from map_processing
		self.initialPDB 	= initialPDB 	# the first dataset pdb code
		self.seriesName 	= seriesName 	# the general series name 
		self.pklSeries		= pklSeries 	# the combined series pkl file from post_processing
		self.doses 			= doses 		# list of increasing doses
		self.plot           = plot 			# (boolian) decide to plot per-residue summary plots per dataset

	def runPipeline(self,map_process,post_process,retrieve_PDBmulti,inputfilename):
		# the following function reads in the above functions one by one in a scripted 
		# pipeline. Takes inputs to specify which parts of above pipeline are to be
		# included. For each function input specify True if this part is to be performed
		# and False otherwise.
		self.titleCaption('ETRACK pipeline')

		# check whether valid inputs to function
		valid = self.checkValidInputs(map_process,post_process,retrieve_PDBmulti,inputfilename)
		if valid is False: return

		# first need to run function above to read in input file containing info
		# on where input files are and where output files should be written
		print 'Reading input file: {}'.format(str(inputfilename))
		success = self.readInputFile(inputfilename)
		if success is False: return

		if map_process is True:
			self.map_processing()
		else:
			print 'Map processing task not chosen...'
		self.fillerLine()

		if post_process is True:
			self.post_processing()
			# save PDBmulti as pkl file
			pklSeries = saveGenericObject(self.combinedAtoms,self.seriesName)
			os.system('mv {} {}{}'.format(pklSeries,self.outputDir,pklSeries))
			self.pklSeries = pklSeries
		else: 
			print 'Post processing job not chosen...'
		self.fillerLine()

		if retrieve_PDBmulti is True:
			self.PDBmulti_retrieve()
		else:
			print 'PDBmulti retrieval from pkl file not chosen...'
		self.fillerLine()

	def checkValidInputs(self,map_process,post_process,retrieve_PDBmulti,inputfilename):
		# check the runPipeline inputs to make sure that they are valid
		for var in [map_process,post_process,retrieve_PDBmulti]:
			if var not in (True,False):
				print 'Only Boolian parameters expected here'
				return False
		try:
			f = open(inputfilename,'r')
			f.close()
		except IOError:
			print 'ETRACK input file "{}" not found'.format(inputfilename)
			return False
		return True

	def readInputFile(self,inputfilename):
		# this function reads in the input file e_Track_input.txt, which 
		# specifies the location of all the relevant input files and where
		# the output files should be written

		pklFiles = [] # preallocate pklFiles list for post_processing pipeline step below
		
		# read in input file and determine location 'where' and name of input files
		inputfile = open(inputfilename,'r')
		for line in inputfile.readlines():
			l = line.split()
			if '#' == line[0]: continue
			elif 'where' in l[0]:
				self.where 		= l[1]
			elif 'damageset_name' in l[0]:
				self.seriesName = l[1]
			elif 'damageset_num' in l[0]:
				datasetNums = l[1]
			elif 'PKLFILE' in l[0]:
				pklFiles.append(l[1])
			elif 'initialPDB' in l[0]:
				self.initialPDB = l[1]
			elif 'PKLMULTIFILE' in l[0]:
				self.pklSeries 	= l[1]
			elif 'doses' in l[0]:
				self.doses = [float(d) for d in l[1].split(',')]
			elif 'plot' in l[0]:
				self.plot = True

		# check that an output directory has been found and make subdirectories if present
		if os.path.isdir(self.where) == False:
			print 'Input file location: {} does not exist. Please select an appropriate directory'.format(self.where)
			return False
		self.outputDir 		= '{}ETRACK_output/'.format(self.where)
		self.outputPlotDir 	= '{}plots/'.format(self.outputDir)

		# add pkl file names as attribute if specified in input file
		if len(pklFiles) != 0: 
			self.pklFiles = [self.outputDir+f for f in pklFiles]
		inputfile.close()

		# next put the dataset name and damage set numbers together in a list
		self.pdbNames = [self.seriesName+num for num in datasetNums.split(',')]
		return True

	def map_processing(self):
		self.titleCaption('Map Processing')

		# create additional subdirectories
		for oDir in (self.outputDir,self.outputPlotDir):
			if not os.path.exists(oDir):
				os.makedirs(oDir)

		pklFileNames = []
		for dataset in self.pdbNames:
			# derive per-atom density metrics from maps
			mapfilname1 	= '{}_atoms.map'.format(dataset)
			mapfilname2 	= '{}_density.map'.format(dataset)
			maps2DensMets 	= maps2DensMetrics(self.where,self.outputDir,dataset,	
											   mapfilname1,'atom_map',
											   mapfilname2,'density_map',
											   self.plot)
   			maps2DensMets.maps2atmdensity()

			# move pkl file to working output directory
			pklFileName = maps2DensMets.pklFileName
			os.system('mv {} {}{}'.format(pklFileName,self.outputDir,pklFileName))
			pklFileNames.append(self.outputDir+pklFileName)

		self.pklFiles = pklFileNames

	def post_processing(self):
		self.titleCaption('Post Processing')
		print 'Input pkl files for post processing chosen from input file:'
		for file in self.pklFiles: 
			print '\t{}'.format(file)

		# next read in the pdb structure file as list of atom objects
		print 'Reading in initial pdb file...'
		initialPDBlist = PDBtoList(self.where+self.initialPDB,[])

		# retrieve object lists of atoms for each damage set
		self.fillerLine()
		print 'Reading in damaged pkl files...'
		dList = [] # dataset list
		for pkl_filename in self.pklFiles:
			print 'Damage file number: {}\n'.format(len(dList)+1)
			PDB_ret = retrieve_objectlist(pkl_filename)

			# add new retrieved damage set list to dList
			dList.append(PDB_ret)

		# create a list of atom objects with attributes as lists varying over 
		# dose range, only including atoms present in ALL damage datasets
		print 'New list of atoms over full dose range calculated...'
		combinedAtoms = combinedAtomList(dList,len(dList),self.doses,initialPDBlist,
										 self.outputDir,True,self.seriesName)
		combinedAtoms.getMultiDoseAtomList()

		# calculate 'standardised' and 'average' variant Dloss metrics
		combinedAtoms.calcStandardisedMetrics('loss')
		combinedAtoms.calcAdditionalMetrics('loss','Standard','average')
		
		# write atom numbers and density metrics to simple text files - one for 
		# each density metric separately
		for densMet in combinedAtoms.getDensMetrics():
			print 'Writing .txt file for per-atom density metric: {}, normalisation: {}'.format(*densMet)
			combinedAtoms.writeMetric2File(self.outputDir,*densMet)
		
		self.combinedAtoms = combinedAtoms
		self.summaryStats()
		#self.sensAtomPlots()

	def PDBmulti_retrieve(self):
		self.titleCaption('Atom List Retrieval')
		print 'Input pkl file for data retrieval chosen from input file:'
		print '\t{}'.format(self.outputDir+self.pklSeries)
		# retrieve the combinedAtoms object from the pkl file
		self.combinedAtoms = retrieveGenericObject(self.outputDir+self.pklSeries)

	def summaryStats(self):
		# produce a selection of per-dataset summary statistics
		summaryFile = open('{}/summaryFile.txt'.format(self.outputDir),'w')
		for i in range(self.combinedAtoms.atomList[0].getNumDatasets()):

			# standard Dloss statistics
			summaryFile.write('\n\nDataset: {}'.format(i))
			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-residue Statistics:\n')
			summaryFile.write('Dloss metric ranked by mean value:\n')
			statsOut = self.combinedAtoms.getPerResidueStats('loss','Standard',i,'mean','all')
			summaryFile.write(statsOut[0])
			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-chain Statistics:\n')
			summaryFile.write('Dloss metric ranked by mean value:\n')
			statsOut = self.combinedAtoms.getPerChainStats('loss','Standard',i,'mean','all')
			summaryFile.write(statsOut[0])
			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-atom-type Statistics:\n')
			summaryFile.write('Dloss metric ranked by mean value:\n')
			statsOut = self.combinedAtoms.getPerAtmtypeStats('loss','Standard',i,'mean',25)
			summaryFile.write(statsOut[0])
			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-atom Statistics:\n')
			summaryFile.write('Top hits ranked by Dloss metric:\n')
			statsOut = self.combinedAtoms.getTopNAtomsString('loss','Standard',i,25)
			summaryFile.write(statsOut)
			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-atom Statistics:\n')
			summaryFile.write('Top Dloss hits grouped by residue type:\n')
			n = self.combinedAtoms.getNumAtoms()*0.1 # take top 10% of atoms here
			statsOut = self.combinedAtoms.breakdownTopNatomsBy('loss','Standard',i,n,['basetype','atomtype'])
			summaryFile.write(statsOut)

			# standard Dmean statistics
			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-residue Statistics:\n')
			summaryFile.write('Dmean metric ranked by mean:\n')
			statsOut = self.combinedAtoms.getPerResidueStats('mean','Standard',i,'mean','all')
			summaryFile.write(statsOut[0])
			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-chain Statistics:\n')
			summaryFile.write('Dmean metric ranked by mean:\n')
			statsOut = self.combinedAtoms.getPerChainStats('mean','Standard',i,'mean','all')
			summaryFile.write(statsOut[0])
			summaryFile.write('\n-----------------------\n')
			summaryFile.write('Per-atom-type Statistics:\n')
			summaryFile.write('Dmean metric ranked by mean value:\n')
			statsOut = self.combinedAtoms.getPerAtmtypeStats('mean','Standard',i,'mean',25)
			summaryFile.write(statsOut[0])
		summaryFile.close()

	def sensAtomPlots(self):
		# set of plots investigating damage progression for sensitive atoms within structure

		# determine ratios of atoms within susceptible residue types and write output txt file and plot
		for met in ('distance','ratio'):
			pairs = [['GLU','CD','CG'],['GLU','CD','OE1'],['ASP','CG','CB'],['ASP','CG','OD1'],
					 ['TYR','OH','CZ'],['TYR','CZ','CE2'],['MET','SD','CG'],['MET','SD','CE'],
					 ['CYS','SG','CB'],['CYS','CA','CB']]
			self.combinedAtoms.findMetricRatioKeyResidues('loss','Standard',met,True,pairs,'susceptRes')

			pairs = [['TYR','CZ','OH'],['TYR','CZ','CE2'],['TYR','CZ','CE1'],['TYR','CZ','CG'],
					 ['PHE','CZ','CE2'],['PHE','CZ','CE1'],['PHE','CZ','CG']]
			self.combinedAtoms.findMetricRatioKeyResidues('loss','Standard',met,True,pairs,'PHEvTYR')

			pairs = [['TYR','CB','CG'],['TYR','CB','CD1'],['TYR','CB','CD2'],['TYR','CB','CE1'],
					 ['TYR','CB','CE2'],['TYR','CB','CZ'],['TYR','CB','OH']]
			self.combinedAtoms.findMetricRatioKeyResidues('loss','Standard',met,True,pairs,'TYR')

			pairs = [['PHE','CB','CG'],['PHE','CB','CD1'],['PHE','CB','CD2'],['PHE','CB','CE1'],
					 ['PHE','CB','CE2'],['PHE','CB','CZ']]
			self.combinedAtoms.findMetricRatioKeyResidues('loss','Standard',met,True,pairs,'PHE')

		self.combinedAtoms.calcAdditionalMetrics('','Standard','netChange') # calculate Dnet metric
		atomtypes = [] 
		atomtypes.append([['GLU','CD'],['ASP','CG'],['TYR','OH'],['CYS','SG'],['MET','SD']])
		atomtypes.append([['GLU','CD'],['ASP','CG'],['TYR','OH'],['CYS','SG'],['MET','SD'],['PHE','CZ'],['TYR','CZ']])
		atomtypes.append([['GLU','CG'],['ASP','CB'],['TYR','OH'],['CYS','CB'],['MET','CG'],['MET','CE'],['PHE','CZ'],['TYR','CZ']])
		atomtypes.append([['TYR','OH'],['PHE','CZ'],['TYR','CZ'],['PHE','CE1'],['TYR','CE1'],['PHE','CE2'],['TYR','CE2']])
		for metric in ('loss','mean','gain','net'):
			for aTypes in atomtypes:
				self.combinedAtoms.plotSusceptibleAtoms(metric,'Standard',True,'y',aTypes)
				self.combinedAtoms.plotSusceptibleAtoms(metric,'Standard',False,'y',aTypes)

	def fillerLine(self):
		print '---------------------------------------------------------------'	

	def titleCaption(self,title):
		print '||========================== {} ==========================||'.format(title)
