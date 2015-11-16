# -*- coding: utf-8 -*-

from readAtomMap import maps2DensMetrics
from savevariables import retrieve_objectlist,save_objectlist
from PDBFileManipulation import multiARRAY_diffatomnumbers
from topDamageHits import topNdamsites_resibarplotter,topNdamsites_chainbarplotter,topNdamsites_printer
import os
from PDBMulti2Txt import objlist2txt
from findMetricChange import find_Bchange
from PDBFileManipulation import PDBtoCLASSARRAY_v2 as pdb2list
from BDamageCalculator import numsurroundatoms_calculate,bdamage_calculate,numsurroundatms_extract
from BDamChangeVsAtomnum import bdamBfac_change_v_atomnum,bdamchange_v_atomnum
from densityAnalysisPlots import numneighbours_scatter
import shutil
 
class eTrack(object):
	# A class for retrieving the eTrack input text file information and running 
	# the eTrack pipeline functions separately or together in full pipeline

	def __init__(self, where = "", pdbname = [], pklfiles = [], initialPDB = "",
				 graph_analysis_topN = 0, graph_analysis_densmet = "",
				 seriesname = "", PDBmultipklname = "",plot = False):

		self.where 					= where
		self.pdbname 				= pdbname
		self.pklfiles 				= pklfiles
		self.initialPDB 			= initialPDB
		self.graph_analysis_topN 	= graph_analysis_topN
		self.graph_analysis_densmet = graph_analysis_densmet
		self.seriesname 			= seriesname
		self.PDBmultipklname 		= PDBmultipklname
		self.plot 					= plot

	def runPipeline(self,map_process,post_process,retrieve_PDBmulti,graphs,inputfilename):
		# the following function reads in the above functions one by one in a scripted 
		# pipeline. Takes inputs to specify which parts of above pipeline are to be
		# included. For each function input specify 'y' if this part is to be performed
		# and 'n' otherwise.
		self.titleCaption('eTrack pipeline')

		# first need to run function above to read in input file containing info
		# on where input files are and where output files should be written
		print 'Reading input file: {}'.format(str(inputfilename))
		self.readInputFile(inputfilename)

		if map_process == 'y':
			self.map_processing()
		else:
			print 'Map processing task not chosen...'
			print 'Input pkl files for post processing chosen from input file:'
			for file in self.pklfiles:
				print '\t{}'.format(file)
		self.fillerLine()

		if post_process == 'y':
			self.post_processing()
			# save PDBmulti as pkl file
			PDBmultipklname = save_objectlist(self.PDBmulti,self.seriesname)
			os.system('mv {} {}{}'.format(PDBmultipklname,self.where +'output/',PDBmultipklname))
			self.PDBmultipklname = self.where +'output/'+ PDBmultipklname
			del self.PDBmulti
		else: 
			print 'Post processing job not chosen...'
		self.fillerLine()

		if retrieve_PDBmulti == 'y':
			self.PDBmulti_retrieve()
			if graphs != 'y':
				return 
		else:
			print 'PDBmulti retrieval from pkl file not chosen...'

		if graphs == 'y':
			self.graph_analysis()
		else:
			print 'Graph analysis job not selected...'
		self.fillerLine()

	def readInputFile(self,inputfilename):
		# this function reads in the input file e_Track_input.txt, which 
		# specifies the location of all the relevant input files and where
		# the output files should be written

		inputfile = open(inputfilename,'r')

		# preallocate pklfiles list for post_processing pipeline step below
		pklfiles = []

		# read in input file and determine location 'where' and name of input files
		for line in inputfile.readlines():
			if '#' == line[0]:
				continue
			elif 'where' in line.split()[0]:
				self.where 					= line.split()[1]
			elif 'damageset_name' in line.split()[0]:
				self.seriesname 			= line.split()[1]
			elif 'damageset_num' in line.split()[0]:
				damageset_nums 				= line.split()[1]
			elif 'PKLFILE' in line.split()[0]:
				pklfiles.append(self.where+'output/'+line.split()[1])
			elif 'initialPDB' in line.split()[0]:
				self.initialPDB 			= line.split()[1]
			elif 'topN' in line.split()[0]:
				self.graph_analysis_topN 	= line.split()[1]
			elif 'densmet' in line.split()[0]:
				self.graph_analysis_densmet = line.split()[1]
			elif 'PKLMULTIFILE' in line.split()[0]:
				self.PDBmultipklname 		= self.where+line.split()[1]
			elif 'plotGraphs' in line.split()[0]:
				self.plot 					= True

		# check that an output directory has been found and make subdirectories if present
		if os.path.isdir(self.where) == False:
			print 'Input file location: {} does not exist. Please select an appropriate directory'.format(self.where)
			return
		self.outputDir 				= '{}output/'.format(self.where)
		self.outputPlotDir 			= '{}plots/'.format(self.outputDir)
		self.outputCombPlotDir 		= '{}combinedDatasets'.format(self.outputPlotDir)

		# add pkl file names as attribute if specified in input file
		if len(pklfiles) != 0:
			self.pklfiles = pklfiles

		inputfile.close()

		# want to find individual damage set numbers from string damageset_nums
		damageset_nums_list = damageset_nums.split(',')

		# next put the dataset name and damage set numbers together in a list
		self.pdbname = [self.seriesname+num for num in damageset_nums_list]

	def map_processing(self):
		self.titleCaption('Map Processing')

		# want to create the following additional subdirectories directories
		where2 = '{}output/'.format(self.where)
		if not os.path.exists(self.outputDir):
		    os.makedirs(self.outputDir)

		if not os.path.exists(self.outputPlotDir):
		    os.makedirs(self.outputPlotDir)

		pklFileNames = []
		for dataset in self.pdbname:
			# derive per-atom density metrics from maps
			mapfilname1 		= '{}_atoms.map'.format(dataset)
			mapfilname2 		= '{}_density.map'.format(dataset)
			maps2DensMets 	= maps2DensMetrics(self.where,dataset,mapfilname1,'atom_map',
											   mapfilname2,'density_map',self.plot)
   			maps2DensMets.maps2atmdensity()

			# move pkl file to working output directory
			pklFileName = maps2DensMets.pklFileName
			os.system('mv {} {}{}'.format(pklFileName,self.outputDir,pklFileName))
			pklFileNames.append(self.outputDir+pklFileName)

		self.pklfiles = pklFileNames

	def post_processing(self):
		self.titleCaption('Post Processing')

		# next read in the pdb structure file as list of atom objects
		print 'Reading in initial pdb file...'
		initialPDBlist = pdb2list(self.where+self.initialPDB,[])

		# retrieve object lists of atoms for each damage set
		self.fillerLine()
		print 'Reading in damaged pkl files...'
		data_list = []
		for pkl_filename in self.pklfiles:
			print 'Damage file number: {}\n'.format(len(data_list)+1)
			PDB_ret = retrieve_objectlist(pkl_filename)

			# add new retrieved damage set list to data_list
			data_list.append(PDB_ret)

		# create a list of atom objects with attributes as lists varying over 
		# dose range, only including atoms present in ALL damage datasets
		print 'New list of atoms over full dose range calculated...'
		self.PDBmulti = multiARRAY_diffatomnumbers(data_list)

		# determine Bfactor change between later datasets and initial 
		find_Bchange(initialPDBlist,self.PDBmulti,'Bfactor')

		# write atom numbers and density metrics to simple text files - one for 
		# each density metric separately
		for densMet in ('mean','min','max','median','std','rsd','mode',
						'min90tile','max90tile','min95tile','max95tile'):
			print 'Writing .txt file for per-atom density metric: {}'.format(densMet)
			objlist2txt(self.PDBmulti,self.outputDir,densMet)

	def PDBmulti_retrieve(self):
		self.titleCaption('Atom List Retrieval')
		# retrieve the PDBmulti list from the pickle list
		self.PDBmulti = retrieve_objectlist(self.PDBmultipklname)

	def graph_analysis(self):
		self.titleCaption('Graph Analysis')
		
		# location of output plots (make folder if doesn't exist)
		if not os.path.exists(self.outputCombPlotDir):
		    os.makedirs(self.outputCombPlotDir)

		# extract number of top damage sites N and density metric
		N 		= self.graph_analysis_topN
		densmet = self.graph_analysis_densmet

		# want to know name of initial pdb here as template file to make new N
		# damage sites file
		initialPDBname = self.where + self.initialPDB

		# determine top N damage sites:
		for normType in ('NOTnormalised','normalised'):
			topNdamsites_resibarplotter(self.PDBmulti,N,self.outputCombPlotDir,densmet,normType)
		topNdamsites_chainbarplotter(self.PDBmulti,N,self.outputCombPlotDir,densmet)
		topNdamsites_printer(self.PDBmulti,N,self.outputCombPlotDir,densmet,initialPDBname)


	def fillerLine(self):
		print '---------------------------------------------------------------'	

	def titleCaption(self,title):
		print '||========================== {} ==========================||'.format(title)


