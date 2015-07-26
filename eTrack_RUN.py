# -*- coding: utf-8 -*-

from convertMaps2PkldObjLists import maps2pkldobjs
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
				 seriesname = "", PDBmultipklname = ""):

		self.where 					= where
		self.pdbname 				= pdbname
		self.pklfiles 				= pklfiles
		self.initialPDB 			= initialPDB
		self.graph_analysis_topN 	= graph_analysis_topN
		self.graph_analysis_densmet = graph_analysis_densmet
		self.seriesname 			= seriesname
		self.PDBmultipklname 		= PDBmultipklname

	def readInputFile(self,inputfilename):
		# this function reads in the input file e_Track_input.txt, which 
		# specifies the location of all the relevant input files and where
		# the output files should be written

		inputfile = open(inputfilename,'r')

		# preallocate pklfiles list for post_processing pipeline step below
		pklfiles = []
		# preallocate initial pdb file name 
		initialPDB = ''

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
				pklfiles.append(line.split()[1])
			elif 'initialPDB' in line.split()[0]:
				self.initialPDB 			= line.split()[1]
			elif 'topN' in line.split()[0]:
				self.graph_analysis_topN 	= line.split()[1]
			elif 'densmet' in line.split()[0]:
				self.graph_analysis_densmet = line.split()[1]
			elif 'PKLMULTIFILE' in line.split()[0]:
				self.PDBmultipklname 		= line.split()[1]

		# add pkl file names as attribute if specified in input file
		if len(pklfiles) != 0:
			self.pklfiles = pklfiles

		inputfile.close()

		# want to find individual damage set numbers from string damageset_nums
		damageset_nums_list = damageset_nums.split(',')

		# next put the dataset name and damage set numbers together in a list
		self.pdbname = [self.seriesname+num for num in damageset_nums_list]

	def map_processing(self):

		# want to create the following additional subdirectories directories
		where2 = self.where + 'output/'
		if not os.path.exists(where2):
		    os.makedirs(where2)

		where3 = where2 + 'plots/'
		if not os.path.exists(where3):
		    os.makedirs(where3)

		pkl_filenames = []
		for dataset in self.pdbname:
			mapfilname1 	= dataset + '_atoms.map'
			mapfilname2 	= dataset + '_density.map'
			pkl_filename 	= maps2pkldobjs(self.where,dataset,mapfilname1,mapfilname2)
			pkl_filenames.append(pkl_filename)

		self.pklfiles = pkl_filenames

	def post_processing(self):

		print '••••••••••••••••••••••••••••••'
		print 'Reading in initial pdb file...'
		# next read in the pdb structure file:
		# run function to fill PDBarray list with atom objects from structure

		initialPDBlist = pdb2list(self.where+self.initialPDB,[])

		# determine the number of surrounding atoms for each atom in structure
		numsurroundatoms_calculate(self.where+self.initialPDB,initialPDBlist,10)

		# plot a scatter plot of # neighbouring atoms vs protons
		numneighbours_scatter(self.where,initialPDBlist,self.seriesname)

		# determine Bdamage metric for initial PDB structure
		bdamage_calculate(initialPDBlist)

		# retrieve object lists of atoms for each damage set
		print '•••••••••••••••••••••••••••••••'
		print 'Reading in damaged pkl files...'
		data_list = []
		for pkl_filename in self.pklfiles:
			print '\nDamage file number {}:'.format(len(data_list)+1)
			PDB_ret = retrieve_objectlist(pkl_filename)

			# extract number of surrounding atoms for each atom in later structure
			# from initial structure
			numsurroundatms_extract(initialPDBlist,PDB_ret)

			# determine Bdamage metric for later PDB structures
			bdamage_calculate(PDB_ret)

			# add new retrieved damage set list to data_list
			data_list.append(PDB_ret)

		# create a list of atom objects with attributes as lists varying over 
		# dose range, only including atoms present in ALL damage datasets
		self.PDBmulti = multiARRAY_diffatomnumbers(data_list)

		# determine Bfactor change between later datasets and initial 
		find_Bchange(initialPDBlist,self.PDBmulti,'Bfactor')

		# determine Bdamage change between later datasets and initial 
		find_Bchange(initialPDBlist,self.PDBmulti,'Bdamage')

		# write atom numbers and density metrics to simple text files - one for 
		# each density metric separately
		where = self.where + 'output/'	    
		objlist2txt(self.PDBmulti,where,'mean')
		objlist2txt(self.PDBmulti,where,'min')
		objlist2txt(self.PDBmulti,where,'max')
		objlist2txt(self.PDBmulti,where,'median')
		objlist2txt(self.PDBmulti,where,'std')
		objlist2txt(self.PDBmulti,where,'rsd')
		objlist2txt(self.PDBmulti,where,'mode')
		objlist2txt(self.PDBmulti,where,'min90tile')
		objlist2txt(self.PDBmulti,where,'max90tile')
		objlist2txt(self.PDBmulti,where,'min95tile')
		objlist2txt(self.PDBmulti,where,'max95tile')

	def PDBmulti_retrieve(self):
		# retrieve the PDBmulti list from the pickle list
		self.PDBmulti = retrieve_objectlist(self.PDBmultipklname)

	def graph_analysis(self):
		
		# location of output plots (make folder if doesn't exist)
		where = self.where + 'output/plots/combineddatasets/'
		if not os.path.exists(where):
		    os.makedirs(where)

		# extract number of top damage sites N and density metric
		N 		= self.graph_analysis_topN
		densmet = self.graph_analysis_densmet

		# want to know name of initial pdb here as template file to make new N
		# damage sites file
		initialPDBname = self.where + self.initialPDB

		# determine top N damage sites:
		topNdamsites_resibarplotter(self.PDBmulti,N,where,densmet,'NOTnormalised')
		topNdamsites_resibarplotter(self.PDBmulti,N,where,densmet,'normalised')
		topNdamsites_chainbarplotter(self.PDBmulti,N,where,densmet)
		topNdamsites_printer(self.PDBmulti,N,where,densmet,initialPDBname)

		# plot graphs of Bfactor change and Bdamage change against atom number
		bdamchange_v_atomnum(where,self.PDBmulti)
		bdamBfac_change_v_atomnum(where,self.PDBmulti)

	def runPipeline(self,map_process,post_process,retrieve_PDBmulti,graphs,inputfilename):
		# the following function reads in the above functions one by one in a scripted 
		# pipeline. Takes inputs to specify which parts of above pipeline are to be
		# included. For each function input specify 'y' if this part is to be performed
		# and 'n' otherwise.

		print '||===========================================================||'
		print '||                  eTrack pipeline                          ||'
		print '||===========================================================||'

		# first need to run function above to read in input file containing info
		# on where input files are and where output files should be written
		print 'Reading input file: {}'.format(str(inputfilename))
		self.readInputFile(inputfilename)

		if map_process == 'y':
			print '\n||===========================================================||'
			print '||                  	  Map Processing                        ||'
			print '||===========================================================||'
			self.map_processing()
		else:
			print 'Map processing task not chosen...'
			print 'Input pkl files for post processing chosen from input file:'
			for file in self.pklfiles:
				print '\t{}'.format(file)
		print '---------------------------------------------------------------'	

		if post_process == 'y':
			print '\n||===========================================================||'
			print '||                  	  Post Processing                       ||'
			print '||===========================================================||'
			self.post_processing()

			# save PDBmulti as pkl file
			self.PDBmultipklname = save_objectlist(self.PDBmulti,self.seriesname)
			del self.PDBmulti
		else: 
			print 'Post processing job not chosen...'
		print '---------------------------------------------------------------'	

		if retrieve_PDBmulti == 'y':
			print '\n||===========================================================||'
			print '||                  	  PDBmulti retrieval                    ||'
			print '||===========================================================||'
			self.PDBmulti_retrieve()
			if graphs != 'y':
				return
		else:
			print 'PDBmulti retrieval from pkl file not chosen...'

		if graphs == 'y':
			print '\n||===========================================================||'
			print '||                  	  Graph Analysis                        ||'
			print '||===========================================================||'
			self.graph_analysis()
		else:
			print 'Graph analysis job not selected...'
		print '---------------------------------------------------------------'	
