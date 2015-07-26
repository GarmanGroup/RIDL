# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 22:23:09 2015

@author: charlie
"""
from convertmaps2PkldObjLists import maps2pkldobjs
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
 
########################################################################
#                      master file for eTrack                          #
########################################################################

class inputfile_info:
	# A class for input text file
    def __init__(self,where="",pdbname=[],pklfiles=[],initialPDB="",
    			 graph_analysis_topN=0,graph_analysis_densmet="",
    			 seriesname="",PDBmultipklname=""):

		self.where = where
		self.pdbname = pdbname
		self.pklfiles = pklfiles
		self.initialPDB = initialPDB
		self.graph_analysis_topN = graph_analysis_topN
		self.graph_analysis_densmet = graph_analysis_densmet
		self.seriesname = seriesname
		self.PDBmultipklname = PDBmultipklname


def read_inputfile(inputfilename):
	# this function reads in the input file e_Track_input.txt, which 
	# specifies the location of all the relevant input files and where
	# the output files should be written

	inputfile = open(inputfilename,'r')

	# use inputfile class above to define inputfile object
	inputfileinfo = inputfile_info()

	# preallocate pklfiles list for post_processing pipeline step below
	pklfiles = []
	# preallocate initial pdb file name 
	initialPDB = ''

	# read in input file and determine location 'where' and name of input files
	for line in inputfile.readlines():
		if '#' == line[0]:
			continue
		elif 'where' in line.split()[0]:
			inputfileinfo.where = line.split()[1]
		elif 'damageset_name' in line.split()[0]:
			inputfileinfo.seriesname = line.split()[1]
		elif 'damageset_num' in line.split()[0]:
			damageset_nums = line.split()[1]
		elif 'PKLFILE' in line.split()[0]:
			pklfiles.append(line.split()[1])
		elif 'initialPDB' in line.split()[0]:
			inputfileinfo.initialPDB = line.split()[1]
		elif 'topN' in line.split()[0]:
			inputfileinfo.graph_analysis_topN = line.split()[1]
		elif 'densmet' in line.split()[0]:
			inputfileinfo.graph_analysis_densmet = line.split()[1]
		elif 'PKLMULTIFILE' in line.split()[0]:
			inputfileinfo.PDBmultipklname = line.split()[1]

	# add pkl file names as attribute if specified in input file
	if len(pklfiles) != 0:
		inputfileinfo.pklfiles = pklfiles

	inputfile.close()

	# want to find individual damage set numbers from string damageset_nums
	damageset_nums_list = damageset_nums.split(',')

	# next put the dataset name and damage set numbers together in a list
	inputfileinfo.pdbname = [inputfileinfo.seriesname+num for num in damageset_nums_list]

	return inputfileinfo

def map_processing(inputfileinfo):

	# want to create the following additional subdirectories directories
	where2 = inputfileinfo.where+'output/'
	if not os.path.exists(where2):
	    os.makedirs(where2)

	where3 = where2+'plots/'
	if not os.path.exists(where3):
	    os.makedirs(where3)

	pkl_filenames = []
	for dataset in inputfileinfo.pdbname:
		mapfilname1 = dataset+'_atoms.map'
		mapfilname2 = dataset+'_density.map'
		pkl_filename = maps2pkldobjs(inputfileinfo.where,dataset,mapfilname1,mapfilname2)
		pkl_filenames.append(pkl_filename)

	inputfileinfo.pklfiles = pkl_filenames


def post_processing(inputfileinfo):

	print '••••••••••••••••••••••••••••••'
	print 'Reading in initial pdb file...'
	# next read in the pdb structure file:
	# run function to fill PDBarray list with atom objects from structure

	initialPDBlist = pdb2list(inputfileinfo.where+inputfileinfo.initialPDB,[])

	# determine the number of surrounding atoms for each atom in structure
	numsurroundatoms_calculate(inputfileinfo.where+inputfileinfo.initialPDB,initialPDBlist,10)

	# plot a scatter plot of # neighbouring atoms vs protons
	numneighbours_scatter(inputfileinfo.where,initialPDBlist,inputfileinfo.seriesname)

	# determine Bdamage metric for initial PDB structure
	bdamage_calculate(initialPDBlist)

	# retrieve object lists of atoms for each damage set
	print '•••••••••••••••••••••••••••••••'
	print 'Reading in damaged pkl files...'
	data_list = []
	for pkl_filename in inputfileinfo.pklfiles:
		print '\nDamage file number %s:' %(len(data_list)+1)
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
	PDBmulti = multiARRAY_diffatomnumbers(data_list)

	# determine Bfactor change between later datasets and initial 
	find_Bchange(initialPDBlist,PDBmulti,'Bfactor')

	# determine Bdamage change between later datasets and initial 
	find_Bchange(initialPDBlist,PDBmulti,'Bdamage')

	# write atom numbers and density metrics to simple text files - one for 
	# each density metric separately
	where = inputfileinfo.where+'output/'	    
	objlist2txt(PDBmulti,where,'mean')
	objlist2txt(PDBmulti,where,'min')
	objlist2txt(PDBmulti,where,'max')
	objlist2txt(PDBmulti,where,'median')
	objlist2txt(PDBmulti,where,'std')
	objlist2txt(PDBmulti,where,'rsd')
	objlist2txt(PDBmulti,where,'mode')
	objlist2txt(PDBmulti,where,'min90tile')
	objlist2txt(PDBmulti,where,'max90tile')
	objlist2txt(PDBmulti,where,'min95tile')
	objlist2txt(PDBmulti,where,'max95tile')
	objlist2txt(PDBmulti,where,'dipstat')
	# return the list of merged objects PDBmulti
	return PDBmulti


def PDBmulti_retrieve(inputfileinfo):
	# retrieve the PDBmulti list from the pickle list
	PDBmulti = retrieve_objectlist(inputfileinfo.PDBmultipklname)
	return PDBmulti


def graph_analysis(inputfileinfo,PDBmulti):

	# location of output plots (make folder if doesn't exist)
	where = inputfileinfo.where+'output/plots/combineddatasets/'
	if not os.path.exists(where):
	    os.makedirs(where)

	# extract number of top damage sites N and density metric
	N = inputfileinfo.graph_analysis_topN
	densmet = inputfileinfo.graph_analysis_densmet

	# want to know name of initial pdb here as template file to make new N
	# damage sites file
	initialPDBname = inputfileinfo.where+inputfileinfo.initialPDB

	# determine top N damage sites:
	topNdamsites_resibarplotter(PDBmulti,N,where,densmet,'NOTnormalised')
	topNdamsites_resibarplotter(PDBmulti,N,where,densmet,'normalised')
	topNdamsites_chainbarplotter(PDBmulti,N,where,densmet)
	topNdamsites_printer(PDBmulti,N,where,densmet,initialPDBname)

	# plot graphs of Bfactor change and Bdamage change against atom number
	bdamchange_v_atomnum(where,PDBmulti)
	bdamBfac_change_v_atomnum(where,PDBmulti)


def e_Track_pipeline(map_process,post_process,retrieve_PDBmulti,graphs,inputfilename):
	# the following function reads in the above functions one by one in a scripted 
	# pipeline. Takes inputs to specify which parts of above pipeline are to be
	# included. For each function input specify 'y' if this part is to be performed
	# and 'n' otherwise.

	print '||===========================================================||'
	print '||                  eTrack pipeline                          ||'
	print '||===========================================================||'

	# first need to run function above to read in input file containing info
	# on where input files are and where output files should be written
	print 'Reading input file: '+ str(inputfilename)
	inputfileinfo = read_inputfile(inputfilename)

	# # want to make a log file for the current eTrack run
	# logfile = open(where+'/logfile.txt','w')
	# logfile.write('Damage series: '+inputfileinfo.seriesname+'\n')


	if map_process == 'y':
		print '\n||===========================================================||'
		print '||                  	  Map Processing                        ||'
		print '||===========================================================||'
		map_processing(inputfileinfo)
	else:
		print 'Map processing task not chosen...'
		print 'Input pkl files for post processing chosen from input file:'
		for file in inputfileinfo.pklfiles:
			print '\t%s' %(file)
	print '---------------------------------------------------------------'	


	if post_process == 'y':
		print '\n||===========================================================||'
		print '||                  	  Post Processing                       ||'
		print '||===========================================================||'
		PDBmulti = post_processing(inputfileinfo)

		# save PDBmulti as pkl file
		inputfileinfo.PDBmultipklname = save_objectlist(PDBmulti,inputfileinfo.seriesname)
		del PDBmulti
	else: 
		print 'Post processing job not chosen...'
	print '---------------------------------------------------------------'	


	if retrieve_PDBmulti == 'y':
		print '\n||===========================================================||'
		print '||                  	  PDBmulti retrieval                    ||'
		print '||===========================================================||'
		PDBmulti = PDBmulti_retrieve(inputfileinfo)
		if graphs != 'y':
			return PDBmulti
	else:
		print 'PDBmulti retrieval from pkl file not chosen...'


	if graphs == 'y':
		print '\n||===========================================================||'
		print '||                  	  Graph Analysis                        ||'
		print '||===========================================================||'
		graph_analysis(inputfileinfo,PDBmulti)
		return PDBmulti
	else:
		print 'Graph analysis job not selected...'
	print '---------------------------------------------------------------'	
