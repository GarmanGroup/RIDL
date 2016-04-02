# -*- coding: utf-8 -*-

from readAtomMap import maps2DensMetrics
from savevariables import retrieve_objectlist,save_objectlist,saveGenericObject,retrieveGenericObject
from PDBFileManipulation import PDBtoList,writePDBline
from combinedAtomList import combinedAtomList
import os
from time import gmtime, strftime
import numpy as np

class eTrack(object):
	# A class for retrieving the eTrack input text file information and running 
	# the eTrack pipeline functions separately or together in full pipeline
	def __init__(self, inDir = "",outDir = "",pdbNames = [], pklFiles = [], initialPDB = "",
				 seriesName = "untitled-series", pklSeries = "",doses=[],plot = False):

		self.inDir 			= inDir 		# the input file directory
		self.outDir  		= outDir 		# the output file directory
		self.pdbNames 		= pdbNames 		# the list of pdb codes for series
		self.pklFiles 		= pklFiles 		# list of pkl files from map_processing
		self.initialPDB 	= initialPDB 	# the first dataset pdb code
		self.seriesName 	= seriesName 	# the general series name 
		self.pklSeries		= pklSeries 	# the combined series pkl file from post_processing
		self.doses 			= doses 		# list of increasing doses
		self.plot           = plot 			# (boolian) decide to plot per-residue summary plots per dataset

	def runPipeline(self,map_process=True,post_process=True,retrieve_PDBmulti=True,inputFileName=''):
		# the following function reads in the above functions one by one in a scripted 
		# pipeline. Takes inputs to specify which parts of above pipeline are to be
		# included. For each function input specify True if this part is to be performed
		# and False otherwise.
		self.titleCaption(title='ETRACK pipeline')
		self.inputFileName = inputFileName

		# check whether valid inputs to function
		valid = self.checkValidInputs(map_process,post_process,retrieve_PDBmulti)
		if valid is False: 
			return

		# first need to run function above to read in input file containing info
		# on where input files are and where output files should be written
		print 'Reading input file: {}'.format(self.inputFileName)
		self.readInputFile()
		success = self.checkInOutDirExist()
		if success is False: 
			return
		self.makeOutputDirs()

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
			# provide summary txt file on Dloss metric metric per-dataset
			self.summaryFile(normType='Standard') 
			self.summaryFile(normType='Calpha normalised') 
			self.writeDamSitesToFile()
		else: 
			print 'Post processing job not chosen...'
		self.fillerLine()
		if retrieve_PDBmulti is True:
			self.PDBmulti_retrieve()
		else:
			print 'PDBmulti retrieval from pkl file not chosen...'
		self.fillerLine()

	def checkValidInputs(self,map_process,post_process,retrieve_PDBmulti):
		# check the runPipeline inputs to make sure that they are valid
		for var in [map_process,post_process,retrieve_PDBmulti]:
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

	def readInputFile(self):
		# read input file e_Track_input.txt to specify location of input files 
		# and where to write output files

		props = {'inDir':'inDir',
				 'outDir':'outDir',
				 'damageset_name':'seriesName',
				 'initialPDB':'initialPDB',
				 'PKLMULTIFILE':'pklSeries'}

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

			# # elif 'inDir' in l[0]:
			# # 	self.inDir 		= l[1]
			# # elif 'outDir' in l[0]:
			# # 	self.outDir 	= l[1]
			# # elif 'damageset_name' in l[0]:
			# # 	self.seriesName = l[1]
			# elif 'damageset_num' in l[0]:
			# 	datasetNums 	= l[1]
			# elif 'PKLFILE' in l[0]:
			# 	self.pklFiles.append(l[1])
			# # elif 'initialPDB' in l[0]:
			# # 	self.initialPDB = l[1]
			# # elif 'PKLMULTIFILE' in l[0]:
			# # 	self.pklSeries 	= l[1]
			# elif 'doses' in l[0]:
			# 	self.doses 		= [float(d) for d in l[1].split(',')]
			# elif 'plot' in l[0]:
			# 	self.plot = True
		inputfile.close()
		# next put the dataset name and damage set numbers together in a list
		self.pdbNames = [self.seriesName+num for num in datasetNums.split(',')]

	def checkInOutDirExist(self):
		# check that an input/output directories have been found and make subdirectories if present
		for dir in ([[self.inDir,'Input'],[self.outDir,'Output']]):
			if os.path.isdir(dir[0]) == False:
				print '{} file location: {} does not exist. Please select an appropriate directory'.format(dir[1],dir[0])
				return False
		return True

	def makeOutputDirs(self):
		self.outputDir 		= '{}ETRACK_output/'.format(self.outDir)
		self.outputPlotDir 	= '{}plots/'.format(self.outputDir)

		# add pkl file names as attribute if specified in input file
		if len(self.pklFiles) != 0: 
			self.pklFiles = [self.outputDir+f for f in self.pklFiles]

	def map_processing(self):
		self.titleCaption(title='Map Processing')

		# create additional subdirectories
		for oDir in (self.outputDir,self.outputPlotDir):
			if not os.path.exists(oDir):
				os.makedirs(oDir)

		# make pklFiles and dir to move all generated per-dataset pkl files to this
		pklFileDir = 'pklFiles-perDataset'
		os.system('mkdir {}{}'.format(self.outputDir,pklFileDir))

		pklFileNames = []
		for dataset in self.pdbNames:
			# derive per-atom density metrics from maps
			mapName1 	= '{}_atoms.map'.format(dataset)
			mapName2 	= '{}_density.map'.format(dataset)
			maps2DensMets 	= maps2DensMetrics(filesIn=self.inDir,filesOut=self.outputDir,
											   pdbName=dataset,mapFileName1=mapName1,
											   mapFileName2=mapName2,plot=self.plot)
   			maps2DensMets.maps2atmdensity()

			# move pkl file to working output directory
			pklFileName = maps2DensMets.pklFileName
			os.system('mv {} {}{}/{}'.format(pklFileName,self.outputDir,pklFileDir,pklFileName))
			pklFileNames.append('{}{}/{}'.format(self.outputDir,pklFileDir,pklFileName))

		self.pklFiles = pklFileNames

	def post_processing(self):
		self.titleCaption(title='Post Processing')
		print 'Input pkl files for post processing chosen from input file:'
		for file in self.pklFiles: 
			print '\t{}'.format(file)

		# next read in the pdb structure file as list of atom objects
		print 'Reading in initial pdb file...'
		initialPDBlist = PDBtoList(self.inDir+self.initialPDB,[])

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

		# calculate 'standardised' and 'average' variant Dloss metrics and Calpha
		# normalised metrics
		combinedAtoms.calcStandardisedMetrics('loss')
		combinedAtoms.calcAdditionalMetrics('loss','Standard','average')
		for m in ('loss','mean','gain','Bfactor'):
			combinedAtoms.calcAdditionalMetrics(m,'Standard','Calpha')

		# write atom numbers and density metrics to simple text files - one for 
		# each density metric separately
		for densMet in combinedAtoms.getDensMetrics():
			print 'Writing .csv file for per-atom density metric: {}, normalisation: {}'.format(*densMet)
			combinedAtoms.writeMetric2File(self.outputDir,'none',*densMet)
		for groupBy in ('residue','atomtype'):
			combinedAtoms.writeMetric2File(self.outputDir,groupBy,'loss','Standard')
			combinedAtoms.writeMetric2File(self.outputDir,groupBy,'loss','Calpha normalised')

		# make csvFiles dir and move all generated csv files to this
		os.system('mkdir {}/csvFiles'.format(self.outputDir))
		for file in os.listdir(self.outputDir):
			if file.endswith(".csv"):
				os.system('mv -f {}{} {}csvFiles/{}'.format(self.outputDir,file,self.outputDir,file))

		self.combinedAtoms = combinedAtoms
		#self.sensAtomPlots()

	def PDBmulti_retrieve(self):
		self.titleCaption(title='Atom List Retrieval')
		print 'Input pkl file for data retrieval chosen from input file:'
		print '\t{}'.format(self.outputDir+self.pklSeries)
		# retrieve the combinedAtoms object from the pkl file
		self.combinedAtoms = retrieveGenericObject(self.outputDir+self.pklSeries)

	def summaryFile(self,fileType='html',metric='loss',normType='Standard'):
		if fileType == 'html':
			self.summaryHTML(metric=metric,normType=normType)
		else:
			self.summaryTxt(metric=metric,normType=normType)

	def summaryTxt(self,metric='loss',normType='Standard'):
		# produce a selection of per-dataset summary statistics
		# by default set metric as 'loss' for Dloss metric
		# by default 'normType' is 'Standard'

		numDsets = self.getNumDatasets()
		summaryFile = open('{}/summaryFile-D{}-{}.txt'.format(self.outputDir,metric,normType.replace(' ','-')),'w')
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
			n = self.combinedAtoms.numAtmsWithMetricAboveLevel(i,metric,normType,0.5,False)
			summaryFile.write('{}\t\t{}\n'.format(0.5,n))
			t = 0
			while n != 0:
				t += 1
				n = self.combinedAtoms.numAtmsWithMetricAboveLevel(i,metric,normType,t,False)
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
			n = self.combinedAtoms.getNumAtoms()*0.1 # take top 10% of atoms here
			statsOut = self.combinedAtoms.breakdownTopNatomsBy(metric,normType,i,n,['basetype','atomtype'])
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

	def summaryHTML(self,metric='loss',normType='Standard'):
		# produce a selection of per-dataset summary statistics

		numDsets = self.getNumDatasets()
		summaryFile = open('{}/summaryFile-D{}-{}.html'.format(self.outputDir,metric,normType.replace(' ','-')),'w')
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
					 '<li>ratio: net ratio of distn values either side of metric distn mode</li>\n</ul>'
		summaryFile.write(bodyString)

		av,std = self.combinedAtoms.getAverageMetricVals(metric,normType) # get structure-wide metric average & std dev
		for i in range(numDsets):
			bodyString = '<h3>Dataset info</h3>\n'+\
						 'Number in series: {}<br>\n'.format(i+1)+\
						 'Dose (MGy): {}<br>\n'.format(self.doses[i])+\
						 '<h3>Structure-wide D{} summary</h3>\n'.format(metric)+\
						 'Average D{}: {}<br>\n'.format(metric,round(av[i],3))+\
						 'Std dev in D{}: {}<br>\n'.format(metric,round(std[i],3))+\
						 '# atoms with {} Dloss metric above N std dev of structure-wide mean:<br><br>\n'.format(normType)
			summaryFile.write(bodyString)

			n = self.combinedAtoms.numAtmsWithMetricAboveLevel(i,metric,normType,0.5,False)
			tableString = '<table border="1" style="width:10%">\n'+\
						  '<tr><th>N</th><th>#atoms</th></tr>\n'+\
						  '<tr><td>{}</td><td>{}</td></tr>\n'.format(0.5,n)
			t = 0
			while n != 0:
				t += 1
				n = self.combinedAtoms.numAtmsWithMetricAboveLevel(i,metric,normType,t,False)
				tableString += '<tr><td>{}</td><td>{}</td></tr>\n'.format(t,n)
			tableString += '</table>\n'
			summaryFile.write(tableString)

			statsOut = self.combinedAtoms.getTopNAtomsString(metric,normType,i,25)
			infoString = '<h3>Per-atom Statistics</h3>\n'+\
						 'Top hits ranked by D{} metric:<br><br>\n'.format(metric)
			infoString += self.convertPlainTxtTable2html(statsOut)
			summaryFile.write(infoString)

			statsOut = self.combinedAtoms.getPerAtmtypeStats(metric,normType,i,'mean',25)
			infoString = '<h3>Per-atom-type Statistics</h3>\n'+\
						 'D{} metric ranked by mean value:<br><br>\n'.format(metric)
			infoString += self.convertPlainTxtTable2html(statsOut[0],width='60%')
			summaryFile.write(infoString)

			n = self.combinedAtoms.getNumAtoms()*0.1 # take top 10% of atoms here
			statsOut = self.combinedAtoms.breakdownTopNatomsBy(metric,normType,i,n,['basetype','atomtype'])
			infoString = '<h3>Per-atom Statistics</h3>\n'+\
						 '{}<br><br>\n'.format(statsOut[0])
			infoString += self.convertPlainTxtTable2html(statsOut[1])
			summaryFile.write(infoString)

			statsOut = self.combinedAtoms.getPerResidueStats(metric,normType,i,'mean',25)
			infoString = '<h3>Per-residue Statistics</h3>\n'+\
						 'D{} metric ranked by mean value:<br><br>\n'.format(metric)
			infoString += self.convertPlainTxtTable2html(statsOut[0],width='60%')
			summaryFile.write(infoString)

			statsOut = self.combinedAtoms.getPerChainStats(metric,normType,i,'mean','all')
			infoString = '<h3>Per-chain Statistics</h3>\n'+\
						 'D{} metric ranked by mean value:<br><br>\n'.format(metric)
			infoString += self.convertPlainTxtTable2html(statsOut[0],width='60%')
			summaryFile.write(infoString)

			self.makeDistnPlots(densMet=metric,normType=normType,set=1)
			plotString = '<h3>Distribution of Dloss for known susceptible residue types</h3>\n'+\
						 '<img src="GLUASPCYSMETTYR_{}D{}-kde-{}.png">'.format(normType.replace(' ',''),metric,i)
			summaryFile.write(plotString)

		summaryFile.write('</body>\n</html>')
		summaryFile.close()

	def convertPlainTxtTable2html(self,plainText,numHeaderLines=1,width='40%',border='1'):
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
		numDsets = self.combinedAtoms.atomList[0].getNumDatasets() # get number of datasets in damage series
		return numDsets

	def writeDamSitesToFile(self,metric='loss',normType='Standard',numDamSites=25):
		# write top damage sites to .pdb file for each dataset
		self.damSitesPDB = []
		numDsets = self.combinedAtoms.atomList[0].getNumDatasets() # get number of datasets in damage series
		for i in range(numDsets): 
			damPDB = self.combinedAtoms.getTopNAtomsPDBfile(metric,normType,i,numDamSites,self.inDir+self.initialPDB)
			self.damSitesPDB.append(damPDB)

	def colorByMetric(self,metric='loss',normType='Standard',dataset=0,singleResidue=''):
		# for the initial pdb file in damage series, convert the Bfactor column 
		# to values for the specified metric.
		# If 'singleResidue' is specified (use 3-letter residue code) then the average 
		# density metric for each atom within this residue is calculated within the structure
		# and only this is output to resulting PDB file, otherwise use ''
		if normType == 'Calpha normalised': 
			self.combinedAtoms.calcAdditionalMetrics(metric,normType,'Calpha')

		pdbIn = open(self.inDir+self.initialPDB,'r')
		fileOut = self.outputDir+self.initialPDB.strip('.pdb')+'_{}D{}_{}.pdb'.format(normType.replace(" ",""),metric,dataset)
		if singleResidue != '':
			fileOut = fileOut.strip('.pdb')+'-{}.pdb'.format(singleResidue)

		pdbOut = open(fileOut,'w')
		pdbOut.write('REMARK\tBfactor column replaced by {} D{} metric values\n'.format(normType,metric))
		for l in pdbIn.readlines():
			if l.split()[0] in ('CRYST1','SCALE1','SCALE2','SCALE3'):
				pdbOut.write(l)
			elif 'ATOM' in l.split()[0]:
				break
		pdbIn.close()

		if singleResidue == '':
			for atm in self.combinedAtoms.atomList:
				dens = atm.densMetric[metric][normType]['values'][dataset]
				if not np.isnan(dens): # don't include atoms for which not calculated properly
					l = writePDBline(atm,dens)
					pdbOut.write(l+'\n')
		else:
			atmDic = self.combinedAtoms.getAvMetricPerAtmInRes(singleResidue,metric,normType,dataset)
			for atm in self.combinedAtoms.atomList:
				if atm.basetype == singleResidue:
					resNum = atm.residuenum # save the residue number now
					break
			for atm in self.combinedAtoms.atomList:
				if atm.residuenum == resNum:
					dens = atmDic[atm.atomtype]
					l = writePDBline(atm,dens)
					pdbOut.write(l+'\n')
		pdbOut.write('END')
		pdbOut.close()

	def visualiseDamSites(self,dataset=0,metric='loss',software='pymol',size=1):
		# open coot/pymol to view top damage sites
		# size is the density metric scale factor used when visualising damage sites as spheres 
		# of vdw = size*{density metric} within pymol
		if software not in ('coot','pymol'):
			print "Damage sites can be visualised in either 'coot' or 'pymol"
			print "Please make sure the paths to these programs are correctly set up before proceeding"
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
			pymolScript.write('load {}\n'.format(self.inDir+self.initialPDB)+\
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
							  'color white, {}\n'.format(structureTag))
			pymolScript.close()
			os.system('pymol {}'.format(scriptName))

	def makeDistnPlots(self,densMet='loss',normType='Standard',plotType='both',set=1):
		# retrieve the distribution for damage metric values for all atoms of 
		# specified residues (see 'resList' below), and then plot a kde plot for each
		if set == 1:
			resList = [['GLU','ASP','CYS','MET','TYR']]
		elif set == 2:
			resList = [['GLU','GLN'],['ASP','ASN'],['ILE','LEU'],['TYR','PHE'],
				   	   ['TYR','ASP','GLU'],['TYR','PHE','GLY'],
				   	   ['GLU','GLY'],['ASP','GLY'],['TYR','GLY'],['PHE','GLY'],
				   	   ['CYS','GLY'],['MET','GLY'],['GLU','ASP','CYS','MET','TYR']]
		for resGroup in resList:
			if len(resGroup) > 4:
				plotType = 'kde'
			else:
				plotType = 'both'
			for i in range(self.getNumDatasets()):
				data = self.combinedAtoms.graphMetricDistn(metric=densMet,normType=normType,
														   valType=i,plotType=plotType,resiType=resGroup)	

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
				self.combinedAtoms.plotSusceptibleAtoms(densMet=metric,susAtms=aTypes)
				self.combinedAtoms.plotSusceptibleAtoms(densMet=metric,errorbars=False,susAtms=aTypes)

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

	def fillerLine(self):
		print '---------------------------------------------------------------'	

	def titleCaption(self,title='unspecified title'):
		print '||========================== {} ==========================||'.format(title)
