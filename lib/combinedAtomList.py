from confidenceIntervalCalculator import mean_confidence_interval
from PDBFileManipulation import writePDBline_DamSite
from findMetricChange import findBchange
from halfDoseCalc import halfDoseApprox
from combinedAtom import combinedAtom
from CalphaWeight import CalphaWeight
import matplotlib.pyplot as plt
from matplotlib import cm,gridspec 
from scipy import stats
import pandas as pd
import string
import scipy.stats
import numpy as np
import operator
import os

from checkSeabornPresent import checkSeabornPresent as checkSns
seabornFound = checkSns()
if seabornFound is True:
	import seaborn as sns

class combinedAtomList(object):

	# class for list of atom objects defined by combinedAtom class

	def __init__(self,
				 datasetList       = [],
				 numLigRegDsets    = 1,
				 doseList          = [],
				 initialPDBList    = [],
				 outputDir         = './',
				 partialDatasets   = True,
				 seriesName        = 'untitled'):

		self.datasetList 		= datasetList
		self.numLigRegDatasets 	= numLigRegDsets    # number of datasets to perform linear regression over
		self.doseList 			= doseList 			# list of (DWD)doses for each dataset
		self.initialPDBList		= initialPDBList 	# list of atom objects for initial (lowest dose) pdb file
		self.outputDir			= outputDir
		self.partialDatasets	= partialDatasets   # Bool, whether atoms only in subset of datasets are included
		self.seriesName			= seriesName
		
		# if number of datasets to perform lin reg not stated then 
		# set to total number of datasets in series
		if self.numLigRegDatasets == 0:
			self.numLigRegDatasets = len(self.datasetList[0].mindensity)

		# hard coded list of amino acids here
		self.aminoAcids = ['ALA',
					  	   'ARG',
					  	   'ASN',
					  	   'ASP',
					  	   'CYS',
					 	   'GLN',
					 	   'GLU',
					 	   'GLY',
					 	   'HIS',
					 	   'ILE',
					  	   'LEU',
					  	   'LYS',
					  	   'MET',
					  	   'PHE',
					  	   'PRO',
					  	   'SER',
					  	   'THR',
					  	   'TRP',
					  	   'TYR',
					  	   'VAL']

	def getMultiDoseAtomList(self,
							 printText = True,
							 logFile   = ''):

		# this function inputs a list of lists of PDB atom objects 
		# (see StructurePDB class) and formats as an object of the 
		# class 'multiPDB'. It is a variant of the function above which 
		# can also cope with structures containing different numbers of 
		# atoms (say if solvent molecules/ligands are included in a subset 
		# of the structures). In this case, the smallest common substructure 
		# between all structures will be used

		# first check that each PDBarray contains the same number of atoms (consistency check)
		if len(self.datasetList) > 1:
			if printText is True:
				ln = 'Multiple datasets detected...'
				self.printOrWriteToLog(logFile = logFile, txt = ln)
			for dataset in self.datasetList:
				if len(dataset) != len(self.datasetList[0]):
					if printText is True:
						txt = 'Not all PDB structures have same number of atoms!\n'+\
							  'Will only include atoms common to ALL structures...'
						self.printOrWriteToLog(logFile = logFile, txt = txt)

					break
		elif len(self.datasetList) == 1:
			if printText is True:
				ln = 'Single dataset detected...'
				self.printOrWriteToLog(logFile = logFile, txt=ln)

		PDBdoses = []
		numNotFoundAtoms = 0

		if printText is True:
			ln = 'Locating common atoms to ALL datasets...:'
			self.printOrWriteToLog(logFile = logFile, txt = ln)

		singDimAttrs = ('atomnum',
						'residuenum',
						'atomtype',
						'basetype',
					    'chaintype',
					    'X_coord',
					    'Y_coord',
					    'Z_coord')

		multiDimAttrs = ('Bfactor',
					     'Occupancy',
					     'meandensity',
					     'maxdensity',
					     'mindensity',
						 'mediandensity',
						 'numvoxels',
						 'stddensity',
						 'min90tile',
						 'max90tile',
						 'min95tile',
						 'max95tile',
						 'reliability')

		for atom in self.datasetList[0]:
			atm_counter = 1
			atomDict = {attr:getattr(atom, attr) for attr in singDimAttrs}.copy()
			atomDict.update({attr:[getattr(atom, attr)] for attr in multiDimAttrs})

			indexindataset = [] # list of index of same atom in each later dataset

			# check whether atom in all datasets:
			atomIndentifier = atom.getAtomID()
			for dataset in self.datasetList[1:]:
				found = False
				k = -1        
				for otheratom in dataset: 
					k += 1
					if atomIndentifier == otheratom.getAtomID():        
						atm_counter += 1 
						for attr in multiDimAttrs:
							atomDict[attr].append(getattr(otheratom, attr))
						indexindataset.append(k)
						found = True
						break
				if found is False: 
					indexindataset.append(-1) # if atom not in dataset
					for attr in multiDimAttrs: 
						atomDict[attr].append(np.nan) # add dummy value if atom not found

			# remove this located atom from the later dataset now that it
			# has been located --> to make loop faster 
			for j in range(1,len(indexindataset)+1):
				if indexindataset[j-1] != -1:
					self.datasetList[j].pop(indexindataset[j-1])

			if atm_counter != len(self.datasetList) and self.partialDatasets is False:
				if printText is True:
					txt = 'Atom "{}" not found in all datasets\n'.format(atomIndentifier)+\
					      '---> not including atom in atom list...'
					self.printOrWriteToLog(logFile = logFile, txt = txt)
				numNotFoundAtoms += 1
				continue

			else:
				newatom = combinedAtom()
				for attr in singDimAttrs+multiDimAttrs:
					if attr[0] != '_':
						if isinstance(atomDict[attr],list) is False:
							setattr(newatom,attr,atomDict[attr])
						else:
							metName = self.findMetricName(attr)
							newatom.getDensMetricInfo(metName,'Standard',atomDict[attr])

				# include reliability value also for Dloss metric
				newatom.getDensMetricInfo('loss','reliability',atomDict['reliability'])

				if atm_counter != len(self.datasetList):
					if printText is True:
						txt = 'Atom "{}" not found in all datasets\n'.format(atomIndentifier)+\
							  '---> including partial information in atom list...'
						self.printOrWriteToLog(logFile = logFile, txt = txt)
					numNotFoundAtoms += 1
				PDBdoses.append(newatom)

		if self.partialDatasets is True:
			if printText is True:
				ln = '# atoms not in all datasets: {}'.format(numNotFoundAtoms)
				self.printOrWriteToLog(logFile = logFile, txt = ln)
		else:
			if printText is True:
				ln = '# atoms removed since not in all datasets: {}'.format(numNotFoundAtoms)
				self.printOrWriteToLog(logFile = logFile, txt = ln)

		if printText is True:
			ln = '---> Finished!'
			self.printOrWriteToLog(logFile = logFile, txt = ln)

		self.atomList = PDBdoses

	def printOrWriteToLog(self,
						  logFile = '',
					      txt     = ''):

		# print to command line or write to log file

		if logFile == '':
			print txt
		else:
			logFile.writeToLog(str = txt)

	def findMetricName(self,metricName):

		# a conversion between the naming convention used by singlePDB class 
		# for each individual dataset and that to be used by dictionary of 
		# density metrics for combinedAtom class

		conversions = (['Bfactor','Bfactor'],
					   ['Occupancy','occupancy'],
					   ['meandensity','mean'],
					   ['maxdensity','gain'],
					   ['mindensity','loss'],
					   ['mediandensity','median'],
					   ['numvoxels','num-voxels'],
					   ['stddensity','standard-deviation'],
	    			   ['min90tile','90tile-loss'],
	    			   ['max90tile','90tile-gain'],
	    			   ['min95tile','95tile-loss'],
	    			   ['max95tile','95tile-gain'])

		for conversion in conversions:
			if metricName == conversion[0]:
				return conversion[1]

	def getNumAtoms(self):
		return len(self.atomList)

	def getAverageMetricVals(self,
							 densMet  = 'loss',
							 normType = 'Standard'):

		# get structure-wide average of selected density metric

		densList = [atom.densMetric[densMet][normType]['values'] for atom in self.atomList]
		densMean = np.nanmean(densList,0)
		densStd  = np.nanstd(densList,0)
		return densMean,densStd

	def getAtom(self,
				chain       = '',
				restype     = '',
				resnum      = '',
				atomtype    = '',
				printOutput = False):

		# get atom(s) matching specified description if found

		foundAtoms = []
		for atom in self.atomList:
			if atom.chaintype == chain or chain == '':
				if atom.residuenum == resnum or resnum == '':
					if atom.basetype == restype or restype == '':
						if atom.atomtype == atomtype or atomtype == '':
							foundAtoms.append(atom)
		if len(foundAtoms) != 0:
			if printOutput is True:
				print 'Found {} atom(s) matching description'.format(len(foundAtoms))
			return foundAtoms
		if printOutput is True:
			print 'Atom matching description not found'
		return False

	def getDensMetrics(self):

		# get a list of all density metrics and normalisations that 
		# currently have been defined with a set of values

		currentMetrics = []
		for metric in self.atomList[0].densMetric.keys():
			for normType in self.atomList[0].densMetric[metric].keys():
				currentMetrics.append([metric,normType])
		return currentMetrics

	def calcBfactorChange(self):

		# calculate Bfactor change between initial 
		# and each later dataset

		findBchange(self.initialPDBList,self.atomList,'Bfactor')

	def calcAdditionalMetrics(self,
							  metric    = 'loss',
							  normType  = 'Standard',
							  newMetric = 'Calpha normalised',
							  printText = False):

		# calculate the Calpha weights for each dataset (see CalphaWeight 
		# class for details) for metric "metric" (loss, gain, mean etc.)
		# 'newMetric' takes values in 'options' list

		options = ['Calpha normalised',
				   'net',
				   'lin reg',
				   'dataset 1 subtracted',
				   'average']

		if newMetric not in options:
			print 'new metric type not recognised.. choose from: {}'.format(options)
			return

		# check if already calculated and skip if so
		if newMetric == ('net','dataset 1 subtracted','Calpha normalised'):
			args = [metric,newMetric,'values']
		if newMetric in ('average','lin reg'):
			args = [metric,normType,newMetric]
		try:
			self.atomList[0].densMetric[args[0]][args[1]][args[2]]
			return
		except (NameError,KeyError):
			pass

		if newMetric == 'Calpha normalised':
			if printText is True:
				print 'Calculating Calpha weights at each dataset '+\
					  'for metric: {}, normalisation: {}'.format(metric,normType)
			CAweights = self.retrieveCalphaWeight(metric=metric)

		# loop over all atoms in list and calculate additional 
		# metrics for each atom in atomList
		for atom in self.atomList:
			if newMetric == 'Calpha normalised':
				atom.CalphaWeightedDensChange(CAweights,metric)
			elif newMetric == 'lin reg':
				atom.calcLinReg(self.numLigRegDatasets,'Standard',metric)
			elif newMetric == 'net':
				atom.calcNetChangeMetric('Standard')
			elif newMetric == 'dataset 1 subtracted':
				atom.calcFirstDatasetSubtractedMetric('Standard',metric)
			elif newMetric == 'average':
				atom.calcAvMetric(normType,metric)

	def retrieveCalphaWeight(self,metric = 'loss'):

		# retrieve the metric average over all Calpha atoms 
		# for current model

		CAweights = CalphaWeight(self.atomList)
		CAweights.calculateWeights(metric)
		return CAweights
	
	def calcVectorWeightedMetric(self,
								 metric   = 'loss',
								 normType = 'Standard',
								 vector   = []):

		# for a specified metric calculate new vector-weighted 
		# values, where the metric over a series of doses
		# is multiplied by a vector of per-dataset scalings

		for atom in self.atomList:
			atom.calcVectorWeightedMetric(metric,normType,vector)

	def calcVectorSubtractedMetric(self,
								   metric   = 'loss',
								   normType = 'Standard',
								   vector   = []):

		# for a specified metric calculate new vector-weighted 
		# values, where the metric over a series of doses
		# is subtracted from a vector of per-dataset scalings

		for atom in self.atomList:
			atom.calcVectorSubtractedMetric(metric,normType,vector)	

	def calcStandardisedMetrics(self,metric = 'loss'):

		# standardise distribution of a given metric to have 
		# (mean=0,std=1) for each dataset within a damage series

		data  = [atm.densMetric[metric]['Standard']['values'] for atm in self.atomList]
		meand = np.nanmean(data,0)
		stdd  = np.nanstd(data,0)
		for atm in self.atomList:
			values = atm.densMetric[metric]['Standard']['values']
			atm.getDensMetricInfo(metric,'Standardised',(values-meand)/stdd)

	def writeMetric2File(self,
						 where    = './',
						 groupBy  = 'none',
						 metric   = 'loss',
						 normType = 'Standard',
						 numDP    = 2):

		# write all metric values to a .csv file to location 'where'
		# 'groupBy' takes values 'none','residue','atomtype'. numDP 
		# is number of decimal points the values in csv should be
		# rounded to (default of 2dp)

		csvName = '{}{}-{}.csv'.format(where,metric,normType.replace(" ",""))
		if groupBy in ('residue','atomtype'):
			csvName = csvName.replace('.csv','_{}.csv'.format(groupBy))
		csvfile = open(csvName,'w')
		csvfile.write('metric: {} normalisation: {}\n'.format(metric,normType))
		if groupBy == 'none':
			csvfile.write('atomnum,atominfo,')
		else:
			csvfile.write('index,atominfo,')
		metricLength = self.atomList[0].getNumDatasets(metric)
		for i in range(0,metricLength):
			if groupBy == 'none':
				csvfile.write('{}'.format(str(i+1)))
			else:
				csvfile.write('{} (mean),{} (std dev)'.format(str(i+1),str(i+1)))
			if i != metricLength-1: csvfile.write(',')
		csvfile.write('\n')

		# retrieve metric stats depending on groupBy
		if groupBy != 'none':
			# create dictionary for stats over datasets
			if groupBy == 'atomtype':
				TypeDic = self.groupByAtmType()
			elif groupBy == 'residue':
				TypeDic = self.groupByResType()
			statDic = {k:[] for k in TypeDic.keys()}

			for d in range(metricLength):
				if groupBy == 'atomtype':
					a,b = self.getPerAtmtypeStats(metric   = metric,
												  normType = normType,
												  dataset  = d,
												  n        = 1)
				elif groupBy == 'residue':
					a,b = self.getPerResidueStats(metric   = metric,
												  normType = normType,
												  dataset  = d,
												  n        = 1)
				for key in b.keys():
					statDic[key].append(b[key]['mean'])
					statDic[key].append(b[key]['std'])
			i = 0
			for k in statDic.keys():
				i+=1
				csvfile.write('{},{},{}\n'.format(i,k,','.join(map(str,statDic[k]))))			
		else:
			self.atomList.sort(key=lambda x: x.atomnum) # sort atom list by atom number
			for atom in self.atomList:
				csvfile.write('{},{},'.format(atom.atomnum,atom.getAtomID()))
				roundedVals = [round(v,numDP) for v in atom.densMetric[metric][normType]['values']]
				csvfile.write(','.join(map(str,roundedVals)))
				csvfile.write('\n')
		csvfile.close()

	def findMetricRatio(self,
						metric    = 'loss',
						normType  = 'Standard',
						resiType  = '',
						atomType1 = '',
						atomType2 = ''):

		# find ratio between metric for two atom 
		# types within a single residue type

		title = '\n{} D{} comparison between {}-{} and {}-{} atoms'.format(normType,metric,
																		   resiType,atomType1,
																		   resiType,atomType2)
		print '\n'+'-'*len(title)+title
		dic,ratioDic = {},{'distance':[],'ratio':[]}
		for atom in self.atomList:
			if atom.basetype == resiType and atom.atomtype == atomType1:
				dic[atom.getAtomID()] = atom.densMetric[metric][normType]['values']
		for atom in self.atomList:
			if atom.basetype == resiType and atom.atomtype == atomType2:
				key   = '-'.join(atom.getAtomID().split('-')[0:3]+[atomType1])
				# check that generated 'key' exists and skip if not (needed for incomplete residues)
				if key not in dic.keys(): continue
				a,b = np.array(dic[key]),atom.densMetric[metric][normType]['values']
				ratioDic['ratio'].append(a/b)
				ratioDic['distance'].append(a-b)
		summary = {}
		for key in ratioDic:
			outputString  = '\nFor metric {} between two atoms'.format(key)	
			outputString += '\n# atom pairs: {}'.format(len(ratioDic[key]))
			outputString += '\nmean\t' + '\t'.join([str(round(val,3)) for val in list(np.mean(ratioDic[key],0))])
			outputString += '\nmedian\t' + '\t'.join([str(round(val,3)) for val in list(np.median(ratioDic[key],0))])
			outputString += '\nstd dev\t' + '\t'.join([str(round(val,3)) for val in list(np.std(ratioDic[key],0))])
			print outputString
			summary[key] = title+outputString

		return ratioDic,summary

	def findMetricRatioKeyResidues(self,
								   metric     = 'loss',
								   normType   = 'Standard',
								   rType      = 'ratio',
								   errorBars  = True,
								   pairs      = [],
								   title      = '',
								   fileType   = '.svg',
								   axesFont   = 18,
								   legendFont = 18,
								   titleFont  = 24):

		# run the above findMetricRatio method for known susceptible 
		# residues and write to file. 'rType' takes values ('ratio','distance')
		# 'pairs' is list of residues and atoms of form 
		# [['GLU','CD','CG'],['GLU','CD','OE1'],..] etc

		if seabornFound is False:
			return

		fileout = open('{}{}_{}-D{}-{}.txt'.format(self.outputDir,title,normType,metric,rType),'w')

		foundPairs = {}
		for pair in pairs:
			found = self.getAtom(restype=pair[0])

			if found == False: 
				continue # don't proceed if no atoms of a residue type

			output = self.findMetricRatio(metric   = metric,
										  normType  = normType,
										  resiType  = pair[0],
										  atomType1 = pair[1],
										  atomType2 = pair[2])

			foundPairs['-'.join(pair)] = output[0]
			fileout.write(output[1][rType]+'\n')
		fileout.close()

		# plot comparison for different atom pairings
		sns.set_palette("deep", desat=.6)
		sns.set_context(rc={"figure.figsize":(10, 10)})
		fig = plt.figure()
		ax = plt.subplot(111)
		colormap = plt.cm.nipy_spectral
		plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(foundPairs.keys()))])
		for key in foundPairs.keys():
			xData = self.doseList
			yData = np.mean(foundPairs[key][rType],0)
			if errorBars is False:
				plt.plot(xData,yData,label=key)
			else:
				yError = np.std(foundPairs[key][rType],0)
				plt.errorbar(xData,
							 yData,
							 yerr     = yError,
							 fmt      = '-o',
							 capthick = 2,
							 label    = key)

		plt.xlabel('Dataset #', fontsize = axesFont)
		plt.ylabel('D{} {}'.format(metric,rType), fontsize = axesFont)

		# place legend outside to right of plot
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax.legend(loc            = 'center left',
				  bbox_to_anchor = (1, 0.5),
				  fontsize       = legendFont)

		fig.suptitle('{} D{} {}'.format(normType,metric,rType),fontsize = titleFont)

		figName = '{}{}-D{}-{}-{}_{}'.format(self.outputDir,normType,metric,rType,title)
		fileName = self.checkUniqueFileName(fileName = figName,
								 			fileType = fileType)
		fig.savefig(fileName)

	def calcMetricDiffFromStructureMean(self,
										metric   = 'loss',
										normType = 'Standard',
										diff     = 'difference'):

		# for a specified metric calculate difference between metric 
		# for atom and from mean of structure. 'diff' takes three 
		# forms: (a) 'difference', (b) 'ratio', (c) 'num-stds'

		av,std = self.getAverageMetricVals(metric,normType)
		for atom in self.atomList:
			if diff == 'difference':
				atom.calcDiffFromMeanMetric(metric,normType,av)
			elif diff == 'ratio':
				atom.calcRatioToMeanMetric(metric,normType,av)
			elif diff == 'std-devs':
				atom.calcNumStdFromMeanMetric(metric,normType,av,std)

	def numAtmsWithMetricAboveLevel(self,
									dataset   = 0,
									metric    = 'loss',
									normType  = 'Standard',
									threshold = 0,
									printStr  = False):

		# determine the number of atoms with metric value above a 
		# specified threshold 'threshold' determined as number of 
		# std dev away from structure-wide average for metric 'printStr' 
		# (Boolian) is True if string printed to command line

		self.calcMetricDiffFromStructureMean(metric   = metric,
											 normType = normType,
											 diff     = 'std-devs')
		count = 0
		for atom in self.atomList:
			if atom.densMetric[metric]['num stds']['values'][dataset] < -threshold:
				count += 1
		if printStr is True:
			print 'For dataset {}:'.format(dataset)
			print '# atoms with {} D{} metric above {}'.format(normType,metric,threshold)+\
				  'std dev of structure-wide mean is {}'.format(count)
		return count

	def findMetricRatioKeyResidues_scatterplot(self,
											   metric     = 'loss',
											   normType   = 'Standard',
											   rType      = 'ratio',
											   pairs      = [['GLU','CD','CG'],['GLU','CD','OE1']],
											   title      = '',
											   fileType   = '.svg',
											   axesFont   = 18,
											   legendFont = 18,
											   titleFont  = 24):

		# run the above findMetricRatio method for known susceptible 
		# residues and write to file. 'rType' takes values ('ratio',
		# 'distance'). 'pairs' is list of residues and atoms of form 
		# [['GLU','CD','CG'],['GLU','CD','OE1'],..] etc.

		if seabornFound is False:
			return

		foundPairs = {}
		for pair in pairs:
			found = self.getAtom(restype=pair[0])
			if found == False: 
				continue # don't proceed if no atoms of a residue type
			output = self.findMetricRatio(metric   = metric,
										  normType  = normType,
										  resiType  = pair[0],
										  atomType1 = pair[1],
										  atomType2 = pair[2])
			foundPairs['-'.join(pair)] = output[0]

		# plot the relationship between xData and yData
		sns.set_palette("deep", desat=.6)
		sns.set_context(rc={"figure.figsize":(10, 10)})
		fig = plt.figure()
		ax = plt.subplot(111)

		total = 0
		for key in foundPairs.keys(): 
			total += len(foundPairs[key][rType])
		xDataTotal = np.linspace(0,1,total)

		i,j = 0,-1
		cmap = cm.Set1
		nsteps = len(foundPairs.keys())

		colors = ['#737474',
				  '#409cd6',
				  '#58bb6b',
				  '#faa71a',
				  '#ff6a6a']

		DataSave = {}
		for key in foundPairs.keys():	
			j += 1
			yData = np.nanmean(foundPairs[key][rType],1)
			xData = xDataTotal[i:i+len(yData)]
			i += len(yData)

			plt.scatter(xData,
						yData,
						marker     = 'o',
						s          = 100,
						c          = colors[j],
						edgecolors = '#FFFFFF',
						label      = key)

			DataSave[key] = {'y' : yData,
							 'x' : xData}

		# place legend outside to right of plot
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax.legend(loc            ='center left',
				  bbox_to_anchor = (1, 0.5),
				  fontsize       = legendFont)

		plt.xlabel(' ', fontsize = axesFont)
		plt.ylabel('D{} {}'.format(metric,rType), fontsize = axesFont)
		figtitle = '{} D{} {}'.format(normType,metric,rType)
		fig.suptitle(figtitle, fontsize = titleFont)

		fileName = self.checkUniqueFileName(fileName = figtitle,
								 			fileType = fileType)
		fig.savefig(fileName)

		return DataSave

	def checkUniqueFileName(self,
							fileName = '',
							fileType = '.svg'):

		# check that a given file name for a figure is unique before 
		# it is saved and append index if not

		strippedName = fileName.replace(' ','') # remove whitespace

		if os.path.isfile(strippedName) is False:
			if strippedName.endswith(fileType):
				return strippedName
			else:
				return strippedName+fileType # this may give wrong name..

		fname = lambda x: strippedName.strip(fileType)+'_{}{}'.format(x,fileType)
		i = 1
		while os.path.isfile(fname(i)): 
			i += 1 
		return fname(i)

	def getAtomListWithoutPartialAtoms(self,dataset = 0):

		# return new list of atoms with partial atoms 
		# removed for that dataset

		newList = []
		for atom in self.atomList:
			if dataset in atom.getPresentDatasets():
				newList.append(atom)
		return newList

	def getTopNAtoms(self,
					 metric   = 'loss',
					 normType = 'Standard',
					 dataset  = 0,
					 n        = 25):

		# for a given metric type, determine top 
		# n damage sites. 'n' is integer or 'all'.

		if dataset != 'all':
			atomList = self.getAtomListWithoutPartialAtoms(dataset = dataset) # remove atoms not in present dataset
			atomList.sort(key=lambda x: x.densMetric[metric][normType]['values'][dataset])

			# for difference maps, more negative metrics indicate 
			# more damage, however if Calpha normalisation used, 
			# then the highest positive metric values are the most 
			# damaged, this must be accounted for here:
			if normType == 'Calpha normalised':
				atomList.reverse()

			if n != 'all': 
				return atomList[:int(n)]
			else: 
				return atomList

		else:
			# TODO: rank with respect to all datasets
			return 

	def getTopNAtomsPDBfile(self,
							metric   = 'loss',
							normType = 'Standard',
							dataset  = 0,
							n        = 'all',
							pdbFile  = 'untitled.pdb'):

		# for a given metric type, determine top n damage sites. 
		# 'n' is integer (not 'all' here). method writes locations 
		# of top damage sites to pdb file

		pdbIn = open(pdbFile,'r')
		pdbOutName = '{}top{}-{}D{}sites-dset{}.pdb'.format(self.outputDir,n,normType,metric,dataset)
		pdbOut = open(pdbOutName,'w')

		headerInfo = 'REMARK Top {} damage sites as indicated by the {} D{} metric.\n'.format(n,normType,metric)
		if dataset == 'all':
			headerInfo += 'REMARK All datasets within damage series are present (increasing chain index).\n'
		else:
			headerInfo += 'REMARK File contains information for dataset {}.\n'.format(dataset)
		pdbOut.write(headerInfo)

		for l in pdbIn.readlines():
			if l.split()[0] in ('CRYST1','SCALE1','SCALE2','SCALE3'):
				pdbOut.write(l)

		if dataset == 'all':
			datasets = range(self.getNumDatasets())
		else:
			datasets = [dataset]

		chains = string.ascii_uppercase
		for d in datasets:
			# next get top n metric sites and write to file
			topAtms = self.getTopNAtoms(metric   = metric,
										normType = normType,
										dataset  = d,
										n        = n)
			count = 0
			for atm in topAtms:
				count += 1
				# generate pdb line - note: absolute value of metric 
				# taken since Bfactor>0 (for displaying in pymol etc.)
				l = writePDBline_DamSite(atom     = atm,
										 damValue = np.abs(atm.densMetric[metric][normType]['values'][d]),
										 index    = count,
										 chain    = chains[d])
				pdbOut.write(l+'\n')

		pdbOut.write('END')
		pdbIn.close()
		pdbOut.close()
		return pdbOutName

	def getTopNAtomsString(self,
						   metric   = 'loss',
						   normType = 'Standard',
						   dataset  = 0,
						   n        = 25):

		# return the top n atoms within structure in terms of 
		# metric 'metric'. 'n' takes values 'all' or an integer

		topN = self.getTopNAtoms(metric   = metric,
								 normType = normType,
								 dataset  = dataset,
								 n        = n)
		atomInfoList = []
		for atom in topN:
			data = atom.getAtomID()+'\t\t'
			data += str(round(atom.densMetric['loss'][normType]['values'][dataset],2))+'\t\t'
			data += str(round(atom.densMetric['loss']['reliability']['values'][dataset],2))+'\t\t'

			for met in ('mean','gain','Bfactor'):
				data += str(round(atom.densMetric[met][normType]['values'][dataset],2))+'\t\t'
			atomInfoList.append(data)
		stringOut = 'Atom-ID\t\t\tDloss\t\tProximity\t\tDmean\t\tDgain\t\tBfactor\n'
		stringOut += '\n'.join(atomInfoList)	
		return stringOut

	def breakdownTopNatomsBy(self,
							 metric   = 'loss',
							 normType = 'Standard',
							 dataset  = 0,
							 n        = 0,
							 sortby   = ['basetype','atomtype']):

		# find the top n atoms and sort them by grouping 'sortby'

		topN = self.getTopNAtoms(metric   = metric,
								 normType = normType,
								 dataset  = dataset,
								 n        = n)
		countDic,totalDic = {},{}
		for atom in topN:
			vals = [getattr(atom,att) for att in sortby]
			atomID = '-'.join(vals)
			if atomID not in countDic.keys():
				countDic[atomID] = 1
				totalDic[atomID] = self.getNumAtomsOfType(atts = sortby,
														  vals = vals)
			else:
				countDic[atomID] += 1
		headerString  = 'Partition of top atom {} D{} values by {}.\n'.format(normType,metric,' & '.join(sortby))
		fracOfTotal   = round((float(n)/self.getNumAtoms())*100,2) # n as fraction of total atoms
		headerString += 'For dataset {}, top {} ({}%) atoms included.\n'.format(dataset,int(n),fracOfTotal)
		outputString  = 'Type\t\tFraction'

		keys = countDic.keys()
		keys.sort()
		for key in keys:
			outputString += '\n{}\t\t{}/{}'.format(key,countDic[key],totalDic[key])
		return headerString,outputString,countDic

	def getNumAtomsOfType(self,
						  atts = ['chaintype','restype'],
						  vals = ['A','GLU']):

		# get number of atoms with attribute 'att' at specified 
		# value 'val' set 'att' and 'val' as lists of strings, 
		# i.e. att = ['chaintype','restype'], val = ['A','GLU']

		i = 0
		for atom in self.atomList:
			if vals == [getattr(atom,att) for att in atts]: i += 1
		return i

	def getTopAtomsStackedBarplot(self,
								  metric     = 'loss',
								  normType   = 'Standard',
								  n          = 40,
								  palette    = 'hls',
								  barWidth   = 0.5,
								  axesFont   = 18,
								  titleFont  = 18,
								  legendFont = 18,
								  saveFig    = True,
								  outputDir  = './',
								  fileType   = '.svg',
								  plotTitle  = ''):

		# get top n atoms per dataset and sort by residue/nucleotide 
		# type. Then create a stacked barplot per dataset to distinguish
		# the breakdown of damage with increasing dose

		plotDic = {}
		numDsets = self.getNumDatasets()

		for d in range(numDsets):
			stats = self.breakdownTopNatomsBy(metric    = metric,
								 	          normType  = normType,
								 			  dataset   = d,
								 			  n         = n,
								 			  sortby    = ['basetype'])

			breakdown = stats[2]

			for key in breakdown.keys():
				if key not in plotDic.keys():
					plotDic[key] = [0]*d + [breakdown[key]]
				else:
					plotDic[key] += [breakdown[key]]
			for key in plotDic.keys():
				if len(plotDic[key]) < d+1:
					plotDic[key] += [0]

		sns.set_style("whitegrid")
		sns.set_context(rc = {"figure.figsize": (10, 6)})
		fig = plt.figure()
		ax = plt.subplot(111)

		accumulatedVals = np.array([n]*numDsets)

		numResis = len(plotDic.keys())
		for i,(k, color) in enumerate( zip( plotDic.keys(), sns.color_palette(palette, n_colors = numResis,desat = 0.9) ) ):

			# plt.barplot(range(numDsets),
			# 		accumulatedVals,
			# 		label = k,
			# 		color = color,
			# 		width = barWidth)

			fd = pd.DataFrame()
			fd['dataset']   = np.array(range(numDsets))+1
			fd['frequency'] = accumulatedVals
			sns.barplot(x     = 'dataset',
						y     = 'frequency',
						data  = fd,
						label = k,
						color = color)

			accumulatedVals -= np.array(plotDic[k])

		# place legend outside to right of plot
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax.legend(loc            = 'center left',
				  bbox_to_anchor = (1, 0.5),
				  fontsize       = legendFont)

		sns.despine()

		# plt.xticks(np.array(range(numDsets)) + barWidth/2.,['{}'.format(i+1) for i in range(numDsets)])

		plt.xlabel('Dataset Number', fontsize = axesFont)
		plt.ylabel('Frequency', fontsize = axesFont)

		if plotTitle == '':
			t = 'Top damage sites by residue/nucleotide type'
		else:
			t = plotTitle
		plt.title(t, fontsize = titleFont)

		if not saveFig: 
			plt.show()
		else:
			saveName = '{}Top{}DamSiteByRestype_Metric-D{}_Normalisation-{}'.format(outputDir,
																				    n,
																					metric,
																					normType.replace(" ",""))

			fileName = self.checkUniqueFileName(fileName = saveName,
									 			fileType = fileType)
			fig.savefig(fileName)

		return fileName

	def getPerAtmtypeStats(self,
						   metric   = 'loss',
						   normType = 'Standard',
						   dataset  = 0,
						   sortby   = 'mean',
						   n        = 25):

		# for a given metric type, determine per-atom-type 
		# statistics on the distribution of damage

		atmtypeDict = self.groupByAtmType(dataset = dataset) # group atoms by atom type

		statsDic 	= self.getStats(metric   = metric,
									normType = normType,
									dataset  = dataset,
									dic      = atmtypeDict)

		statsString	= self.reportStats(stats    = statsDic,
									   name     = 'Type',
									   sortby   = sortby,
									   n        = n,
									   normType = normType)
		return statsString,statsDic

	def getPerResidueStats(self,
						   metric   = 'loss',
						   normType = 'Standard',
						   dataset  = 0,
						   sortby   = 'mean',
						   n        = 25):

		# for a given metric type, determine per-residue 
		# statistics on the distribution of damage

		resDict     = self.groupByResType(dataset = dataset) # group atoms by residue type

		statsDic 	= self.getStats(metric   = metric,
									normType = normType,
									dataset  = dataset,
									dic      = resDict)

		statsString	= self.reportStats(stats    = statsDic,
									   name     = 'Residue',
									   sortby   = sortby,
									   n        = n,
									   normType = normType)

		return statsString,statsDic

	def getPerChainStats(self,
						   metric   = 'loss',
						   normType = 'Standard',
						   dataset  = 0,
						   sortby   = 'mean',
						   n        = 25):

		# for a given metric type, determine per-chain 
		# statistics on the distribution of damage

		chainDict 	= self.groupByChainType(dataset = dataset) # group atoms by chain type

		statsDic 	= self.getStats(metric   = metric,
									normType = normType,
									dataset  = dataset,
									dic      = chainDict)

		statsString	= self.reportStats(stats    = statsDic,
									   name     = 'Chain',
									   sortby   = sortby,
									   n        = n,
									   normType = normType)

		return statsString,statsDic

	def getStats(self,
				 metric   = 'loss',
				 normType = 'Standard',
				 dataset  = 0,
				 dic      = {}):

		# output distribution stats for each element 
		# in dictionary dic as a new dictionary

		getStatsPerKey = {}
		for key in dic.keys():
			metricList = [atom.densMetric[metric][normType]['values'][dataset] for atom in dic[key]]
			getStatsPerKey[key] = self.getStatsForList(metricList = metricList,
													   metric     = metric,
													   normType   = normType)
		return getStatsPerKey

	def getStatsForList(self,
					    metricList = [],
					    metric     = 'loss',
					    normType   = 'Standard'):

		# calculate measures of the distribution 
		# of values in a list metricList

		statsDic = {}
		statsDic['#atoms']	 	= len(metricList)
		statsDic['mean'] 		= np.mean(metricList)
		statsDic['std'] 		= np.std(metricList)
		statsDic['skew'] 		= self.calculateSkew(metricList = metricList)
		statsDic['ratio'] 		= self.calcNetRatio(metricList = metricList)
		statsDic['outliers'] 	= self.calcNumOutliers(metricList = metricList,
													   metric     = metric,
													   normType   = normType)
		statsDic['normality']   = self.testForNormality(metricList = metricList)
		statsDic['returnOrder'] = ['mean',
								   'std',
								   '#atoms',
								   'outliers',
								   'skew',
								   'normality']
		return statsDic

	def reportStats(self,
					stats    = {},
					name     = 'Residue',
					sortby   = 'mean',
					n        = 20,
					normType = 'Standard',
					format   = 'txt'):

		# format stats defined in 'stats' dictionary into an 
		# output string. 'sortby' specifies how the output should 
		# be ranked. 'n' is number of entries to print (specify 
		# an integer or 'all') format takes values in ('txt','html')

		sortbyDict   = self.getStatsForList(metricList=[1,1])
		headerOrder  = sortbyDict['returnOrder']
		headerString = '{}\t\t{}\t\t{}\t\t{}\t\t{}\t{}\t\t{}\n'.format(name,*headerOrder)
		if sortby not in sortbyDict.keys(): return 'Unexpected ranking specified'
		if n != 'all' and not isinstance(n,int): return 'input "n" must be integer or "all", unable to complete'
		
		list1,list2 = [],[]
		for key in stats.keys():
			statsDic = stats[key]
			statsFmtd = []		
			for key2 in headerOrder: 
				if isinstance(statsDic[key2],float): statsFmtd.append('{0:.3f}'.format(statsDic[key2]))
				else: statsFmtd.append(str(statsDic[key2]))
			string = '\t\t'.join([str(key)]+statsFmtd)
			list1.append(string)
			list2.append(statsDic[sortby])

		sortedList1 = [x for (y,x) in sorted(zip(list2,list1))] # sort by chosen 'sortby' parameter

		# for difference map analysis, more negative metrics indicate 
		# more damage, however if Calpha normalisation used, then the 
		# highest positive metric values are the most damaged, this 
		# must be accounted for here:
		if normType == 'Calpha normalised': 
			sortedList1.reverse()

		if n != 'all':
			stringOut =	'\n'.join(sortedList1[:n])
		else:
			stringOut =	'\n'.join(sortedList1)	
		return headerString+stringOut

	def calcDiscreteDistMode(self,metricList = []):

		# calculate mode of discrete histogram with 100 bins

		h = np.histogram(metricList,bins=100)
		mode_ind = np.argmax(h[0])
		mode_val = np.mean(h[1][mode_ind:mode_ind+2])
		return mode_val

	def calcNetRatio(self,metricList = []):

		# calculate net ratio of distn values either side of distn mode

		switchVal = self.calcDiscreteDistMode(metricList = metricList)
		group = {'above':[],'below':[]}
		for val in metricList:
			if val >= switchVal: 
				group['above'].append(val-switchVal)
			else:
				group['below'].append(switchVal-val)
		netAbove = np.sum(group['above'])
		netBelow = np.sum(group['below'])
		netRatio = float(netBelow)/netAbove
		return netRatio

	def calcNumOutliers(self,
						metricList = [],
						metric    = 'loss',
						normType  = 'Standard'):

		# by assuming a symmetric distn around the mode, flag 
		# specific atoms that fall outside this range

		distMode = self.calcDiscreteDistMode(metricList = metricList)
		if normType == 'Standard':
			distMax  = np.max(metricList)
			sudoMin  = distMode - np.linalg.norm(distMax-distMode)
			count 	 = 0
			for val in metricList:
				if val < sudoMin: count+=1
		elif normType == 'Calpha normalised':
			distMin  = np.min(metricList)
			sudoMax  = distMode + np.linalg.norm(distMin-distMode)
			count 	 = 0
			for val in metricList:
				if val > sudoMax: count+=1
		return count

	def calculateSkew(self,metricList = []):

		# calculate skewness for a distribution of 
		# metric values for an input list of atoms

		skew = scipy.stats.skew(metricList,
							    axis = 0,
							    bias = True)
		return skew

	def testForNormality(self,
						 metricList   = [],
						 suppressText = True):

		# test whether the metric values for a list of 
		# atoms differs from a normal distribution

		if len(metricList) < 8:
			if suppressText is False:
				print 'Skew-test not valid when less than 8 values provided - skipping'
			return 'n/a'
		elif len(metricList) < 20:
			if suppressText is False:
				print 'kurtosis-test not valid when less than 20 values provided - skipping'
			return 'n/a'
		try:
			(k2,pvalue) = scipy.stats.mstats.normaltest(metricList, axis=0)
		except np.ma.core.MaskError:
			return 'n/a'
		return pvalue

	def groupByAtmType(self,dataset = 0):

		# group atoms in a dictionary by atom type

		atmtypeDict = {}
		for atom in self.atomList:
			if dataset not in atom.getPresentDatasets(): continue
			atmtype = atom.basetype+'-'+atom.atomtype
			if atmtype not in atmtypeDict.keys():
				atmtypeDict[atmtype] = [atom]
			else:
				atmtypeDict[atmtype].append(atom)
		return atmtypeDict

	def groupByResType(self,dataset = 0):

		# group atoms in a dictionary by residue type

		resDict = {}
		for atom in self.atomList:
			if dataset not in atom.getPresentDatasets(): continue
			if atom.basetype not in resDict.keys():
				resDict[atom.basetype] = [atom]
			else:
				resDict[atom.basetype].append(atom)
		return resDict

	def groupByChainType(self,dataset = 0):

		# group atoms in a dictionary by chain type

		chainDict = {}
		for atom in self.atomList:
			if dataset not in atom.getPresentDatasets(): continue
			if atom.chaintype not in chainDict.keys():
				chainDict[atom.chaintype] = [atom]
			else:
				chainDict[atom.chaintype].append(atom)
		return chainDict

	def checkCalphaAtomsExist(self,printText = True):

		# check that Calpha backbone protein atoms 
		# actually exist within structure

		count = 0
		for res in self.aminoAcids:
			atom = self.getAtom(restype     = res,
								atomtype    = 'CA',
								printOutput = False)
			if atom is False:
				count += 1
		if count == len(self.aminoAcids):
			if printText is True:
				print 'Warning: no Calpha backbone atoms found in structure!'
				return False
		else:
			return True

	def densityMetricHeatMap(self,
							 metric       = 'loss',
							 normType     = 'Standard',
							 saveFig      = False,
							 fileType     = '.svg',
							 titleFont    = 24,
							 outputDir    = './',
							 forceSquares = True):

		# plot a heatmap of density metric values per atom per dataset.
		# depending on the metric and normalisation type chosen, a different
		# colour scheme is chosen
		chains = self.getChains()
		data = {c : {'Atom'    : [],
					 'Dataset' : [],
					 'metric'  : []} for c in chains}

		# data = {'Atom'    : [],
		# 		'Dataset' : [],
		# 		'metric'  : []}

		for atom in self.atomList:
			c = atom.chaintype
			data[c]['metric']  += atom.densMetric[metric][normType]['values']
			data[c]['Atom']    += [atom.getAtomID()]*self.getNumDatasets()
			data[c]['Dataset'] += range(1, self.getNumDatasets() + 1)

		# stretch the plot in the vertical direction, depending 
		# on how many atoms are present within the structure
		fontsize_pt = plt.rcParams['ytick.labelsize']
		dpi = 60
		chainLengths = [len(data[c]['Atom']) for c in chains]
		maxLength = float(max(chainLengths)) # get largest chain length
		matrix_height_pt = fontsize_pt * maxLength
		matrix_height_in = matrix_height_pt / dpi
		top_margin    = 0.01   
		bottom_margin = 0.04
		figure_height = matrix_height_in / (1 - top_margin - bottom_margin)
		figure_width = 5*self.getNumDatasets()

		# fig, ax = plt.subplots(1,len(chains),figsize = (figure_width,figure_height))
		# 					# gridspec_kw = dict(width_ratios = [1,1,1]))

		fig = plt.figure(figsize = (figure_width,figure_height))

		# gs = gridspec.GridSpec(1, 3, height_ratios = np.array(chainLengths)/maxLength,width_ratios =[1,1,1] )
		
		gs = gridspec.GridSpec(2, 3)

		ax = []
		for i in range(3):
			ax.append(plt.subplot(gs[i]))
		# cbar_ax = fig.add_axes()

			# , 
		 #        			   gridspec_kw = dict(height_ratios = (.001, 1),width_ratios = (0.5,1)))

		# top    = 1 - top_margin,
		#         				 				  bottom = bottom_margin,

		if metric == 'loss' and normType == 'Standard':
			pal = 'Reds_r'
		else:
			pal = sns.diverging_palette(220, 
										10, 
										sep     = 80,
										n       = 7,
										as_cmap = True)


		for i,c in enumerate(chains):
			d = pd.DataFrame(data = data[c])
			d = d.pivot("Atom",
						"Dataset",
						"metric")
						
			f = sns.heatmap(d,
							ax       = ax[i],
							cbar     = False,
							cmap     = pal,
							robust   = True,
							square   = forceSquares,
							linewidths = 1,
							linecolor = 'white',
							cbar_kws = {"orientation" : "horizontal"})


		# cbar=i == 0,
	 #               			cbar_ax=None if i else cbar_ax,

		for i,a in enumerate(ax):
			plt.sca(a)
			plt.yticks(rotation = 0) 
			# a.set_aspect(chainLengths[i]/maxLength)

		figTitle = plt.suptitle('Metric: D{},\nNormalisation: {}\n'.format(metric,normType),
				  fontsize = titleFont, y = 1.02)

		fig.tight_layout()

		if saveFig is False:
			plt.show()

		else:
			saveName = '{}heatmap_metric-{}_normalisation-{}'.format(outputDir,
																	 metric,
																	 normType)
			fileName = self.checkUniqueFileName(fileName = saveName,
									 			fileType = fileType)
			fig.savefig(fileName, bbox_extra_artists = (figTitle,), bbox_inches = "tight")

	def getChains(self):

		# find different chains in list of atoms

		chains = []
		for atom in self.atomList:
			chain = atom.chaintype
			if chain not in chains:
				chains.append(chain)
		return chains

	def detectSuspiciousAtoms(self,
							  dataset      = 0,
							  metric       = 'loss',
							  normType     = 'Standard',
							  threshold    = 5,
							  suppressText = True):

		# detect any atoms within a speicific type (e.g. LYS-NZ) 
		# that do not behaviour like rest of that type. Returns
		# a list of atom identifiers and a list of states (higher
		# or lower than expected for each atom detected).

		atmDic = self.groupByAtmType(dataset = dataset)
		suspAtoms = []
		highOrLow = [] 
		for k in atmDic.keys():
			vals = [atm.densMetric[metric][normType]['values'][dataset] for atm in atmDic[k]]
			meanVal = np.mean(vals)
			stdVal = np.std(vals)
			for atm in atmDic[k]:
				val = atm.densMetric[metric][normType]['values'][dataset]
				if np.linalg.norm(val-meanVal) > stdVal*threshold:
					suspAtoms.append(atm.getAtomID())
					if val < meanVal:
						if suppressText is False:
							print 'Unusually low compared to mean value ({} < {})'.format(round(val,3),
																						  round(meanVal,3))
						highOrLow.append('low')
					else:
						if suppressText is False:
							print 'Unusually high compared to mean value ({} > {})'.format(round(val,3),
																						   round(meanVal,3))
						highOrLow.append('high')
		if suppressText is False:
			print '{} atoms found with suspiciously high/low damage relative to average of that atom type'
		return suspAtoms,highOrLow

	def getNumDatasets(self):

		# get the number of datasets within the damage series

		return self.atomList[0].getNumDatasets()

	def graphMetricDistn(self,
						 metric    = 'loss',
						 normType  = 'Standard',
						 valType   = 'average',
						 plotType  = 'both',
						 resiType  = 'all',
						 save      = True,
						 fileType  = '.svg',
						 outputDir = './',
						 printText = True,
						 axesFont  = 18,
						 titleFont = 20,
						 plotTitle = ''):

		# histogram/kde plot of density metric per atom.
		# plotType is 'histogram' or 'kde'.
		# resiType is 'all' or list of residue types.
		# 'valType' takes values 'average' (compute average),
		# 'all' (plot all datasets), or 'n' for int (not str)
		#  n (plot dataset n).

		if seabornFound is False:
			return

		# attempt to find atoms of type 'resiType' and 
		# flag if not all found
		if resiType != 'all':
			count = 0
			for res in resiType:
				atms = self.getAtom(restype = res)
				if atms is False: 
					count += 1
			if count > 0:
				if printText is True:
					print 'Warning: not all selected residue types found!'
			if count == len(resiType):
				if printText is True:
					print 'Warning: no residues found for current plot'
				return {}

		if plotType not in ('hist','kde','both'): 
			return 'Unknown plotting type selected.. cannot plot..'

		if normType == 'Calpha normalised':
			self.calcAdditionalMetrics(metric = metric)
		if self.checkMetricPresent(metric = metric,normType = normType) is False: 
			return # check metric valid

		sns.set_style("whitegrid")
		sns.set_context(rc = {"figure.figsize": (10, 6)})
		fig = plt.figure()
		ax = plt.subplot(111)

		plotData = {}
		if valType == 'all':
			numDsets = self.getNumDatasets()
			for j,(i, color) in enumerate(zip(range(numDsets), sns.color_palette('Blues', numDsets))):
				presentAtms = self.getAtomListWithoutPartialAtoms(dataset = i)

				if resiType == 'all':
						datax = [atm.densMetric[metric][normType]['values'][i] for atm in presentAtms]
						self.plotHist(plotType = plotType,
									  datax    = datax,
									  lbl      = 'Dataset {}'.format(i),
									  color    = color)
						plotData[resiType] = datax

				else:
					for res in resiType:
						datax = [atm.densMetric[metric][normType]['values'][i] for atm in presentAtms if atm.basetype == res]
						if len(datax) > 0:
							self.plotHist(plotType = plotType,
										  datax    = datax,
										  lbl      = 'Dataset {},{}'.format(i,res),
										  color    = color)
						plotData[res] = datax	

		else:
			self.calcAdditionalMetrics(metric    = metric,
									   normType  = normType,
									   newMetric = 'average')
			if resiType == 'all':

				if valType == 'average':
					datax = [atm.densMetric[metric][normType]['average'] for atm in self.atomList]
					lbl = 'average'
				else:
					presentAtms = self.getAtomListWithoutPartialAtoms(dataset = valType)
					datax = [atm.densMetric[metric][normType]['values'][valType] for atm in presentAtms]
					lbl = 'Dataset '+str(valType)

				self.plotHist(plotType = plotType,
							  datax    = datax,
							  lbl      = lbl,
							  color    = 'r')
				plotData[resiType] = datax

			else:
				for j,(res, color) in enumerate(zip(resiType, sns.color_palette('hls', len(resiType)))):

					if valType == 'average':
						datax = [atm.densMetric[metric][normType]['average'] for atm in self.atomList if atm.basetype == res]
						lbl = 'average, {}'.format(res)

					else:
						presentAtms = self.getAtomListWithoutPartialAtoms(dataset = valType)
						datax = [atm.densMetric[metric][normType]['values'][valType] for atm in presentAtms if atm.basetype == res]
						lbl = 'Dataset {}, {}'.format(valType,res)

					if len(datax) > 0:
						self.plotHist(plotType = plotType,
									  datax    = datax,
									  lbl      = lbl,
									  color    = color)
					plotData[res] = datax	

		plt.legend()
		plt.xlabel('{} D{} per atom'.format(normType,metric), fontsize = axesFont)
		plt.ylabel('Normed-frequency', fontsize = axesFont)
		sns.despine()
		if plotTitle == '':
			t = '{} D{} per atom, residues: {}'.format(normType,metric,resiType)
		else:
			t = plotTitle
		plt.title(t,fontsize = titleFont)

		if not save: 
			plt.show()
		else:
			saveName = '{}DistnPlot_Residues-{}_Metric-D{}_Normalisation-{}'.format(outputDir,
																					''.join(resiType),
																					metric,
																					normType.replace(" ",""))
			if valType != 'all':
				if valType == 'average':
					saveName += '_{}'.format(valType)
				else:
					saveName += '_Dataset-{}'.format(valType)

			fileName = self.checkUniqueFileName(fileName = saveName,
									 			fileType = fileType)
			fig.savefig(fileName)

		return plotData

	def plotHist(self,
				 plotType = 'both',
				 nBins    = 300,
				 datax    = [],
				 lbl      = '',
				 color    = 'b'):

		# plot histogram or kde plot for datax and give current label
		# 'nBins' is number of bins (only used if plotType is 'hist' or 'both')

		if seabornFound is False:
			return

		if plotType == 'hist':
			plt.hist(datax,
					 nBins,
					 histtype = "stepfilled",
					 alpha    = .7,
					 label    = lbl,
					 color=color)
		elif plotType == 'kde':
			sns.kdeplot(np.array(datax), shade=True, label=lbl, color=color)
		elif plotType == 'both':
			sns.distplot(np.array(datax), label=lbl, color=color)

	def graphMetric(self,
					densMet   = 'loss',
					normType  = 'Standard',
					atomType  = '',
					restype   = '',
					resiNum   = '',
					errorBars = 'NONE',
					saveFig   = False,
					outputDir = './',
					fileType  = '.svg',
					axesFont  = 18,
					titleFont = 24,
					palette   = 'hls'):
	
		# produce a graph of selected metric against dataset number 
		# for a specified atom. 'errorBars' takes values in ('NONE',
		# 'RESNUM','ATOMTYPE')

		if seabornFound is False:
			return

		errorOptions = ('NONE',
					    'RESNUM',
					    'ATOMTYPE')

		if errorBars not in errorOptions: 
			print 'invalid error parameter..\n choose from: {}'.format(errorOptions)
			return

		# find atoms of interest 
		foundAtoms = self.getAtom(restype  = restype,
								  resnum   = resiNum,
								  atomtype = atomType)
		if len(foundAtoms) == 0:
			print 'No atoms found..'
			return
			
		# define x range here (number in damage series)
		x = range(self.getNumDatasets())

		# determine y values here dependent on density metric type specified 
		if errorBars == 'NONE':

			sns.set_palette(palette  = palette,
							n_colors = len(foundAtoms),
							desat    = .6)
			sns.set_context(rc={"figure.figsize":(10, 10)})
			f = plt.figure()	

			for atom in foundAtoms:

				if self.checkMetricPresent(atom=atom,metric=densMet,normType=normType) is False: 
					return # check metric valid

				y = atom.densMetric[densMet][normType]['values']
				plt.plot(x,y,label=atom.getAtomID())

		else: # error bars will be plotted
			errorBarKey = {'ATOMTYPE' : 'atomtype',
						   'RESNUM'   : 'residuenum'}
			errorAttr 	= errorBarKey[errorBars]
			yDict 		= {}
			for atom in foundAtoms:

				if self.checkMetricPresent(atom = atom,metric = densMet,normType = normType) is False: 
					return # check metric valid

				vals = atom.densMetric[densMet][normType]['values']

				if getattr(atom,errorAttr) not in yDict.keys():
					yDict[getattr(atom,errorAttr)] = [vals]
				else:
					yDict[getattr(atom,errorAttr)].append(vals)

			sns.set_palette(palette  = palette,
							n_colors = len(yDict.keys()),
							desat    = .6)
			sns.set_context(rc={"figure.figsize":(10, 10)})
			f = plt.figure()

			for key in yDict.keys():
				plt.errorbar(x,
							 np.mean(yDict[key],0),
							 yerr     = np.std(yDict[key],0),
							 fmt      = '-o',
							 capthick = 2,
							 label    = key)

		plt.xlabel('Dataset', fontsize = axesFont)
		plt.ylabel('{} D{}'.format(normType,densMet), fontsize = axesFont)
		plt.legend()

		args = [normType,
				densMet,
				restype,
				resiNum,
				atomType]

		f.suptitle('{} D{}: {} {} {}'.format(*args), fontsize = titleFont)
		if saveFig is False:
			plt.show()
		else:
			name = '{}D{}-{}-{}-{}'.format(*args)
			if errorBars != 'NONE':
				name += 'witherrorbars'

			fileName = self.checkUniqueFileName(fileName = outputDir+name,
								 				fileType = fileType)
			fig.savefig(fileName)

	def plotSusceptibleAtoms(self,
							 densMet   = 'loss',
							 errorbars = True,
							 saveFig   = False,
							 susAtms   = [],
							 fileType  = '.svg',
							 axesFont  = 18,
							 titleFont = 24):

		# create a line plot of metric values for susceptible atoms within
		# structure. Susceptible atoms are defined as below

		if seabornFound is False:
			return

		if susAtms == []: 
			susAtms = [['GLU','CD'],
					   ['ASP','CG'],
					   ['TYR','OH'],
					   ['CYS','SG'],
					   ['MET','SD']]
		sDic = {}
		for s in susAtms:
			findAtms = self.getAtom(restype  = s[0],
									atomtype = s[1])
			if findAtms is False: 
				continue
			if not isinstance(findAtms,list): 
				findAtms = [findAtms]
			sDic['-'.join(s)] = findAtms
		
		# define x range here (number in damage series)
		x = range(self.getNumDatasets())
		sns.set(style="white", context="talk")
		f, axes = plt.subplots(1, 1, figsize=(12, 12), sharex=True)
		x = range(self.getNumDatasets())
		cmap = cm.Set1
		i = -1
		nsteps = len(sDic.keys())
		for key in sDic.keys():
			i += 1
			if errorbars is False:
				for atom in sDic[key]:
					y = atom.densMetric[densMet][normType]['values']
					plt.plot(x,y,label=key,color=cmap(i/float(nsteps)))
			else:
				yData = np.mean([atom.densMetric[densMet][normType]['values'] for atom in sDic[key]],0)
				yErr  = np.std([atom.densMetric[densMet][normType]['values'] for atom in sDic[key]],0)
				plt.errorbar(x,
							 yData,
							 yerr     = yErr,
							 fmt      = '-o',
							 capthick = 2,
							 label    = key,
							 color    = cmap(i/float(nsteps)))

		av,std = self.getAverageMetricVals(densMet,normType) # get structure-wide average and std dev
		plt.errorbar(x,
					 av,
					 yerr     = std,
					 fmt      = ':o',
					 capthick = 2,
					 color    = 'k') 

		plt.xlabel('Dataset', fontsize = axesFont)
		plt.ylabel('{} D{}'.format(normType,densMet), fontsize = axesFont)
		plt.legend()
		f.suptitle('{} D{}: susceptible residues'.format(normType,densMet),fontsize = titleFont)
		if saveFig is False:
			plt.show()
		else:
			i = 0
			figName = '{}{}_D{}_susceptResis'.format(self.outputDir,
													 normType,
													 densMet)

			fileName = self.checkUniqueFileName(fileName = figName,
								 				fileType = fileType)
			fig.savefig(fileName)

	def plotNumAtomsWithMetricAboveStructureWideMean(self,
													 startThreshold = 0.5,
													 thresholdDiv   = 0.05,
													 metric         = 'loss',
													 normType       = 'Standard',
													 dataset        = 0,
													 axesFont       = 18,
													 saveFig        = True,
													 fileType       = '.svg',
													 outputLoc      = './',
													 titleFont      = 18):

		# plot the number of atoms present within the structure
		# with metric value of n standard deviations above the
		# the structure wide mean (for increasing n until no 
		# atoms remain)

		n = self.numAtmsWithMetricAboveLevel(dataset   = dataset,
											 metric    = metric,
											 normType  = normType,
											 threshold = startThreshold)
		numAtomsFound = [n]
		numStds       = [startThreshold]
		t = startThreshold
		while n != 0:
			t += thresholdDiv
			n = self.numAtmsWithMetricAboveLevel(dataset   = dataset,
										         metric    = metric,
										         normType  = normType,
										         threshold = t)
			numAtomsFound.append(n)
			numStds.append(t)

		sns.set_context(rc={"figure.figsize":(10, 10)})
		f = plt.figure()
		plt.plot(numStds,numAtomsFound)

		plt.xlabel('# standard deviations from structure mean', fontsize = axesFont)
		plt.ylabel('Number of atoms located', fontsize = axesFont)
		f.suptitle('# atoms with metric value of n standard deviations from structure-wide mean',fontsize = titleFont)
		if saveFig is False:
			plt.show()
		else:
			figName = '{}NumAtomsVsStdsFromStructureMean_metric-{}_norm-{}_dataset-{}'.format(outputLoc,
																 						      metric,
																 							  normType,
																 							  dataset)

			fileName = self.checkUniqueFileName(fileName = figName,
								 				fileType = fileType)
			f.savefig(fileName)

		return fileName

	def checkMetricPresent(self,
						   atom     = '',
						   metric   = 'loss',
						   normType = 'Standard'):

		# method to check whether a method is valid 
		# for selected atom object

		if atom == '':
			atom = self.atomList[0]
		try:
			atom.densMetric[metric][normType]['values']
		except NameError:
			print 'Unexpected density metric name.."{} {}"'.format(normType,metric)
			return False
		return True

	def calcHalfDose(self,
					 chain    = '',
					 restype  = '',
					 resnum   = '',
					 atomtype = '',
					 densMet  = 'loss',
					 normType = 'Standard',
					 shift    = True,
					 offset   = True):

		# call the external halfDoseCalc class to calculate a rough half 
		# dose decay value for a specified atom. 'shift' (Bool) determines 
		# whether half dose is dose for density to reach half of initial 
		# density (=False) or half dose is dose for density to reach 
		# av(initial density,end density limit (=True). 'offset' (Bool)
		# determines whether exponential decay function is allowed 
		# a non-zero end density value

		atom = self.getAtom(chain    = chain,
							restype  = restype,
							resnum   = resnum,
							atomtype = atomtype)

		atom = atom[0] # convert from list
		if len(self.doseList) != self.getNumDatasets():
			return 'need to specify doses list as class attribute before this can be calculated'
		halfDoseApprox(atom            = atom,
					   atoms           = self.atomList,
					   doses           = self.doseList,
					   densityMetric   = densMet,
					   normType        = normType,
					   plotDir         = self.outputDir,
					   shiftedHalfDose = shift,
					   zeroOffset      = offset)

	def calcHalfDoseForAtomtype(self,
								restype  = '',
								atomtype = '',
								densMet  = 'loss',
								normType = 'Standard',
								shift    = True,
								offset   = True,
								n        = 1,
								fileType = '.svg'):

		# call the above calcHalfDose method for instances 
		# of restype and atomtype

		atoms = self.getAtom(restype  = restype,
							 atomtype = atomtype)
		for atom in atoms:
			self.calcHalfDose(chain    = atom.chaintype,
							  restype  = restype,
							  resnum   = atom.residuenum,
							  atomtype = atomtype,
							  densMet  = densMet,
							  normType = normType,
							  shift    = shift,
							  offset   = offset)

		print '--------------------'
		print 'Summary here of run:'

		residuals = [atom.densMetric[densMet][normType]['Half-dose']['Residuals'] for atom in atoms]
		certainties = [atom.densMetric[densMet][normType]['Half-dose']['Certainty'] for atom in atoms]
		residualThreshold = np.mean(residuals) + n*np.std(residuals)
		certaintyThreshold = np.mean(certainties) - n*np.std(certainties)
		print 'Will remove atoms with residual score > {} or certainty score < {}'.format(residualThreshold,certaintyThreshold)

		xData,yData = [],[]
		count = 0
		for atom in atoms:
			halfDoseStatDic = atom.densMetric[densMet][normType]['Half-dose']
			print '{} --> {}'.format(halfDoseStatDic['Initial density'],halfDoseStatDic['Half-dose'])

			if (halfDoseStatDic['Residuals'] <= residualThreshold and 
				halfDoseStatDic['Certainty'] >= certaintyThreshold):
				xData.append(halfDoseStatDic['Initial density'])
				yData.append(halfDoseStatDic['Half-dose'])
			else: 
				count += 1
				print 'atom {} not included in final plot:'.format(atom.getAtomID())
				print '\tResidual: {}, Certainty: {}'.format(halfDoseStatDic['Residuals'],halfDoseStatDic['Certainty'])
		print '{}/{} atoms not included in final scatter plot'.format(count,len(atoms))

		# plot the relationship between intitial density and half-dose
		self.plotScatterPlot(xData     = xData,
							 yData     = yData,
							 xLabel    = 'Initial density',
							 ylabel    = 'Half-dose',
							 figtitle  = 'D{}: {}-{}'.format(densMet,restype,atomtype),
							 fileType  = fileType,
							 saveTitle = 'D{}Scatter_{}-{}_thres{}'.format(densMet,restype,atomtype,n))

	def compareAtomHalfDoses(self,
						     restype   = '',
						     atomtype1 = '',
						     atomtype2 = '',
						     densMet   = 'loss',
						     normType  = 'Standard',
						     n         = 1,
						     fileType  = '.svg'):

		# after above method calcHalfDoseForAtomtype is run for 
		# two atom types atomtype1 and atomtype2 can compare 
		# half-doses between two atoms within same side-chain

		atoms1 = self.getAtom(restype  = restype,
							  atomtype = atomtype1)
		atoms2 = self.getAtom(restype  = restype,
							  atomtype = atomtype2)

		residuals = [atom.densMetric[densMet][normType]['Half-dose']['Residuals'] for atom in atoms1]
		certainties = [atom.densMetric[densMet][normType]['Half-dose']['Certainty'] for atom in atoms2]
		residualThreshold = np.mean(residuals) + n*np.std(residuals)
		certaintyThreshold = np.mean(certainties) - n*np.std(certainties)

		xData,yData = [],[]
		count = 0
		for atom1 in atoms1:
			for atom2 in atoms2:
				if (atom1.chaintype == atom2.chaintype and
					atom1.residuenum == atom2.residuenum):
					if (atom1.densMetric[densMet][normType]['Half-dose']['Residuals'] <= residualThreshold and 
						atom1.densMetric[densMet][normType]['Half-dose']['Certainty'] >= certaintyThreshold and
						atom2.densMetric[densMet][normType]['Half-dose']['Residuals'] <= residualThreshold and 
						atom2.densMetric[densMet][normType]['Half-dose']['Certainty'] >= certaintyThreshold):
						xData.append(atom1.densMetric[densMet][normType]['Half-dose']['Half-dose'])
						yData.append(atom2.densMetric[densMet][normType]['Half-dose']['Half-dose'])
					else:
						count += 1
					break
		print '{}/{} atoms not included in final scatter plot'.format(count,len(atoms1))
		
		# plot this data as scatter plot
		self.plotScatterPlot(xData     = xData,
							 yData     = yData,
							 xLabel    = '{}-{} Half-dose'.format(restype,atomtype1),
							 ylabel    = '{}-{} Half-dose'.format(restype,atomtype2),
							 figtitle  = 'D{} half-dose for {} atoms'.format(densMet,restype),
							 fileType  = fileType,
							 saveTitle = 'D{}_halfDose{}_{}-{}'.format(densMet,restype,atomtype1,atomtype2))

	def HalfDoseDistn(self,
					  metric   = 'loss',
					  normType = 'Standard',
					  plotType = 'kde',
					  restype  = '',
					  atomtype = '',
					  save     = False,
					  nBins    = 300,
					  colour   = 'b',
					  fileType = '.svg'):

		# temporary method to do a specific plot.
		# histogram/kde plot of density metric per atom.
		# plotType is 'histogram' or 'kde'

		if seabornFound is False:
			return

		if plotType not in ('hist','kde','both'): 
			return 'Unknown plotting type selected.. cannot plot..'
		if self.checkMetricPresent(metric=metric,normType=normType) is False: 
			return # check metric valid

		sns.set_palette("deep", desat=.6)
		sns.set_context(rc={"figure.figsize": (10, 6)})
		fig = plt.figure()

		atoms = self.getAtom(restype  = restype,
							 atomtype = atomtype)

		datax = [atm.densMetric[metric][normType]['Half-dose']['Half-dose'] for atm in atoms if atm.densMetric[metric][normType]['Half-dose']['Half-dose'] < 100 ]
		self.plotHist(plotType=plotType,nBins=nBins,datax=datax,
					  lbl='',colour=colour)
		plt.xlabel('Half-dose (MGy)')
		plt.ylabel('Normed frequency')
		plt.title('Half-dose distribution: {}-{}'.format(restype,atomtype))
		if not save: 
			plt.show()
		else:
			figName = '{}halfDoseDistn-{}{}'.format(self.outputDir,restype,atomtype)
			fileName = self.checkUniqueFileName(fileName = figName,
								 				fileType = fileType)
			fig.savefig(fileName)

	def compareSensAtoms(self,
						 densMet  = 'loss',
						 normType = 'Standard'):

		rSquaredDic = {}
		numPairsDic = {}
		dataDic 	= {}
		sensAtoms = [['GLU','CD','CG'],
					 ['GLU','CD','OE1'],
					 ['ASP','CG','CB'],
					 ['ASP','CG','OD1'],
					 ['TYR','OH','CZ'],
					 ['TYR','CZ','CE2']]
		for i in sensAtoms:
			rSquared,numPairs,data = self.compareMetricsBetweenAtoms(restype   = i[0],
																	 atomtype1 = i[1],
																	 atomtype2 = i[2],
																	 densMet   = densMet,
																	 normType  = normType,
																	 dataset   = 'av')
			if rSquared != False: 
				rSquaredDic['-'.join(i)] = rSquared
				numPairsDic['-'.join(i)] = numPairs
				dataDic['_'.join(i)] 	 = data
		return rSquaredDic,numPairsDic,dataDic

	def compareMetricsBetweenAtoms(self,
								   restype   = '',
								   atomtype1 = '',
								   atomtype2 = '',
								   densMet   = 'loss',
								   normType  = 'Standard',
								   fileType  = '.svg',
								   dataset   = 0):

		# for two atom types atomtype1 and atomtype2, compare 
		# density metric between two atoms within same side-chain

		atoms1 = self.getAtom(restype  = restype,
							  atomtype = atomtype1)

		atoms2 = self.getAtom(restype  = restype,
							  atomtype = atomtype2)

		if atoms1 is False or atoms2 is False: 
			return False,0,{'x':[],'y':[]}

		if dataset == 'av': 
			self.calcAdditionalMetrics(metric    = densMet,
									   normType  = normType,
									   newMetric = 'average')

		xData,yData = [],[]
		for atom1 in atoms1:
			for atom2 in atoms2:
				if (atom1.chaintype == atom2.chaintype and
					atom1.residuenum == atom2.residuenum):
					if dataset != 'av':
						xData.append(atom1.densMetric[densMet][normType]['values'][dataset])
						yData.append(atom2.densMetric[densMet][normType]['values'][dataset])
					else:
						xData.append(atom1.densMetric[densMet][normType]['average'])
						yData.append(atom2.densMetric[densMet][normType]['average'])
					break
			
		# plot this data as scatter plot
		rSquared = self.plotScatterPlot(xData       = xData,
										yData       = yData,
										xLabel      = '{}-{} D{}'.format(restype,atomtype1,densMet),
							 			ylabel      = '{}-{} D{}'.format(restype,atomtype2,densMet),
							 			figtitle    = 'D{} for {} atoms'.format(densMet,restype),
							 			fileType    = fileType,
							 			saveTitle   = 'D{}_{}-{}-{}_scatterplot_{}'.format(densMet,restype,atomtype1,atomtype2,dataset),
							 			lineBestFit = True)

		data = {'x':xData,'y':yData}
		return rSquared,len(xData),data

	def compareMetrics(self,
					   restype   = '',
					   atomtype  = '',
					   metric1   = 'loss',
					   metric2   = 'loss',
					   normType1 = 'Standard',
					   normType2 = 'Standard',
					   dSet      = 0,
					   fileType  = '.svg'):

		# for all atoms of type restype and atomtype plot scatter plot 
		# comparing metric1 and 2 for chosen dataset dSet set restype = '' 
		# and atomtype = '' to include all atoms

		atoms = self.getAtom(restype  = restype,
							 atomtype = atomtype)

		xData,yData = [],[]
		for atom in atoms:
			xData.append(atom.densMetric[metric1][normType1]['values'][dSet])
			yData.append(atom.densMetric[metric2][normType2]['values'][dSet])

		rSquared = self.plotScatterPlot(xData       = xData,
										yData       = yData,
										xLabel      = '{} D{}'.format(normType1,metric1),
										ylabel      = '{} D{}'.format(normType2,metric2),
										fileType    = fileType,
										figtitle    = '{} D{} vs {} D{} for res:{} atoms:{}'.format(normType1,metric1,normType2,metric2,restype,atomtype),
										saveTitle   = '{}D{}-vs-{}D{}-{}-{}-{}'.format(normType1,metric1,normType2,metric2,restype,atomtype,dSet),
										lineBestFit = True)

	def plotScatterPlot(self,
						xData        = [],
						yData        = [],
						xLabel       = '',
						ylabel       = '',
						figtitle     = '',
						saveTitle    = 'untitled',
						fileType     = '.svg',
						lineBestFit  = False,
						yequalsxLine = False,
						colors       = '#50527A',
						axesFont     = 18,
						titleFont    = 24):

		# plot the relationship between xData and yData.
		# if 'lineBestFit' is True then linear line of best fit 
		# is calculated and plotted. if 'yequalsxLine' is 
		# True then line y=x plotted

		if seabornFound is False:
			return

		sns.set_context("talk")
		fig = plt.figure(figsize=(10, 10))
		ax = plt.subplot(111)

		plt.scatter(xData,
					yData,
					marker     = 'o',
					s          = 100,
					c          = colors,
					edgecolors = '#FFFFFF',
					cmap       = 'copper')

		if lineBestFit is True: # plot line of best if specified
			slope, intercept, r_value, p_value, std_err = stats.linregress(xData,yData)
			x = np.linspace(min(xData), max(xData),10)
			y = slope*np.array(x) + np.array(intercept)
			plt.plot(x,
					 y,
					 '-',
					 color     = '#d11141',
					 linewidth = 3)
			figtitle += ', R^2:{}'.format(round(r_value**2,2))

		if yequalsxLine is True:
			x = np.linspace(min(xData), max(xData),10)
			y = x
			plt.plot(x, 
					 y,
					 ':',
					 color     = '#9890af',
					 linewidth = 3)

		plt.xlabel(xLabel, fontsize = axesFont)
		plt.ylabel(ylabel, fontsize = axesFont)
		fig.suptitle(figtitle, fontsize = titleFont)
		sns.despine()

		figName = self.outputDir+saveTitle
		fileName = self.checkUniqueFileName(fileName = figName,
								 			fileType = fileType)
		fig.savefig(fileName)

		if lineBestFit is True: 
			return r_value**2

	def damageVdistPlot(self,
						type1    = '',
						type2    = '',
						dataset  = 0,
						densMet  = 'loss',
						normType = 'Standard',
						fileType = '.svg'):

		# for atom of type1 determine a correlation between the closest 
		# distance to atom type2 and the damage to atom type1.
		# type1,type2 of form [['GLU','CD'],['ASP','CG'],..] depending 
		# on how many atom types required

		atoms1,atoms2 = [],[]
		for ind in type1: atoms1 += self.getAtom(restype  = ind[0],
												 atomtype = ind[1])
		for ind in type2: atoms2 += self.getAtom(restype  = ind[0],
												 atomtype = ind[1])

		densVals,minDists = [],[]
		for atom1 in atoms1:
			densVals.append(atom1.densMetric[densMet][normType]['values'][dataset])
			distList = []
			for atom2 in atoms2:
				distList.append(self.getDistanceBetweenAtoms(atom1 = atom1,
															 atom2 = atom2))
			minDists.append(min(distList))

		yLabel = '{} D{}'.format(normType,densMet)
		xLabel = '{} {} distance'.format(type1,type2)
		figTitle = '{} D{} against {} distance'.format(type1,densMet,type2)
		saveTitle = figTitle.strip()

		self.plotScatterPlot(xData     = minDists,
							 yData     = densVals,
							 xLabel    = xLabel,
							 yLabel    = yLabel,
							 figtitle  = figTitle,
							 fileType  = fileType,
							 saveTitle = saveTitle)

	def getDistanceBetweenAtoms(self,
							    atom1='',
							    atom2=''):

		# get distance between two atoms

		xyz1 = np.array([atom1.X_coord,atom1.Y_coord,atom1.Z_coord])
		xyz2 = np.array([atom2.X_coord,atom2.Y_coord,atom2.Z_coord])
		return np.linalg.norm(xyz1-xyz2)

	def getAtomsWithinDist(self,atom,distLim):

		# find all atoms within 'distLim' angstroms 
		# of selected atomtype

		nearAtmDic = {'atoms':[],'distances':[]}
		for otheratom in self.atomList:
			if (otheratom.getAtomID()).split('-')[:-1] == (atom.getAtomID()).split('-')[:-1]: 
				continue
			dist = self.getDistanceBetweenAtoms(atom1 = atom,
												atom2 = otheratom)
			if dist < distLim:
				nearAtmDic['atoms'].append(otheratom)
				nearAtmDic['distances'].append(dist)
		print '{} in total within {} Angstrom of {}'.format(len(nearAtmDic['atoms']),distLim,atom.getAtomID())
		return nearAtmDic

	def densMetSurroundAtmsCorrel(self,
								  restype       = '',
								  atomtype      = '',
								  distLim       = 4,
								  normType      = 'Standard',
								  densMet       = 'loss',
								  crystConts    = True,
								  pdbName       = '',
								  symmetrygroup = '',
								  fileType      = '.svg'):

		# determine whether there is a correlation between density metric 
		# and types of surrounding atoms if 'crystConts' is True then crystal
		# contacts will also be located here

		if crystConts is True:
			from NCONTjob import NCONTjob

			ncont = NCONTjob(pdbName,self.outputDir,symmetrygroup,self.seriesName)
			success = ncont.run()
			if success is False: return
			logFile = self.outputDir+'/'+ncont.outputLogfile

		atoms = self.getAtom(restype  = restype,
							 atomtype = atomtype)
		groupA,groupB,groupAContacts = [],[],[]
		for atom in atoms:
			nearCarboxyl = False

			if crystConts is True:
				atomFound = False
				readLog = open(logFile,'r')
				for l in readLog.readlines():
					if atomFound is True and len(l[0:5].strip()) != 0: 
						atomFound = False # reset if new atom info reached

					if (atom.basetype == l[11:14] and
						str(atom.residuenum) == l[7:10].strip() and
						atom.chaintype == l[4] and
						atom.atomtype == l[19:21]):
						atomFound = True
						print 'found atom: {}'.format(atom.getAtomID())

					if atomFound is True:
						if (l[39:42] in ('GLU','ASP') and
							l[47:50] in ['OE1','OE2','OD1','OD2'] and
							float(l[58:62]) < distLim):
							resNum = int(l[35:38].strip())
							chain = l[32]
							print 'contact found! - {}-{}-{}-{}'.format(chain,resNum,l[39:42],l[47:50])
							contactAtom = self.getAtom(chain    = chain,
													   restype  = l[39:42],
													   resnum   = resNum,
													   atomtype = l[47:50])
							nearCarboxyl = True
							break
				readLog.close()
			else:
				nearAtms = self.getAtomsWithinDist(atom,distLim)
				nearCarboxyl = False
				for atom2 in nearAtms['atoms']:
					if atom2.atomtype in ['OE1','OE2','OD1','OD2']:
						if atom2.basetype in ['GLU','ASP']:
							contactAtom = atom2
							nearCarboxyl = True
							break

			if nearCarboxyl is True:
				groupA.append(atom)
				groupAContacts.append(contactAtom[0])
			else:
				groupB.append(atom)

		densValsA 		  = [atom.densMetric[densMet][normType]['average'] for atom in groupA]
		densValsAContacts = [atom.densMetric[densMet][normType]['average'] for atom in groupAContacts]
		densValsB 		  = [atom.densMetric[densMet][normType]['average'] for atom in groupB]
		strA = 'For the {} atoms near carboxyl groups, mean avDloss: {}, std: {}'.format(len(densValsA),np.mean(densValsA),np.std(densValsA))
		strB = 'For the {} atoms not near carboxyl groups, mean avDloss: {}, std: {}'.format(len(densValsB),np.mean(densValsB),np.std(densValsB))
		print strA+'\n'+strB

		# determine how solvent accessibility changes for each of these atoms upon Glu/Asp decarboxylation.
		# this will dictate the colour scheme for the resulting scatter plot
		scatterColor = []
		solvAccDic,plotData = self.compareSolvAccessWithAndWithoutGluAspGroups(pdbName       = pdbName,
																			   symmetrygroup = symmetrygroup,
																			   resType       = restype,
																			   atomType      = atomtype,
																			   specialAtoms  = groupA,
																			   densMet       = densMet,
																			   normType      = normType)
	
		count = 0
		for atm in groupA: 
			solvAcc = solvAccDic[atm.getAtomID()]
			diff = float(solvAcc[1]) - float(solvAcc[0])

			scatterColor.append(float(solvAcc[1]))

			if diff == 0:
				print 'No change in solvent accessibility'
				# scatterColor.append('r')
			else:
				print 'Solvent accessibility reduced by {} upon nearby decarboxylation'.format(diff)
				# scatterColor.append('b')
				count += 1
		if count < len(groupA):
			print 'Not all {}-{} atoms experience change in solvent accessibility upon nearby decarboxylation'.format(restype,atomtype)

		# plot the relationship between TYR-OH density and nearby carboxyl atom density
		self.plotScatterPlot(xData    = densValsA,
							 yData    = densValsAContacts,
							 xLabel   = '{}-{} D{}'.format(restype,atomtype,densMet),
							 ylabel   = 'Carboxyl-contact D{}'.format(densMet),
							 figtitle  = '{}-{} vs carboxyl-contact D{}'.format(restype,atomtype,densMet),
							 saveTitle = 'D{}Scatter_{}-{}_carboxylContacts'.format(densMet,restype,atomtype),
							 fileType  = fileType,
							 colors    = scatterColor)

		return strA+'\n'+strB,densValsA,densValsAContacts,scatterColor

	def susceptAtmComparisonBarplot(self,
									metric   = 'loss',
									normType = 'Standard',
									dataset  = 0,
									set      = 1,
									box      = 'Box',
									fileType = '.svg',
									axesFont = 18):

		# produce a barplot to compare the damage metric of 
		# susceptible atom types at a given dataset

		if seabornFound is False:
			return

		if set == 1:
			keyAtomTypes = [['GLU','CD'],
							['ASP','CG'],
							['TYR','OH'],
							['CYS','SG'],
							['MET','SD'],
							['LYS','NZ'],
							['ARG','NE'],
							['ASN','OD1'],
							['GLN','OE1'],
							['ALA','CB']]
		if set == 2:
			keyAtomTypes = [['GLU','OE1'],
							['GLU','OE2'],
							['ASP','OD1'],
							['ASP','OD2'],
							['TYR','OH'],
							['ASN','OD1'],
							['GLN','OE1'],
							['SER','OG'],
							['THR','OG1'],
							['','O']]
	
		plotData = {'x'       : [],
					'y'       : [],
					'y-error' : []}

		for atm in keyAtomTypes:
			foundAtms = self.getAtom(restype  = atm[0],
									 atomtype = atm[1])
			if foundAtms == False: 
				continue 
			for foundAtm in foundAtms:
				atmId 	= '-'.join(atm)
				plotData['x'].append(atmId)
				plotData['y'].append(foundAtm.densMetric[metric][normType]['values'][dataset])

		fig = plt.figure()
		sns.set_style("whitegrid")
		if box == 'Box':
			ax = sns.boxplot(x    = "x",
							 y    = "y",
							 data = plotData)
		elif box == 'Bar':
			ax = sns.barplot(x    = "x",
							 y    = "y",
							 data = plotData)
		else:
			print 'Unknown graph type - specify "Box" or "Bar"'
			return

		plt.xlabel('Atom type', fontsize = axesFont)
		plt.ylabel('{} D{}'.format(normType,metric), fontsize = axesFont)

		figName = '{}SusceptAtms{}plot-{}-D{}_{}-set{}'.format(self.outputDir,
															   box,
															   normType,
															   metric,
															   dataset,
															   set)
		fileName = self.checkUniqueFileName(fileName = figName,
								 			fileType = fileType)
		fig.savefig(fileName)

	def compareSolvAccessWithAndWithoutGluAspGroups(self,
													pdbName       = '',
													symmetrygroup = '',
													resType       = '',
													atomType      = '',
													specialAtoms  = '',
													densMet       = 'loss',
													normType      = 'Standard',
													fileType      = '.svg'):

		# determine the change in solvent accessibility for 'resType'-'atomType' 
		# atoms before & after Glu/Asp decarboxylation events 'specialAtoms' 
		# is an optional parameter, which specifies a subset of atoms objects 
		# which should be coloured differently on the resulting scatter plot 
		# below - the idea is to run .densMetSurroundAtmsCorrel() above first 
		# to get the list of atoms exhibiting contacts to Glu/Asp CO2 groups 
		# and then colour these differently below

		from ccp4Job import ccp4Job

		pdbIn = open(pdbName,'r')
		strippedPdb = pdbName.strip('.pdb')+'_noGluAspCO2groups.pdb'
		pdbOut = open(strippedPdb,'w')
		for l in pdbIn.readlines():
			if l.split()[0] not in ('ATOM','ANISOU'):
				pdbOut.write(l)
			else:
				if l[13:21].replace(' ','') not in ('CDGLU','OE1GLU','OE2GLU','CGASP','OD1ASP','OD2ASP'):
					pdbOut.write(l)
		pdbIn.close()
		pdbOut.close()

		#  run AREAIMOL in CCP4 to calculate solvent accessibility for each pdb file
		count = 0
		for pdb in (pdbName,strippedPdb):
			count += 1
			commandInput1 = 'areaimol XYZIN "{}" XYZOUT "{}"'.format(pdb,pdb.strip('.pdb')+'_areaimol.pdb')
			commandInput2 = 'DIFFMODE OFF\nMODE -\nNOHOH\nSMODE IMOL\nREPORT -\nCONTACT -\n'+\
							'YES -\nRESAREA -\nYES\nPNTDEN 10\nPROBE 1.4\nOUTPUT\nEND'
			run = ccp4Job('areaimol',commandInput1,commandInput2,self.outputDir,'areaimolRun{}.txt'.format(count),'')

		# these are names of resulting pdb files
		fullPDB 	= pdbName.strip('.pdb')+'_areaimol.pdb'
		partialPDB 	= strippedPdb.strip('.pdb')+'_areaimol.pdb'

		# get Tyr-OH atoms in structure
		TyrOHatms = self.getAtom(restype  = resType,
								 atomtype = atomType)

		solvAccDic = {}
		for atm in TyrOHatms:
			solvAccDic[atm.getAtomID()] = []

			# get other Tyr-ring atoms for each Tyr-OH atom
			TyrRingAtms = []
			if resType == 'TYR' and atomType == 'OH':
				for atmtype in ['CG','CD1','CD2','CE1','CE2','CZ']:
					TyrRingAtms.append(self.getAtom(chain    = atm.chaintype,
													restype  = resType,
													resnum   = atm.residuenum,
													atomtype = atmtype)[0])

			for pdb in (fullPDB,partialPDB):
				print pdb
				solvAcc = []
				atms = [atm]+TyrRingAtms
				for a in atms:
					solvAcc.append(float(a.findSolventAccessibility(pdb)))
				solvAccDic[atm.getAtomID()].append(np.mean(solvAcc))

		# remove pdb files since no longer needed after job run
		for f in [fullPDB,partialPDB]:
			os.system('rm {}'.format(f))
			os.remove(f)

		plotData = {'x':[],'y':[],'colours':[]}
		for atm in TyrOHatms:
			if atm in specialAtoms:
				plotData['colours'].append('r')
			else:
				plotData['colours'].append('b')
			key = atm.getAtomID()

			plotData['y'].append(atm.densMetric[densMet][normType]['average'])
			plotData['x'].append(float(solvAccDic[key][1]))
			print '{}: {} --> {} ... Dloss: {}'.format(key,round(solvAccDic[key][0],2),round(solvAccDic[key][1],2),plotData['y'][-1])

		# plot the relationship between TYR-OH density and nearby carboxyl atom density
		self.plotScatterPlot(xData     = plotData['x'],
			                 yData     = plotData['y'],
							 xLabel    = 'Solvent Accessibility',
							 ylabel    = '{}-{} D{}'.format(resType,atomType,densMet),
							 fileType  = fileType,
							 saveTitle = 'D{}Scatter_{}-{}_vsSolvAccChange'.format(densMet,resType,atomType),
							 colors    = plotData['colours'])

		return solvAccDic,plotData

	def calculateLocalDloss(self,
							resType  = '',
							atomType = '',
							distance = 4,
							densMet  = 'loss',
							normType = 'Standard',
							weighted = False,
							fileType = '.svg'):

		# for a given residue, calculate the average density metric 
		# (for a specified metric) within a given region of space 
		# around that particular residue. if 'weighted' is True 
		# then a distance-weighted mean is instead calculated

		plotData = {'x':[],'y':[]}
		# find specified atoms within structure
		atms = self.getAtom(restype  = resType,
							atomtype = atomType)
		if atms is False: 
			return plotData # don't continue if no atoms exist
		self.calcAdditionalMetrics(metric    = densMet,
								   normtype  = normType,
								   newMetric = 'average')

		for atm in atms:
			nearAtms = self.getAtomsWithinDist(atm,distance)

			densList,distList = [],[]
			i = -1
			for nearAtm in nearAtms['atoms']:
				i += 1
				dist = nearAtms['distances'][i]
				dens = nearAtm.densMetric[densMet][normType]['average']
				densList.append(dens)
				distList.append(dist)
			if len(densList) == 0: 
				continue
			if weighted is True:
				meanDens = np.average(densList, weights=np.square(np.reciprocal(distList)))
			else:
				meanDens = np.mean(densList)

			dens = atm.densMetric[densMet][normType]['average']
			print '{}: D{}: {} --> local D{}: {}'.format(atm.getAtomID(),densMet,dens,densMet,meanDens)
			plotData['x'].append(dens)
			plotData['y'].append(meanDens)

		# plot the relationship between atom Dloss and local environment average Dloss
		self.plotScatterPlot(xData     = plotData['x'],
			                 yData     = plotData['y'],
							 xLabel    = '{}-{} D{}'.format(resType,atomType,densMet),
							 ylabel    = 'local environment D{}'.format(densMet),
							 fileType  = fileType,
							 figtitle  = 'Average D{} within {} Angstrom of {}-{}'.format(densMet,distance,resType,atomType),
							 saveTitle = 'D{}Scatter_{}-{}_localEnvironmentComparison'.format(densMet,resType,atomType))

		return plotData

	def getAvMetricPerAtmInRes(self,
							   resType  = '',
							   densMet  = 'loss',
							   normType = 'Standard',
							   dataset  = 0):

		# for specified residue type, calculate the average 
		# metric for each atom type within the residue

		atms = self.getAtom(restype=resType)
		if atms is False: 
			return
		atmDic,meanDic = {},{}
		for atm in atms:
			atmType = atm.atomtype
			dens = atm.densMetric[densMet][normType]['values'][dataset]
			if atmType not in atmDic.keys():
				atmDic[atmType] = [dens]
			else:
				atmDic[atmType] += [dens]
		for k in atmDic.keys():
			mean = np.mean(atmDic[k])
			print '{} --> {}'.format(k,mean)
			meanDic[k] = mean
		return meanDic


















	################################################################################
	# SOME METHODS WHICH DEPEND ON TOM'S BDAMAGE PYTHON SCRIPT

	def parseBDamageFile(self,
						 BDamageFile = '',
						 atomtype    = '',
						 restype     = ''):

		# parse the output file to a separately generated BDamage calculation job
		# (see https://github.com/td93/B_Damage/tree/master/Python for code to do this)

		fileIn = open(BDamageFile,'r')
		BdamDic = {}
		for l in fileIn.readlines():
			if 'ATOM' in l and len(l.split())>10:
				if atomtype in l or atomtype == '':
					if restype in l or restype == '':
						atomID = '-'.join([l.split()[2:6][i] for i in [2,3,1,0]])
						Bdam = l.split()[-1]
						print '{} --> {}'.format(atomID,Bdam)
						BdamDic[atomID] = Bdam
		return BdamDic

	def getBdamageStats(self,BDamageFile = ''):

		# retrieve BDamage values for each atom type and
		#  calculate mean & std dev

		csvFile = open('{}bDamageStats.csv'.format(self.outputDir),'w')
		csvFile.write('id,mean,std dev\n')
		atmTypeDic = self.groupByAtmType()
		for aType in atmTypeDic.keys():
			aInfo 	= aType.split('-')
			BdamDic = self.parseBDamageFile(BDamageFile = BDamageFile,
											atomtype    = aInfo[1],
											restype     = aInfo[0])
			mean 	= np.mean(map(float,BdamDic.values()))
			std 	= np.std(map(float,BdamDic.values()))
			csvFile.write('{},{},{}\n'.format(aType,mean,std))
		csvFile.close()

	def compareDensityToBdamage(self,
								BDamageFile = '',
								resType     = '',
								atomType    = '',
								densMet     = 'loss',
								normType    = 'Standard',
								fileType    = '.svg'):

		# compare calculated density metric values to 
		# calculated Bdamage values

		# get required atoms in structure
		atms = self.getAtom(restype  = resType,
							atomtype = atomType)

		self.calcAdditionalMetrics(metric    = densMet,
								   normType  = normType,
								   newMetric = 'average')

		# retrieve Bdamage values from Bdamage run log file
		BdamDic = self.parseBDamageFile(BDamageFile = BDamageFile,
										atomType    = atomType,
										resType     = resType)

		# get density metric values and Bdamage values
		plotData = {'x':[],
					'y':[]}

		for atm in atms:
			Bdam = BdamDic[atm.getAtomID()]
			dens = atm.densMetric[densMet][normType]['average']
			plotData['x'].append(Bdam)
			plotData['y'].append(dens)

		# plot the relationship between Bdamage and specified density metric
		self.plotScatterPlot(xData     = plotData['x'],
						     yData     = plotData['y'],
						     xLabel    = 'Bdamage',
							 ylabel    = '{} D{}'.format(normType,densMet),
							 figtitle  = '{}-{}'.format(resType,atomType),
							 fileType  = fileType,
							 saveTitle = 'D{}Scatter_{}-{}_BdamageVsDensMetric'.format(densMet,resType,atomType))























