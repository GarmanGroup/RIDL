import matplotlib.pyplot as plt
from PDBFileManipulation import writePDBline_DamSite
from matplotlib.gridspec import GridSpec
from findMetricChange import findBchange
from combinedAtom import combinedAtom
from CalphaWeight import CalphaWeight
from pandas import DataFrame
import string
from scipy.stats import ttest_ind,skew,kurtosis,mstats,ks_2samp,linregress,f,anderson_ksamp
import numpy as np
import operator
import os
import warnings
warnings.filterwarnings('ignore')
from bioInfo import bioInfo
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
				 seriesName        = 'untitled',
				 inclFCmetrics     = False):

		self.datasetList 		= datasetList
		self.numLigRegDatasets 	= numLigRegDsets    # number of datasets to perform linear regression over
		self.doseList 			= doseList 			# list of (DWD)doses for each dataset
		self.initialPDBList		= initialPDBList 	# list of atom objects for initial (lowest dose) pdb file
		self.outputDir			= outputDir
		self.partialDatasets	= partialDatasets   # Bool, whether atoms only in subset of datasets are included
		self.seriesName			= seriesName

		# only include density-weighted metrics if they exist. 
		# This is only in case where FC maps have been created 
		# within the RIDL map generation stage of th pipeline
		self.inclFCderivedMetrics = inclFCmetrics
		
		# if number of datasets to perform lin reg not stated then 
		# set to total number of datasets in series
		if self.numLigRegDatasets == 0:
			self.numLigRegDatasets = len(self.datasetList[0].mindensity)

		bInfo = bioInfo()
		self.aminoAcids = bInfo.getAminoAcids()

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

		# first check that each PDBarray contains the 
		# same number of atoms (consistency check)
		if len(self.datasetList) > 1:
			if printText:
				ln = 'Multiple datasets detected...'
				self.printOrWriteToLog(logFile = logFile,
									   txt     = ln)
			for dataset in self.datasetList:
				if len(dataset) != len(self.datasetList[0]):
					if printText:
						txt = 'Not all PDB structures have same number of atoms!\n'+\
							  'Will only include atoms common to ALL structures...'
						self.printOrWriteToLog(logFile = logFile, txt = txt)

					break
		elif len(self.datasetList) == 1:
			if printText:
				ln = 'Single dataset detected...'
				self.printOrWriteToLog(logFile = logFile, txt=ln)

		PDBdoses = []
		numNotFoundAtoms = 0

		if printText:
			ln = 'Locating common atoms to ALL datasets...:'
			self.printOrWriteToLog(logFile = logFile, txt = ln)

		singDimAttrs = ['atomnum',
						'residuenum',
						'atomtype',
						'basetype',
					    'chaintype',
					    'X_coord',
					    'Y_coord',
					    'Z_coord']

		multiDimAttrs = ['Bfactor',
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
						 'meanNegOnly',
						 'meanPosOnly']

		if self.inclFCderivedMetrics:
			extraMetrics = ['fracOfMaxAtomDensAtMin',
							'densityWeightedMean',
							'densityWeightedMin',
							'densityWeightedMax',
							'densityWeightedMeanNegOnly',
							'densityWeightedMeanPosOnly']
			for w in extraMetrics:
				multiDimAttrs.append(w)

		for atom in self.datasetList[0]:
			atm_counter = 1
			atomDict = { attr : getattr(atom, attr) for attr in singDimAttrs }.copy()
			atomDict.update( { attr : [getattr(atom, attr)] for attr in multiDimAttrs } )

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
				if not found: 
					# if atom not in dataset, add dummy value 
					indexindataset.append(-1)
					for attr in multiDimAttrs: 
						atomDict[attr].append(np.nan)

			# remove this located atom from the later dataset now 
			# that it has been located --> to make loop faster 
			for j in range(1,len(indexindataset)+1):
				if indexindataset[j-1] != -1:
					self.datasetList[j].pop( indexindataset[j-1] )

			if atm_counter != len(self.datasetList) and self.partialDatasets is False:
				if printText:
					txt = 'Atom "{}" not found in all datasets\n'.format(atomIndentifier)+\
					      '---> not including atom in atom list...'
					self.printOrWriteToLog(logFile = logFile, txt = txt)
				numNotFoundAtoms += 1
				continue

			else:
				newatom = combinedAtom()
				for attr in singDimAttrs+multiDimAttrs:
					if attr[0] != '_':
						if not isinstance(atomDict[attr],list):
							setattr(newatom,attr,atomDict[attr])
						else:
							metName = self.findMetricName(attr)
							newatom.getDensMetricInfo(metric   = metName,
													  normType = 'Standard',
													  values   = atomDict[attr])

				if atm_counter != len(self.datasetList):
					if printText:
						txt = 'Atom "{}" not found in all datasets\n'.format(atomIndentifier)+\
							  '---> including partial information in atom list...'
						self.printOrWriteToLog(logFile = logFile, txt = txt)
					numNotFoundAtoms += 1
				PDBdoses.append(newatom)

		if self.partialDatasets:
			if printText:
				ln = '# atoms not in all datasets: {}'.format(numNotFoundAtoms)
				self.printOrWriteToLog(logFile = logFile,
									   txt     = ln)
		else:
			if printText:
				ln = '# atoms removed since not in all datasets: {}'.format(numNotFoundAtoms)
				self.printOrWriteToLog(logFile = logFile, 
									   txt     = ln)

		if printText:
			ln = '---> Finished!'
			self.printOrWriteToLog(logFile = logFile, 
								   txt     = ln)

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

		conversions = [['Bfactor','Bfactor'],
					   ['Occupancy','occupancy'],
					   ['meandensity','mean'],
					   ['maxdensity','gain'],
					   ['mindensity','loss'],
					   ['mediandensity','median'],
					   ['numvoxels','num_voxels'],
					   ['stddensity','standard_deviation'],
	    			   ['min90tile','90tile_loss'],
	    			   ['max90tile','90tile_gain'],
	    			   ['min95tile','95tile_loss'],
	    			   ['max95tile','95tile_gain'],
	    			   ['meanNegOnly','mean_negOnly'],
	    			   ['meanPosOnly','mean_posOnly']]

		if self.inclFCderivedMetrics:

			conversions += [['fracOfMaxAtomDensAtMin','atomic_density_fraction_at_loss'],
							['densityWeightedMean','density_weighted_mean'],
							['densityWeightedMeanNegOnly','density_weighted_mean_negOnly'],
							['densityWeightedMeanPosOnly','density_weighted_mean_posOnly'],
							['densityWeightedMin','density_weighted_loss'],
							['densityWeightedMax','density_weighted_gain']]

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
			if printOutput:
				print 'Found {} atom(s) matching description'.format(len(foundAtoms))
			return foundAtoms
		if printOutput:
			print 'Atom matching description not found'
		return []

	def getDensMetrics(self):

		# get a list of all density metrics and normalisations that 
		# currently have been defined with a set of values

		currentMetrics = []
		for metric in self.atomList[0].densMetric.keys():
			for normType in self.atomList[0].densMetric[metric].keys():
				currentMetrics.append([metric,normType])
		return currentMetrics

	def getFormattedmetricName(self,
		                       metric='loss', norm='Standard', form='HTML'):

		# return a formatted version of the metric name

		if form == 'HTML':
			if norm == 'Standard':
				if metric == 'loss':
					name = r"\(D_\text{loss}\text{(atom)}\)"
				elif metric == 'density_weighted_loss':
					name = r"\(D_\text{loss}^{\rho}\text{(atom)}\)"
				elif metric == 'density_weighted_mean_negOnly':
					name = r"\(D_\text{mean}^{-,\,\rho}\text{(atom)}\)"
				else:
					name = metric
			elif norm == 'Calpha normalised':
				if metric == 'loss':
					name = r"\(C_\alpha\text{-normalised}\ D_\text{loss}\text{(atom)}\)"
				elif metric == 'density_weighted_loss':
					name = r"\(C_\alpha\text{-normalised}\ D_\text{loss}^{\rho}\text{(atom)}\)"
				elif metric == 'density_weighted_mean_negOnly':
					name = r"\(C_\alpha\text{-normalised}\ D_\text{mean}^{-,\,\rho}\text{(atom)}\)"
				else:
					name = 'Calpha-normalised ' + metric
			else:
				name = norm + ' ' + metric

		elif form == 'TEX':
			if norm == 'Standard':
				if metric == 'loss':
					name = r"$D_\mathrm{loss}\mathrm{(atom)}$"
				elif metric == 'density_weighted_loss':
					name = r"$D_\mathrm{loss}^{\rho}\mathrm{(atom)}$"
				elif metric == 'density_weighted_mean_negOnly':
					name = r"$D_\mathrm{mean}^{-,\,\rho}\mathrm{(atom)}$"
				else:
					name = metric
			elif norm == 'Calpha normalised':
				if metric == 'loss':
					name = r"$C_\alpha\mathrm{-normalised}\ D_\mathrm{loss}\mathrm{(atom)}$"
				elif metric == 'density_weighted_loss':
					name = r"$C_\alpha\mathrm{-normalised}\ D_\mathrm{loss}^{\rho}\mathrm{(atom)}$"
				elif metric == 'density_weighted_mean_negOnly':
					name = r"$C_\alpha\mathrm{-normalised}\ D_\mathrm{mean}^{-,\,\rho}\mathrm{(atom)}$"
				else:
					name = 'Calpha-normalised ' + metric
			else:
				name = norm + ' ' + metric

		return name

	def calcBfactorChange(self):

		# calculate Bfactor change between initial 
		# and each later dataset

		BmetDic = findBchange(self.initialPDBList,self.atomList,'Bfactor')
		for atom in self.atomList:
			atom.getDensMetricInfo(metric = 'BfactorChange',
					   			   values = BmetDic[atom.getAtomID()])

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
			if printText:
				print 'Calculating Calpha weights at each dataset '+\
					  'for metric: {}, normalisation: {}'.format(metric,normType)
			CAweights = self.retrieveCalphaWeight(metric = metric)

		# loop over all atoms in list and calculate additional 
		# metrics for each atom in atomList
		for atom in self.atomList:
			if newMetric == 'Calpha normalised':
				atom.CalphaWeightedDensChange(CalphaWeights = CAweights,
											  metric        = metric)
			elif newMetric == 'lin reg':
				atom.calcLinReg(numLinRegDsets = self.numLigRegDatasets,
								normType       = 'Standard',
								metric         = metric)
			elif newMetric == 'net':
				atom.calcNetChangeMetric(normType = 'Standard')
			elif newMetric == 'dataset 1 subtracted':
				atom.calcFirstDatasetSubtractedMetric(normType = 'Standard',
													  metric   = metric)
			elif newMetric == 'average':
				atom.calcAvMetric(normType,metric)

	def retrieveCalphaWeight(self,
							 metric = 'loss'):

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
			atom.calcVectorWeightedMetric(metric   = metric,
										  normType = normType,
										  vector   = vector)

	def calcVectorSubtractedMetric(self,
								   metric   = 'loss',
								   normType = 'Standard',
								   vector   = []):

		# for a specified metric calculate new vector-weighted 
		# values, where the metric over a series of doses
		# is subtracted from a vector of per-dataset scalings

		for atom in self.atomList:
			atom.calcVectorSubtractedMetric(metric   = metric,
											normType = normType,
											vector   = vector)	

	def calcStandardisedMetrics(self,
								metric = 'loss'):

		# standardise distribution of a given metric to have 
		# (mean=0,std=1) for each dataset within a damage series

		data  = [atm.densMetric[metric]['Standard']['values'] for atm in self.atomList]
		meand = np.nanmean(data,0)
		stdd  = np.nanstd(data,0)
		for atm in self.atomList:
			values = atm.densMetric[metric]['Standard']['values']
			atm.getDensMetricInfo(metric   = metric,
								  normType = 'Standardised',
								  values   = (values-meand)/stdd)

	def writeMetric2File(self,
						 where    = './',
						 groupBy  = 'none',
						 metric   = 'loss',
						 normType = 'Standard',
						 sortby   = 'atomnum',
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

			for i,k in enumerate(statDic.keys()):
				roundedVals = [round(v,numDP) for v in statDic[k]]
				csvfile.write('{},{},{}\n'.format(i,k,','.join(map(str,roundedVals))))			
		else:

			# sort atom list by atom number or damage metric
			if sortby == 'atomnum':
				self.atomList.sort(key=lambda x: x.atomnum)
			else:
				self.atomList.sort(key=lambda x: x.densMetric[metric][normType]['values'][0], reverse=True)
			for atom in self.atomList:
				csvfile.write('{},{},'.format(atom.atomnum,atom.getAtomID()))
				roundedVals = [round(v,numDP) for v in atom.densMetric[metric][normType]['values']]
				csvfile.write(','.join(map(str,roundedVals)))
				csvfile.write('\n')
		csvfile.close()

		return csvName

	def reorderByAtomNumber(self):

		# reorder the list of atoms by atom number.
		# this is to put the atoms back into a
		# known consistent order

		self.atomList.sort(key = lambda x: x.atomnum)

	def findProbAboveAvDam(self,
						   metric    = 'loss',
						   normType  = 'Calpha normalised',
						   threshold = 0):

		# a function that only will work for the calpha
		# normalised metric. Assuming zero is the base
		# Dloss level, determine overall probability 
		# that a randomly picked atom will have above
		# average Dloss. Could be used as an overall 
		# indicator for damage to the structure

		for d in self.getDsetList():
			count = 0
			for atm in self.atomList:
				if d in atm.getPresentDatasets():
					if atm.densMetric[metric][normType]['values'][d] > threshold:
						count += 1

			prob = float(count)/self.getNumAtoms()

			print 'Dataset {}: damage probability: {}'.format(d,round(prob,2))

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

		atms1 = self.getAtom(restype  = resiType,
							 atomtype = atomType1)
		if atms1 == []:
			return

		for atm in atms1:
			dic[atm.getAtomID()] = atm.densMetric[metric][normType]['values']

		for atom in self.atomList:
			if atom.basetype == resiType and atom.atomtype == atomType2:
				key   = '-'.join(atom.getAtomID().split('-')[0:3]+[atomType1])
				# check that generated 'key' exists and skip if not (needed for incomplete residues)
				if key not in dic.keys(): 
					continue
				a,b = np.array(dic[key]),atom.densMetric[metric][normType]['values']
				ratioDic['ratio'].append(a/b)
				ratioDic['distance'].append(a-b)
		summary = {}
		for key in ratioDic:
			print np.mean(ratioDic[key],0)
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

		fileout = open('{}{}_{}-D{}-{}.txt'.format(self.outputDir,title,normType,metric,rType),'w')

		foundPairs = {}
		for pair in pairs:
			found = self.getAtom(restype = pair[0])

			if not found: 
				continue # don't proceed if no atoms of a residue type

			output = self.findMetricRatio(metric    = metric,
										  normType  = normType,
										  resiType  = pair[0],
										  atomType1 = pair[1],
										  atomType2 = pair[2])
			if output == None:
				continue

			foundPairs['-'.join(pair)] = output[0]
			fileout.write(output[1][rType]+'\n')
		fileout.close()

		# plot comparison for different atom pairings
		sns.set_palette("deep", desat=.8)
		sns.set_context(rc={"figure.figsize":(10, 10)})
		fig = plt.figure()
		ax  = plt.subplot(111)
		colormap = plt.cm.nipy_spectral
		plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(foundPairs.keys()))])
		for key in foundPairs.keys():
			xData = self.doseList
			yData = np.mean(foundPairs[key][rType],0)
			if not errorBars:
				plt.plot(xData,
					     yData,
					     label = key)
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

		figName = '{}metric-D{}_normalisation-{}_{}_{}'.format(self.outputDir,metric,normType,rType,title)
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
				atom.calcRatioToMeanMetric(metric   = metric,
										   normType = normType,
										   avMetric = av)
			elif diff == 'std-devs':
				atom.calcNumStdFromMeanMetric(metric    = metric,
											  normType  = normType,
											  avMetric  = av,
											  stdMetric = std)

	def numAtmsWithMetricAboveLevel(self,
									dataset   = 0,
									metric    = 'loss',
									normType  = 'Standard',
									threshold = 0,
									printStr  = False,
									firstTime = True,
									atomType  = '',
									resType   = ''):

		# determine the number of atoms with metric value above a 
		# specified threshold 'threshold' determined as number of 
		# std dev away from structure-wide average for metric 'printStr' 
		# (Boolian) is True if string printed to command line.
		# if 'dataset' is 'all', take average over datasets

		if firstTime:
			self.calcMetricDiffFromStructureMean(metric   = metric,
												 normType = normType,
												 diff     = 'std-devs')
		if atomType == '' and resType == '':
			aList = self.atomList
		else:
			aList = self.getAtom(restype  = resType,
							     atomtype = atomType)

		if dataset == 'all':
			dList = self.getDsetList()
		else:
			dList = [dataset]

		countList = []
		for d in dList:
			count = 0
			for atom in aList:
				if atom.densMetric[metric]['num stds']['values'][d] > threshold:
					count += 1
			if printStr:
				print 'For dataset {}:'.format(d)
				print '# atoms with {} D{} metric above {}'.format(normType,metric,threshold)+\
					  'std dev of structure-wide mean is {}'.format(count)
			countList.append(count)

		if dataset != 'all':
			return count
		else:
			meanCount = np.mean(countList)
			try:
				fracCount = round(float(meanCount)/len(aList),2)
			except ZeroDivisionError:
				fracCount = np.nan
			return meanCount,fracCount

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

		foundPairs = {}
		for pair in pairs:
			found = self.getAtom(restype = pair[0])
			if len(found) == 0: 
				continue # don't proceed if no atoms of a residue type
			output = self.findMetricRatio(metric   = metric,
										  normType  = normType,
										  resiType  = pair[0],
										  atomType1 = pair[1],
										  atomType2 = pair[2])
			foundPairs['-'.join(pair)] = output[0]

		# plot the relationship between xData and yData
		sns.set_palette("deep", desat=.8)
		sns.set_context(rc={"figure.figsize":(10, 10)})
		fig = plt.figure()
		ax  = plt.subplot(111)

		total = 0
		for key in foundPairs.keys(): 
			total += len(foundPairs[key][rType])
		xDataTotal = np.linspace(0,1,total)

		i,j = 0,-1
		cmap = plt.cm.Set1
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

		if not os.path.isfile(strippedName):
			if strippedName.endswith(fileType):
				return strippedName
			else:
				return strippedName+fileType # this may give wrong name..

		fname = lambda x: strippedName.replace(fileType,'_{}{}'.format(x,fileType))
		i = 1
		while os.path.isfile(fname(i)): 
			i += 1 
		return fname(i)

	def getAtomListWithoutPartialAtoms(self,
									   dataset   = 0,
									   atmList   = []):

		# return new list of atoms with partial atoms 
		# removed for that dataset

		if atmList == []:
			atmList = self.atomList

		newList = []
		for atom in atmList:
			if dataset in atom.getPresentDatasets():
				newList.append(atom)
		return newList

	def findTopAtomOfType(self,
						  metric    = 'loss',
						  normType  = 'Standard',
						  dataset   = 'all',
						  atomType  = 'OH',
						  resType   = 'TYR',
						  printText = True):

		# find the top atom of a particular type for a
		# specified dataset for the currently set 
		# damage metric

		if dataset != 'all':
			dList = [dataset]
		else:
			dList = self.getDsetList()

		for d in dList:
			if printText:
				print 'For dataset: {}'.format(d)
			orderedAtms = self.getTopNAtoms(metric   = metric,
											normType = normType,
											dataset  = d,
											n        = 'all')

			for i,atm in enumerate(orderedAtms):
				if atm.basetype == resType:
					if atm.atomtype == atomType:
						foundAtm = atm
						break
			try:
				foundAtm
				if printText:
					print 'Top {}-{} atom is "{}" at position {}'.format(resType,
																		 atomType,
																		 foundAtm.getAtomID(),
																		 i)

			except NameError:
				if printText:
					print 'No atom of type: "{}-{}"'.format(resType,atomType)

	def getTopNAtoms(self,
					 metric   = 'loss',
					 normType = 'Standard',
					 dataset  = 0,
					 n        = 25,
					 topOrBot = 'top'):

		# for a given metric type, determine top 
		# n damage sites. 'n' is integer or 'all'.

		# decide whether to take top or least damaged atoms
		if topOrBot == 'top':
			r = True
		else:
			r = False

		if dataset != 'all':
			atomList = self.getAtomListWithoutPartialAtoms(dataset = dataset) # remove atoms not in present dataset
			atomList.sort(key=lambda x: x.densMetric[metric][normType]['values'][dataset], reverse = r)

			if n != 'all': 
				return atomList[:int(n)]
			else: 
				return atomList

		else:
			# rank with respect to all datasets
			ranks = {}
			for d in self.getDsetList():

				atomList = self.getAtomListWithoutPartialAtoms(dataset = d)
				atomList.sort(key = lambda x: x.densMetric[metric][normType]['values'][d], reverse = False)

				for i,atom in enumerate(atomList):
					if d == 0:
						ranks[atom.getAtomID()] = [i+1]
					else:
						ranks[atom.getAtomID()].append(i+1)

			overallRank = {}
			for k in ranks.keys():
				overallRank[k] = np.product(map(float,ranks[k]))

			top = sorted(overallRank.iteritems(), key = operator.itemgetter(1), reverse = r)[:n]
			atomList = [t[0] for t in top]
			return atomList

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
			headerInfo += 'REMARK All datasets within damage series '+\
						  'are present (increasing chain index).\n'
		else:
			headerInfo += 'REMARK File contains information for dataset {}.\n'.format(dataset)
		pdbOut.write(headerInfo)

		for l in pdbIn.readlines():
			if l.split()[0] in ('CRYST1','SCALE1','SCALE2','SCALE3'):
				pdbOut.write(l)

		if dataset == 'all':
			datasets = self.getDsetList()
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
						   metric          = 'loss',
						   normType        = 'Standard',
						   incOtherMetrics = True,
						   dataset  	   = 0,
						   useInHTML       = True,
						   splitBy         = ',',
						   n        	   = 25):

		# return the top n atoms within structure in terms of 
		# metric 'metric'. 'n' takes values 'all' or an integer

		# if normalisation type not recognised, set to default
		if normType not in ('Standard','Calpha normalised'):
			normType = 'Standard'

		topN = self.getTopNAtoms(metric   = metric,
								 normType = normType,
								 dataset  = dataset,
								 n        = n)
		atomInfoList = []
		for atom in topN:
			data = atom.getAtomID()+splitBy
			for norm in ('Standard','Calpha normalised'):
				data += str(round(atom.densMetric['loss'][norm]['values'][dataset],2))+splitBy


			if metric == 'loss':
				if self.inclFCderivedMetrics:
					# this gives an indication of how far away the Dloss-valued voxel is from the atomic centre
					data += str(round(atom.densMetric['loss-reliability']['Standard']['values'][dataset],2))+splitBy

			if incOtherMetrics:
				extra = []

				if metric == 'loss':
					metricsToInclude = ['mean',
										'gain',
										'Bfactor']

				if metric == 'density_weighted_mean_negOnly':
					metricsToInclude = ['density_weighted_loss',
										'loss',
										'Bfactor']

				else:
					metricsToInclude = ['Bfactor']

				for met in metricsToInclude:
					extra.append(str(round(atom.densMetric[met]['Standard']['values'][dataset],2)))
				data += splitBy.join(extra)
			atomInfoList.append(data)

		if useInHTML:
			calphaStr = 'C<sub>&#945</sub>'
		else:
			calphaStr = 'Calpha'

		if metric == 'loss':
			stringOut = 'Atom-ID{}Dloss{}Dloss ({}-normalised){}'.format(splitBy,splitBy,calphaStr,splitBy)
			if self.inclFCderivedMetrics:
				stringOut += 'Peak Proximity{}'.format(splitBy)
			stringOut += 'Dmean{}Dgain{}Bfactor\n'.format(splitBy,splitBy)+\
					 	 '\n'.join(atomInfoList)

		elif metric == 'density_weighted_mean_negOnly':
			stringOut = 'Atom-ID{}DW-mean (-ve only){}DW-mean (-ve only) ({}-normalised){}'.format(splitBy,splitBy,calphaStr,splitBy)	
			stringOut += 'DW-loss{}Dloss{}Bfactor\n'.format(splitBy,splitBy)+\
					 	 '\n'.join(atomInfoList)

		else:
			stringOut = 'Atom-ID{}{}{}Bfactor'.format(splitBy,metric,splitBy)
		return stringOut

	def getTopNAtomsDotPlot(self,
						    metrics   = ['loss', 'loss'],
						    normTypes = ['Standard', 'Calpha normalised'],
						    dataset   = 0,
						    numHits   = 25,
						    savePlot  = True,
						    fileType  = '.png'):

		# make a dot plot for top n atoms in structure in terms of
		# metric 'metric'. 'numHits' takes values 'all' or an integer
		# Top n atoms are ranked based of the primary metric as
		# specified by the lists 'metrics' and 'normTypes'

		if len(metrics) != len(normTypes):
			print 'getTopNAtomsDotPlot() requires "metrics" and ' +\
				  '"normTypes" inputs to be lists of same length'
			return

		primaryMetric = self.getFormattedmetricName(metrics[0], normTypes[0], form='TEX')

		topN = self.getTopNAtoms(metric   = metrics[0],
								 normType = normTypes[0],
								 dataset  = dataset,
								 n        = numHits)

		metricsToPlot = [self.getFormattedmetricName(m, n, form='TEX') for m, n in zip(metrics, normTypes)]
		plotData = {m: [] for m in metricsToPlot}
		plotData['Atom ID'] = []

		for atom in topN:
			plotData['Atom ID'].append(atom.getAtomID())
			for m, n in zip(metrics, normTypes):
				metName = self.getFormattedmetricName(m, n, form='TEX')
				val = atom.densMetric[m][n]['values'][dataset]
				plotData[metName].append(val)

		plotdf = DataFrame(plotData)

		# Make the PairGrid		
		g = sns.PairGrid(plotdf.sort_values(primaryMetric,
						 ascending = False),
		                 x_vars    = metricsToPlot,
		                 y_vars    = ['Atom ID'],
		                 size      = 14,
		                 aspect    = .25)

		# Draw a dot plot using the stripplot function
		g.map(sns.stripplot,
			  size      = 20,
			  orient    = "h",
			  palette   = sns.color_palette("Blues_r",numHits),
		      edgecolor = "gray")

		g.set(ylabel = "")

		titles = metricsToPlot # this might need changing to give nicer labels
		for ax, title in zip(g.axes.flat, titles):
			ax.xaxis.grid(False)
			ax.yaxis.grid(True)

		sns.despine(left   = True,
					bottom = True)

		if savePlot:
			plotDir = 'topNatoms/'
			if plotDir.replace('/','') not in os.listdir(self.outputDir+'plots/'):
				os.mkdir(self.outputDir+'plots/'+plotDir)

			saveName = '{}plots/{}dotPlot_top-{}-DamSites_dset-{}{}'.format(self.outputDir,plotDir,numHits,dataset,fileType)
			saveName = self.checkUniqueFileName(fileName = saveName,
									 			fileType = fileType)
			g.savefig(saveName)

		return saveName

	def plotAtomtypeRankingWithDataset(self,
		 							   metric    = 'loss',
						   			   normType  = 'Standard',
						   			   residue   = ['TYR','TYR','TYR','TYR','GLU','GLU','GLU','GLU'],
						   			   atomtype  = ['OH','CZ','CE1','CE2','CG','CD','OE1','OE2'],
						   			   figTitle  = '',
						   			   axesFont  = 18,
						   			   legFont   = 12,
						   			   titleFont = 20,
						   			   saveFig   = True,
						   			   percent   = True,
						   			   outputDir = '',
						   			   fileType  = '.svg'):

		# plot atom type ranking (see method below)
		# for series of datasets in damage series
		# for the selected atom type of interest

		rankList = {}
		for res,atm in zip(residue,atomtype):
			aTag = res + '-' + atm
			rankList[aTag] = []
			for d in self.getDsetList():
				rank = self.getAtomtypeRanking(metric   = metric,
											   normType = normType,
											   dataset  = d,
											   residue  = [res],
											   atomtype = [atm])
				if percent:
					rankList[aTag].append(rank[aTag]['Ratio'])
				else:
					rankList[aTag].append(rank[aTag]['Rank'])

		args = [metric,
				normType,
				residue,
				atomtype]

		if figTitle == '':
			title = 'metric: D{}, normalisation:{}\n{}-{} ranking'.format(*args)
		else:
			title = figTitle

		ylabel = 'Ranking'
		if percent:
			ylabel += ' (%)'

		name = 'Lineplot_Metric-D{}_Normalisation-{}-{}_ranking'.format(*args)
		fileName = self.makeLinePlot(vals2Plot   = rankList,
								 	 xVals       = 'doses',
								 	 yLabel      = ylabel,
								 	 axesFont    = axesFont,
								 	 legFont     = legFont,
								 	 titleFont   = titleFont,
								 	 title       = title,
								 	 fileType    = fileType,
					 				 saveFig     = saveFig,
					 				 fileName    = name,
					 				 outputDir   = outputDir)

		if saveFig:
			return fileName

	def getAtomtypeRanking(self,
						   metric    = 'loss',
						   normType  = 'Standard',
						   dataset   = 0,
						   residue   = ['TYR'],
						   atomtype  = ['OH'],
						   printText = False):

		# for specified dataset and metric type, determine
		# the rank of a particular atomtype in terms of 
		# the currently specified metric

		if printText:
			info = '---------------------------------'+\
				   'Determining ranking of {}-{} atoms'.format(residue,atomtype)+\
				   'Dataset: {}'.format(dataset)+\
				   'metric: {}; normalisation: {}'.format(metric,normType)

		if dataset != 'all':
			dList = [dataset]
		else:
			dList = self.getDsetList()

		atomRanks = {}
		for d in dList:

			a,b = self.getPerAtmtypeStats(metric   = metric,
								  		  normType = normType,
								  		  dataset  = d,
								  		  n        = 1)
			
			atomVals = {}
			for key in b.keys():
				atomVals[key] = [b[key]['mean']]

			sortedVals = sorted(atomVals.items(), key=operator.itemgetter(1))
			for i,k in enumerate(sortedVals):
				try:
					atomRanks[k[0]].append(float(i))
				except KeyError:
					atomRanks[k[0]] = [float(i)]

		atomRanksOverall = {}
		for key in atomRanks.keys():
			atomRanksOverall[key] = np.product(atomRanks[key])

		sortedRanks = sorted(atomRanksOverall.items(), key=operator.itemgetter(1))
		sortedRanks2 = [s[0] for s in sortedRanks]
		sortedRanks2.reverse()

		ranks = {}
		for res,atm in zip(residue,atomtype):
			keyAtm     = res + '-' + atm
			try:
				keyAtmPos  =  sortedRanks2.index(keyAtm)
			except ValueError:
				keyAtmPos  = np.nan
			keyAtmRank = round(float(keyAtmPos)/len(sortedRanks2),2)
			ranks[keyAtm] = {'Rank'  : keyAtmPos,
							 'Ratio' : keyAtmRank}

		return ranks

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
		fracOfTotal   = round( (float(n)/self.getNumAtoms())*100, 2 ) # n as fraction of total atoms
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
								  n          = 25,
								  palette    = 'hls',
								  barWidth   = 0.5,
								  axesFont   = 18,
								  titleFont  = 18,
								  legendFont = 14,
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

		sns.set_style("ticks")
		sns.set_context(rc = {"figure.figsize": (10, 6)})
		fig = plt.figure()
		ax  = plt.subplot(111)

		accumulatedVals = np.array([n]*numDsets)

		numResis = len(plotDic.keys())
		for i,(k, color) in enumerate( zip( plotDic.keys(), sns.color_palette(palette, n_colors = numResis,desat = 0.9) ) ):

			fd = DataFrame()
			fd['Dataset']   = np.array(range(numDsets))+1
			fd['Frequency'] = accumulatedVals
			sns.barplot(x     = 'Dataset',
						y     = 'Frequency',
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

		plt.xlabel('Dataset Number',
			   	   fontsize = axesFont)
		plt.ylabel('Frequency', 
			       fontsize = axesFont)

		if plotTitle == '':
			t = 'Top damage sites by residue/nucleotide type\n'
		else:
			t = plotTitle
		plt.title(t,
			 	  fontsize = titleFont)

		if not saveFig: 
			plt.show()
		else:

			args = [outputDir,
					n,
					metric,
					normType.replace(" ","")]

			saveName = '{}Top{}DamSiteByRestype_Metric-D{}_Normalisation-{}'.format(*args)

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

	def getStructureStats(self,
						  metric   = 'loss',
						  normType = 'Standard',
						  dataset  = 0):

		# for a given metric type, determine per-chain 
		# statistics on the distribution of damage

		statsDic 	= self.getStats(metric   = metric,
									normType = normType,
									dataset  = dataset,
									dic      = {'Structure':self.atomList})

		statsString	= self.reportStats(stats    = statsDic,
									   name     = 'Structure',
									   normType = normType,
									   n        = 1)

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
		m = metricList
		statsDic['#atoms']	 	= len(m)
		statsDic['mean'] 		= np.mean(m)
		statsDic['std'] 		= np.std(m)
		statsDic['skew'] 		= self.calculateSkew(metricList = m)
		statsDic['kurtosis']    = self.calculateKurtosis(metricList = m)
		statsDic['ratio'] 		= self.calcNetRatio(metricList = m)
		statsDic['outliers'] 	= self.calcNumOutliers(metricList = m,
													   metric     = metric,
													   normType   = normType)
		statsDic['normality']   = self.testForNormality(metricList = m)
		statsDic['returnOrder'] = ['mean',
								   'std',
								   '#atoms',
								   'outliers',
								   'skew',
								   'kurtosis']
		return statsDic

	def reportStats(self,
					stats    = {},
					name     = 'Residue',
					sortby   = 'mean',
					n        = 20,
					numDp    = 2,
					normType = 'Standard',
					format   = 'txt'):

		# format stats defined in 'stats' dictionary into an 
		# output string. 'sortby' specifies how the output should 
		# be ranked. 'n' is number of entries to print (specify 
		# an integer or 'all') format takes values in ('txt','html')

		sortbyDict   = self.getStatsForList(metricList = [1,1])
		headerOrder  = sortbyDict['returnOrder']
		headerString = '{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}\n'.format(name,*headerOrder)
		if sortby not in sortbyDict.keys(): 
			return 'Unexpected ranking specified'
		if n != 'all' and not isinstance(n,int): 
			return 'input "n" must be integer or "all", unable to complete'
		
		list1,list2 = [],[]
		for key in stats.keys():
			statsDic = stats[key]
			statsFmtd = []		
			for key2 in headerOrder: 
				if isinstance(statsDic[key2],float): statsFmtd.append('{}'.format(round(statsDic[key2],numDp)))
				else: statsFmtd.append(str(statsDic[key2]))
			string = '\t\t'.join([str(key)]+statsFmtd)
			list1.append(string)
			list2.append(statsDic[sortby])

		sortedList1 = [ x for (y,x) in sorted( zip(list2,list1) ) ] # sort by chosen 'sortby' parameter
		sortedList1.reverse()

		if n != 'all':
			stringOut =	'\n'.join(sortedList1[:n])
		else:
			stringOut =	'\n'.join(sortedList1)	
		return headerString + stringOut

	def getStatPerDataset(self,
					      metric          = 'loss',
					      normType        = 'Standard',
					      stat            = ['mean','skew'],
					      inclPartialAtms = True):

		# return a dictionary of stat per dataset 
		# for each residue in structure
		# 'stat' is either a single input (e.g. 'skew')
		# or a list of inputs (e.g. ['skew','mean'])

		if not isinstance(stat,list):
			stat = [stat]

		statPerDset = {s : {} for s in stat}
		for i in self.getDsetList():
			stats = self.getPerResidueStats(metric   = metric,
							   				normType = normType,
							   				dataset  = i,
							   				n        = 'all')
			if i == 0:
				for k in stats[1].keys():
					for s in stat:
						statPerDset[s][k] = [stats[1][k][s]]

			else:
				for k in stats[1].keys():
					for s in stat:
						statPerDset[s][k].append(stats[1][k][s])

		if not inclPartialAtms:
			# only include residue types that are present within 
			# all datasets in series. This is used to remove water
			# and ligands that only appear in subset of datasets
			for s in statPerDset.keys():
				for k in statPerDset[s].keys():
					if len(statPerDset[s][k]) != self.getNumDatasets():
						del statPerDset[s][k]

		return statPerDset

	def plotStatVsDataset(self,
		 				  metric    = 'loss',
						  normType  = 'Standard',
						  stat      = ['skew','kurtosis','mean','std'],
						  saveFig   = True,
						  axesFont  = 18,
						  titleFont = 16,
						  legFont   = 10,
						  figTitle  = '',
						  outputDir = './',
						  fileType  = '.svg',
						  palette   = 'hls',
						  printTxt  = False):

		# produce a plot of skew versus dataset number for
		# each residue type within the structure

		if not isinstance(stat,list):
			stat = [stat]

		statPerDset = self.getStatPerDataset(stat            = stat,
							  				 metric          = metric,
							   				 normType        = normType,
							   				 inclPartialAtms = False)

		fileNames = []
		for s in stat:

			if printTxt:
				print 'Plotting stat "{}" versus dataset'.format(s)

			args = [metric,
					normType,
					s]

			if figTitle == '':
				title = 'metric: D{}, normalisation:{}\n {} per residue '.format(*args)
			else:
				title = figTitle

			name = 'Lineplot_Metric-D{}_Normalisation-{}-{}-per-residue'.format(*args)

			fileName = self.makeLinePlot(palette     = palette,
									 	 vals2Plot   = statPerDset[s],
									 	 xVals       = 'doses',
									 	 yLabel      = s,
									 	 axesFont    = axesFont,
									 	 legFont     = legFont,
									 	 titleFont   = titleFont,
									 	 title       = title,
									 	 fileType    = fileType,
					 					 saveFig     = saveFig,
					 					 fileName    = name,
					 					 outputDir   = outputDir)
			fileNames.append(fileName)

		if saveFig:
			return fileNames

	def plotStatVsStat(self,
					   dataset   = 'all',
	 				   metric    = 'loss',
					   normType  = 'Standard',
					   stat1     = ['mean','mean', 'mean','skew'],
					   stat2     = ['skew','kurtosis','std','kurtosis'],
					   saveFig   = True,
					   axesFont  = 18,
					   titleFont = 20,
					   legFont   = 14,
					   figTitle  = '',
					   outputDir = './',
					   fileType  = '.svg',
					   palette   = 'husl',
					   markSize  = 40,
					   printTxt  = False):

		# produce a plot of stat1 versus stat2 for
		# each residue type within the structure
		# created for the specified dataset number

		if not isinstance(stat1,list):
			stat1 = [stat1]
		if not isinstance(stat2,list):
			stat2 = [stat2]

		statPerDset = {}
		stats = list(set(stat1 + stat2))
		statPerDset = self.getStatPerDataset(stat            = stats,
								  		     metric          = metric,
								   		     normType        = normType,
								   		     inclPartialAtms = False)

		fileNames = []
		for (s1,s2) in zip(stat1,stat2):

			if printTxt:
				print 'plotting stat "{}" vs "{}" for each dataset'.format(s1,s2)

			sns.set_style("whitegrid")
			sns.set_context(rc = {"figure.figsize": (14, 12)})
			num_cols = int( np.ceil( np.sqrt( self.getNumDatasets() ) ) )
			num_rows = int( np.floor( np.sqrt( self.getNumDatasets() ) ) )
			f, axes = plt.subplots(num_cols, 
							   	   num_rows+1, 
							       sharex = True, 
							       sharey = True)

			col = -1
			row = 0
			uniqPoints = []
			uniqLabels = []
			for d in self.getDsetList():
				col += 1
				if col >= num_cols:
					col  = 0
					row += 1 
				numResis = len(statPerDset[s1].keys())

				for (k,color) in zip( statPerDset[s1].keys(), sns.color_palette(palette,numResis) ):
					p = axes[col,row].scatter(statPerDset[s1][k][d],
										 	  statPerDset[s2][k][d],
										  	  label = k,
										  	  c     = color,
										  	  s     = markSize)
					if col == 0 and row == 0:
						uniqPoints.append(p)
						uniqLabels.append(k)
				sns.despine()

			f.text(0.5,
				   0.04,
				   s1,
				   ha       = 'center',
				   fontsize = axesFont)

			f.text(0.04,
				   0.5,
				   s2,
				   va       = 'center', 
				   rotation = 'vertical',
				   fontsize = axesFont)		

			# place legend outside to right of plot
			box = axes[col,row].get_position()
			axes[num_cols-1,num_rows/2].legend(loc            = 'top right',
					  			               bbox_to_anchor = (1.5, 1.3),
					  			               fontsize       = legFont)

			args = [metric,
					normType,
					s1,
					s2]

			if figTitle == '':
				title = 'metric: D{}, normalisation:{}\n {} vs {}'.format(*args)
			else:
				title = figTitle

			f.suptitle(title,
					   fontsize = titleFont)

			if not saveFig:
				plt.show()
			else:
				name = 'scatterplot_Metric-D{}_Normalisation-{}_{}vs{}'.format(*args)

				fileName = self.checkUniqueFileName(fileName = outputDir+name,
									 				fileType = fileType)
				f.savefig(fileName)
				fileNames.append(fileName)
		if saveFig:
			return fileNames

	def twoAtomTypeTtest(self,
						 metric    = 'loss',
						 normType  = 'Standard',
						 dataset   = 0,
						 atomType1 = ['OH'],
						 atomType2 = ['OE1','OE2','OD1','OD2'],
						 resType1  = ['TYR'],
						 resType2  = ['GLU','GLU','ASP','ASP'],
						 alpha     = 0.01):

		# perform Student t-test to compare metric for
		# two distinct sets of atoms within structure

		atmList1 = []
		for atm,res in zip(atomType1,resType1):
			atmList1 += self.getAtom(restype  = res,
									 atomtype = atm)
		atmList2 = []
		for atm,res in zip(atomType2,resType2):
			atmList2 += self.getAtom(restype  = res,
									 atomtype = atm)

		vals1 = [atm.densMetric[metric][normType]['values'][dataset] for atm in atmList1]
		vals2 = [atm.densMetric[metric][normType]['values'][dataset] for atm in atmList2]

		(tstat,pval) = ttest_ind(vals1,vals2, equal_var = False)

		if pval < alpha:
			reject = True
		else:
			reject = False
		return (tstat,pval,reject)

	def calcDiscreteDistMode(self,
							 metricList = []):

		# calculate mode of discrete histogram with 100 bins

		h = np.histogram(metricList,
						 bins = 100)
		mode_ind = np.argmax(h[0])
		mode_val = np.mean(h[1][mode_ind:mode_ind+2])
		return mode_val

	def calcNetRatio(self,
					 metricList = []):

		# calculate net ratio of distn values either side of distn mode

		switchVal = self.calcDiscreteDistMode(metricList = metricList)

		group = {'above' : [],
				 'below' : []}

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
		distMin  = np.min(metricList)
		sudoMax  = distMode + np.linalg.norm(distMin - distMode)
		count 	 = 0
		for val in metricList:
			if val > sudoMax: count += 1
		return count

	def calculateSkew(self,
					  metricList = []):

		# calculate skewness for a distribution of 
		# metric values for an input list of atoms

		skewVal = skew(metricList,
					   axis = 0,
					   bias = True)
		return skewVal

	def calculateKurtosis(self,
					      metricList = []):

		# calculate kurtosis for a distribution of 
		# metric values for an input list of atoms

		kurt = kurtosis(metricList,
				        axis = 0,
				        bias = True)
		return kurt

	def testForNormality(self,
						 metricList   = [],
						 suppressText = True):

		# test whether the metric values for a list of 
		# atoms differs from a normal distribution

		if len(metricList) < 8:
			if not suppressText:
				print 'Skew-test not valid when less '+\
					  'than 8 values provided - skipping'
			return 'n/a'
		elif len(metricList) < 20:
			if not suppressText:
				print 'kurtosis-test not valid when less '+\
					  'than 20 values provided - skipping'
			return 'n/a'
		try:
			(k2,pvalue) = mstats.normaltest(metricList,
											axis = 0)
		except np.ma.core.MaskError:
			return 'n/a'
		return pvalue

	def groupByAtmType(self,
					   dataset = 0):

		# group atoms in a dictionary by atom type

		atmtypeDict = {}
		for atom in self.atomList:
			if dataset not in atom.getPresentDatasets(): 
				continue
			atmtype = '-'.join([atom.basetype,atom.atomtype])
			if atmtype not in atmtypeDict.keys():
				atmtypeDict[atmtype] = [atom]
			else:
				atmtypeDict[atmtype].append(atom)
		return atmtypeDict

	def groupByResType(self,
					   dataset = 0):

		# group atoms in a dictionary by residue type

		resDict = {}
		for atom in self.atomList:
			if dataset not in atom.getPresentDatasets(): 
				continue
			if atom.basetype not in resDict.keys():
				resDict[atom.basetype] = [atom]
			else:
				resDict[atom.basetype].append(atom)
		return resDict

	def groupByChainType(self,
		 				 dataset = 0):

		# group atoms in a dictionary by chain type

		chainDict = {}
		for atom in self.atomList:
			if dataset not in atom.getPresentDatasets(): 
				continue
			if atom.chaintype not in chainDict.keys():
				chainDict[atom.chaintype] = [atom]
			else:
				chainDict[atom.chaintype].append(atom)
		return chainDict

	def checkCalphaAtomsExist(self,
							  printText = True):

		# check that Calpha backbone protein atoms 
		# actually exist within structure

		count = 0
		for res in self.aminoAcids:
			atom = self.getAtom(restype     = res,
								atomtype    = 'CA',
								printOutput = False)
			if len(atom) == 0:
				count += 1
		if count == len(self.aminoAcids):
			if printText:
				print 'Warning: no Calpha backbone atoms found in structure!'
				return False
		else:
			return True

	def densityMetricHeatMap(self,
							 metric         = 'loss',
							 normType       = 'Standard',
							 saveFig        = False,
							 fileType       = '.svg',
							 titleFont      = 24,
							 outputDir      = './',
							 forceSquares   = True,
							 columnPerChain = False):

		# plot a heatmap of density metric values per atom 
		# per dataset. Depending on the metric and normalisation
		# type chosen, a different colour scheme is chosen

		if columnPerChain:
			chains       = self.getChains()
			data = { c : { 'Atom'    : [],
					 	   'Dataset' : [],
						   'metric'  : [] } for c in chains }

			for atom in self.atomList:
				c = atom.chaintype
				data[c]['metric']  += atom.densMetric[metric][normType]['values']
				data[c]['Atom']    += [atom.getAtomID()]*self.getNumDatasets()
				data[c]['Dataset'] += range(1, self.getNumDatasets() + 1)
		else:
			chains      = self.getChains()
			chainCounts = {c:0 for c in chains}
			newChains   = []
			data        = {}
			for atom in self.atomList:

				c = atom.chaintype
				chainCounts[c] += 1
				i = chainCounts[c]
				j = i/200

				cj = c + str(j)
				if cj not in data.keys():

					data[cj] = {'Atom'    : [],
						  	    'Dataset' : [],
						  	    'metric'  : []}
					newChains.append(cj)

				data[cj]['metric']  += atom.densMetric[metric][normType]['values']
				data[cj]['Atom']    += [atom.getAtomID()]*self.getNumDatasets()
				data[cj]['Dataset'] += range(1, self.getNumDatasets() + 1)
			chains = newChains

		numChains = len(data.keys())

		sns.set_context("notebook")
		chainLengths = np.array([len( data[c]['Atom'] ) for c in chains])/self.getNumDatasets()
		maxLength    = int(2.0*max(chainLengths)) # get largest chain length

		# stretch the plot in the vertical direction, depending 
		# on how many atoms are present within the structure
		fontsize_pt = plt.rcParams['ytick.labelsize']
		dpi = 80
		matrix_height_pt = fontsize_pt * (maxLength+1)
		matrix_height_in = matrix_height_pt / dpi

		matrix_width_pt  = fontsize_pt * (20+self.getNumDatasets())*numChains
		matrix_width_in  = matrix_width_pt / dpi

		top_margin       = 0.01   
		bottom_margin    = 0.04
		figure_height    = matrix_height_in / (1 - top_margin - bottom_margin)
		figure_width     = matrix_width_in

		fig = plt.figure(figsize = (figure_width,figure_height))
		gs = GridSpec(maxLength+1,
					  numChains+1,
					  width_ratios = [0.1]+[1]*numChains)
		ax = []
		ax.append( plt.subplot( gs[:20,0] ) )
		for i in range(numChains):
			ax.append( plt.subplot( gs[1:chainLengths[i]+1,i+1] ) )

		if metric == 'loss' and normType == 'Standard':
			pal = 'Reds'
		else:
			pal = sns.diverging_palette(220, 
										10, 
										sep     = 80,
										n       = 7,
										as_cmap = True)

		# find overall max and min metric values
		vals    = [data[k]['metric'] for k in data.keys()]
		vals1d  = [x for sublist in vals for x in sublist]

		quantiles = mstats.mquantiles(vals1d, prob = [0.05, 0.95])
		highVal = quantiles[1]
		lowVal  = quantiles[0]

		for i,c in enumerate(chains):

			d = DataFrame(data = data[c])
			d = d.pivot("Atom",
						"Dataset",
						"metric")
			
			f = sns.heatmap(d,
							ax       = ax[i+1],
							cbar     = i == 0,
							cbar_ax  = None if i else ax[0],
							cmap     = pal,
							vmin     = lowVal, 
							vmax     = highVal,
							square   = forceSquares,
							linewidths = 1,
							linecolor = 'white',
							cbar_kws = {"orientation" : "vertical"})

		for i,a in enumerate(ax):
			if i == 0:
				continue
			plt.sca(a)
			plt.yticks(rotation = 0) 
			a.xaxis.tick_top()
			a.xaxis.set_label_position("top")

		figTitle = plt.suptitle('Metric: $D_{' + str(metric) + '}$'+'\nNormalisation: {}'.format(normType),
				  				fontsize = titleFont,
				  				y        = 1.02)

		fig.tight_layout()

		if not saveFig:
			plt.show()

		else:
			saveName = '{}heatmap_metric-{}_normalisation-{}'.format(outputDir,
																	 metric,
																	 normType)
			fileName = self.checkUniqueFileName(fileName = saveName,
									 			fileType = fileType)
			fig.savefig(fileName,
						bbox_extra_artists = (figTitle,),
						bbox_inches        = "tight")

	def getChains(self):

		# find different chains in list of atoms

		chains = self.getListOf(prop = 'chain')
		return chains

	def getResidues(self):

		# find different residues in list of atoms

		residues = self.getListOf(prop = 'residue')
		return residues

	def getListOf(self,
				  prop = 'chain'):

		# return list of different 'prop'
		#  values in list of atoms

		if prop == 'chain':
			attr = 'chaintype'
		elif prop == 'residue':
			attr = 'basetype'
		elif prop == 'nucleotide':
			attr = 'basetype'
		elif prop == 'atom':
			attr = 'atomtype'

		vals = []
		for atom in self.atomList:
			val = getattr(atom,attr)
			if val not in vals:
				vals.append(val)
		return vals

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
						if not suppressText:
							print 'Unusually low compared to mean '+\
								  'value ({} < {})'.format(round(val,3),
														   round(meanVal,3))
						highOrLow.append('low')
					else:
						if not suppressText:
							print 'Unusually high compared to mean '+\
								  'value ({} > {})'.format(round(val,3),
														   round(meanVal,3))
						highOrLow.append('high')
		if not suppressText:
			print '{} atoms found with suspiciously high/low damage '+\
				  'relative to average of that atom type'
		return suspAtoms,highOrLow

	def getNumDatasets(self):

		# get the number of datasets within the damage series

		return self.atomList[0].getNumDatasets()

	def getDsetList(self):

		# get a list of dataset numbers

		return range(self.getNumDatasets())

	def getSidechainAtoms(self):

		# get list of side chain atoms

		sideList = []
		for atm in self.atomList:
			if atm.side_or_main() != 'mainchain':
				sideList.append(atm)

		return sideList

	def graphMetricDistn(self,
						 metric     = 'loss',
						 normType   = 'Standard',
						 valType    = 'average',
						 plotType   = 'both',
						 resiType   = 'all',
						 resSplit   = False,
						 sideOnly   = False,
						 requireAll = False,
						 save       = True,
						 xLabel     = '',
						 fileType   = '.svg',
						 outputDir  = './',
						 printText  = True,
						 axesFont   = 18,
						 legFont    = 14,
						 titleFont  = 20,
						 plotTitle  = '',
						 inclTitle  = True,
						 calcKSstat = False,
						 calcADstat = False):

		# histogram/kde plot of density metric per atom.
		# plotType is 'histogram' or 'kde'.
		# resiType is 'all' or list of residue types.
		# 'valType' takes values 'average' (compute average),
		# 'all' (plot all datasets), or 'n' for int (not str)
		# n (plot dataset n).
		# If 'calcKSstat' is True, perform the 2-sample 
		# Kolmogorov-Smirnov test (will only calculate when 
		# two histograms plotted on 1 axis)
		# 'requireAll' set to True will only plot if ALL 
		# specified residue types are found
		# 'resSplit' means that separate distributions are
		# plotted for side and main chain (proteins) and 
		# base, phosphate, sugar (DNA/RNA)

		# attempt to find atoms of type 'resiType'
		# and flag if not all found
		if resiType != 'all':
			count = 0
			for res in resiType:
				atms = self.getAtom(restype = res)
				if len(atms) == 0: 
					count += 1
			if count > 0:
				if printText:
					print 'Warning: not all selected residue/nucleotide types found!'
				if requireAll:
					return {},''
			if count == len(resiType):
				if printText:
					print 'Warning: no residues/nucleotides found for current plot'
				return {},''

		if plotType not in ('hist','kde','both'): 
			print 'Unknown plotting type selected.. cannot plot..'
			return {},''

		if normType == 'Calpha normalised':
			self.calcAdditionalMetrics(metric = metric)
		if not self.checkMetricPresent(metric = metric,normType = normType): 
			print 'Specified metric to plot not present..'
			return {},''

		if sideOnly:
			atmList = self.getSidechainAtoms()
		else:
			atmList = self.atomList

		sns.set_style("ticks")
		sns.set_context(rc = {"figure.figsize": (10, 7)})
		fig = plt.figure()
		ax  = plt.subplot(111)

		datax      = {}
		dataxOrder = []
		stats      = {}

		if valType == 'all':
			numDsets = self.getNumDatasets()
			for (i, color) in zip( self.getDsetList(), sns.color_palette('rainbow', numDsets) ):
				presentAtms = self.getAtomListWithoutPartialAtoms(dataset = i,
																  atmList = atmList)

				if resiType == 'all':
					data = [atm.densMetric[metric][normType]['values'][i] for atm in presentAtms]
					lbl = 'Dataset {}'.format(i+1)
					datax[lbl] = {'values' : data,
								  'colors' : color}
					dataxOrder.append(lbl)

				else:
					for res in resiType:
						data = [atm.densMetric[metric][normType]['values'][i] for atm in presentAtms if atm.basetype == res]

						if len(data) > 0:
							lbl = 'Dataset {},{}'.format(i+1,res)
							datax[lbl] = {'values' : data,
									      'colors' : color}
							dataxOrder.append(lbl)

		else:
			self.calcAverageMetricOverDoses(metric   = metric,
									   	    normType = normType)

			if resiType == 'all':

				if valType == 'average':

					data = [atm.densMetric[metric][normType]['average'] for atm in atmList]

					datax['average'] = {'values' : data,
								  		'colors' : 'r'}

				else:
					presentAtms = self.getAtomListWithoutPartialAtoms(dataset = valType,
																	  atmList = atmList)
					lbl  = 'Dataset {} '.format(str(valType + 1))

					if resSplit:

						for atm in presentAtms:

							lblFull = lbl + atm.categorise()
							vals = atm.densMetric[metric][normType]['values'][valType]

							try:
								datax[lblFull]['values'].append(vals)
							except KeyError:
								datax[lblFull] = {'values' : [vals],
												  'colors' : ''}
						
							if lblFull not in dataxOrder:
								dataxOrder.append(lblFull)

						for k,col in zip( dataxOrder,sns.color_palette( 'hls', len( dataxOrder ) ) ):
							datax[k]['colors'] = col

					else: 
						data = [atm.densMetric[metric][normType]['values'][valType] for atm in presentAtms]
						datax[lbl] = {'values' : data,
									  'colors' : 'r'}
						dataxOrder.append(lbl)
			else:
				for (res, color) in zip( resiType, sns.color_palette('hls', len(resiType) ) ):
			
					if valType == 'average':
					
						data = [atm.densMetric[metric][normType]['average'] for atm in atmList if atm.basetype == res]
						lbl  = 'average, {}'.format(res)

						datax[lbl] = {'values' : data,
									  'colors' : color}
						dataxOrder.append(lbl)

					else:
						presentAtms = self.getAtomListWithoutPartialAtoms(dataset = valType,
																		  atmList = atmList)
						lbl  = 'Dataset {}, {}'.format(valType+1,res)

						if resSplit:

							# currently should only be performed if 1 single residue type plotted 
							# (otherwise colour scheme will loop around...)
							if len(resiType) != 1:
								print 'resSplit = True only supported here when 1 single '+\
									  'residue type is to be plotted..'
								return {},''

							for atm in presentAtms:
								if atm.basetype == res:
									lblFull = lbl + '; ' + atm.categorise()
									vals = atm.densMetric[metric][normType]['values'][valType]

									try:
										datax[lblFull]['values'].append(vals)
									except KeyError:
										datax[lblFull] = {'values' : [vals],
														  'colors' : ''}

									if lblFull not in dataxOrder:
										dataxOrder.append(lblFull)

							for k,col in zip( dataxOrder,sns.color_palette( 'hls', len( dataxOrder ) ) ):
								datax[k]['colors'] = col

						else: 
							data = [atm.densMetric[metric][normType]['values'][valType] for atm in presentAtms if atm.basetype == res]					

							if data != []:
								datax[lbl] = {'values' : data,
											  'colors' : color}
								dataxOrder.append(lbl)

				# calculate Kolmogrov-Smirnov statistic in case where strictly 2 residue types plotted
				vals = [datax[k]['values'] for k in datax.keys()] 

				if calcKSstat and len( datax.keys() ) == 2:
					(KSstat,KSpval) = ks_2samp(*vals)
					stats['KS'] = {'statistic' : KSstat,
								   'p-value'   : KSpval}

				if calcADstat and len (datax.keys() ) > 1:
					try: 
						(ADstat,critVals,ADpval) = anderson_ksamp(vals)
						stats['AD'] = {'statistic' : ADstat,
									   'p-value'   : ADpval}
					except OverflowError:
						pass # anderson_ksamp can crash unexpectedly

		if len(datax.keys()) == 0:
			return {},''

		self.plotHist(plotType=plotType, datax=datax, order=dataxOrder)

		if xLabel == '':
			xLabel = '{} D{} per atom'.format(normType,metric)

		plt.legend(loc='best', fontsize=legFont)
		plt.xlabel(xLabel, fontsize=axesFont)
		plt.ylabel('Normed-frequency', fontsize=axesFont)
		sns.despine()

		if plotTitle == '':
			t = '{} D{} per atom'.format(normType,metric)
			if resiType != 'all':
				t += ': {}'.format(', '.join(resiType))
			if sideOnly:
				t+= '; side chain only'

			if stats != {}:
				for k in stats.keys():
					t += '\n{} stat: {}, p-value: {}'.format(k,
															 round(stats[k]['statistic'],3),
															 round(stats[k]['p-value'],3))
		else:
			t = plotTitle

		if inclTitle:
			plt.title(t,fontsize = titleFont)

		if not save: 
			plt.show()
		else:

			args = [outputDir,
					''.join(resiType),
					metric,
					normType.replace(" ","")]

			saveName = '{}DistnPlot_Residues-{}_Metric-D{}_Normalisation-{}'.format(*args)

			if sideOnly:
				saveName += '_SideChainOnly'

			if valType != 'all':
				if valType == 'average':
					saveName += '_{}'.format(valType)
				else:
					saveName += '_Dataset-{}'.format(valType)

			fileName = self.checkUniqueFileName(fileName = saveName,
									 			fileType = fileType)
			fig.savefig(fileName)

		return datax,fileName

	def plotHist(self,
				 plotType = 'both',
				 nBins    = 300,
				 datax    = {},
				 order    = []):

		# plot histogram or kde plot for datax and give current label
		# 'nBins' is number of bins (only used if plotType is 'hist' or 'both'

		if order == []:
			order = datax.keys()

		for k in order:
			if plotType == 'hist':
				plt.hist(datax[k]['values'],
						 nBins,
						 histtype = "stepfilled",
						 alpha    = .7,
						 label    = k,
						 color    = datax[k]['colors'])

			elif plotType == 'kde':
				sns.kdeplot(np.array(datax[k]['values']),
							shade = True,
							label = k,
							color = datax[k]['colors'])

			elif plotType == 'both':
				sns.distplot(np.array(datax[k]['values']),
							 label = k,
							 color = datax[k]['colors'])	

	def getKSstatBetween2Residues(self,
								  resType1 = 'Glu',
								  resType2 = 'Gln',
								  sideOnly = False,
								  metric   = 'loss',
								  normType = 'Standard'):

		# calculate the Kolmogrov-Smirnov statistic at each dataset
		# in a series to compare two specified residue types

		if sideOnly:
			atmList = self.getSidechainAtoms()
		else:
			atmList = self.atomList

		KSstatList = []
		for d in self.getDsetList():

			presentAtms = self.getAtomListWithoutPartialAtoms(dataset  = d,
															  atmList  = atmList)

			res1Vals = [atm.densMetric[metric][normType]['values'][d] for atm in presentAtms if atm.basetype == resType1]
			res2Vals = [atm.densMetric[metric][normType]['values'][d] for atm in presentAtms if atm.basetype == resType2]

			if res1Vals == [] or res2Vals == []:
				continue

			(KSstat,pval) = ks_2samp(res1Vals,res2Vals)
			KSstatList.append(KSstat)

		return KSstatList

	def plotKSstatVsDataset(self,
			 				metric    = 'loss',
							normType  = 'Standard',
							residues  = 'all',
							reference = 'GLY',
							sideOnly  = False,
							inclLine  = True,
							plotType  = 'bar',
							saveFig   = True,
							axesFont  = 18,
							titleFont = 16,
							legFont   = 10,
							figTitle  = '',
							saveName  = '',
							outputDir = './',
							fileType  = '.svg',
							palette   = 'hls'):

		# produce a plot of Kolmogrov-Smirnov statistic versus 
		# dataset number for each residue type within the structure, 
		# relative to a reference residue type (typically Gly or Ala)

		if residues == 'all':
			resList = self.getResidues()
		elif isinstance(residues,list):
			resList = residues

		if not isinstance(reference,list):
			reference = [reference]*len(resList)

		if len(reference) != len(resList):
			print 'List of residues and list of reference '+\
				  'residues must be of same length'
			return

		KSperDset = {}
		for res,ref in zip(resList,reference):
			if res != reference:
				KSperDset[res+' relative to '+ref] = self.getKSstatBetween2Residues(resType1 = res,
																			        resType2 = ref,
																			        sideOnly = sideOnly,
																			        metric   = metric,
																			        normType = normType)
		args = [plotType,
				metric,
				normType]

		if figTitle == '':
			title = '{}-plot; metric: D{}, normalisation:{}\n KS-stat'.format(*args)
		else:
			title = figTitle

		# only include residue types present in all datasets
		for k in KSperDset.keys():
			if len(KSperDset[k]) != self.getNumDatasets():
				del KSperDset[k]

		if saveName == '':
			name = '{}plot_Metric-D{}_Normalisation-{}-KSstat'.format(*args)
			if sideOnly:
				name += '_sideChainsOnly'
		else:
			name = saveName

		if plotType == 'line':
			fileName = self.makeLinePlot(palette     = palette,
									 	 vals2Plot   = KSperDset,
									 	 inclLine    = inclLine,
									 	 xVals       = 'doses',
									 	 yLabel      = 'Kolmogrov-Smirnov stat',
									 	 axesFont    = axesFont,
									 	 legFont     = legFont,
									 	 titleFont   = titleFont,
									 	 title       = title,
										 fileType    = fileType,
						 				 saveFig     = saveFig,
						 				 fileName    = name,
						 				 outputDir   = outputDir)

		elif plotType in ('bar','box'):
			pltDic = {'comparison' : [],
					  'dataset'    : []}

			for k in KSperDset.keys():
				for d in KSperDset[k]:
					pltDic['comparison'].append(k)
					pltDic['dataset'].append(round(d,3))
		
			print pltDic
			dFrame = DataFrame(pltDic)
			fileName = self.makeBarPlot(plotData  = dFrame,
									    xName     = "comparison",
									    yName     = "dataset",
									    yLabel    = 'Kolmogrov-Smirnov stat',
									    plotType  = plotType,
									    figName   = name,
									    saveFig   = saveFig,
									    outputDir = outputDir)

		if saveFig:
			return fileName

	def makeLinePlot(self,
					 figSize     = [10,8],
					 palette     = 'hls',
					 vals2Plot   = {},
					 xVals       = 'doses',
					 errors2Plot = {},
					 inclErrors  = False,
					 inclLine    = True,
					 yLabel      = '',
					 axesFont    = 14,
					 legFont     = 10,
					 titleFont   = 16,
					 title       = '',
					 fileType    = '.svg',
					 saveFig     = True,
					 fileName    = 'untitled',
					 outputDir   = '',
					 fullDespine = False):

		# make a generical line plot. vals2Plot and
		# errors2Plot are dictionaries, to allow
		# multiple values to be plotted on same axes

		if palette != 'hls':
			colors = sns.color_palette(palette  = palette,
									   n_colors = len(vals2Plot),
									   desat    = .9)
		else:
			colors = sns.hls_palette(n_colors = len(vals2Plot), 
								 	 l        = .4,
									 s        = .8)

		if inclLine:
			lStyle = '-'
		else:
			lStyle = 'None'

		sns.set_style("ticks")
		sns.set_context("talk",rc = {"figure.figsize":(figSize[0], figSize[1])})
		f  = plt.figure()
		ax = plt.subplot(111)

		# determine y values here dependent 
		# on density metric type specified 
		if xVals == 'doses':
			x = self.doseList
			xLabel = 'Dose (MGy)'
		elif xVals == 'dataset':
			x = self.getDsetList()
			xLabel = 'Dataset'
		else:
			# custom x-values as dictionary
			x      = xVals.values()[0]
			xLabel = xVals.keys()[0]

		numColors = len(vals2Plot.keys())
		for i,k in enumerate(vals2Plot.keys()):
			col = colors[i]
			if not inclErrors:

				if isinstance(vals2Plot[k],list):
					plt.plot(x,
							 vals2Plot[k],
							 label      = k,
							 linestyle  = lStyle,
							 marker     = 'o',
							 markersize = 4,
							 color      = col)
				else:
					for i,key in enumerate(vals2Plot[k].keys()):
						if i == 0:
							plt.plot(x,
								 	 vals2Plot[k][key],
								 	 label      = k,
								 	 linestyle  = lStyle,
								 	 marker     = 'o',
								 	 markersize = 4,
								 	 color      = col)
						else:
							plt.plot(x,
								 	 vals2Plot[k][key],
								 	 linestyle  = lStyle,
								 	 marker     = 'o',
								 	 markersize = 4,
								 	 color      = col)
			else:
				plt.errorbar(x,
							 vals2Plot[k],
							 yerr     = errors2Plot[k],
							 fmt      = '-o',
							 capthick = 2,
							 label    = k,
							 color    = col)

		if fullDespine:
			sns.despine(offset = 0, 
						trim   = True)
		else:
			sns.despine()

		ax.set_xlim([0, x[-1] + 1])

		# place legend outside to right of plot
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax.legend(loc            = 'center left',
				  bbox_to_anchor = (1, 0.5),
				  fontsize       = legFont)
		plt.xlabel(xLabel,
			       fontsize = axesFont)
		plt.ylabel(yLabel,
			       fontsize = axesFont)
		f.suptitle(title,
				   fontsize = titleFont)

		fileName = self.savePlot(outputDir = outputDir,
				                 name      = fileName,
				                 fileType  = fileType,
				         		 save      = saveFig,
				         		 f         = f)

		if saveFig:
			return fileName

	def savePlot(self,
				 outputDir = '',
				 name      = 'untitled',
				 fileType  = '.svg',
				 save      = True,
				 f         = []):

		# save a line plot

		if not save:
			plt.show()
		else:
			fileName = self.checkUniqueFileName(fileName = outputDir + name,
								 				fileType = fileType)
			f.savefig(fileName)

			return fileName

	def graphMetric(self,
					densMet   = 'loss',
					normType  = 'Standard',
					atomType  = '',
					resType   = '',
					resiNum   = '',
					chainType = '',
					errorBars = 'NONE',
					yLabel    = '',
					saveFig   = False,
					outputDir = './',
					fileType  = '.svg',
					axesFont  = 18,
					legFont   = 14,
					titleFont = 20,
					figTitle  = '',
					saveName  = '',
					figSize   = [14,8],
					useDoses  = True,
					palette   = 'hls'):
	
		# produce a graph of selected metric against dataset number 
		# for a specified atom. 'errorBars' takes values in ('NONE',
		# 'RESNUM','ATOMTYPE')

		errorOptions = ('NONE',
					    'RESNUM',
					    'ATOMTYPE')

		if errorBars not in errorOptions: 
			print 'invalid error parameter..\n choose from: {}'.format(errorOptions)
			return

		# find atoms of interest. Can now also handle special case
		# where lists of atom identities are supplied, provided all
		# info (res num, chain, atom type, res type) are supplied.
		# this can only be done for non-error bar plots currently
		if not isinstance(atomType,list):		
			foundAtoms = self.getAtom(restype  = resType,
									  resnum   = resiNum,
									  atomtype = atomType,
									  chain    = chainType)
		else:
			errorBars = 'NONE'
			if (isinstance(resType,list) is False or
				isinstance(chainType,list) is False or
				isinstance(resiNum,list) is False):

				print 'Non-compatible input format. All of "chainType", '+\
					  '"resiNum", "resType" and "atomType" inputs must be '+\
					  'lists if one of them is specified to be a list here.'
				return

			l = len(atomType)
			if (len(resType)   != l or
				len(chainType) != l or
				len(resiNum)   != l):

				print 'Lengths of input lists must be identical..'
				return

			foundAtoms = []
			for i in range(l):
				foundAtom = self.getAtom(restype  = resType[i],
										 resnum   = resiNum[i],
										 chain    = chainType[i],
										 atomtype = atomType[i])

				foundAtoms.append(foundAtom[0])

		if len(foundAtoms) == 0:
			print 'No atoms found..'
			return
		
		if yLabel == '':
			if normType == 'Calpha normalised':
				yLabel = r'$C_\alpha$-normalised D{}'.format(densMet)
			else:
				yLabel = r'{} D{}'.format(normType,densMet)

		# determine y values here dependent 
		# on density metric type specified 
		yValues = {}
		yErrors = {}
		if errorBars == 'NONE':
			inclErrors = False

			for atom in foundAtoms:

				if self.checkMetricPresent(atom = atom, metric = densMet, normType = normType) is False: 
					return # check metric valid

				yValues[atom.getAtomID()] = atom.densMetric[densMet][normType]['values']

		else: 
			# error bars will be plotted
			inclErrors = True
			errorBarKey = {'ATOMTYPE' : 'atomtype',
						   'RESNUM'   : 'residuenum'}
			errorAttr 	= errorBarKey[errorBars]
			yDict 		= {}

			for atom in foundAtoms:

				if self.checkMetricPresent(atom = atom, metric = densMet, normType = normType) is False: 
					return # check metric valid

				vals = atom.densMetric[densMet][normType]['values']

				if getattr(atom,errorAttr) not in yDict.keys():
					yDict[getattr(atom,errorAttr)] = [vals]
				else:
					yDict[getattr(atom,errorAttr)].append(vals)

			for key in yDict.keys():
				yValues[key] = np.mean(yDict[key],0)
				yErrors[key] = np.std(yDict[key],0)

		sns.despine()

		args = [densMet,
				normType,
				resType,
				resiNum,
				atomType]

		if figTitle == '':
			title = 'metric: D{}, normalisation:{}, atom info: {} {} {}'.format(*args)
		else:
			title = figTitle

		if saveName == '':
			name = 'Lineplot_Metric-D{}_Normalisation-{}_{}-{}-{}'.format(*args)
			if errorBars != 'NONE':
				name += '_witherrorbars'
		else:
			name = saveName

		if useDoses:
			xVals = 'doses'
		else:
			xVals = 'dataset'

		fileName = self.makeLinePlot(palette     = palette,
								 	 vals2Plot   = yValues,
								 	 xVals       = xVals,
								 	 errors2Plot = yErrors,
								 	 inclErrors  = inclErrors,
								 	 yLabel      = yLabel,
								 	 axesFont    = axesFont,
								 	 legFont     = legFont,
								 	 titleFont   = titleFont,
								 	 title       = title,
								 	 fileType    = fileType,
					 				 saveFig     = saveFig,
					 				 fileName    = name,
					 				 figSize     = figSize,
					 				 outputDir   = outputDir)

		if saveFig:
			return fileName

	def plotSusceptibleAtoms(self,
							 densMet   = 'loss',
							 normType  = 'Standard',
							 errorbars = True,
							 saveFig   = True,
							 susAtms   = [],
							 fileType  = '.svg',
							 axesFont  = 18,
							 titleFont = 24,
							 outputDir = './'):

		# create a line plot of metric values for susceptible atoms within
		# structure. Susceptible atoms are defined as below

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
			if len(findAtms) == 0: 
				continue
			if not isinstance(findAtms,list): 
				findAtms = [findAtms]
			sDic['-'.join(s)] = findAtms
		
		# define x range here (number in damage series)
		x = self.getDsetList()
		sns.set(style   = "white",
			    context = "talk")
		f, axes = plt.subplots(1, 1, figsize = (12, 12), sharex = True)
		cmap = plt.cm.Set1
		i = -1
		nsteps = len(sDic.keys())
		for key in sDic.keys():
			i += 1
			if not errorbars:
				for atom in sDic[key]:
					y = atom.densMetric[densMet][normType]['values']
					plt.plot(x,
							 y,
							 label = key,
							 color = cmap(i/float(nsteps)))
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

		plt.xlabel('Dataset',
			       fontsize = axesFont)
		plt.ylabel('{} D{}'.format(normType,densMet),
			       fontsize = axesFont)
		plt.legend()
		f.suptitle('{} D{}: susceptible residues'.format(normType,densMet),
			       fontsize = titleFont)
		if not saveFig:
			plt.show()
		else:
			i = 0
			figName = '{}{}_D{}_susceptResis'.format(outputDir,
													 normType,
													 densMet)

			fileName = self.checkUniqueFileName(fileName = figName,
								 				fileType = fileType)
			f.savefig(fileName)

	def plotNumAtomsWithMetricAboveStructureWideMean(self,
													 startThreshold = 0.5,
													 thresholdDiv   = 0.1,
													 metric         = 'loss',
													 normType       = 'Standard',
													 dataset        = 0,
													 atomType       = '',
													 resType        = '',
													 axesFont       = 18,
													 saveFig        = True,
													 fileType       = '.svg',
													 outputLoc      = './',
													 titleFont      = 18,
													 markerType     = 'o',
													 markerColor    = '#4872ae'):

		# plot the number of atoms present within the structure
		# with metric value of n standard deviations above the
		# the structure wide mean (for increasing n until no 
		# atoms remain). Can use 'atomType' and 'resType' to 
		# restrict to atoms of these types only
		
		if dataset == 'all':
			dList = self.getDsetList()
			colours = sns.color_palette('hls', len(dList))
		else:
			dList = [dataset]
			colours = [markerColor]

		sns.set_context(rc = {"figure.figsize":(10, 10)})
		f = plt.figure()

		for d,colour in zip(dList,colours):

			n = self.numAtmsWithMetricAboveLevel(dataset   = d,
												 metric    = metric,
												 normType  = normType,
												 threshold = startThreshold,
												 firstTime = True,
												 atomType  = atomType,
												 resType   = resType)
			numAtomsFound = [n]
			numStds       = [startThreshold]
			t = startThreshold
			while n != 0:
				t += thresholdDiv
				n = self.numAtmsWithMetricAboveLevel(dataset   = d,
											         metric    = metric,
											         normType  = normType,
											         threshold = t,
											         firstTime = False,
											         atomType  = atomType,
												     resType   = resType)
				numAtomsFound.append(n)
				numStds.append(t)

			plt.plot(numStds,
					 numAtomsFound,
					 markerType,
					 markerfacecolor = 'w',
					 markeredgecolor = colour,
					 markeredgewidth = 2,
					 label           = 'Dataset {}'.format(d))

		sns.despine()
		plt.legend(loc = 'best')
		plt.xlabel('# standard deviations from structure mean',
			       fontsize = axesFont)
		plt.ylabel('Number of atoms located',
			       fontsize = axesFont)
		f.suptitle('# atoms with metric value of n standard deviations from structure-wide mean',
				   fontsize = titleFont)

		if not saveFig:
			plt.show()
		else:
			args = [outputLoc,
					metric,
					normType,
					dataset]
			figName = '{}zScorePlot_metric-{}_norm-{}_dataset-{}'.format(*args)

			if atomType != '' and resType != '':
				figName += '_res-{}_atoms-{}'.format(resType,atomType)

			fileName = self.checkUniqueFileName(fileName = figName,
								 				fileType = fileType)
			f.savefig(fileName)

		return fileName

	def checkMetricPresent(self,
						   atom     = '',
						   metric   = 'loss',
						   normType = 'Standard'):

		# method to check whether a method is 
		# valid for selected atom object

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

		try:
			halfDoseApprox
		except NameError:
			from halfDoseCalc import halfDoseApprox

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

		if plotType not in ('hist','kde','both'): 
			return 'Unknown plotting type selected.. cannot plot..'
		if self.checkMetricPresent(metric = metric, normType = normType) is False: 
			return # check metric valid

		sns.set_palette("deep", desat=.8)
		sns.set_context(rc={"figure.figsize": (10, 6)})
		fig = plt.figure()

		atoms = self.getAtom(restype  = restype,
							 atomtype = atomtype)

		datax = {'' : [atm.densMetric[metric][normType]['Half-dose']['Half-dose'] for atm in atoms if atm.densMetric[metric][normType]['Half-dose']['Half-dose'] < 100 ]}
		
		self.plotHist(plotType = plotType,
					  nBins    = nBins,
					  datax    = datax,
					  colour   = {'':colour})

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
						 densMet   = 'loss',
						 normType  = 'Standard',
						 dataset   = 'av',
						 outputDir = './',):

		# plot scatter plots of density metric for 1 atom type against
		# another (both appearing in same residue type). Can be plotted
		# for a single dataset or an average set of values over all 
		# datasets is calculated (dataset = 'av')

		rSquaredDic = {}
		numPairsDic = {}
		dataDic 	= {}
		sensAtoms = [['GLU','CD','CG'],
					 ['GLU','CD','OE1'],
					 ['ASP','CG','CB'],
					 ['ASP','CG','OD1'],
					 ['TYR','OH','CZ']]

		for i in sensAtoms:
			rSquared,numPairs,data = self.compareMetricsBetweenAtoms(restype   = i[0],
																	 atomtype1 = i[1],
																	 atomtype2 = i[2],
																	 densMet   = densMet,
																	 normType  = normType,
																	 dataset   = dataset,
																	 outputDir = outputDir)
			if rSquared != False: 
				rSquaredDic['-'.join(i)] = rSquared
				numPairsDic['-'.join(i)] = numPairs
				dataDic['_'.join(i)] 	 = data
		return rSquaredDic,numPairsDic,dataDic

	def compareSensAtoms_RsquaredLinePlotWithDataset(self,
													 densMet     = 'loss',
						 							 normType    = 'Standard',
						 							 outputDir   = './',
						 							 saveFig     = True,
						 							 fileType    = '.svg',
						 							 plotType    = 'line'):

		# plot a line plot of Rsquared for each scatter plot created by 
		# method .compareSensAtoms() above, over the series of datasets

		yValDic     = {}
		for d in self.getDsetList():
			out = self.compareSensAtoms(densMet   = densMet,
							 		    normType  = normType,
							 		    dataset   = d,
							 		    outputDir = outputDir)
			RsquaredDic = out[0]
			for k in RsquaredDic.keys():
				if k not in yValDic.keys():
					yValDic[k] = [RsquaredDic[k]]
				else:
					yValDic[k].append(RsquaredDic[k])

		title = 'Metric: {}, Normalisation: {}\nCoefficient Of Determination'.format(densMet,normType)
		name  = '{}plot_Metric-D{}_Normalisation-{}_RsquaredValues'.format(plotType,densMet,normType)

		if plotType == 'line':
			fileName = self.makeLinePlot(vals2Plot   = yValDic,
									 	 xVals       = 'doses',
									 	 yLabel      = 'Rsquared',
									 	 title       = title,
									 	 fileType    = fileType,
						 				 saveFig     = saveFig,
						 				 fileName    = name,
						 				 outputDir   = outputDir)

		elif plotType in ('bar','box'):
				pltDic = {'comparison' : [],
						  'dataset'    : []}

				for k in yValDic.keys():
					for d in yValDic[k]:
						pltDic['comparison'].append(k)
						pltDic['dataset'].append(round(d,3))
			
				print pltDic
				dFrame   = DataFrame(pltDic)
				fileName = self.makeBarPlot(plotData  = dFrame,
										    xName     = "comparison",
										    yName     = "dataset",
										    yLabel    = 'R-squared stat',
										    plotType  = plotType,
										    figName   = name,
										    saveFig   = saveFig,
										    outputDir = outputDir)

		if saveFig:
			return fileName

	def compareMetricsBetweenAtoms(self,
								   restype   = '',
								   atomtype1 = '',
								   atomtype2 = '',
								   densMet   = 'loss',
								   normType  = 'Standard',
								   fileType  = '.svg',
								   dataset   = 0,
								   outputDir = ''):

		# for two atom types atomtype1 and atomtype2, compare 
		# density metric between two atoms within same side-chain

		atoms1 = self.getAtom(restype  = restype,
							  atomtype = atomtype1)

		atoms2 = self.getAtom(restype  = restype,
							  atomtype = atomtype2)

		if len(atoms1) == 0 or len(atoms2) == 0: 
			return False,0,{'x':[],'y':[]}

		if dataset == 'av': 
			self.calcAverageMetricOverDoses(metric   = densMet,
									    	normType = normType)

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
		figTitle  = 'Metric: {}, Normalisation: {}\n for {} atoms, Dataset: {}'.format(densMet,
																					   normType,
																					   restype,
																					   dataset)

		saveTitle = 'scatterplot_metric-{}_normalisation-{}_{}-{}-{}_Dset-{}'.format(densMet,
																				     normType.replace(' ',''),
																				     restype,
																				     atomtype1,
																				     atomtype2,
																				     dataset)
		rSquared = self.plotScatterPlot(xData       = xData,
										yData       = yData,
										xLabel      = '{}-{} D{}'.format(restype,atomtype1,densMet),
							 			ylabel      = '{}-{} D{}'.format(restype,atomtype2,densMet),
							 			figtitle    = figTitle,
							 			fileType    = fileType,
							 			saveTitle   = saveTitle,
							 			lineBestFit = True,
							 			outputDir   = outputDir)

		data = {'x' : xData,
				'y' : yData}

		return rSquared,len(xData),data

	def compareMetrics(self,
					   restype   = '',
					   atomtype  = '',
					   metric1   = 'loss',
					   metric2   = 'loss',
					   normType1 = 'Standard',
					   normType2 = 'Standard',
					   dSet      = 0,
					   fileType  = '.svg',
					   outputDir = './'):

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
										lineBestFit = True,
										outputDir   = outputDir)

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
						titleFont    = 24,
						outputDir    = ''):

		# plot the relationship between xData and yData.
		# if 'lineBestFit' is True then linear line of best fit 
		# is calculated and plotted. if 'yequalsxLine' is 
		# True then line y=x plotted

		sns.set_context("talk")
		fig = plt.figure(figsize = (10, 10))
		ax  = plt.subplot(111)

		plt.scatter(xData,
					yData,
					marker     = 'o',
					s          = 100,
					c          = colors,
					edgecolors = '#FFFFFF',
					cmap       = 'copper')

		if lineBestFit: 
			# plot line of best fit if specified
			slope, intercept, r_value, p_value, std_err = linregress(xData,yData)
			x = np.linspace(min(xData), max(xData),10)
			y = slope*np.array(x) + np.array(intercept)
			plt.plot(x,
					 y,
					 '-',
					 color     = '#d11141',
					 linewidth = 3)
			figtitle += ', R^2: {}'.format(round(r_value**2,2))

		if yequalsxLine:
			# plot y=x line if specified
			x = np.linspace(min(xData), max(xData),10)
			y = x
			plt.plot(x, 
					 y,
					 ':',
					 color     = '#9890af',
					 linewidth = 3)

		plt.xlabel(xLabel,
			       fontsize = axesFont)
		plt.ylabel(ylabel,
			       fontsize = axesFont)
		fig.suptitle(figtitle,
					 fontsize = titleFont)
		sns.despine()

		if outputDir == '':
			outDir = self.outputDir
		else:
			outDir = outputDir

		figName = outDir + saveTitle
		fileName = self.checkUniqueFileName(fileName = figName,
								 			fileType = fileType)
		fig.savefig(fileName)

		if lineBestFit: 
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
							    atom1 = '',
							    atom2 = ''):

		# get distance between two atoms

		xyz1 = np.array([atom1.X_coord,atom1.Y_coord,atom1.Z_coord])
		xyz2 = np.array([atom2.X_coord,atom2.Y_coord,atom2.Z_coord])
		return np.linalg.norm(xyz1-xyz2)

	def getMinDistToAtomType(self,
							 atomType1  = 'OH',
							 resType1   = 'TYR',
							 resNum1    = '330',
							 chainType1 = 'M',
							 atomType2  = 'SG',
							 resType2   = 'CYS'):

		# get min distance from specified atom 1 
		# to atom of type resType2 & atomType2

		atm1 = self.getAtom(chain    = chainType1,
						    restype  = resType1,
						    resnum   = resNum1,
						    atomtype = atomType1)

		atms2 = self.getAtom(restype  = resType2,
						     atomtype = atomType2)

		minAtm  = []
		minDist = 1e6
		for atm2 in atms2:
			dist = self.getDistanceBetweenAtoms(atom1 = atm1[0],
										 		atom2 = atm2)
			try:
				minDist
			except NameError:
				minDist = dist
				minAtm  = atm2
				continue

			if minDist > dist:
				minDist = dist
				minAtm  = atm2

		return minAtm,minDist

	def scatterAtmsByDistToOtherAtms(self,
								     atomType1 = 'OH',
								     resType1  = 'TYR',
								     atomType2 = 'CB',
								     resType2  = 'CYS',
								     metric    = 'loss',
								     normType  = 'Standard',
								     outputDir = './'):

		# rank a particular type of atom (e.g. TYR-OH)
		# by shortest distance to another atom type
		# (e.g. CYS-SG) and make scatter plot

		atms1 = self.getAtom(restype   = resType1,
						     atomtype  = atomType1)

		xVals = []
		yVals = []
		for atm1 in atms1:
			info = self.getMinDistToAtomType(atomType1  = atm1.atomtype,
										     resType1   = atm1.basetype,
										     resNum1    = atm1.residuenum,
										     chainType1 = atm1.chaintype,
										     atomType2  = atomType2,
										     resType2   = resType2)
			xVals.append(info[1])
			yVals.append(atm1.densMetric[metric][normType]['values'])

		for d in self.getDsetList():
			Rsquared = self.plotScatterPlot(xData       = xVals,
										    yData       = [y[d] for y in yVals],
										    lineBestFit = True,
										    xLabel      = 'Distance to {}-{}'.format(resType2,atomType2),
										    ylabel      = '{} D{}'.format(normType,metric),
										    figtitle    = 'Min {}-{} distance vs {}-{} D{}'.format(resType2,atomType2,resType1,atomType1,metric),
										    saveTitle   = 'metric-{}_normalisation-{}_scatter_{}-{}_minDistTo_{}-{}_Dset-{}'.format(metric,normType,resType1,atomType1,resType2,atomType2,d),
										    outputDir   = outputDir)
			
			print 'Dataset {} --> R^2 = {}'.format(d+1,round(Rsquared,4))

	def getAtomsWithinDist(self,
						   atom       = '',
						   distLimMax = 4,
						   distLimMin = 0,
						   printText  = False,
						   ignoreSame = 'res'):

		# find all atoms within 'distLim' angstroms 
		# of selected atomtype. Disregards atoms of 
		# same exact residue by default

		nearAtmDic = {'atoms'     : [],
					  'distances' : []}

		for otheratom in self.atomList:
			if ignoreSame == 'res':
				if (otheratom.getAtomID()).split('-')[:-1] == (atom.getAtomID()).split('-')[:-1]: 
					continue
			elif ignoreSame == 'atom':
				if otheratom.getAtomID() == atom.getAtomID():
					continue
			dist = self.getDistanceBetweenAtoms(atom1 = atom,
												atom2 = otheratom)
			if dist < distLimMax:
				if dist > distLimMin:
					nearAtmDic['atoms'].append(otheratom)
					nearAtmDic['distances'].append(dist)
		if printText:
			print '{} in total within {}-{} Angstrom of {}'.format(len(nearAtmDic['atoms']),
																   round(distLimMin,2),
																   round(distLimMax,2),
																   atom.getAtomID())
		return nearAtmDic

	def groupAtomsByDist(self,
						 atomType  = 'OH',
						 resType   = 'TYR',
						 resNum    = '0',
						 chainType = 'A'):

		# group atoms by their distance 
		# away from a specified atom

		atm = self.getAtom(chain    = chainType,
						   restype  = resType,
						   resnum   = resNum,
						   atomtype = atomType)

		if atm == []:
			print 'Error: no atom found matching description'
			return

		# find max distance away from atom
		keyAtm = atm[0]
		maxDist = 0
		for othAtm in self.atomList:
			dist = self.getDistanceBetweenAtoms(atom1 = keyAtm,
												atom2 = othAtm)
			if dist > maxDist:
				maxDist = dist

		foundAtms = self.getAtomsWithinDist(atom       = keyAtm,
											distLimMax = maxDist + 1)

		for i,atm in enumerate(foundAtms['atoms']):
			dist = foundAtms['distances'][i]

		distArr  = np.array(foundAtms['distances'])
		groupedDic = {}
		order      = []
		for d in range(1,int(maxDist)+1):
			inds = np.where((distArr>(d-1)) & (distArr<d))

			try:
				atmsAtDist = [foundAtms['atoms'][i] for i in inds[0]]
			except TypeError:
				atmsAtDist = []
			print '{} atoms {}-{} Angstroms away'.format(len(atmsAtDist),d-1,d)
			dicID = '{}-{}'.format(d-1,d)
			groupedDic[dicID] = atmsAtDist
			order.append(dicID)
		return groupedDic,order

	def metricValsAtDistFromAtom(self,
								 atomType  = 'OH',
								 resType   = 'TYR',
								 resNum    = '0',
								 chainType = 'A',
								 metric    = 'loss',
								 normType  = 'Standard',
								 saveFig   = True,
								 outputDir = './'):

		# calculate metric statistics for groups of 
		# atoms are certain distances from a specified
		# fixed atom positon in structure

		atmDic,order = self.groupAtomsByDist(atomType  = atomType,
										     resType   = resType,
										     resNum    = resNum,
										     chainType = chainType)

		for stat in ('mean','max','min'):

			xDic  = {'Distance': []}
			yDic  = {'Dataset {}'.format(d) : [] for d in self.getDsetList()}

			for k in order:
				atmVals = [atm.densMetric[metric][normType]['values'] for atm in atmDic[k]]
				if atmVals == []:
					# no atoms found at this distance
					statVals = [np.nan]*self.getNumDatasets()
				else:
					if stat == 'mean':
						statVals = np.mean(atmVals,0)
					elif stat == 'max':
						statVals = np.max(atmVals,0)
					elif stat == 'min':
						statVals = np.min(atmVals,0)

				for d in self.getDsetList():
					yDic['Dataset {}'.format(d)].append(statVals[d])
				xDic['Distance'].append(float(k.split('-')[-1]))

			args = [normType,
					metric,
					chainType,
					resType,
					resNum,
					atomType,
					stat]

			plotTitle = '{} D{}: distance from {}-{}-{}-{}; {}'.format(*args)
			fileName = 'normalisation-{}_metric-{}_distFrom-{}-{}-{}-{}_{}'.format(*args)

			fileName = self.makeLinePlot(vals2Plot   = yDic,
									 	 xVals       = xDic,
									 	 yLabel      = 'D{}'.format(metric),
									 	 title       = plotTitle,
						 				 saveFig     = saveFig,
						 				 fileName    = fileName,
						 				 outputDir   = outputDir)

	def hotellingTsquareTest(self,
							 atomType1  = 'OH',
							 atomType2  = 'CD',
							 resType1   = 'TYR',
							 resType2   = 'GLU',
							 resNum1    = '',
							 resNum2    = '',
							 metric     = 'loss',
							 normType   = 'Standard',
							 atmVals1   = [],
							 atmVals2   = [],
							 alpha      = 0.05,
							 printText  = True,
							 roundTo    = 2,
							 roundMVals = 8):

		# this function is designed to calculate Hotelling's
		# T-square test statistic in order to compare two sets 
		# of atoms over multiple datasets simultaneously 
		# (each dose level is a different dimension to the data)
		# If atmVals1 and atmVals2 given, then these override
		# any new atoms found by method

		if printText:
			print '-------------------------\n'+\
				  'Hotelling T-squared Test:'
		atmDens1 = []
		atmDens2 = []

		if atmVals1 == [] or atmVals2 == []:

			if printText:
				args = [resType1,
						atmType1,
						resNum1,
						resType2,
						atmType2,
						resNum2]
				print 'Comparing atoms {}-{}-{} and {}-{}-{}'.format(*args)

			# get atoms to be compared
			atoms1 = self.getAtom(restype  = resType1,
								  resnum   = resNum1,
								  atomtype = atomType1)

			atoms2 = self.getAtom(restype  = resType2,
								  resnum   = resNum2,
								  atomtype = atomType2)

			# get densities for atoms
			for atms,atmDens in zip([atoms1,atoms2],[atmDens1,atmDens2]):
				for atm in atms:
					col = []
					for dens in atm.densMetric[metric][normType]['values']:
						col.append([dens])
					atmDens.append(col)
		else:
			if printText:
				print 'Comparing pre-determined values'
			for atmsVals,atmDens in zip([atmVals1,atmVals2],[atmDens1,atmDens2]):
				for val in atmsVals:
					vals = [[v] for v in val]
					atmDens.append(vals)

		# calculate means for unbound and bound vectors
		nx = len(atmDens1)
		ny = len(atmDens2)
		mean1 = np.sum( np.array(atmDens1), axis = 0 )/nx
		mean2 = np.sum( np.array(atmDens2), axis = 0 )/ny

		# calculate the unbiased pooled covariance matrix estimate
		W = 0
		for i in range(0,nx):
			part1 = (np.array(atmDens1)[i] - mean1)*np.transpose((np.array(atmDens1)[i] - mean1))
			part2 = (np.array(atmDens2)[i] - mean2)*np.transpose((np.array(atmDens2)[i] - mean2))
			combinedParts = part1 + part2
			W = W + combinedParts
		W = W/(nx+ny-2)

		if roundMVals != 0:
			W = np.around(W,roundMVals)

		# now calculate the Hotelling's two-sample T-squared statistic
		print W
		tSquared = (float(nx*ny)/(nx+ny))*np.dot(np.transpose(mean1-mean2),np.dot(np.linalg.inv(W),(mean1-mean2)))
		if printText:
			print 'tSquared is as follows: {}'.format(round(tSquared[0][0],roundTo))

		# now calculate the related F test statistic from F distribution
		p = self.getNumDatasets()
		F = (float(nx+ny-p-1)/((nx+ny-2)*p))*tSquared
		if printText:
			print 'F-value is as follows: {}'.format(round(F[0][0],roundTo))

		# generate a p-value for the given F test statistic
		df1 = p
		df2 = nx+ny-1-p
		p_value = 1-f.cdf(F, df1, df2)
		if printText:
			print 'p-value is as follows: {}'.format(round(p_value[0][0],roundTo))

		reject = 'notsureyet'
		if p_value < alpha:
			if printText:
				print 'p_value < {}'.format(alpha)
				print '---> reject null hypothesis that means are equal'
			reject = 'YES'
		else:
			if printText:
				print 'p_value > {}'.format(alpha)
				print '---> cannot reject null hypothesis that means are equal'
			reject = 'NO'

		return F,p_value,reject

	def plot_densMetSurroundAtmsCorrel(self,
									   restype       = 'TYR',
								  	   atomtype      = 'OH',
								  	   distLim       = 4,
								  	   errorbars     = True,
								  	   normType      = 'Standard',
								  	   densMet       = 'loss',
								  	   crystConts    = True,
								  	   pdbName       = '',
								  	   symmetrygroup = '',
								  	   outputDir     = './',
								  	   fileType      = '.svg',
								  	   saveFig       = True,
								  	   hotelling     = True,
								  	   printText     = True):

		# run .densMetSurroundAtmsCorrel() method for 
		# each dataset and plot versus dose

		if errorbars:
			plotDic = {'H-bond' : [],
					   'Other'  : []}
		else:
			plotDic = {'H-bond' : {},
					   'Other'  : {}}

		errBars = {'H-bond' : [],
				   'Other'  : []}

		for i,d in enumerate(self.getDsetList()):
			output = self.densMetSurroundAtmsCorrel(restype       = restype,
												    atomtype      = atomtype,
												    distLim       = distLim,
												    dataset       = d,
												    normType      = normType,
												    densMet       = densMet,
												    crystConts    = crystConts,
												    pdbName       = pdbName,
												    symmetrygroup = symmetrygroup,
												    plotScatter   = True,
												    printText     = printText,
												    outputDir     = outputDir)


			if hotelling:
				if i == 0:
					HTgroup1 = [[v] for v in output[1]]
					HTgroup2 = [[v] for v in output[3]]
				else:
					for n,group in zip([1,3],[HTgroup1,HTgroup2]):
						for j,v in enumerate(output[n]):
							group[j].append(v)

			if errorbars:
				for n,t in zip([1,3],['H-bond','Other']):
					plotDic[t].append(np.mean(output[n]))
					errBars[t].append(np.std(output[n]))
			else:
				for n,t in zip([1,3],['H-bond','Other']):
					for i,j in enumerate(output[n]):
						try:
							plotDic[t][i].append(j)
						except KeyError:
							plotDic[t][i] = [j]

		ylabel = '{} D{}'.format(normType,densMet)
		title  = 'Effect of carboxylate H-bonds on {}-{} damage'.format(restype,atomtype)
		name   = 'carboxylateHbond_vs_{}-{}_damage-normalisation-{}_metric-{}'.format(restype,atomtype,normType,densMet)

		fileName = self.makeLinePlot(vals2Plot   = plotDic,
									 figSize     = [12,8],
								 	 xVals       = 'doses',
								 	 errors2Plot = errBars,
					 				 inclErrors  = errorbars,
								 	 yLabel      = ylabel,
								 	 title       = title,
					 				 saveFig     = saveFig,
					 				 fileType    = fileType,
					 				 fileName    = name,
					 				 outputDir   = outputDir,
					 				 palette     = 'deep')

		if hotelling:
			self.hotellingTsquareTest(atmVals1  = HTgroup1,
							 		  atmVals2  = HTgroup2)

		if saveFig:
			return fileName

	def densMetSurroundAtmsCorrel(self,
								  restype       = 'TYR',
								  atomtype      = 'OH',
								  distLim       = 4,
								  dataset       = '0',
								  normType      = 'Standard',
								  densMet       = 'loss',
								  crystConts    = True,
								  pdbName       = '',
								  symmetrygroup = '',
								  fileType      = '.svg',
								  plotScatter   = True,
								  printText     = True,
								  outputDir     = './'):

		# determine whether there is a correlation between density metric 
		# and types of surrounding atoms if 'crystConts' is True then crystal
		# contacts will also be located here

		if crystConts:
			from NCONTjob import NCONTjob

			ncont = NCONTjob(inputPDBfile = pdbName,
							 outputDir    = self.outputDir,
							 symGroup     = symmetrygroup,
							 seriesName   = self.seriesName,
							 printText    = printText)
			success = ncont.run()

			if not success: 
				return
			logFile = self.outputDir+'/'+ncont.outputLogfile

		atoms = self.getAtom(restype  = restype,
							 atomtype = atomtype)

		groupA,groupB,groupAContacts = [],[],[]
		for atom in atoms:
			nearCarboxyl = False

			if crystConts:
				atomFound = False
				readLog = open(logFile,'r')
				for l in readLog.readlines():
					if atomFound and len(l[0:5].strip()) != 0: 
						atomFound = False # reset if new atom info reached

					if (atom.basetype == l[11:14] and
						str(atom.residuenum) == l[7:10].strip() and
						atom.chaintype == l[4] and
						atom.atomtype == l[19:21]):
						atomFound = True
						if printText:
							print 'found atom: {}'.format(atom.getAtomID())

					if atomFound:
						if (l[39:42] in ('GLU','ASP') and
							l[47:50] in ['OE1','OE2','OD1','OD2'] and
							float(l[58:62]) < distLim):
							resNum = l[35:38].strip()
							chain = l[32]
							if printText:
								print 'contact found! - {}-{}-{}-{}'.format(chain,resNum,l[39:42],l[47:50])
							contactAtom = self.getAtom(chain    = chain,
													   restype  = l[39:42],
													   resnum   = resNum,
													   atomtype = l[47:50])
							nearCarboxyl = True
							break
				readLog.close()

			else:
				nearAtms = self.getAtomsWithinDist(atom       = atom,
												   distLimMax = distLim)
				nearCarboxyl = False
				for atom2 in nearAtms['atoms']:
					if atom2.atomtype in ['OE1','OE2','OD1','OD2']:
						if atom2.basetype in ['GLU','ASP']:
							contactAtom = atom2
							nearCarboxyl = True
							break

			if nearCarboxyl:
				groupA.append(atom)
				groupAContacts.append(contactAtom[0])
			else:
				groupB.append(atom)

		if dataset == 'average':
			self.calcAverageMetricOverDoses(metric   = densMet,
										    normType = normType)
			densValsA 		  = [atom.densMetric[densMet][normType]['average'] for atom in groupA]
			densValsAContacts = [atom.densMetric[densMet][normType]['average'] for atom in groupAContacts]
			densValsB 		  = [atom.densMetric[densMet][normType]['average'] for atom in groupB]
		else: 
			if printText:
				print '\nFor Dataset {}:'.format(dataset)
			densValsA 		  = [atom.densMetric[densMet][normType]['values'][dataset] for atom in groupA]
			densValsAContacts = [atom.densMetric[densMet][normType]['values'][dataset] for atom in groupAContacts]
			densValsB 		  = [atom.densMetric[densMet][normType]['values'][dataset] for atom in groupB]
		strA = 'For the {} atoms near carboxyl groups, mean avDloss: {}, std: {}'.format(len(densValsA),round(np.mean(densValsA),3),round(np.std(densValsA),2))
		strB = 'For the {} atoms not near carboxyl groups, mean avDloss: {}, std: {}'.format(len(densValsB),round(np.mean(densValsB),3),round(np.std(densValsB),2))
		
		if printText:
			print strA+'\n'+strB

		# determine how solvent accessibility changes for each of these atoms upon Glu/Asp decarboxylation.
		# this will dictate the colour scheme for the resulting scatter plot
		scatterColor = []
		solvAccDic,rSqd = self.compareSolvAccessWithAndWithoutGluAspGroups(pdbName       = pdbName,
																	       symmetrygroup = symmetrygroup,
																	       dataset       = dataset,
																		   resType       = restype,
																		   atomType      = atomtype,
																		   specialAtoms  = groupA,
																		   densMet       = densMet,
																		   normType      = normType,
																		   plotScatter   = plotScatter,
																		   printText     = printText,
																		   outputDir     = outputDir)
		print 'Dataset {} --> R^2 = {}'.format(dataset,round(rSqd,4))
		count = 0
		for atm in groupA: 
			atmID = atm.getAtomID()
			solvAcc = solvAccDic[atmID]
			diff = float(solvAcc[1]) - float(solvAcc[0])

			scatterColor.append(float(solvAcc[1]))

			if diff == 0:
				if printText:
					print '{}: No change in solvent accessibility'.format(atmID)
			else:
				if printText:
					print '{}: Solvent accessibility increased by {} upon nearby decarboxylation'.format(atmID,round(diff,2))
				count += 1
		if count < len(groupA):
			if printText:
				print 'Not all {}-{} atoms experience change in solvent accessibility upon nearby decarboxylation'.format(restype,atomtype)

		if plotScatter:
			# plot the relationship between TYR-OH density and nearby carboxyl atom density
			self.plotScatterPlot(xData    = densValsA,
								 yData    = densValsAContacts,
								 xLabel   = '{}-{} D{}'.format(restype,atomtype,densMet),
								 ylabel   = 'Carboxyl-contact D{}'.format(densMet),
								 figtitle  = '{}-{} vs carboxyl-contact D{}'.format(restype,atomtype,densMet),
								 saveTitle = 'D{}Scatter_{}-{}_carboxylContacts'.format(densMet,restype,atomtype),
								 fileType  = fileType,
								 colors    = scatterColor,
								 outputDir = outputDir)

		return strA+'\n'+strB,densValsA,densValsAContacts,densValsB,scatterColor

	def susceptAtmComparisonBarplot(self,
									metric    = 'loss',
									normType  = 'Standard',
									dataset   = 0,
									set       = 1,
									box       = 'Box',
									fileType  = '.svg',
									axesFont  = 18,
									outputDir = './',
									saveFig   = True):

		# produce a barplot to compare the damage metric of 
		# susceptible atom types at a given dataset

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
			if len(foundAtms) == 0: 
				continue 
			for foundAtm in foundAtms:
				atmId 	= '-'.join(atm)
				plotData['x'].append(atmId)
				plotData['y'].append(foundAtm.densMetric[metric][normType]['values'][dataset])

		args = [box,
			    normType,
			    metric,
			    dataset,
			    set]

		figName = 'SusceptAtms{}plot-{}-D{}_{}-set{}'.format(*args)

		if box == 'Box':
			ax = sns.boxplot(x    = "x",
							 y    = "y",
							 data = plotData)
		elif box == 'Bar':
			self.makeBarPlot(plotData  = plotData,
							 xName     = "x",
							 yName     = "y",
							 xLabel    = 'Atom type',
							 yLabel    = '{} D{}'.format(normType,metric),
							 figName   = figName,
							 axesFont  = axesFont,
							 saveFig   = saveFig,
							 outputDir = outputDir)

		else:
			print 'Unknown graph type - specify "Box" or "Bar"'
			return

	def makeBarPlot(self,
					plotData  = {},
					figSize   = [10,8],
					color     = '#0095b6',
					xName     = 'x-data',
					yName     = 'y-data',
					plotType  = 'bar',
					xLabel    = '',
					yLabel    = '',
					xRot      = 90,
					yRot      = 0,
					axesFont  = 14,
					figName   = '',
					fileType  = '.svg',
					saveFig   = True,
					outputDir = './'):

		# make a general bar plot 

		sns.set_style("ticks")
		sns.set_context("talk", rc = {"figure.figsize":(figSize[0], figSize[1])})
		f  = plt.figure()
		ax = plt.subplot(111)

		if plotType == 'bar':
			ax = sns.barplot(x          = xName,
							 y          = yName,
							 data       = plotData,
							 color      = color,
							 saturation = .5)

		elif plotType == 'box':
			ax = sns.boxplot(x          = xName,
							 y          = yName,
							 data       = plotData,
							 color      = color,
							 saturation = .5)
		else:
			print 'Unknown plot type "{}"'.format(plotType)
			return

		if xLabel == '':
			xLabel = xName
		if yLabel == '':
			yLabel = yName

		plt.xlabel(xLabel,
				   fontsize = axesFont)
		plt.ylabel(yLabel,
				   fontsize = axesFont)
		plt.xticks(rotation = xRot) 
		plt.yticks(rotation = yRot) 
		sns.despine()

		fileName = self.savePlot(outputDir = outputDir,
				                 name      = figName,
				                 fileType  = fileType,
				         		 save      = saveFig,
				         		 f         = f)

		if saveFig:
			return fileName

	def compareSolvAccessWithAndWithoutGluAspGroups(self,
													pdbName       = '',
													symmetrygroup = '',
													dataset       = 'average',
													resType       = '',
													atomType      = '',
													specialAtoms  = '',
													densMet       = 'loss',
													normType      = 'Standard',
													fileType      = '.svg',
													plotScatter   = True,
													printText     = True,
													outputDir     = './'):

		# determine the change in solvent accessibility for 'resType'-'atomType' 
		# atoms before & after Glu/Asp decarboxylation events 'specialAtoms' 
		# is an optional parameter, which specifies a subset of atoms objects 
		# which should be coloured differently on the resulting scatter plot 
		# below - the idea is to run .densMetSurroundAtmsCorrel() above first 
		# to get the list of atoms exhibiting contacts to Glu/Asp CO2 groups 
		# and then colour these differently below

		from ccp4Job import ccp4Job

		pdbIn = open(pdbName,'r')
		strippedPdb = pdbName.replace('.pdb','_noGluAspCO2groups.pdb')
		pdbOut = open(strippedPdb,'w')
		for l in pdbIn.readlines():
			if l[0:6].strip() not in ('ATOM','ANISOU'):
				pdbOut.write(l)
			elif l[0:6] != 'ANISOU':
				if l[13:21].replace(' ','') not in ('CDGLU','OE1GLU','OE2GLU','CGASP','OD1ASP','OD2ASP'):
					pdbOut.write(l)
		pdbIn.close()
		pdbOut.close()

		# run AREAIMOL in CCP4 to calculate solvent accessibility for each pdb file
		count = 0
		for pdb in (pdbName,strippedPdb):
			count += 1
			commandInput1 = 'areaimol XYZIN "{}" XYZOUT "{}"'.format(pdb,pdb.replace('.pdb','_areaimol.pdb'))
			commandInput2 = 'DIFFMODE OFF\nMODE -\nNOHOH\nSMODE IMOL\nREPORT -\nCONTACT -\n'+\
							'YES -\nRESAREA -\nYES\nPNTDEN 10\nPROBE 1.4\nOUTPUT\nEND'
			run = ccp4Job('areaimol',commandInput1,commandInput2,self.outputDir,'areaimolRun{}.txt'.format(count),'')

		# these are names of resulting pdb files
		fullPDB 	= pdbName.replace('.pdb','_areaimol.pdb')
		partialPDB 	= strippedPdb.replace('.pdb','_areaimol.pdb')

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
				solvAcc = []
				atms = [atm]+TyrRingAtms
				for a in atms:
					solvAcc.append(float(a.findSolventAccessibility(inputPDBfile = pdb)))
				solvAccDic[atm.getAtomID()].append(np.mean(solvAcc))

		# remove pdb files since no longer needed after job run
		for f in [fullPDB,partialPDB]:
			os.remove(f)

		plotData = {'x'       : [],
					'y'       : [],
					'colours' : []}

		for atm in TyrOHatms:
			if atm in specialAtoms:
				plotData['colours'].append('r')
			else:
				plotData['colours'].append('b')
			key = atm.getAtomID()

			plotData['x'].append(float(solvAccDic[key][1]))

			if dataset == 'average':
				plotData['y'].append(atm.densMetric[densMet][normType]['average'])
			else:
				plotData['y'].append(atm.densMetric[densMet][normType]['values'][dataset])

			args = [key,
					round(solvAccDic[key][0],2),
					round(solvAccDic[key][1],2),
					round(plotData['y'][-1],2)]

			if printText:
				print '{}: {} --> {} ... Dloss: {}'.format(*args)

		rSquared = np.nan
		if plotScatter:
			# plot the relationship between TYR-OH density and nearby carboxyl atom density
			rSquared = self.plotScatterPlot(xData      = plotData['x'],
						                    yData       = plotData['y'],
										    xLabel      = 'Solvent Accessibility',
										    lineBestFit = True,
										    ylabel      = '{}-{} D{}'.format(resType,atomType,densMet),
										    fileType    = fileType,
										    saveTitle   = 'D{}Scatter_{}-{}_vsSolvAccChange_Dset-{}'.format(densMet,resType,atomType,dataset),
										    colors      = plotData['colours'],
										    outputDir   = outputDir)

		return solvAccDic,rSquared

	def calcAverageMetricOverDoses(self,
								   metric   = 'loss',
								   normType = 'Standard'):

		# for a specified metric calculate the average
		# over dose range if not already calculated

		try:
			self.atomList[0].densMetric[metric][normType]['average']
		except KeyError:
			self.calcAdditionalMetrics(metric    = metric,
								       normType  = normType,
								       newMetric = 'average')

	def findProbHighNeighbourGivenHighAtom(self,
										   distance = 5,
										   densMet  = 'loss',
										   normType = 'Standard',
										   dataset  = 0,
										   criteria = '1std',
										   sign     = 'above'):

		# attempts to determine the probability that given that a 
		# randomly picked atom has a defined 'high' metric value
		# then there exists a neighbouring atom that is also mean.
		# Here 'high' is defined by 'criteria' as either 'mean',
		# meaning that any value above/below the average structure 
		# value for the metric is considered as 'high' (depending on 
		# the specification of 'sign'), or can be set as a specific 
		# float threshold

		if criteria == 'mean' or 'std' in criteria:
			av,std = self.getAverageMetricVals(densMet  = densMet,
				                          	   normType = normType)
			if criteria == 'mean':
				thres = av[dataset]
			elif 'std' in criteria:
				numStd = int(criteria.replace('std',''))
				if sign == 'above':
					thres = av[dataset] + numStd*std[dataset]
				else:
					thres = av[dataset] - numStd*std[dataset]
		else:
			try: 
				thres = float(criteria)
			except ValueError:
				print 'Unexpected assignment of "criteria" parameter'

		numHighAtoms     = 0
		numHighNearAtoms = 0
		for atm in self.atomList:
			atmMetric = atm.densMetric[densMet][normType]['values'][dataset]
			if sign == 'above':
				if atmMetric <= thres:
					continue
			elif sign == 'below':
				if atmMetric >= thres:
					continue
			else:
				print '"sign" parameter must take either "above" or "below"'
				return
			numHighAtoms += 1
			nearAtms = self.getAtomsWithinDist(atom       = atm,
										   	   distLimMax = distance)
			for nearAtm in nearAtms['atoms']:
				nearAtmMetric = nearAtm.densMetric[densMet][normType]['values'][dataset]
				if sign == 'above':
					if nearAtmMetric > thres:
						numHighNearAtoms += 1
						break
				else:
					if nearAtmMetric < thres:
						numHighNearAtoms += 1
						break
		probHighAtom = float(numHighAtoms)/self.getNumAtoms()
		probHighNeighbourAndHighAtom = float(numHighNearAtoms)/self.getNumAtoms()
		probHighNeighGivenHighAtom = probHighNeighbourAndHighAtom/probHighAtom

		print 'Probability of high ({} {}) neighbour atom GIVEN high atom: {}'.format(sign,thres,probHighNeighGivenHighAtom)
		
		return probHighNeighGivenHighAtom

	def calculateLocalDloss(self,
							resType   = '',
							atomType  = '',
							distance  = 4,
							densMet   = 'loss',
							normType  = 'Standard',
							weighted  = False,
							fileType  = '.svg',
							printText = True):

		# for a given residue, calculate the average density metric 
		# (for a specified metric) within a given region of space 
		# around that particular residue. if 'weighted' is True 
		# then a density-weighted mean is instead calculated

		plotData = {'x' : [],
					'y' : []}

		# find specified atoms within structure
		atms = self.getAtom(restype  = resType,
							atomtype = atomType)
		if len(atms) == 0: 
			return plotData # don't continue if no atoms exist

		self.calcAverageMetricOverDoses(metric   = densMet,
									    normType = normType)

		for atm in atms:
			nearAtms = self.getAtomsWithinDist(atom       = atm,
											   distLimMax = distance)

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
			if weighted:
				meanDens = np.average(densList, weights = np.square(np.reciprocal(distList)))
			else:
				meanDens = np.mean(densList)

			dens = atm.densMetric[densMet][normType]['average']

			args = [atm.getAtomID(),
					densMet,
					dens,
					densMet,
					meanDens]

			if printText:
				print '{}: D{}: {} --> local D{}: {}'.format(*args)
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
							   resType   = '',
							   densMet   = 'loss',
							   normType  = 'Standard',
							   dataset   = 0,
							   printText = True):

		# for specified residue type, calculate the average 
		# metric for each atom type within the residue

		atms = self.getAtom(restype = resType)
		if len(atms) == 0: 
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
			if printText:
				print '{} --> {}'.format(k,mean)
			meanDic[k] = mean
		return meanDic


















	################################################################################
	# SOME METHODS WHICH DEPEND ON TOM'S BDAMAGE PYTHON SCRIPT

	def parseBDamageFile(self,
						 BDamageFile = 'untitled.txt',
						 atomtype    = 'OH',
						 restype     = 'TYR',
						 Bmetric     = 'Bdamage'):

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
						Bfac = l.split()[-12]
						if Bfac[0:4] == '1.00':
							Bfac = Bfac[4:]
						if Bmetric == 'Bdamage':
							BdamDic[atomID] = Bdam
						elif Bmetric == 'Bfactor':
							BdamDic[atomID] = Bfac
						print '{} --> {}'.format(atomID,BdamDic[atomID] )
		return BdamDic

	def getBdamageStats(self,
						BDamageFile = '',
						atomtype    = 'OH',
						restype     = 'TYR'):

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

	def getBdamageChange(self,
						 BDamageFile1  = 'untitled1.txt',
						 BDamageFile2  = 'untitled2.txt',
						 dataset       = 0,
						 atomtype      = 'OH',
						 restype       = 'TYR',
						 Bmetric       = 'Bdamage',
						 percentChange = True,
						 saveFig       = True,
						 outputDir     = './'):

		# parse two Bdamage output files and calculate 
		# change in Bdamage for each atom

		output1 = self.parseBDamageFile(BDamageFile = BDamageFile1,
										atomtype    = atomtype,
										restype     = restype,
										Bmetric     = Bmetric)

		output2 = self.parseBDamageFile(BDamageFile = BDamageFile2,
										atomtype    = atomtype,
										restype     = restype,
										Bmetric     = Bmetric)

		if not percentChange:
			BdamChange = [round(float(output2[k]) - float(output1[k]),2) for k in output1.keys()]
		else:
			BdamChange = [round((float(output2[k]) - float(output1[k]))/float(output1[k]),2) for k in output1.keys()]

		atomLbls   = output1.keys()
		BdamChange, atomLbls = (list(t) for t in zip(*sorted(zip(BdamChange, atomLbls))))

		dic = {'Atom'     : ['-'.join(lbl.split('-')[0:2]) for lbl in atomLbls],
			   'B-change' : BdamChange}

		dFrame = DataFrame(dic)

		figName = '{}Change_{}-{}_Dset-{}'.format(Bmetric,restype,atomtype,dataset)
		yLabel  = 'Change in {}'.format(Bmetric)
		if percentChange:
			figName = 'Relative-' + figName
			yLabel  = 'Relative '+ yLabel

		self.makeBarPlot(plotData  = dFrame,
						 xName     = "Atom",
						 yName     = "B-change",
						 yLabel    = yLabel,
						 figName   = figName,
						 saveFig   = saveFig,
						 outputDir = outputDir,
						 axesFont  = 14)
		
		outDic = {lbl:val for lbl,val in zip(atomLbls,BdamChange)}

		return outDic

	def compareDensToBdamChange(self,
								BDamageFile1  = 'untitled1.txt',
								BDamageFile2  = 'untitled2.txt',
								BdamChange    = True,
								resType       = 'TYR',
								atomType      = 'OH',
								Bmetric       = 'Bdamage',
								percentChange = True,
								densMet       = 'loss',
								normType      = 'Standard',
								fileType      = '.svg',
								saveFig       = True,
								outputDir     = './',
								dataset       = 0):

		# compare calculated density metric values to 
		# calculated Bdamage values

		# get required atoms in structure
		atms = self.getAtom(restype  = resType,
							atomtype = atomType)

		# retrieve Bdamage values from Bdamage run log file
		if not BdamChange:
			BdamDic = self.parseBDamageFile(BDamageFile = BDamageFile2,
											atomtype    = atomType,
											restype     = resType,
											Bmetric     = Bmetric)
			xLabel = Bmetric
		else:
			BdamDic = self.getBdamageChange(BDamageFile1  = BDamageFile1,
									 	    BDamageFile2  = BDamageFile2,
									 	    atomtype      = atomType,
									 	    restype       = resType,
									 	    Bmetric       = Bmetric,
									 	    percentChange = percentChange,
									 	    saveFig       = saveFig,
									 	    dataset       = dataset,
									 	    outputDir     = outputDir)
			if not percentChange:
				xLabel = '{} Change'.format(Bmetric)
			else:
				xLabel = '{} Relative Change'.format(Bmetric)

		# get density metric values and Bdamage values
		plotData = {'x' : [],
					'y' : []}

		for atm in atms:
			Bdam = BdamDic[atm.getAtomID()]
			dens = atm.densMetric[densMet][normType]['values'][dataset]
			plotData['x'].append(float(Bdam))
			plotData['y'].append(dens)

		# plot the relationship between Bdamage and specified density metric
		self.plotScatterPlot(xData       = plotData['x'],
						     yData       = plotData['y'],
						     xLabel      = xLabel,
							 ylabel      = '{} D{}'.format(normType,densMet),
							 lineBestFit = True,
							 figtitle    = '{}-{}; Dataset {}'.format(resType,atomType,dataset),
							 fileType    = fileType,
							 saveTitle   = 'D{}Scatter_{}-{}_{}VsDensMetric_Dset-{}'.format(densMet,resType,atomType,xLabel.replace(' ',''),dataset),
							 outputDir   = outputDir)

























