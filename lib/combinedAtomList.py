from combinedAtom import combinedAtom
from CalphaWeight import CalphaWeight
from progbar import progress
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import operator
import os
from confidenceIntervalCalculator import mean_confidence_interval
from findMetricChange import findBchange
import scipy.stats
from halfDoseCalc import halfDoseApprox
from PDBFileManipulation import writePDBline_DamSite

class combinedAtomList(object):
	# class for list of atom objects defined by combinedAtom class

	def __init__(self,datasetList,numLigRegDatasets,doseList,initialPDBList,outputDir,partialDatasets,seriesName):
		self.datasetList 		= datasetList
		self.numLigRegDatasets 	= numLigRegDatasets 	# number of datasets to perform linear regression over
		self.doseList 			= doseList 				# list of (DWD)doses for each dataset
		self.initialPDBList		= initialPDBList 		# list of atom objects for initial (lowest dose) pdb file
		# if number of datasets to perform lin reg not stated then set to total number of datasets in series
		self.outputDir			= outputDir
		if self.numLigRegDatasets == 0:
			self.numLigRegDatasets = len(self.datasetList[0].mindensity)
		self.partialDatasets	= partialDatasets # Boolian, whether atoms only in subset of datasets are included
		self.seriesName			= seriesName

	def getMultiDoseAtomList(self):
		# this function inputs a list of lists of PDB atom objects (see StructurePDB class)
		# and formats as an object of the class 'multiPDB'. It is a variant of the function 
		# above which can also cope with structures containing different numbers of atoms 
		# (say if solvent molecules/ligands are included in a subset of the structures). 
		# In this case, the smallest common substructure between all structures will be used

		# first check that each PDBarray contains the same number of atoms (consistency check)
		if len(self.datasetList) > 1:
			print 'Multiple datasets detected...'
			for dataset in self.datasetList:
				if len(dataset) != len(self.datasetList[0]):
					print 'Not all PDB structures have same number of atoms!\n'+\
						  'Will only include atoms common to ALL structures...'
					break
		elif len(self.datasetList) == 1:
			print 'Single dataset detected...'

		PDBdoses = []
		numNotFoundAtoms = 0
		print '------------------------------------------------'
		print 'Locating common atoms to ALL datasets...:'
		singDimAttrs = ('atomnum','residuenum','atomtype','basetype',
					    'chaintype','X_coord','Y_coord','Z_coord')
		multiDimAttrs = ('Bfactor','Occupancy','meandensity','maxdensity','mindensity',
						 'mediandensity','numvoxels','stddensity','min90tile','max90tile',
						 'min95tile','max95tile','rsddensity','rangedensity')
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
				print 'Atom not found in all datasets!'
				print 'Atom: {}'.format(atomIndentifier)
				print '---> not including atom in atom list...'
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

				if atm_counter != len(self.datasetList):
					print 'Atom not found in all datasets!'
					print 'Atom: {}'.format(atomIndentifier)
					print '---> including partial information in atom list...'
					numNotFoundAtoms += 1
				PDBdoses.append(newatom)

		print '\n------------------------------------------------'   
		if self.partialDatasets is True:
			print 'Number of atoms not in all datasets: {}'.format(numNotFoundAtoms)
		else:
			print 'Number of atoms removed since not in all datasets: {}'.format(numNotFoundAtoms)
		print '---> Finished!'
		self.atomList = PDBdoses

	def findMetricName(self,metricName):
		# a conversion between the naming convention used by singlePDB class for each individual
		# dataset and that to be used by dictionary of density metrics for combinedAtom class
		conversions = (['Bfactor','Bfactor'],['Occupancy','occupancy'],
					   ['meandensity','mean'],['maxdensity','gain'],
					   ['mindensity','loss'],['mediandensity','median'],
					   ['numvoxels','num-voxels'],['stddensity','standard-deviation'],
	    			   ['min90tile','90tile-loss'],['max90tile','90tile-gain'],
	    			   ['min95tile','95tile-loss'],['max95tile','95tile-gain'],
	    			   ['rsddensity','relative-std'],['rangedensity','range'])

		for conversion in conversions:
			if metricName == conversion[0]:
				return conversion[1]

	def getNumAtoms(self):
		return len(self.atomList)

	def getAverageMetricVals(self,densMet,normType):
		# get structure-wide average of selected density metric
		densList = [atom.densMetric[densMet][normType]['values'] for atom in self.atomList]
		densMean = np.nanmean(densList,0)
		densStd = np.nanstd(densList,0)
		return densMean,densStd

	def getAtom(self,chain,restype,resnum,atomtype):
		# get atom(s) matching specified description if found
		foundAtoms = []
		for atom in self.atomList:
			if atom.chaintype == chain or chain == '':
				if atom.residuenum == resnum or resnum == '':
					if atom.basetype == restype or restype == '':
						if atom.atomtype == atomtype or atomtype == '':
							foundAtoms.append(atom)
		if len(foundAtoms) != 0:
			print 'Found {} atom(s) matching description'.format(len(foundAtoms))
			return foundAtoms
		print 'Atom matching description not found'
		return False

	def getDensMetrics(self):
		# get a list of all density metrics and normalisations that currently have been defined
		# with a set of values
		currentMetrics = []
		for metric in self.atomList[0].densMetric.keys():
			for normType in self.atomList[0].densMetric[metric].keys():
				currentMetrics.append([metric,normType])
		return currentMetrics

	def calcBfactorChange(self):
		# calculate Bfactor change between initial and each later dataset
		findBchange(self.initialPDBList,self.atomList,'Bfactor')

	def calcAdditionalMetrics(self,metric,normType,newMetric):
		# calculate the Calpha weights for each dataset (see CalphaWeight class for details)
		# for metric "metric" (loss, gain, mean etc.)
		# 'newMetric' takes values ('Calpha','netChange','linreg','subtract1')
		options = ['Calpha','netChange','linreg','subtract1','average']
		if newMetric == 'Calpha':
			print 'Calculating Calpha weights at each dataset...'
			CAweights = CalphaWeight(self.atomList)
			CAweights.calculateWeights(metric)

		# loop over all atoms in list and calculate additional metrics for each atom in atomList
		counter = 0
		numAtoms = self.getNumAtoms()
		for atom in self.atomList:
			counter += 1
			progress(counter, numAtoms, suffix='') # unessential loading bar add-in
			if newMetric == 'Calpha':
				atom.CalphaWeightedDensChange(CAweights,metric)
			elif newMetric == 'linreg':
				atom.calcLinReg(self.numLigRegDatasets,'Standard',metric)
			elif newMetric == 'netChange':
				atom.calcNetChangeMetric('Standard')
			elif newMetric == 'subtract1':
				atom.calcFirstDatasetSubtractedMetric('Standard',metric)
			elif newMetric == 'average':
				atom.calcAvMetric(normType,metric)
			else:
				print 'new metric type not recognised.. choose from: {}'.format(options)
				return

	def calcVectorWeightedMetric(self,metric,normType,vector):
		# for a specified metric calculate new vector-weighted values, where the metric over a series of doses
		# is multiplied by a vector of per-dataset scalings
		for atom in self.atomList:
			atom.calcVectorWeightedMetric(metric,normType,vector)

	def calcVectorSubtractedMetric(self,metric,normType,vector):
		# for a specified metric calculate new vector-weighted values, where the metric over a series of doses
		# is subtracted from a vector of per-dataset scalings
		for atom in self.atomList:
			atom.calcVectorSubtractedMetric(metric,normType,vector)	

	def calcStandardisedMetrics(self,metric):
		# standardise distribution of a given metric to have (mean=0,std=1) for each 
		# dataset within a damage series
		data  = [atm.densMetric[metric]['Standard']['values'] for atm in self.atomList]
		meand = np.nanmean(data,0)
		stdd  = np.nanstd(data,0)
		for atm in self.atomList:
			values = atm.densMetric[metric]['Standard']['values']
			atm.getDensMetricInfo(metric,'Standardised',(values-meand)/stdd)

	def writeMetric2File(self,where,groupBy,metric,normType):
		# write all metric values to a .csv file to location 'where'
		# 'groupBy' takes values 'none','residue','atomtype'

		csvName = '{}{}-{}.csv'.format(where,metric,normType)
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
				TypeDic = self.groupByAtmType(0)
			elif groupBy == 'residue':
				TypeDic = self.groupByResType(0)
			statDic = {k:[] for k in TypeDic.keys()}

			for d in range(metricLength):
				if groupBy == 'atomtype':
					a,b = self.getPerAtmtypeStats(metric,normType,d,'mean',1)
				elif groupBy == 'residue':
					a,b = self.getPerResidueStats(metric,normType,d,'mean',1)
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
				csvfile.write(','.join(map(str,atom.densMetric[metric][normType]['values'])))
				csvfile.write('\n')
		csvfile.close()

	def findMetricRatio(self,metric,normType,resiType,atomType1,atomType2):
		# find ratio between metric for two atom types within a single residue type
		title = '\n{} D{} comparison between {}-{} and {}-{} atoms'.format(normType,metric,resiType,
																  		  atomType1,resiType,atomType2)
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

	def findMetricRatioKeyResidues(self,metric,normType,rType,errorBars,pairs,title):
		# run the above findMetricRatio method for known susceptible residues and write to file
		# 'rType' takes values ('ratio','distance')
		# 'pairs' is list of residues and atoms of form [['GLU','CD','CG'],['GLU','CD','OE1'],..] etc
		fileout = open('{}{}_{}-D{}-{}.txt'.format(self.outputDir,title,normType,metric,rType),'w')

		foundPairs = {}
		for pair in pairs:
			found = self.getAtom('',pair[0],'','')
			if found == False: continue # don't proceed if no atoms of a residue type
			output = self.findMetricRatio(metric,normType,*pair)
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
				plt.errorbar(xData,yData,yerr=yError,fmt='-o',capthick=2,label=key)
		plt.xlabel('Dataset #', fontsize=18)
		plt.ylabel('D{} {}'.format(metric,rType), fontsize=18)
		# place legend outside to right of plot
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=18)
		fig.suptitle('{} D{} {}'.format(normType,metric,rType),fontsize=24)
		fig.savefig('{}{}-D{}-{}-{}_{}.png'.format(self.outputDir,normType,metric,rType,title))

	def calcMetricDiffFromStructureMean(self,metric,normType,diff):
		# for a specified metric calculate difference between metric for atom and from mean of structure
		# 'diff' takes three forms: (a) 'difference', (b) 'ratio', (c) 'num-stds'
		av,std = self.getAverageMetricVals(metric,normType)
		for atom in self.atomList:
			if diff == 'difference':
				atom.calcDiffFromMeanMetric(metric,normType,av)
			elif diff == 'ratio':
				atom.calcRatioToMeanMetric(metric,normType,av)
			elif diff == 'std-devs':
				atom.calcNumStdFromMeanMetric(metric,normType,av,std)

	def numAtmsWithMetricAboveLevel(self,dataset,metric,normType,threshold,printStr):
		# determine the number of atoms with metric value above a specified threshold
		# 'threshold' determined as number of std dev away from structure-wide average for metric
		# 'printStr' (Boolian) is True if string printed to command line
		self.calcMetricDiffFromStructureMean(metric,normType,'std-devs')
		count = 0
		for atom in self.atomList:
			if atom.densMetric[metric]['num stds']['values'][dataset] < -threshold:
				count += 1
		if printStr is True:
			print 'For dataset {}:'.format(dataset)
			print '# atoms with {} D{} metric above {} std dev of structure-wide mean is {}'.format(normType,metric,threshold,count)
		return count

	def findMetricRatioKeyResidues_scatterplot(self,metric,normType,rType,pairs,title):
		# run the above findMetricRatio method for known susceptible residues and write to file
		# 'rType' takes values ('ratio','distance')
		# 'pairs' is list of residues and atoms of form [['GLU','CD','CG'],['GLU','CD','OE1'],..] etc
		from matplotlib import cm
		foundPairs = {}
		for pair in pairs:
			found = self.getAtom('',pair[0],'','')
			if found == False: continue # don't proceed if no atoms of a residue type
			output = self.findMetricRatio(metric,normType,*pair)
			foundPairs['-'.join(pair)] = output[0]

		# plot the relationship between xData and yData
		sns.set_palette("deep", desat=.6)
		sns.set_context(rc={"figure.figsize":(10, 10)})
		fig = plt.figure()
		ax = plt.subplot(111)

		total = 0
		for key in foundPairs.keys(): total += len(foundPairs[key][rType])
		xDataTotal = np.linspace(0,1,total)

		i,j = 0,-1
		cmap = cm.Set1
		nsteps = len(foundPairs.keys())

		colors = ['#737474','#409cd6','#58bb6b','#faa71a','#ff6a6a']

		DataSave = {}
		for key in foundPairs.keys():	
			j += 1
			yData = np.nanmean(foundPairs[key][rType],1)
			xData = xDataTotal[i:i+len(yData)]
			i += len(yData)
			plt.scatter(xData,yData,marker='o',s=100,c=colors[j],edgecolors='#FFFFFF',label=key)
			DataSave[key] = {'y':yData,'x':xData}

		# place legend outside to right of plot
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=18)

		plt.xlabel(' ', fontsize=18)
		plt.ylabel('D{} {}'.format(metric,rType), fontsize=18)
		figtitle = '{} D{} {}'.format(normType,metric,rType)
		fig.suptitle(figtitle,fontsize=24)
		saveTitle = figtitle.replace(' ','_')
		fname = lambda x: saveTitle.strip('.png')+'_{}.png'.format(x)
		i = 0
		while os.path.isfile(fname(i)): i += 1 
		fig.savefig(fname(i))
		return DataSave

	def getAtomListWithoutPartialAtoms(self,dataset):
		# return new list of atoms with partial atoms removed for that dataset
		newList = []
		for atom in self.atomList:
			if dataset in atom.getPresentDatasets():
				newList.append(atom)
		return newList

	def getTopNAtoms(self,metric,normType,dataset,n):
		# for a given metric type, determine top n damage sites. 'n' is integer or 'all'
		atomList = self.getAtomListWithoutPartialAtoms(dataset) # remove atoms not in present dataset
		atomList.sort(key=lambda x: x.densMetric[metric][normType]['values'][dataset])

		# for difference maps, more negative metrics indicate more damage, however if Calpha normalisation 
		# used, then the highest positive metric values are the most damaged, this must be accounted for here:
		if normType == 'Calpha normalised':
			atomList.reverse()

		if n != 'all': 
			return atomList[:int(n)]
		else: 
			return atomList

	def getTopNAtomsPDBfile(self,metric,normType,dataset,n,pdbFile):
		# for a given metric type, determine top n damage sites. 'n' is integer (not 'all' here)
		# method writes locations of top damage sites to pdb file
		pdbIn = open(pdbFile,'r')
		pdbOutName = '{}top{}-{}D{}sites-dset{}.pdb'.format(self.outputDir,n,normType,metric,dataset)
		pdbOut = open(pdbOutName,'w')
		for l in pdbIn.readlines():
			if l.split()[0] in ('CRYST1','SCALE1','SCALE2','SCALE3'):
				pdbOut.write(l)
		# next get top n metric sites and write to file
		topAtms = self.getTopNAtoms(metric,normType,dataset,n)
		count = 0
		for atm in topAtms:
			count += 1
			# generate pdb line - note: absolute value of metric taken since Bfactor>0 (for displaying in pymol etc.)
			l = writePDBline_DamSite(atm,np.abs(atm.densMetric[metric][normType]['values'][dataset]),count)
			pdbOut.write(l+'\n')
		pdbOut.write('END')
		pdbIn.close()
		pdbOut.close()
		return pdbOutName

	def getTopNAtomsString(self,metric,normType,dataset,n):
		# return the top n atoms within structure in terms of metric 'metric'
		# 'n' takes values 'all' or an integer
		topN = self.getTopNAtoms(metric,normType,dataset,n)
		atomInfoList = []
		for atom in topN:
			data = atom.getAtomID()+'\t\t'
			for met in ('loss','mean','gain','Bfactor'):
				data += str(round(atom.densMetric[met][normType]['values'][dataset],2))+'\t\t'
			atomInfoList.append(data)
		stringOut = 'Atom-ID\t\t\tDloss\t\tDmean\t\tDgain\t\tBfactor\n'
		stringOut += '\n'.join(atomInfoList)	
		return stringOut

	def breakdownTopNatomsBy(self,metric,normType,dataset,n,by):
		# 'by' is grouping, i.e ['chaintype'] or ['basetype','atomtype']
		topN = self.getTopNAtoms(metric,normType,dataset,n)
		countDic,totalDic = {},{}
		for atom in topN:
			atts = [getattr(atom,att) for att in by]
			atomID = '-'.join(atts)
			if atomID not in countDic.keys():
				countDic[atomID] = 1
				totalDic[atomID] = self.getNumAtomsOfType(by,atts)
			else:
				countDic[atomID] += 1
		headerString  = 'Partition of atom {} D{} values by {}.\n'.format(normType,metric,' & '.join(by))
		fracOfTotal   = round((float(n)/self.getNumAtoms())*100,2) # n as fraction of total atoms
		headerString += 'For dataset {}, top {} ({}%) atoms included.\n'.format(dataset,n,fracOfTotal)
		outputString  = 'Type\t\tFraction\tPercent(%)'

		keys = countDic.keys()
		keys.sort()
		for key in keys:
			outputString += '\n{}\t\t{}/{}\t\t{}%'.format(key,countDic[key],totalDic[key],
											  				   round((float(countDic[key])/totalDic[key])*100,2))
		return headerString,outputString

	def getNumAtomsOfType(self,atts,vals):
		# get number of atoms with attribute 'att' at specified value 'val'
		# set 'att' and 'val' as lists of strings, i.e. att = ['chaintype','restype'], val = ['A','GLU']
		i = 0
		for atom in self.atomList:
			if vals == [getattr(atom,att) for att in atts]: i += 1
		return i

	def getPerAtmtypeStats(self,metric,normType,dataset,sortby,n):
		# for a given metric type, determine per-atom-type statistics on the distribution of damage
		atmtypeDict 		= self.groupByAtmType(dataset) # group atoms by atom type
		statsDic 			= self.getStats(metric,normType,dataset,atmtypeDict)
		statsString		 	= self.reportStats(statsDic,'Type',sortby,n,normType)
		return statsString,statsDic

	def getPerResidueStats(self,metric,normType,dataset,sortby,n):
		# for a given metric type, determine per-residue statistics on the distribution of damage
		resDict 			= self.groupByResType(dataset) # group atoms by residue type
		statsDic 			= self.getStats(metric,normType,dataset,resDict)
		statsString		 	= self.reportStats(statsDic,'Residue',sortby,n,normType)
		return statsString,statsDic

	def getPerChainStats(self,metric,normType,dataset,sortby,n):
		# for a given metric type, determine per-chain statistics on the distribution of damage
		chainDict 			= self.groupByChainType(dataset) # group atoms by chain type
		statsDic 			= self.getStats(metric,normType,dataset,chainDict)
		statsString		 	= self.reportStats(statsDic,'Chain',sortby,n,normType)
		return statsString,statsDic

	def getStats(self,metric,normType,dataset,dic):
		# output distribution stats for each element in dictionary dic as a new dictionary
		getStatsPerKey = {}
		for key in dic.keys():
			metricList = [atom.densMetric[metric][normType]['values'][dataset] for atom in dic[key]]
			getStatsPerKey[key] = self.getStatsForList(metricList)
		return getStatsPerKey

	def getStatsForList(self,metricList):
		# calculate measures of the distribution of values in a list metricList
		statsDic = {}
		statsDic['#atoms']	 	= len(metricList)
		statsDic['mean'] 		= np.mean(metricList)
		statsDic['std'] 		= np.std(metricList)
		statsDic['skew'] 		= scipy.stats.skew(metricList, axis=0, bias=True)
		statsDic['ratio'] 		= self.calcNetRatio(metricList)
		statsDic['outliers'] 	= self.calcNumOutliers(metricList)
		statsDic['returnOrder'] = ['mean','std','#atoms','outliers','skew','ratio']
		return statsDic

	def reportStats(self,stats,name,sortby,n,normType,format='txt'):
		# format stats defined in 'stats' dictionary into an output string. 'sortby' specifies how 
		# the output should be ranked. 'n' is number of entries to print (specify an integer or 'all')
		# format takes values in ('txt','html')
		sortbyDict   = self.getStatsForList([1,1])
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

		# for difference map analysis, more negative metrics indicate more damage, however if Calpha normalisation 
		# used, then the highest positive metric values are the most damaged, this must be accounted for here:
		if normType == 'Calpha normalised': 
			sortedList1.reverse()

		if n != 'all':
			stringOut =	'\n'.join(sortedList1[:n])
		else:
			stringOut =	'\n'.join(sortedList1)	
		return headerString+stringOut

	def calcDiscreteDistMode(self,metricList):
		# calculate mode of discrete histogram with 100 bins
		h = np.histogram(metricList,bins=100)
		mode_ind = np.argmax(h[0])
		mode_val = np.mean(h[1][mode_ind:mode_ind+2])
		return mode_val

	def calcNetRatio(self,metricList):
		# calculate net ratio of distn values either side of distn mode
		switchVal = self.calcDiscreteDistMode(metricList)
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

	def calcNumOutliers(self,metricList):
		# by assuming a symmetric distn around the mode, flag specific atoms that fall outside
		# this range
		distMode = self.calcDiscreteDistMode(metricList)
		distMax  = np.max(metricList)
		sudoMin  = distMode-(distMax - distMode)
		count 	 = 0
		for val in metricList:
			if val < sudoMin: count+=1
		return count

	def groupByAtmType(self,dataset):
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

	def groupByResType(self,dataset):
		# group atoms in a dictionary by residue type
		resDict = {}
		for atom in self.atomList:
			if dataset not in atom.getPresentDatasets(): continue
			if atom.basetype not in resDict.keys():
				resDict[atom.basetype] = [atom]
			else:
				resDict[atom.basetype].append(atom)
		return resDict

	def groupByChainType(self,dataset):
		# group atoms in a dictionary by chain type
		chainDict = {}
		for atom in self.atomList:
			if dataset not in atom.getPresentDatasets(): continue
			if atom.chaintype not in chainDict.keys():
				chainDict[atom.chaintype] = [atom]
			else:
				chainDict[atom.chaintype].append(atom)
		return chainDict

	def getNumDatasets(self):
		return self.atomList[0].getNumDatasets()

	def graphMetricDistn(self,metric='loss',normType='Standard',valType='average',
						 plotType='both',resiType='all',save=True,fileType='png'):
		# histogram/kde plot of density metric per atom
		# plotType is 'histogram' or 'kde'
		# resiType is 'all' or list of residue types
		# save is Boolian to save or not
		# 'valType' takes values 'average' (compute average),'all' (plot all datasets),
		# or 'n' for int (not str) n (plot dataset n)

		# attempt to find atoms of type 'resiType' and flag if not all found
		if resiType != 'all':
			count = 0
			for res in resiType:
				atms = self.getAtom('',res,'','')
				if atms is False: 
					count += 1
			if count > 0:
				print 'Warning: not all selected residue types found!'

		if plotType not in ('hist','kde','both'): return 'Unknown plotting type selected.. cannot plot..'

		if normType == 'Calpha normalised':
			self.calcAdditionalMetrics(metric,normType,'Calpha')
		if self.checkMetricPresent(self.atomList[0],metric,normType) is False: return # check metric valid

		# sns.set_palette("deep", desat=.6)
		sns.set_style("whitegrid")
		sns.set_context(rc={"figure.figsize": (10, 6)})
		fig = plt.figure()
		ax = plt.subplot(111)

		plotData = {}
		if valType == 'all':
			numDsets = self.getNumDatasets()
			for j,(i, color) in enumerate(zip(range(numDsets), sns.color_palette('Blues', numDsets))):
				if resiType == 'all':
						datax = [atm.densMetric[metric][normType]['values'][i] for atm in self.atomList]
						self.plotHist(plotType=plotType,datax=datax,
									  lbl='Dataset {}'.format(i),color=color)
						plotData[resiType] = datax
				else:
					for res in resiType:
						datax = [atm.densMetric[metric][normType]['values'][i] for atm in self.atomList if atm.basetype == res]
						if len(datax) > 0:
							self.plotHist(plotType=plotType,datax=datax,
										  lbl='Dataset {},{}'.format(i,res),color=color)
						plotData[res] = datax					
		else:
			self.calcAdditionalMetrics(metric,normType,'average')
			if resiType == 'all':
				if valType == 'average':
					datax = [atm.densMetric[metric][normType]['average'] for atm in self.atomList]
					lbl = 'average'
				else:
					datax = [atm.densMetric[metric][normType]['values'][valType] for atm in self.atomList]
					lbl = 'Dataset '+valType
				self.plotHist(plotType=plotType,datax=datax,lbl=lbl,color='r')
				plotData[resiType] = datax
			else:
				for j,(res, color) in enumerate(zip(resiType, sns.color_palette('hls', len(resiType)))):
					if valType == 'average':
						datax = [atm.densMetric[metric][normType]['average'] for atm in self.atomList if atm.basetype == res]
						lbl = 'average, {}'.format(res)
					else:
						datax = [atm.densMetric[metric][normType]['values'][valType] for atm in self.atomList if atm.basetype == res]
						lbl = 'Dataset {}, {}'.format(valType,res)
					if len(datax) > 0:
						self.plotHist(plotType=plotType,datax=datax,lbl=lbl,color=color)
					plotData[res] = datax	

		plt.legend()
		plt.xlabel('{} D{} per atom'.format(normType,metric),fontsize=18)
		plt.ylabel('Frequency',fontsize=18)
		plt.title('{} D{} per atom, residues: {}'.format(normType,metric,resiType))

		if not save: 
			plt.show()
		else:
			saveName = '{}{}_{}D{}-{}.{}'.format(self.outputDir,''.join(resiType),
												 normType.replace(" ",""),metric,
												 plotType,fileType)
			if valType != 'all':
				saveName = saveName.replace('.'+fileType,'-{}.'.format(valType,fileType))
			fig.savefig(saveName)
		return plotData

	def plotHist(self,plotType='both',nBins=300,datax=[],lbl='',color='b'):
		# plot histogram or kde plot for datax and give current label
		# 'nBins' is number of bins (only used if plotType is 'hist' or 'both')
		if plotType == 'hist':
			plt.hist(datax, nBins, histtype="stepfilled", alpha=.7,label=lbl,color=color)
		elif plotType == 'kde':
			sns.kdeplot(np.array(datax), shade=True,label=lbl,color=color)
		elif plotType == 'both':
			sns.distplot(np.array(datax),label=lbl,color=color)

	def graphMetric(self,*params):
		# produce a graph of selected metric against dataset number for a specified atom
		# get equivalent atoms of specified type (command line input to specify)
		from matplotlib import cm

		if len(params) == 0:
			densMet 	= raw_input("Density metric type: ")
			normType  	= raw_input("Normalisation type: ")
			atomType 	= raw_input("Atom type: ")
			baseType 	= raw_input("Residue/nucleotide type: ")
			residueNum 	= raw_input("Residue number: ")
			errorBars 	= raw_input("Include error bars? Options are 'NONE','RESNUM','ATOMTYPE':")
			save 		= raw_input("Save plot as .png file? (y/n)")
		else:
			if len(params) != 7:
				print 'error in number of plot parameters specified'
				return
			else:
				[densMet,normType,atomType,baseType,residueNum,errorBars,save] = params

		# check valid inputs 
		errorOptions = ('NONE','RESNUM','ATOMTYPE')
		if errorBars not in errorOptions: 
			print 'invalid error parameter..\n choose from: {}'.format(errorOptions)
			return

		# find atoms of interest 
		foundAtoms = []
		for atom in self.atomList:
			if atom.basetype == baseType:
				if residueNum == 'all' or str(atom.residuenum) == residueNum:
					if atomType == 'all' or atom.atomtype == atomType:
						foundAtoms.append(atom)

		# check at least one atom has been found
		if len(foundAtoms) == 0:
			print 'No atoms found..'
			return
			
		# define x range here (number in damage series)
		x = range(self.getNumDatasets())
		sns.set(style="white", context="talk")
		f, axes = plt.subplots(1, 1, figsize=(12, 12), sharex=True)

		# determine y values here dependent on density metric type specified 
		if errorBars == 'NONE':
			colormap = plt.cm.nipy_spectral
			plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(foundAtoms))])			
			for atom in foundAtoms:
				if self.checkMetricPresent(atom,densMet,normType) is False: return # check metric valid
				y = atom.densMetric[densMet][normType]['values']
				plt.plot(x,y,label=atom.getAtomID())

		else: # error bars will be plotted
			errorBarKey = {'ATOMTYPE':'atomtype','RESNUM':'residuenum'}
			errorAttr 	= errorBarKey[errorBars]
			yDict 		= {}
			for atom in foundAtoms:
				if self.checkMetricPresent(atom,densMet,normType) is False: return # check metric valid
				vals = atom.densMetric[densMet][normType]['values']
				if getattr(atom,errorAttr) not in yDict.keys():
					yDict[getattr(atom,errorAttr)] = [vals]
				else:
					yDict[getattr(atom,errorAttr)].append(vals)
			for key in yDict.keys():
				plt.errorbar(x,np.mean(yDict[key],0),yerr=np.std(yDict[key],0),fmt='-o',capthick=2,label=key)

		plt.xlabel('Dataset', fontsize=18)
		plt.ylabel('{} D{}'.format(normType,densMet), fontsize=18)
		plt.legend()
		f.suptitle('{} D{}: {} {} {}'.format(normType,densMet,baseType,residueNum,atomType),fontsize=24)
		if save != 'y':
			plt.show()
		else:
			f.savefig('{}{}_D{}_{}_{}_{}.png'.format(self.outputDir,normType,densMet,baseType,residueNum,atomType))

	def plotSusceptibleAtoms(self,densMet='loss',normType='Standard',
							 errorbars=True,save='y',susAtms=[],fileType='png'):
		from matplotlib import cm
		if susAtms == []: 
			susAtms = [['GLU','CD'],['ASP','CG'],['TYR','OH'],['CYS','SG'],['MET','SD']]
		sDic = {}
		for s in susAtms:
			findAtms = self.getAtom('',s[0],'',s[1])
			if findAtms is False: continue
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
				plt.errorbar(x,yData,yerr=yErr,fmt='-o',capthick=2,label=key,color=cmap(i/float(nsteps)))

		av,std = self.getAverageMetricVals(densMet,normType) # get structure-wide average and std dev
		plt.errorbar(x,av,yerr=std,fmt=':o',capthick=2,color='k') 

		plt.xlabel('Dataset', fontsize=18)
		plt.ylabel('{} D{}'.format(normType,densMet), fontsize=18)
		plt.legend()
		f.suptitle('{} D{}: susceptible residues'.format(normType,densMet),fontsize=24)
		if save != 'y':
			plt.show()
		else:
			i = 0
			fname = lambda x: '{}{}_D{}_susceptResis_{}.{}'.format(self.outputDir,normType,densMet,x,fileType)
			while os.path.isfile(fname(i)):
				i += 1 
			f.savefig(fname(i))

	def checkMetricPresent(self,atom,densMet,normType):
		# method to check whether a method is valid for selected atom object
		try:
			atom.densMetric[densMet][normType]['values']
		except NameError:
			print 'Unexpected density metric name..'
			return False
		return True

	def calcHalfDose(self,chain,restype,resnum,atomtype,densMet,normType,shift,offset):
		# call the external halfDoseCalc class to calculate a rough half dose decay value for a specified atom
		# 'shift' (Boolian) determines whether half dose is dose for density to reach half of initial density (=False)
		# or half dose is dose for density to reach av(initial density,end density limit (=True)
		# 'offset' (Boolian) determines whether exponential decay function is allowed a non-zero end density value
		atom = self.getAtom(chain,restype,resnum,atomtype)
		atom = atom[0] # convert from list
		if len(self.doseList) != self.getNumDatasets():
			return 'need to specify doses list as class attribute before this can be calculated'
		halfDoseApprox(atom,self.atomList,True,self.doseList,False,0.5,densMet,normType,self.outputDir,shift,offset)

	def calcHalfDoseForAtomtype(self,restype,atomtype,densMet,normType,shift,offset,n):
		# call the above calcHalfDose method for instances of of restype and atomtype
		atoms = self.getAtom('',restype,'',atomtype)
		for atom in atoms:
			self.calcHalfDose(atom.chaintype,restype,atom.residuenum,atomtype,densMet,normType,shift,offset)
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
		self.plotScatterPlot(xData,yData,'Initial density','Half-dose',
							'D{}: {}-{}'.format(densMet,restype,atomtype),
							'D{}Scatter_{}-{}_thres{}.png'.format(densMet,restype,atomtype,n),False,False,'')

	def compareAtomHalfDoses(self,restype,atomtype1,atomtype2,densMet,normType,n):
		# after above method calcHalfDoseForAtomtype is run for two atom types atomtype1 and atomtype2
		# can compare half-doses between two atoms within same side-chain
		atoms1 = self.getAtom('',restype,'',atomtype1)
		atoms2 = self.getAtom('',restype,'',atomtype2)

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
		self.plotScatterPlot(xData,yData,'{}-{} Half-dose'.format(restype,atomtype1),
							 '{}-{} Half-dose'.format(restype,atomtype2),
							 'D{} half-dose for {} atoms'.format(densMet,restype),
							 'D{}_halfDose{}_{}-{}.png'.format(densMet,restype,atomtype1,atomtype2),False,False,'')

	# temporary method to do a specific plot:
	def HalfDoseDistn(self,metric,normType,plotType,resiType,atomtype,save,nBins,colour):
		# histogram/kde plot of density metric per atom
		# plotType is 'histogram' or 'kde'
		# save is Boolian to save or not
		if plotType not in ('hist','kde','both'): return 'Unknown plotting type selected.. cannot plot..'
		if self.checkMetricPresent(self.atomList[0],metric,normType) is False: return # check metric valid

		sns.set_palette("deep", desat=.6)
		sns.set_context(rc={"figure.figsize": (10, 6)})
		fig = plt.figure()

		atoms = self.getAtom('',resiType,'',atomtype)
		datax = [atm.densMetric[metric][normType]['Half-dose']['Half-dose'] for atm in atoms if atm.densMetric[metric][normType]['Half-dose']['Half-dose'] < 100 ]
		self.plotHist(plotType=plotType,nBins=nBins,datax=datax,
					  lbl='',colour=colour)
		plt.xlabel('Half-dose (MGy)')
		plt.ylabel('Normed frequency')
		plt.title('Half-dose distribution: {}-{}'.format(resiType,atomtype))
		if not save: 
			plt.show()
		else:
			fig.savefig('{}halfDoseDistn-{}{}.svg'.format(self.outputDir,resiType,atomtype))

	def compareSensAtoms(self,densMet,normType):
		rSquaredDic = {}
		numPairsDic = {}
		dataDic 	= {}
		for i in [['GLU','CD','CG'],['GLU','CD','OE1'],['ASP','CG','CB'],['ASP','CG','OD1'],['TYR','OH','CZ'],['TYR','CZ','CE2']]:
			rSquared,numPairs,data = self.compareMetricsBetweenAtoms(i[0],i[1],i[2],densMet,normType,'av')
			if rSquared != False: 
				rSquaredDic['-'.join(i)] = rSquared
				numPairsDic['-'.join(i)] = numPairs
				dataDic['_'.join(i)] 	 = data
		return rSquaredDic,numPairsDic,dataDic

	def compareMetricsBetweenAtoms(self,restype,atomtype1,atomtype2,densMet,normType,dataset):
		# for two atom types atomtype1 and atomtype2, compare density metric 
		# between two atoms within same side-chain
		atoms1 = self.getAtom('',restype,'',atomtype1)
		atoms2 = self.getAtom('',restype,'',atomtype2)
		if atoms1 is False or atoms2 is False: return False,0,{'x':[],'y':[]}

		if dataset == 'av': self.calcAdditionalMetrics(densMet,normType,'average')

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
		rSquared = self.plotScatterPlot(xData,yData,'{}-{} D{}'.format(restype,atomtype1,densMet),
							 '{}-{} D{}'.format(restype,atomtype2,densMet),
							 'D{} for {} atoms'.format(densMet,restype),
							 'D{}_{}-{}-{}_scatterplot_{}.svg'.format(densMet,restype,atomtype1,atomtype2,dataset),
							 True,False,'')
		data = {'x':xData,'y':yData}
		return rSquared,len(xData),data

	def compareMetrics(self,restype,atomtype,metric1,metric2,normType1,normType2,dSet):
		# for all atoms of type restype and atomtype plot scatter plot comparing metric1 and 2 for chosen dataset dSet
		# set restype = '' and atomtype = '' to include all atoms
		atoms = self.getAtom('',restype,'',atomtype)
		xData,yData = [],[]
		for atom in atoms:
			xData.append(atom.densMetric[metric1][normType1]['values'][dSet])
			yData.append(atom.densMetric[metric2][normType2]['values'][dSet])
		rSquared = self.plotScatterPlot(xData,yData,'{} D{}'.format(normType1,metric1),
										 '{} D{}'.format(normType2,metric2),
										 '{} D{} vs {} D{} for res:{} atoms:{}'.format(normType1,metric1,normType2,metric2,restype,atomtype),
										 '{}D{}-vs-{}D{}-{}-{}-{}.png'.format(normType1,metric1,normType2,metric2,restype,atomtype,dSet),
										 True,False,'')

	def plotScatterPlot(self,xData,yData,xLabel,ylabel,figtitle,saveTitle,lineBestFit,yequalsxLine,colors):
		# plot the relationship between xData and yData
		# if 'lineBestFit' is True then linear line of best fit calculated and plotted
		# if 'yequalsxLine' is True then line y=x plotted
		# sns.set_palette("deep", desat=.6)
		# sns.set_context(rc={"figure.figsize":(10, 10)})
		sns.set_context("talk")

		# sns.set(font_scale=1)
		fig = plt.figure(figsize=(10, 10))
		ax = plt.subplot(111)

		# choose default color if none specified
		if colors == '': colors = '#50527A'

		plt.scatter(xData,yData,marker='o',s=100,c=colors,edgecolors='#FFFFFF',cmap='copper')

		if lineBestFit is True: # plot line of best if specified
			slope, intercept, r_value, p_value, std_err = stats.linregress(xData,yData)
			x = np.linspace(min(xData), max(xData),10)
			y = slope*np.array(x) + np.array(intercept)
			plt.plot(x, y, '-',color='#d11141',linewidth=3)
			figtitle += ', R^2:{}'.format(round(r_value**2,2))

		if yequalsxLine is True:
			x = np.linspace(min(xData), max(xData),10)
			y = x
			plt.plot(x, y, ':',color='#9890af',linewidth=3)

		plt.xlabel(xLabel, fontsize=24)
		plt.ylabel(ylabel, fontsize=24)
		fig.suptitle(figtitle,fontsize=24)
		sns.despine()
		if saveTitle.split('.')[-1] == 'png':
			fname = lambda x: self.outputDir+saveTitle.strip('.png')+'_{}.png'.format(x)
		elif saveTitle.split('.')[-1] == 'svg':
			fname = lambda x: self.outputDir+saveTitle.strip('.svg')+'_{}.svg'.format(x)
		i = 0
		while os.path.isfile(fname(i)): i += 1 
		fig.savefig(fname(i))
		if lineBestFit is True: return r_value**2

	def damageVdistPlot(self,type1,type2,dataset,densMet,normType):
		# for atom of type1 determine a correlation between the closest distance to atom type2
		# and the damage to atom type1
		# type1,type2 of form [['GLU','CD'],['ASP','CG'],..] depending on how many atom types required

		atoms1,atoms2 = [],[]
		for ind in type1: atoms1 += self.getAtom('',ind[0],'',ind[1])
		for ind in type2: atoms2 += self.getAtom('',ind[0],'',ind[1])

		densVals,minDists = [],[]
		for atom1 in atoms1:
			densVals.append(atom1.densMetric[densMet][normType]['values'][dataset])
			distList = []
			for atom2 in atoms2:
				distList.append(self.getDistanceBetweenAtoms(atom1,atom2))
			minDists.append(min(distList))

		yLabel = '{} D{}'.format(normType,densMet)
		xLabel = '{} {} distance'.format(type1,type2)
		figTitle = '{} D{} against {} distance'.format(type1,densMet,type2)
		saveTitle = figTitle.strip()+'.png'
		self.plotScatterPlot(minDists,densVals,xLabel,yLabel,figTitle,saveTitle,False,False,'')

	def getDistanceBetweenAtoms(self,atom1,atom2):
		# get distance between two atoms
		xyz1 = np.array([atom1.X_coord,atom1.Y_coord,atom1.Z_coord])
		xyz2 = np.array([atom2.X_coord,atom2.Y_coord,atom2.Z_coord])
		return np.linalg.norm(xyz1-xyz2)

	def getAtomsWithinDist(self,atom,distLim):
		# find all atoms within 'distLim' angstroms of selected atomtype
		nearAtmDic = {'atoms':[],'distances':[]}
		for otheratom in self.atomList:
			if (otheratom.getAtomID()).split('-')[:-1] == (atom.getAtomID()).split('-')[:-1]: 
				continue
			dist = self.getDistanceBetweenAtoms(atom,otheratom)
			if dist < distLim:
				nearAtmDic['atoms'].append(otheratom)
				nearAtmDic['distances'].append(dist)
		print '{} in total within {} Angstrom of {}'.format(len(nearAtmDic['atoms']),distLim,atom.getAtomID())
		return nearAtmDic

	def densMetSurroundAtmsCorrel(self,restype,atomtype,distLim,normType,densMet,crystConts,pdbName,symmetrygroup):
		# determine whether there is a correlation between density metric and types of surrounding atoms 
		# if 'crystConts' is True then crystal contacts will also be located her

		if crystConts is True:
			from NCONTjob import NCONTjob

			ncont = NCONTjob(pdbName,self.outputDir,symmetrygroup,self.seriesName)
			success = ncont.run()
			if success is False: return
			logFile = self.outputDir+'/'+ncont.outputLogfile

		atoms = self.getAtom('',restype,'',atomtype)
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
							contactAtom = self.getAtom(chain,l[39:42],resNum,l[47:50])
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
		solvAccDic,plotData = self.compareSolvAccessWithAndWithoutGluAspGroups(pdbName,symmetrygroup,restype,atomtype,groupA,densMet,normType)
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
		self.plotScatterPlot(densValsA,densValsAContacts,'{}-{} D{}'.format(restype,atomtype,densMet),'Carboxyl-contact D{}'.format(densMet),
							'{}-{} vs carboxyl-contact D{}'.format(restype,atomtype,densMet),
							'D{}Scatter_{}-{}_carboxylContacts.png'.format(densMet,restype,atomtype),False,False,scatterColor)

		return strA+'\n'+strB,densValsA,densValsAContacts,scatterColor

	def susceptAtmComparisonBarplot(self,metric,normType,dataset,set,box):
		# produce a barplot to compare the damage metric of susceptible atom types at a given dataset

		if set == 1:
			keyAtomTypes = [['GLU','CD'],['ASP','CG'],['TYR','OH'],['CYS','SG'],['MET','SD'],
							['LYS','NZ'],['ARG','NE'],['ASN','OD1'],['GLN','OE1'],['ALA','CB']]
		if set == 2:
			keyAtomTypes = [['GLU','OE1'],['GLU','OE2'],['ASP','OD1'],['ASP','OD2'],['TYR','OH'],
							['ASN','OD1'],['GLN','OE1'],['SER','OG'],['THR','OG1'],['','O']]
	
		plotData = {'x':[],'y':[],'y-error':[]}
		for atm in keyAtomTypes:
			foundAtms = self.getAtom('',atm[0],'',atm[1])
			if foundAtms == False: continue 
			for foundAtm in foundAtms:
				atmId 	= '-'.join(atm)
				plotData['x'].append(atmId)
				plotData['y'].append(foundAtm.densMetric[metric][normType]['values'][dataset])

		fig = plt.figure()
		sns.set_style("whitegrid")
		if box == 'Box':
			ax = sns.boxplot(x="x", y="y", data=plotData)
		elif box == 'Bar':
			ax = sns.barplot(x="x", y="y", data=plotData)
		else:
			print 'Unknown graph type - specify "Box" or "Bar"'
			return

		plt.xlabel('Atom type', fontsize=18)
		plt.ylabel('{} D{}'.format(normType,metric), fontsize=18)

		fig.savefig('{}SusceptAtms{}plot-{}-D{}_{}-set{}.png'.format(self.outputDir,box,normType,metric,dataset,set))

	def compareSolvAccessWithAndWithoutGluAspGroups(self,pdbName,symmetrygroup,resType,atomType,specialAtoms,densMet,normType):
		# determine the change in solvent accessibility for 'resType'-'atomType' atoms before & after Glu/Asp decarboxylation events
		# 'specialAtoms' is an optional parameter, which specifies a subset of atoms objects which should be coloured differently
		# on the resulting scatter plot below - the idea is to run .densMetSurroundAtmsCorrel() above first to get the list of atoms
		# exhibiting contacts to Glu/Asp CO2 groups and then colour these differently below

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
		TyrOHatms = self.getAtom('',resType,'',atomType)

		solvAccDic = {}
		for atm in TyrOHatms:
			solvAccDic[atm.getAtomID()] = []

			# get other Tyr-ring atoms for each Tyr-OH atom
			TyrRingAtms = []
			if resType == 'TYR' and atomType == 'OH':
				for atmtype in ['CG','CD1','CD2','CE1','CE2','CZ']:
					TyrRingAtms.append(self.getAtom(atm.chaintype,resType,atm.residuenum,atmtype)[0])

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
		self.plotScatterPlot(plotData['x'],plotData['y'],
							 'Solvent Accessibility',
							 '{}-{} D{}'.format(resType,atomType,densMet),'',
							 'D{}Scatter_{}-{}_vsSolvAccChange.png'.format(densMet,resType,atomType),
							 False,False,plotData['colours'])

		return solvAccDic,plotData

	def calculateLocalDloss(self,resType,atomType,distance,densMet,normType,weighted):
		# for a given residue, calculate the average density metric (for a specified metric)
		# within a given region of space around that particular residue.
		# if 'weighted' is True then a distance-weighted mean is instead calculated

		plotData = {'x':[],'y':[]}
		# find specified atoms within structure
		atms = self.getAtom('',resType,'',atomType)
		if atms is False: return plotData # don't continue if no atoms exist
		self.calcAdditionalMetrics(densMet,normType,'average')

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
		self.plotScatterPlot(plotData['x'],plotData['y'],
							 '{}-{} D{}'.format(resType,atomType,densMet),
							 'local environment D{}'.format(densMet),
							 'Average D{} within {} Angstrom of {}-{}'.format(densMet,distance,resType,atomType),
							 'D{}Scatter_{}-{}_localEnvironmentComparison.png'.format(densMet,resType,atomType),
							  False,False,'')
		return plotData

	def getAvMetricPerAtmInRes(self,resType,densMet,normType,dataset):
		# for specified residue type, calculate the average metric for each atom type within the residue
		atms = self.getAtom('',resType,'','')
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
	def parseBDamageFile(self,BDamageFile,atomtype,restype):
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

	def getBdamageStats(self,BDamageFile):
		# retrieve BDamage values for each atom type and calculate mean & std dev
		csvFile = open('{}bDamageStats.csv'.format(self.outputDir),'w')
		csvFile.write('id,mean,std dev\n')
		atmTypeDic = self.groupByAtmType(0)
		for aType in atmTypeDic.keys():
			aInfo 	= aType.split('-')
			BdamDic = self.parseBDamageFile(BDamageFile,aInfo[1],aInfo[0])
			mean 	= np.mean(map(float,BdamDic.values()))
			std 	= np.std(map(float,BdamDic.values()))
			csvFile.write('{},{},{}\n'.format(aType,mean,std))
		csvFile.close()

	def compareDensityToBdamage(self,BDamageFile,resType,atomType,densMet,normType):
		# compare calculated density metric values to calculated Bdamage values

		# get required atoms in structure
		atms = self.getAtom('',resType,'',atomType)
		self.calcAdditionalMetrics(densMet,normType,'average')

		# retrieve Bdamage values from Bdamage run log file
		BdamDic = self.parseBDamageFile(BDamageFile,atomType,resType)

		# get density metric values and Bdamage values
		plotData = {'x':[],'y':[]}
		for atm in atms:
			Bdam = BdamDic[atm.getAtomID()]
			dens = atm.densMetric[densMet][normType]['average']
			plotData['x'].append(Bdam)
			plotData['y'].append(dens)

		# plot the relationship between Bdamage and specified density metric
		self.plotScatterPlot(plotData['x'],plotData['y'],'Bdamage',
							 '{} D{}'.format(normType,densMet),
							 '{}-{}'.format(resType,atomType),
							 'D{}Scatter_{}-{}_BdamageVsDensMetric.png'.format(densMet,resType,atomType),
							 False,False,'')

























