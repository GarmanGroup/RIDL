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

class combinedAtomList(object):
	# class for list of atom objects defined by combinedAtom class

	def __init__(self,datasetList,numLigRegDatasets,doseList,initialPDBList):
		self.datasetList 		= datasetList
		self.numLigRegDatasets 	= numLigRegDatasets 	# number of datasets to perform linear regression over
		self.doseList 			= doseList 				# list of (DWD)doses for each dataset
		self.initialPDBList		= initialPDBList 		# list of atom objects for initial (lowest dose) pdb file
		# if number of datasets to perform lin reg not stated then set to total number of datasets in series
		if self.numLigRegDatasets == 0:
			self.numLigRegDatasets = len(self.datasetList[0].mindensity)

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

	        # list of index of same atom in each later dataset
	        indexindataset = []

	        # check whether atom in all datasets:
	        Inds = ('residuenum','atomtype','basetype','chaintype')
	        atomIndentifier = [getattr(atom, attr) for attr in Inds]
	        for dataset in self.datasetList[1:]:
	            k = -1        
	            for otheratom in dataset: 
	                k += 1
	                if atomIndentifier == [getattr(otheratom, att) for att in Inds]:        
	                    atm_counter += 1 
	                    for attr in multiDimAttrs:
	                        atomDict[attr].append(getattr(otheratom, attr))
	                    indexindataset.append(k)
	                    break

	        # remove this located atom from the later dataset now that it
	        # has been located --> to make loop faster 
	        for j in range(1,len(indexindataset)+1):
	            if indexindataset[j-1] != -1:
	                self.datasetList[j].pop(indexindataset[j-1])
	                        
	        if atm_counter != len(self.datasetList):
	            print 'Atom not found in all datasets!'
	            print 'Atom: {} {}{} {}'.format(*atomIndentifier)
	            print '---> not including atom in atom list...'
	            numNotFoundAtoms += 1
	            continue

	        else:                 
				atom = combinedAtom()
				for attr in singDimAttrs+multiDimAttrs:
					if attr[0] != '_':
						if isinstance(atomDict[attr],list) is False:
							setattr(atom,attr,atomDict[attr])
						else:
							metName = self.findMetricName(attr)
							atom.getDensMetricInfo(metName,'Standard',atomDict[attr])
				PDBdoses.append(atom)

	    print '\n------------------------------------------------'    
	    print 'Number of atoms removed since not in all datasets: {}'.format(numNotFoundAtoms)
	    print '---> Finished!'
	    self.atomList = PDBdoses

	def findMetricName(self,metricName):
		# a conversion between the naming convention used by singlePDB class for each individual
		# dataset and that to be used by dictionary of density metrics for combinedAtom class
		conversions = (['Bfactor','Bfactor'],['Occupancy','occupancy'],
					   ['meandensity','mean'],['maxdensity','gain'],
					   ['mindensity','loss'],['mediandensity','median'],
					   ['numvoxels','num voxels'],['stddensity','standard deviation'],
	    			   ['min90tile','90tile loss'],['max90tile','90tile gain'],
	    			   ['min95tile','95tile loss'],['max95tile','95tile gain'],
	    			   ['rsddensity','relative-std'],['rangedensity','range'])

		for conversion in conversions:
			if metricName == conversion[0]:
				return conversion[1]

	def getNumAtoms(self):
		return len(self.atomList)

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

	def calcAdditionalMetrics(self,metric):
		# calculate the Calpha weights for each dataset (see CalphaWeight class for details)
		# for metric "metric" (loss, gain, mean etc.)
		print 'Calculating Calpha weights at each dataset...'
		CAweights = CalphaWeight(self.atomList)
		CAweights.calculateWeights(metric)

		# loop over all atoms in list and calculate additional metrics for each atom in atomList
		counter = 0
		numAtoms = self.getNumAtoms()
		for atom in self.atomList:
			counter += 1
			progress(counter, numAtoms, suffix='') # unessential loading bar add-in
			atom.CalphaWeightedDensChange(CAweights,metric)
			atom.calculateLinReg(self.numLigRegDatasets,'Standard',metric)
			atom.calculateNetChangeMetric('Standard')

	def writeMetric2File(self,where,metric,normType):

		# make a new txt file in location where
		csvfile = open('{}{}-{}.csv'.format(where,metric,normType),'w')
		csvfile.write('metric: {} normalisation: {}\n'.format(metric,normType))
		csvfile.write('atomnum,atominfo,')
		metricLength = self.atomList[0].getNumDatasets(metric)
		for i in range(0,metricLength):
			csvfile.write('{}'.format(str(i+1)))
			if i != metricLength-1: csvfile.write(',')
		csvfile.write('\n')

		# sort atom list by atom number
		self.atomList.sort(key=lambda x: x.atomnum)

		# for each atom in atomList
		for atom in self.atomList:
			csvfile.write('{},'.format(atom.atomnum))
			csvfile.write('{},'.format(atom.getAtomID()))
			for i in range(0,metricLength):
				csvfile.write('{}'.format(atom.densMetric[metric][normType]['values'][i]))
				if i != metricLength-1: csvfile.write(',')
			csvfile.write('\n')
		csvfile.close()

	def getTopDamSites(self,metric,normType,dataset,n):
		# for a given metric type, determine top n damage sites
		self.atomList.sort(key=lambda x: x.densMetric[metric][normType]['values'][dataset])
		i = -1
		while i < n:
			i += 1
			print self.atomList[i].getAtomID()

	def getPerAtmtypeStats(self,metric,normType,dataset,sortby,n):
		# for a given metric type, determine per-atom-type statistics on the distribution of damage
		atmtypeDict = self.groupByAtmType() # group atoms by atom type
		statsOut = self.reportStats(metric,normType,dataset,atmtypeDict,'Atom type',sortby,n)
		return statsOut

	def getPerResidueStats(self,metric,normType,dataset,sortby,n):
		# for a given metric type, determine per-residue statistics on the distribution of damage
		resDict = self.groupByResType() # group atoms by residue type
		statsOut = self.reportStats(metric,normType,dataset,resDict,'Residue',sortby,n)
		return statsOut

	def getPerChainStats(self,metric,normType,dataset,sortby,n):
		# for a given metric type, determine per-chain statistics on the distribution of damage
		chainDict = self.groupByChainType() # group atoms by chain type
		statsOut = self.reportStats(metric,normType,dataset,chainDict,'Chain',sortby,n)
		return statsOut

	def getDistnStats(self,metricList):
		# calculate some measures of the distribution of values in a list metricList
		numAtoms 	= len(metricList)
		mean 		= np.mean(metricList)
		std 		= np.std(metricList)
		skew 		= scipy.stats.skew(metricList, axis=0, bias=True)
		return numAtoms,mean,std,skew

	def reportStats(self,metric,normType,dataset,dic,name,sortby,n):
		# output distribution stats for each element in dictionary. 'sortby' specifies how 
		# the output should be ranked. 'n' is number of entries to print (specify an integer or 'all')
		list1,list2 = [],[]
		for key in dic.keys():
			metricList = [atom.densMetric[metric][normType]['values'][dataset] for atom in dic[key]]
			stats = self.getDistnStats(metricList)
			string = '{} type: {} # atoms: {} mean: {} std: {} skewness: {}'.format(name,key,*stats)
			list1.append(string)
			if sortby == 'mean':
				list2.append(stats[1])
			elif sortby == 'skew':
				list2.append(stats[3])
			elif sortby == 'std':
				list2.append(stats[2])
			else:
				print 'Unexpected ranking specified'
				return
		sortedList1 = [x for (y,x) in sorted(zip(list2,list1))] # sort by chosen 'sortby' parameter
		i = 0
		for string in sortedList1:
			if n != 'all':
				if i > n:
					break
			print string
			i += 1
		if n != 'all':
			stringOut =	'\n'.join(sortedList1[:n])	
		else:
			stringOut =	'\n'.join(sortedList1)	
		return stringOut

	def groupByAtmType(self):
		# group atoms in a dictionary by atom type
		atmtypeDict = {}
		for atom in self.atomList:
			atmtype = atom.basetype+'-'+atom.atomtype
			if atmtype not in atmtypeDict.keys():
				atmtypeDict[atmtype] = [atom]
			else:
				atmtypeDict[atmtype].append(atom)
		return atmtypeDict

	def groupByResType(self):
		# group atoms in a dictionary by residue type
		resDict = {}
		for atom in self.atomList:
			if atom.basetype not in resDict.keys():
				resDict[atom.basetype] = [atom]
			else:
				resDict[atom.basetype].append(atom)
		return resDict

	def groupByChainType(self):
		# group atoms in a dictionary by chain type
		chainDict = {}
		for atom in self.atomList:
			if atom.chaintype not in chainDict.keys():
				chainDict[atom.chaintype] = [atom]
			else:
				chainDict[atom.chaintype].append(atom)
		return chainDict

	def graphMetric(self):
		# produce a graph of selected metric against dataset number for a specified atom
		# get equivalent atoms of specified type (command line input to specify)
		densMet 	= raw_input("Density metric type: ")
		normType  	= raw_input("Normalisation type: ")
		atomType 	= raw_input("Atom type: ")
		baseType 	= raw_input("Residue/nucleotide type: ")
		residueNum 	= int(raw_input("Residue number: "))

		# find atoms of interest 
		foundAtoms = []
		for atom in self.atomList:
			if (atom.atomtype 	== atomType and
				atom.basetype 	== baseType and
				atom.residuenum == residueNum):
				foundAtoms.append(atom)

		# check at least one atom has been found
		if len(foundAtoms) == 0:
			print 'No atoms found..'
			return

		# define x range here (damage set numbers)
		x = range(2,len(foundAtoms[0].densMetric['loss'])+2)

		sns.set(style="white", context="talk")
		f, axes = plt.subplots(1, 1, figsize=(12, 12), sharex=True)

		# determine y values here dependent on density metric type specified 
		for atom in foundAtoms:
			y = atom.densMetric[self.densityMetric][self.normalise]['values']
			try:
				y
			except NameError:
				print 'Unexpected density metric name.. please try again'
				return

			# plot the graph here
			plt.plot(x,y)

		plt.xlabel('Dataset', fontsize=18)
		plt.ylabel('{} D{}'.format(normType,densMet), fontsize=18)
		f.suptitle('{} D{}: {} {} {}'.format(normType,densMet,baseType,residueNum,atomType),fontsize=24)
		
		plt.show()

