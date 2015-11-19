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
	                print 'Not all PDB structures have same number of atoms!'\
	                ' Will only include atoms common to ALL structures...'
	    elif len(self.datasetList) == 1:
	        print 'Single dataset detected...'

	    PDBdoses = []
	    numNotFoundAtoms = 0
	   
	    print '------------------------------------------------'
	    print 'Locating common atoms to ALL datasets...:'
	    singDimAttrs = ('atomnum','residuenum','atomtype','basetype',
	                    'chaintype','X_coord','Y_coord','Z_coord')
	    multiDimAttrs = ('Bfactor','Occupancy','meandensity','maxdensity',
	                     'mindensity','mediandensity','Bfactorchange',
	                     'numvoxels','stddensity','min90tile','max90tile',
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
				y = combinedAtom()
				for attr in singDimAttrs+multiDimAttrs:
					if attr[0] != '_':
						if len(atomDict[attr]) == 1:
							setattr(y, attr,atomDict[attr])
				else:
					y.densMetric[attr]['Standard']['values'] = atomDict[attr]
				PDBdoses.append(y)

	    print '\n------------------------------------------------'    
	    print 'Number of atoms removed since not in all datasets: {}'.format(numNotFoundAtoms)
	    print '---> Finished!'
	    self.atomList = PDBdoses

	def calcBfactorChange(self):
		# calculate Bfactor change between initial and each later dataset
		findBchange(self.initialPDBList,self.atomList,'Bfactor')

	def calcAdditionalMetrics(self):
		# calculate the Calpha weights for each dataset (see CalphaWeight class for details)
		print 'Calculating Calpha weights at each dataset...'
		CAweights = CalphaWeight(self.atomList)
		CAweights.calculateWeights()

		# loop over all atoms in list and calculate additional metrics for each atom in atomList
		counter = 0
		num_atoms = len(self.atomList)
		for atom in self.atomList:
			counter += 1
			progress(counter, num_atoms, suffix='') # unessential loading bar add-in
			atom.CalphaWeightedDensChange(CAweights)
			atom.calculateAdditionalMetrics()
			atom.calculateLinReg(self.numLigRegDatasets,'Standard')

	def graphMetric(self):
		# produce a graph of selected metric against dataset number for a specified atom
		# get equivalent atoms of specified type (command line input to specify)
		self.densityMetric = raw_input("Density metric type: ")
		self.normalise  = raw_input("Normalisation type: ")
		self.atomType 	= raw_input("Atom type: ")
		self.baseType 	= raw_input("Residue/nucleotide type: ")
		self.residueNum = int(raw_input("Residue number: "))

		# find atoms of interest 
		foundAtoms = []
		for atom in self.atomList:
			if (atom.atomtype 	== self.atomType and
				atom.basetype 	== self.baseType and
				atom.residuenum == self.residueNum):
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
		plt.ylabel('{} D{}'.format(self.normalise,self.densityMetric), fontsize=18)
		f.suptitle('{} D{}: {} {} {}'.format(self.normalise,self.densityMetric,
											 self.baseType,self.residueNum,
											 self.atomType),fontsize=24)
		
		plt.show()

