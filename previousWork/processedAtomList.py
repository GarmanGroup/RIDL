# script to convert multiPDB list of atom objects to list of processedAtom class objects
# (which have additional attributes and methods over those of multiPDB)
from processedAtom import processedAtom
from CalphaWeight import CalphaWeight
from progbar import progress
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import operator
import os
from confidenceIntervalCalculator import mean_confidence_interval

class processedAtomList(object):
	# class for list of atom objects defined by processedAtom class

	def __init__(self,unprocessedAtomList,numDatasets,doseList):
		self.unprocessedAtomList = unprocessedAtomList
		self.numDatasets = numDatasets # number of datasets to perform linear regression over
		self.doseList = doseList # list of (DWD)doses for each dataset

		# if number of datasets to perform lin reg not stated then default
		# to number of datasets present within unprocessedAtomList
		if self.numDatasets == 0:
			self.numDatasets = len(self.unprocessedAtomList[0].mindensity)

	def processAtomList(self):
		# process the input multiPDB list of atom objects to create new
		# list of atom objects from processedAtom class
		processedList = []

		# calculate the Calpha weights for each dataset (see CalphaWeight class for details)
		print 'Calculating Calpha weights at each dataset...'
		CAweights = CalphaWeight(self.unprocessedAtomList)
		CAweights.calculateWeights()

		# loop over all atoms in list and determine new atom info (defined by processedAtom class)
		print 'Creating new list of atom objects within class processedAtom...'
		counter = 0
		num_atoms = len(self.unprocessedAtomList)
		for oldAtom in self.unprocessedAtomList:
			counter += 1
			progress(counter, num_atoms, suffix='') #unessential loading bar add-in
			newAtom = processedAtom()
			newAtom.cloneInfo(oldAtom)
			newAtom.CalphaWeightedDensChange(CAweights)
			newAtom.calculateAdditionalMetrics()
			newAtom.calculateLinReg(self.numDatasets,'Standard')
			# newAtom.calculateLinReg(self.numDatasets,'Calpha normalised')

			processedList.append(newAtom)
		self.processedAtomList = processedList

	def getEquivalentAtoms(self,auto):
		# user inputs atom indentifier on command line and method finds all 
		# equivalent atoms in different TRAP RNA chains

		# only run this function if atom list has been processed already
		try:
			self.processedAtomList
		except AttributeError:
			print 'Need to process atom list before this function can be run'
			return

		if auto == False:
			self.atomType 	= raw_input("Atom type: ")
			self.baseType 	= raw_input("Residue/nucleotide type: ")
			self.residueNum = int(raw_input("Residue number: "))

		# find the correct atom specified by above properties
		equivAtoms = []

		# if protein atoms to be found
		if self.baseType not in ('A','C','G','T','U'):
			for atom in self.processedAtomList:
				if (atom.atomtype == self.atomType and atom.basetype == self.baseType 
					and atom.residuenum == self.residueNum):
					equivAtoms.append(atom)

		# if RNA/DNA atoms to be found
		else:
			for atom in self.processedAtomList:
				if (atom.atomtype == self.atomType and atom.basetype == self.baseType):
					equivAtoms.append(atom)

		# if an atom with specified properties not found, don't proceed
		if len(equivAtoms) == 0:
			print 'Atom not found within structure.. make sure atom properties correctly specified'
			return
		else:
			self.equivAtoms = equivAtoms

	def getMetricValForEquivAtoms(self,densMet,normType):
		# this method uses the above method to find a set of equivalent atoms within TRAP
		# then writes the density values specified by densMet and normType to a 
		# csv file for each dose
		self.atomType 	= raw_input("Atom type: ")
		self.baseType 	= raw_input("Residue/nucleotide type: ")
		self.residueNum = int(raw_input("Residue number: "))
		fileName = '{}{}{}_{}{}.csv'.format(self.baseType,self.residueNum,self.atomType,densMet,normType)
		csvOut = open(fileName,'w')
		self.getEquivalentAtoms(True)
		csvOut.write(','.join(range(self.equivAtoms[0].densMetric[densMet][normType]['values'])))
		csvOut.write(',Slope\n')
		for atom in self.equivAtoms:
			densVals = atom.densMetric[densMet][normType]['values']
			densSlope = atom.densMetric[densMet][normType]['lin reg']
			csvOut.write(','.join(densVals))
			csvOut.write(',{}'.format(densSlope))
			csvOut.write('\n')
		csvOut.close()

	def calculateSigDiffsForList(self,densMet,normType,valType):
		# this method calculated the Hotelling T-squared test statistic for all
		# atoms within the TRAP-RNA structure. valType takes values of 'values' or 'slope'

		TsquareDict = {}
		RNAdistDict = {}
		for atom in self.processedAtomList:

			if atom.boundOrUnbound() in ('rna','unbound protein'):
				continue # only consider protein chains here. Find unique atoms from bound protein only

			dictKey = '{} {} {}'.format(atom.basetype,atom.residuenum,atom.atomtype)
			try:
				TsquareDict[dictKey]
				continue # if an equivalent atom has already been found, skip
			except KeyError:
				pass
				
			self.atomType 	= atom.atomtype
			self.baseType 	= atom.basetype
			self.residueNum = atom.residuenum
			self.getEquivalentAtoms(True)

			if len(self.equivAtoms) != 22:
				continue # don't include residues only in a subset of protein chains

			TsquareDict[dictKey] = {}
			if valType == 'values':
				F,p_value,reject = self.hotellingTsquareTest(densMet,normType) # run Hotelling's T squared test to distinguish between bound and unbound TRAP rings
				TsquareDict[dictKey]['p value'] = p_value 
			elif valType == 'lin reg':
				tStat,p_value = self.studentTsquareTest(densMet,normType) # calculate Student T-test if only lin reg slopes tested
				TsquareDict[dictKey]['p value'] = p_value 

			# determine whether RNA has stabilising effect (bound lin reg slope lower than unbound)
			TsquareDict[dictKey]['stabilising?'] = self.findIfStabilising(densMet,normType)

			# calculate the min distance of the protein atom from the RNA
			if atom.boundOrUnbound() == 'bound protein':
				distList = []
				for otheratom in self.processedAtomList:
					if otheratom.boundOrUnbound() == 'rna':
						atomPosition = np.array([atom.X_coord,atom.Y_coord,atom.Z_coord])
						otheratomPosition = np.array([otheratom.X_coord,otheratom.Y_coord,otheratom.Z_coord])
						distList.append(np.linalg.norm(atomPosition-otheratomPosition))
				RNAdistDict[dictKey] = min(distList)

		# sort dictionary by p value 
		sorted_HotTsquareDict = sorted(TsquareDict.items(), key=operator.itemgetter(1))
		return TsquareDict,RNAdistDict

	def findIfStabilising(self,densMet,normType):
		# determine whether RNA is having a stabilising effect on the gradients of density metrics or not
		# Uses linear regression slopes for bound and unbound groups

		boundMetVals,unboundMetVals = [],[]
		for atom in self.equivAtoms:
			if atom.boundOrUnbound() == 'unbound protein':
				unboundMetVals.append(atom.densMetric[densMet][normType]['lin reg']['slope'])
			elif atom.boundOrUnbound() == 'bound protein':
				boundMetVals.append(atom.densMetric[densMet][normType]['lin reg']['slope'])

		# determine here whether RNA binding is reducing the density metric gradient, or increasing it
		if np.abs(np.mean(unboundMetVals)) > np.abs(np.mean(boundMetVals)):
			return True
		else:
			return False	

	def graphMetric(self):
		# produce a graph of selected metric against dataset number for a specified atom

		# get equivalent atoms of specified type (command line input to specify)
		self.getEquivalentAtoms(False)

		self.densityMetric = raw_input("Density metric type: ")
		self.normalise  = raw_input("Normalisation type: ")

		# define x range here (damage set numbers)
		x = range(2,len(self.equivAtoms[0].mindensity)+2)

		sns.set(style="white", context="talk")
		f, axes = plt.subplots(1, 1, figsize=(12, 12), sharex=True)

		# determine y values here dependent on density metric type specified 
		for foundAtom in self.equivAtoms:
			y = foundAtom.densMetric[self.densityMetric][self.normalise]['values']
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

	def densMetricErrorbarGraphs(self,auto,where,metricTypes,confInt):
		# function to plot density change as function of dataset number
		# for a specific atom in the structure, with the mean value over 
		# all protein chains plotted, along with error bars for the 22 
		# equivalent atoms present. specify auto=False to specify atom type
		# on the command line
		# 'where' specifies where to plot, if doesn't exist, makes directory in
		# current directory
		# metricTypes takes values 1 or 2 if auto == True
		# If 'confInt' is True then error bars are 95% confidence intervals, 
		# otherwise, 1 SD used at each dose

		# get equivalent atoms of specified type (command line input to specify)
		self.getEquivalentAtoms(auto)

		# determine whether dealing with protein or RNA atoms 
		if self.equivAtoms[0].boundOrUnbound() in ('unbound protein','bound protein'):
			protein = True
		else:
			protein = False

		sns.set(style="white", context="talk")
		f = plt.figure(figsize=(16, 8))

		# define x range here (damage set numbers or doses if specified)
		if self.doseList == []:
			x = range(2,len(self.equivAtoms[0].meandensity)+2)[0:10]
			x_label = 'Damage set'
		else:
			x = self.doseList
			x_label = "Dose (MGy)"

		# Determine density metric set to plot here
		# Currently two distinct options given below
		densMets1 = ['loss','net','mean','gain']
		densMets2 = ['loss','net','mean','gain','bfactor','bdamage']
		densMets3 = ['|loss|','loss','net','mean','gain','bfactor']
		densMets4 = ['max-simple','median-simple','gain','median','mean','loss']

		if auto == False:
			print 'Which metrics would you like to plot...'
			userInput = raw_input("1 or 2?: ")
			if userInput == str(1):
				densMets = densMets1
				normTypes = ['Standard','Calpha normalised']
			else:
				densMets = densMets2
				normTypes = ['Standard']
		if auto == True:
			if metricTypes == 1:
				densMets = densMets1
				normTypes = ['Standard','Calpha normalised']
			elif metricTypes == 2:
				densMets = densMets2
				normTypes = ['Standard']
			elif metricTypes == 3:
				densMets = densMets3
				normTypes = ['Standard']
			elif metricTypes == 4:
				densMets = densMets4
				normTypes = ['Standard']
		i = 0
		HotellingTsquareDict = {}
		for densMet in densMets:
			for normType in normTypes:
				if densMet in ('mean','gain') and normType in ('Calpha normalised'):
					continue
				i+=1
				yValue = {}
				# for protein atoms, group by bound and unbound chains
				if protein == True:
					for boundType in ('unbound','bound'):
						yValue[boundType] = {}
						for valType in ('mean','std','95ConfInt'):
							yValue[boundType][valType] = []
						for j in range(0,len(x)):
							yValue[boundType]['mean'].append(np.mean([atom.densMetric[densMet][normType]['values'][j] for atom in self.equivAtoms if atom.boundOrUnbound() == '{} protein'.format(boundType)]))
							yValue[boundType]['std'].append(np.std([atom.densMetric[densMet][normType]['values'][j] for atom in self.equivAtoms if atom.boundOrUnbound() == '{} protein'.format(boundType)]))
							yValue[boundType]['95ConfInt'].append(mean_confidence_interval([atom.densMetric[densMet][normType]['values'][j] for atom in self.equivAtoms if atom.boundOrUnbound() == '{} protein'.format(boundType)]))

				# for RNA atoms, just create a list of density values
				if protein == False:
					yValue['RNA 1'] = {}
					yValue['RNA 2'] = {}
					for valType in ('mean','std'):
						yValue['RNA 1'][valType] = []
						yValue['RNA 2'][valType] = []
					for j in range(0,len(x)):	
						yValue['RNA 1']['mean'].append(np.mean([atom.densMetric[densMet][normType]['values'][j] for atom in self.equivAtoms if atom.residuenum%5 == self.residueNum]))
						yValue['RNA 1']['std'].append(np.std([atom.densMetric[densMet][normType]['values'][j] for atom in self.equivAtoms if atom.residuenum%5 == self.residueNum]))
						yValue['RNA 2']['mean'].append(np.mean([atom.densMetric[densMet][normType]['values'][j] for atom in self.equivAtoms if atom.residuenum%5 != self.residueNum]))
						yValue['RNA 2']['std'].append(np.std([atom.densMetric[densMet][normType]['values'][j] for atom in self.equivAtoms if atom.residuenum%5 != self.residueNum]))

				ax = plt.subplot(2,3,i)
				ax.set_xlim([0, 29])
				if protein == True:
					if confInt == True:
						plt.errorbar(x,yValue['unbound']['mean'],yerr=yValue['unbound']['95ConfInt'], fmt='-o',capthick=2,color='#99ccff',label='Non-bound')
						plt.errorbar(x,yValue['bound']['mean'],yerr=yValue['bound']['95ConfInt'],fmt='-o',capthick=2,color='#f47835',label='Bound')
					else:
						plt.errorbar(x,yValue['unbound']['mean'],yerr=yValue['unbound']['std'], fmt='-o',capthick=2,color='#99ccff',label='Non-bound')
						plt.errorbar(x,yValue['bound']['mean'],yerr=yValue['bound']['std'],fmt='-o',capthick=2,color='#f47835',label='Bound')
				else:
					try:
						plt.errorbar(x,yValue['RNA 1']['mean'],yerr=yValue['RNA 1']['std'], fmt='-o',capthick=2,color='r',label='G1')
						plt.errorbar(x,yValue['RNA 2']['mean'],yerr=yValue['RNA 2']['std'], fmt='-o',capthick=2,color='g',label='G3')
					except KeyError:
						plt.errorbar(x,yValue['RNA 1']['mean'],yerr=yValue['RNA 1']['std'], fmt='-o',capthick=2,color='r')

				# ax.set_xlim([1, 11])
				ax.legend(loc='best')
				plt.xlabel(x_label)
				if normType == 'Calpha normalised':
					plt.ylabel('Normalised D{} change'.format(densMet))
				else:
					plt.ylabel('{} D{} change'.format(normType,densMet))

				# # perform Hotelling T-squared test if protein atoms
				# if protein == True:
				# 	keyVal = '{} D{}'.format(normType,densMet)
				# 	HotellingTsquareDict[keyVal] = {}
				# 	F,p_value,reject = self.hotellingTsquareTest(densMet,normType) # run Hotelling's T squared test to distinguish between bound and unbound TRAP rings
				# 	HotellingTsquareDict[keyVal]['F value'] = F 
				# 	HotellingTsquareDict[keyVal]['p value'] = p_value 
				# 	HotellingTsquareDict[keyVal]['reject?'] = reject 

		plt.subplots_adjust(top=0.90)
		f.subplots_adjust(hspace=0.4)
		f.subplots_adjust(wspace=0.5)

		f.suptitle('damage metrics vs damage set: {} {} {}'.format(atom.basetype,atom.residuenum,atom.atomtype),fontsize=20)

		# check if directory exists to save graphs to and make if not:
		if not os.path.exists(where):
			os.makedirs(where)

    	# save graphs to directly specified by 'where'
		if confInt == True:
			f.savefig('{}/6DamageSubplots_{}_{}_{}_95confIntErrorbars.png'.format(where,atom.basetype,atom.residuenum,atom.atomtype))
		else:
			f.savefig('{}/6DamageSubplots_{}_{}_{}_SDerrorbars.png'.format(where,atom.basetype,atom.residuenum,atom.atomtype))

		# # print the results of Hotelling T-square test
		# if protein == True:
		# 	for key in HotellingTsquareDict.keys():
		# 		print '{}: F value: {} p value: {} --> reject?: {}'.format(key,HotellingTsquareDict[key]['F value'],HotellingTsquareDict[key]['p value'],HotellingTsquareDict[key]['reject?'])

	def studentTsquareTest(self,densMet,normType):
		# calculates Student T-test if only lin reg slopes tested between unbound and bound atoms

		boundMetVals,unboundMetVals = [],[]
		for atom in self.equivAtoms:
			if atom.boundOrUnbound() == 'unbound protein':
				unboundMetVals.append(atom.densMetric[densMet][normType]['lin reg']['slope'])
			elif atom.boundOrUnbound() == 'bound protein':
				boundMetVals.append(atom.densMetric[densMet][normType]['lin reg']['slope'])

		tStat,p_value = stats.ttest_ind(unboundMetVals,boundMetVals) 
		return tStat,p_value

	def hotellingTsquareTest(self,densMet,normType):
		# this function is designed to calculate Hotelling's T-square test statistic
		# in order to compare bound vs bound sets of atoms over multiple datasets 
		# simultaneously (each dose level is a different dimension to the data)

		p = len(self.doseList)

		atomDensities_unbound,atomDensities_bound = [],[]
		for atom in self.equivAtoms:
			col = []
			for dens in atom.densMetric[densMet][normType]['values'][0:p]:
				col.append([dens])
			if atom.boundOrUnbound() == 'unbound protein':
				atomDensities_unbound.append(col)
			elif atom.boundOrUnbound() == 'bound protein':
				atomDensities_bound.append(col)

		# calculate means for unbound and bound vectors
		nx,ny = 11,11
		unbound_mean = np.sum(np.array(atomDensities_unbound),axis=0)/nx
		bound_mean = np.sum(np.array(atomDensities_bound),axis=0)/ny

		# calculate the unbiased pooled covariance matrix estimate:
		W = 0
		for i in range(0,nx):
			unboundPart = (np.array(atomDensities_unbound)[i] - unbound_mean)*np.transpose((np.array(atomDensities_unbound)[i] - unbound_mean))
			boundPart = (np.array(atomDensities_bound)[i] - bound_mean)*np.transpose((np.array(atomDensities_bound)[i] - bound_mean))
			combinedParts = unboundPart + boundPart
			W = W + combinedParts
		W = W/(nx+ny-2)

		# now calculate the Hotelling's two-sample T-squared statistic:
		tSquared = (float(nx*ny)/(nx+ny))*np.dot(np.transpose(unbound_mean-bound_mean),np.dot(np.linalg.inv(W),(unbound_mean-bound_mean)))

		# now calculate the related F test statistic from F distribution:
		F = (float(nx+ny-p-1)/((nx+ny-2)*p))*tSquared

		# generate a p-value for the given F test statistic
		alpha = 0.05 
		df1 = p
		df2 = nx+ny-1-p
		p_value = 1-stats.f.cdf(F, df1, df2)
		reject = 'notsureyet'
		if p_value < alpha:
			reject = 'YES'
		else:
			reject = 'NO'
		return F,p_value,reject
				
