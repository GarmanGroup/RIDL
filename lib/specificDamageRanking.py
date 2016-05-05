import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy import stats

from checkSeabornPresent import checkSeabornPresent as checkSns
seabornFound = checkSns()
if seabornFound is True:
	import seaborn as sns

class specificDamageRank(object):

	def __init__(self,
				 residueName    = "",
				 atomType       = "",
				 densMetric     = "",
				 boundOrUnbound = "",
				 damageRank     = 0 ,
				 densSlope      = 0,
				 slopeError     = 0):

		self.residueName 		= residueName
		self.atomType			= atomType
		self.densMetric			= densMetric
		self.boundOrUnbound		= boundOrUnbound # this property is specific to TRAP-RNA binding
		self.damageRank 		= damageRank
		self.densSlope 			= densSlope
		self.slopeError			= slopeError

class specificDamageRanking(object):

	# class to determine the ordering of specific damage within a MX structure
	# using multiple difference maps collected over increasing doses from the 
	# same crystal

	def __init__(self,
				 atomList   = [],
				 densMetric = "Loss",
				 numLines   = 0,
				 normalised = False):

		self.atomList 	= atomList
		self.densMetric = densMetric # the density metric for which the ranking will be based
		self.numLines 	= numLines   # specify the number of lines of output to display to command line
		self.normalised = normalised # (Boolian) specifies whether the 'Standard' or 'Calpha normalised' metric values are used
	
	def normaliseType(self):

		# determine the normalisation type

		if self.normalised == False:
			# no normalisation (default)
			return 'Standard'
		else:
			return 'Calpha normalised'

	def calculateRanks(self):

		# for each atom within the atomList list of atom objects, group into 
		# residue/nucleotide types and determine the average slope for straight
		# lines of best fit for each residue/nucleotide type (for plots of 
		# dataset number vs max density loss)

		# don't run if density metric not specified
		if self.densMetric == "":
			print 'Need to specify density metric before this function can be run'
			return
		
		# determine the normalisation type for density metrics
		normType = self.normaliseType()

		residueDict = {}
		for atom in self.atomList:
			identifier = "{} {} {}".format(atom.basetype,atom.atomtype,atom.boundOrUnbound())
			if identifier not in residueDict.keys():
				residueDict[identifier] = {}
				residueDict[identifier]['slope'] = []
				residueDict[identifier]['std_err'] = []

			residueDict[identifier]['slope'].append(atom.densMetric[self.densMetric.lower()][normType]['lin reg']['slope'])
			residueDict[identifier]['std_err'].append(atom.densMetric[self.densMetric.lower()][normType]['lin reg']['std_err'])
	
		# calculate the average straight line slope for each residue/nucleotide type
		specificDamageRanks = []
		for key in residueDict.keys():
			damageStats = specificDamageRank(key.split()[0],key.split()[1],self.densMetric,key.split()[2])
			damageStats.densSlope = np.mean(residueDict[key]['slope'])
			damageStats.slopeError = np.std(residueDict[key]['slope'])
			specificDamageRanks.append(damageStats)

		# rank residues/nucleotides by slope value 
		if normType == 'Standard' and self.densMetric.lower() in ('loss','mean'):
			specificDamageRanks.sort(key = lambda x: x.densSlope)
		elif normType == 'Standard' and self.densMetric.lower() in ('gain','net'):
			specificDamageRanks.sort(key = lambda x: x.densSlope, reverse = True)
		elif normType == 'Calpha normalised' and self.densMetric.lower() in ('loss','mean','net'):
			specificDamageRanks.sort(key = lambda x: x.densSlope,reverse = True)
		elif normType == 'Calpha normalised' and self.densMetric.lower() in ('gain'):
			specificDamageRanks.sort(key = lambda x: x.densSlope)

		counter = 0
		for residue in specificDamageRanks:
			residue.damageRank = counter
			counter += 1

		self.specificDamageRanks = specificDamageRanks

	def getDamageRanks(self):

		# following on from the above function, output the calculated 
		# ordering of damage within the current structure

		# don't run if damage ranks not yet determined above
		try:
			self.specificDamageRanks
		except AttributeError:
			print 'Need to calculate damage rankings before this function can be run'
			return

		print '--------------------------------------------------------------------'
		print 'Ordering of damage with {} D{} metric as follows:'.format(self.normaliseType(),str(self.densMetric).lower())
		
		# if the number of output lines is specified, use this value, otherwise print all lines
		if self.numLines != 0:
			numLines = self.numLines
		else:
			numLines = len(self.specificDamageRanks)
		for res in self.specificDamageRanks[0:numLines+1]:
			print '{}\t{} {} {}\tSlope: {}\tStd Dev: {}'.format(str(res.damageRank),res.residueName, res.atomType,
																res.boundOrUnbound,str(res.densSlope)[:6],
																str(res.slopeError)[:6])

	def printComparisonPlots(self):

		# run the compareDamageRanks method multiple times to 
		# compare multiple pairs of density metrics

		metrics = ['loss',
				   'gain',
				   'mean',
				   'net']

		norms = [0,1]
		linFitSummary = ""
		for i in range(0,len(metrics)):
			for j in range(i,len(metrics)):
				for k in range(0,len(norms)):
					for l in range(k,len(norms)):
						if not (i == j and k == l):
							linFitSummary += self.compareDamageRanks(metrics[i],norms[k],metrics[j],norms[l])
		print linFitSummary

	def compareDamageRanks(self,met1,norm1,met2,norm2):

		# for two different density change metrics, compare 
		# the rankings for correlations

		self.densMetric = met1
		self.normalised = bool(norm1)
		label1 = '{} D{}'.format(self.normaliseType(),met1)
		self.calculateRanks()
		ranking1 = self.specificDamageRanks

		self.densMetric = met2
		self.normalised = bool(norm2)
		label2 = '{} D{}'.format(self.normaliseType(),met2)
		self.calculateRanks()
		ranking2 = self.specificDamageRanks

		ranks1,ranks2 = [],[]
		for atom in ranking1:
			for otheratom in ranking2:
				if (atom.residueName == otheratom.residueName and atom.atomType == otheratom.atomType
					and atom.boundOrUnbound == otheratom.boundOrUnbound):
					ranks1.append(atom.damageRank)
					ranks2.append(otheratom.damageRank)

		# Create a figure instance
		if seabornFound is True:
			sns.set_palette("deep", desat=.6)
			sns.set_context(rc={"figure.figsize": (12, 12)})

		f = plt.figure()
		ax = plt.subplot(1,1,1)
		plt.scatter(ranks1, ranks2, s=100, c='#d64d4d')
		plt.xlabel('{} ranking'.format(label1), fontsize=18)
		plt.ylabel('{} ranking'.format(label2), fontsize=18)
		f.suptitle('{} vs {} rankings'.format(label1,label2),fontsize=30)
		plt.setp(f.axes)
		f.savefig('{}VS{}Rankings.png'.format(label1.replace(" ", ""),label2.replace(" ", "")))

		# calculate linear fitting for correlation scatter plot
		slope, intercept, r_value, p_value, std_err = stats.linregress(ranks1, ranks2)
		linFitSummary = '{} vs {} rankings'.format(label1,label2)+ ' R2: '+"{0:.2f}".format(r_value**2) +' p-value: '+"{0:.2f}\n".format(p_value)
		return linFitSummary 


