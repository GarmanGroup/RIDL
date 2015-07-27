# class to determine the ordering of specific damage within a MX structure
# using multiple difference maps collected over increasing doses from the 
# same crystal

import numpy as np
from scipy import stats

class specificDamageRank(object):
	def __init__(self, residueName = "", atomType = "", damageRank = 0 ,
				 maxDensLossSlope = 0, slopeError = 0):

		self.residueName 		= residueName
		self.atomType			= atomType
		self.damageRank 		= damageRank
		self.maxDensLossSlope 	= maxDensLossSlope
		self.slopeError			= slopeError

class specificDamageRanking(object):
	def __init__(self,atomList = []):
		self.atomList = atomList

	def calculateRanks(self):
		# for each atom within the atomList list of atom objects, group into 
		# residue/nucleotide types and determine the average slope for straight
		# lines of best fit for each residue/nucleotide type (for plots of 
		# dataset number vs max density loss)

		residueDict = {}
		for atom in self.atomList:
			identifier = "{} {}".format(atom.basetype,atom.atomtype)
			if identifier not in residueDict.keys():
				residueDict[identifier] = {}
				residueDict[identifier]['slope'] = []
				residueDict[identifier]['std_err'] = []

			residueDict[identifier]['slope'].append(atom.maxDensLoss['lin reg']['slope'])
			residueDict[identifier]['std_err'].append(atom.maxDensLoss['lin reg']['std_err'])

		# calculate the average straight line slope for each residue/nucleotide type
		specificDamageRanks = []
		for key in residueDict.keys():
			damageStats = specificDamageRank(key.split()[0],key.split()[1])
			damageStats.maxDensLossSlope = np.mean(residueDict[key]['slope'])
			damageStats.slopeError = np.std(residueDict[key]['slope'])
			specificDamageRanks.append(damageStats)

		# rank residues/nucleotides by slope value
		specificDamageRanks.sort(key = lambda x: x.maxDensLossSlope)
		counter = 0
		for residue in specificDamageRanks:
			residue.damageRank = counter
			counter += 1

		self.specificDamageRanks = specificDamageRanks

	def getDamageRanks(self):
		# following on from the above function, output the calculated ordering of
		# damage within the current structure
		print '--------------------------------------------------------------------'
		print 'Ordering of damage as follows:'
		for res in self.specificDamageRanks:
			print '{}\t{} {}\tSlope: {}\tStd Dev: {}'.format(str(res.damageRank),res.residueName, res.atomType,
														  	 str(res.maxDensLossSlope)[:6], str(res.slopeError)[:6])











