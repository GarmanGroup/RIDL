# class to determine the ordering of specific damage within a MX structure
# using multiple difference maps collected over increasing doses from the 
# same crystal

import numpy as np
from scipy import stats

class specificDamageRank(object):
	def __init__(self, residueName = "", damageRank = 0 , maxDensLossSlope = 0,
				 slopeError = 0):

		self.residueName 		= residueName
		self.damageRank 		= damageRank
		self.maxDensLossSlope 	= maxDensLossSlope
		self.slopeError			= slopeError

class specificDamageRanking(object):
	def __init__(self,PDBmulti = []):
		self.PDBmulti = PDBmulti

	def calculateRanks(self):
		# for each atom within the PDBmulti list of atom objects, group into 
		# residue/nucleotide types and determine the average slope for straight
		# lines of best fit for each residue/nucleotide type (for plots of 
		# dataset number vs max density loss)

		p = 10 # number of damage datasets to include in slope calculation

		residueDict = {}
		for atom in self.PDBmulti:
			if atom.basetype not in residueDict.keys():
				residueDict[atom.basetype]['slope'] = []
				residueDict[atom.basetype]['std_err'] = []

			y = np.array(atom.mindensity[0:p-1])
			x = np.array(range(2,p+1))
			slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
			residueDict[atom.basetype]['slope'].append(slope)
			residueDict[atom.basetype]['std_err'].append(std_err)

		# calculate the average straight line slope for each residue/nucleotide type
		specificDamageRanks = []
		for key in residueDict.keys():
			damageStats = specificDamageRank(key)
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
			print '{}\t{}\tSlope: {}\tStd Dev: {}'.format(str(res.damageRank), res.residueName,
														  str(res.maxDensLossSlope), str(res.slopeError))












