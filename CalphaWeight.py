# small script to determine the Calpha weighting per dataset, to monitor the 
# overall degregation of the electron density map with increasing dose
# (due to global effects)
import numpy as np 

class CalphaWeight(object):
	def __init__(self, atomList = []):

		self.atomList = atomList

	def calculateWeights(self):
		# collect the Calpha atoms and for each dataset number determine
		# the Calpha weight 

		# don't run function if no list of atoms included
		if len(self.atomList) == 0:
			print 'Need to add list of atoms first'
			return

		# collect the Calpha atoms here
		CalphaAtoms = []
		for atom in self.atomList:
			if atom.atomtype == 'CA':
				CalphaAtoms.append(atom)

		# calculate the weighting for each dataset and for each density
		# metric individually here
		self.weight_MaxDensLoss = np.mean([atom.mindensity for atom in CalphaAtoms],0)
		self.weight_MaxDensGain = np.mean([atom.maxdensity for atom in CalphaAtoms],0)
		self.weight_MeanDensChange = np.mean([atom.meandensity for atom in CalphaAtoms],0)

	def printWeights(self):
		# print the Calpha weights per dataset to command line

		# don't run function if weights have not been calculated yet
		try:
			self.weight_MaxDensLoss
		except AttributeError:
			print 'Need to calculate the weights first (use calculateWeights method)'
			return 

		print '-------------------------------------------------------'
		print 'Calpha weights as follows:'
		print '\tDataset\tDloss\tDmean\tDgain'
		for counter in range(0,len(self.atomList[0].mindensity)):
			print '\t{}\t{}\t{}\t{}'.format(counter,str(self.weight_MaxDensLoss[counter])[:5],
													str(self.weight_MeanDensChange[counter])[:5],
													str(self.weight_MaxDensGain[counter])[:5])




