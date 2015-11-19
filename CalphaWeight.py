# small script to determine the Calpha weighting per dataset, to monitor the 
# overall degregation of the electron density map with increasing dose
# (due to global effects)
import numpy as np 

class CalphaWeight(object):
	def __init__(self, atomList = []):

		self.atomList 	= atomList  # list of atoms over multiple doses
		self.weight 	= {}  		# dictionary of weightings per metric type (mean, gain etc.)

	def calculateWeights(self,metric):
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
		self.weight[metric] = np.mean([atom.densMetric[metric]['Standard']['values'] for atom in CalphaAtoms],0)

	def printWeights(self,metric):
		# print the Calpha weights per dataset to command line

		# don't run function if weights have not been calculated yet
		try:
			self.weight[metric]
		except AttributeError:
			print 'Need to calculate the weights for metric "{}" first (use calculateWeights method)'.format(metric)
			return 

		print '-------------------------------------------------------'
		print 'Calpha weights as follows:'
		print '\tDataset\tD{}'.format(metric)
		for i in range(0,len(self.atomList[0].densMetric[metric]['Standard']['values'])):
			print '\t{}:\t{}'.format(i+1,str(self.weight[metric][i])[:5])




