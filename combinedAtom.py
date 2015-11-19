# processedAtom class to refine the PDBmulti list of objects with additional 
# attributes and methods

from classHolder import StructurePDB
from scipy import stats
import numpy as np
import string

class combinedAtom(StructurePDB):
	# A subclass extension for a collection of multiple different dose pdb file structures as defined by the StructurePDB class. 
	# This class adds additonal attributes and methods
	def __init__(self,atomnum=0,residuenum=0,atomtype="",basetype="",chaintype="",
				 X_coord=0,Y_coord=0,Z_coord=0,Bfactor=[],Occupancy=[],atomidentifier="",
				 numsurroundatoms=0,numsurroundprotons=0,densMetric={}):

		super(combinedAtom, self).__init__(atomnum,residuenum,atomtype,basetype,chaintype,
										   X_coord,Y_coord,Z_coord,atomidentifier,numsurroundatoms,
										   numsurroundprotons)
		self.getInfo()   

	def getInfo(self):
		# these attributes are dictionaries and will contain values for multiple
		# variations of the density change metrics
		for metricType in ('loss','gain','mean','bfactor','bdamage','|loss|','median'):
			self.densMetric[metricType] = {}
			self.densMetric[metricType]['Standard'] = {}
			self.densMetric[attr]['Standard']['values'] = {}

	def calculateLinReg(self,numLinRegDatasets,type,densityMetric):
		# Calculates linear regression for the density metric 'densityMetric' and
		# determines the linear slope of density change
		# 'numLinRegDatasets' in the number of difference map datasets across which 
		# linear regression will be preformed. 'type' specifies whether 'Standard'
		# or 'Calpha normalised' metric values should be used
		x = np.array(range(2,numLinRegDatasets+1))

		# this ensures no calculation performed if metric values are empty list
		if len(self.densMetric[metricType][type]['values']) == 0: return

		self.densMetric[metricType][type]['lin reg'] = {}
		y = np.array(self.densMetric[metricType][type]['values'][0:len(x)])
		linRegRslts = stats.linregress(x,y)

		for v1,v2 in zip(['slope','intercept','r_squared','p_value','std_err'],linRegRslts):
			if v1 == 'r_value': v2 = v2**2
			self.densMetric[densityMetric][type]['lin reg'][v1] = v2

	def CalphaWeightedDensChange(self,CalphaWeights,metric):
		# calculates the Calpha weighted density metrics for metric "metric",
		# compensating for the overall background degregation of the electron 
		# density maps with increasing dose/dataset number. Since the same 
		# CalphaWeights (per dataset) are applied to all atoms within the structure,
		# these are an input (and thus not recalculated for every single atom)
		# - see the CalphaWeights class for details on how to calculate these

		try:
			CalphaWeights.weight[metric]
		except AttributeError:
			print 'Calpha weights not yet calculated for metric "{}"\n'.format(metric)+\
				  'Need to calculate first this weight first --> see CalphaWeight class'
			return

		self.densMetric[metric]['Calpha normalised'] = {}
		weight = CalphaWeights.weight[metric]
		metricVals = self.densMetric[metric]['Standard']['values']
		self.densMetric[metric]['Calpha normalised']['values'] 	= list(np.divide(metricVals-weight,weight))

	def getNumDatasets(self):
		return len(self.densMetric['loss'][normType]['values'])

	def calculateNetChangeMetric(self,normType):
		# the following metric attempts to determine whether there is a net loss, gain or disordering
		# of density associated with a specific atom
		self.densMetric['net'] = {}
		self.densMetric['net'][normType] = {}
		self.densMetric['net'][normType]['values'] = []
		for dataset in range(0,self.getNumDatasets()):
			abs_maxDensLoss = np.abs(self.densMetric['loss'][normType]['values'][dataset]) 
			abs_maxDensGain = np.abs(self.densMetric['gain'][normType]['values'][dataset]) 
			metricVal = (abs_maxDensLoss - abs_maxDensGain)
			self.densMetric['net'][normType]['values'].append(metricVal)

	def findSolventAccessibility(self,inputPDBfile):
		# read in a pdb file output by ccp4 program 'areaimol' which calculates solvent 
		# accessibility for each atom within a pdb file and writes value in Bfactor column

		# check if atom already has solvent accessibility calculated
		try:
			self.solventAccess
			return self.solventAccess
		except AttributeError:
			# if not already calculated then..
			openFile =  open(str(inputPDBfile), "r")
			for line in openFile.readlines():
				if (self.atomtype == str(line[12:16].strip())
					and self.residuenum == int(line[22:26].strip())
					and self.chaintype == str(line[21])  
					and self.basetype == str(line[17:20].strip())):
					self.solventAccess = str(line[60:66].strip()) 
					return self.solventAccess
		openFile.close()
