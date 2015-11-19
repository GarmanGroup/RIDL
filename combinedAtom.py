# processedAtom class to refine the PDBmulti list of objects with additional 
# attributes and methods

from classHolder import singlePDB
from scipy import stats
import numpy as np
import string

class combinedAtom(singlePDB):
	# A subclass extension for a collection of multiple different dose pdb file structures as defined by the singlePDB class. 
	# This class adds additonal attributes and methods
	def __init__(self,atomnum=0,residuenum=0,atomtype="",basetype="",chaintype="",
				 X_coord=0,Y_coord=0,Z_coord=0,Bfactor=[],Occupancy=[],meandensity=[],
				 maxdensity=[],mindensity=[],mediandensity=[],atomidentifier="",
				 numsurroundatoms=0,numsurroundprotons=0,bdam=[],bdamchange=[],Bfactorchange=[],
				 meandensity_norm=[],maxdensity_norm=[],mindensity_norm=[],
				 mediandensity_norm=[],numvoxels=[],stddensity=[],min90tile=[],max90tile=[],
				 min95tile=[],max95tile=[],rsddensity=[],rangedensity=[]):

		super(combinedAtom, self).__init__(atomnum,residuenum,atomtype,basetype,chaintype,
										   X_coord,Y_coord,Z_coord,atomidentifier,numsurroundatoms,
										   numsurroundprotons,bdam,bdamchange,Bfactorchange,
									       numvoxels,stddensity,min90tile,max90tile,min95tile,
									       max95tile,rsddensity,rangedensity)
		self.getInfo()   

	def getInfo(self):
		# these attributes are dictionaries and will contain values for multiple
		# variations of the density change metrics
		self.densMetric	= {}
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
		for vals in (['slope',slope],['intercept',intercept],['r_value',r_value**2],
					 ['p_value',p_value],['std_err',std_err]):	

		for v1,v2 in zip(['slope','intercept','r_squared','p_value','std_err'],linRegRslts)
			if v1 == 'r_value': v2 = v2**2
			self.densMetric[densityMetric][type]['lin reg'][v1] = v2

	def CalphaWeightedDensChange(self,CalphaWeights):
		# calculates the Calpha weighted density metrics, compensating for 
		# the overall background degregation of the electron density maps 
		# with increasing dose/dataset number. Since the same CalphaWeights
		# (per dataset) are applied to all atoms within the structure,
		# these are an input (and thus not recalculated for every single atom)
		# - see the CalphaWeights class for details on how to calculate these

		try:
			CalphaWeights.weight['loss']
		except AttributeError:
			print 'Calpha weights not yet calculated.. need to calculate first, see CalphaWeight class'
			return

		for metricType in ('loss','gain','mean'):
			self.densMetric[metricType]['Calpha normalised'] = {}
			weight = CalphaWeights.weight[metricType]
			metric = self.densMetric[metricType]['Standard']['values']
			self.densMetric[metricType]['Calpha normalised']['values'] 	= list(np.divide(metric-weight,weight))

	def calculateAdditionalMetrics(self):
		# calculates addition metrics for each atom. These secondary metrics are 
		# algebraic expressions of the primary metrics (Dmean, Dloss, Dgain, 
		# B-factor etc..)

		self.densMetric['weighted loss'] = {}
		for normaliseType in ('Standard','Calpha normalised'):
			self.densMetric['weighted loss'][normaliseType] = {}
			self.densMetric['weighted loss'][normaliseType]['values'] = []
			for dataset in range(0,len(self.mindensity)):
				absMaxDensLoss = np.abs(self.densMetric['loss'][normaliseType]['values'][dataset])
				absMaxDensGain = np.abs(self.densMetric['gain'][normaliseType]['values'][dataset])
				metricVal = absMaxDensLoss/(absMaxDensLoss + absMaxDensGain)
				self.densMetric['weighted loss'][normaliseType]['values'].append(metricVal)

		# the following metric attempts to determine whether there is a net loss, gain or disordering
		# of density associated with a specific atom
		self.densMetric['net'] = {}
		for normaliseType in ('Standard','Calpha normalised'):
			self.densMetric['net'][normaliseType] = {}
			self.densMetric['net'][normaliseType]['values'] = []
			for dataset in range(0,len(self.mindensity)):
				abs_maxDensLoss = np.abs(self.densMetric['loss'][normaliseType]['values'][dataset]) 
				abs_maxDensGain = np.abs(self.densMetric['gain'][normaliseType]['values'][dataset]) 
				metricVal = (abs_maxDensLoss - abs_maxDensGain)
				self.densMetric['net'][normaliseType]['values'].append(metricVal)

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
