# processedAtom class to refine the PDBmulti list of objects with additional 
# attributes and methods

from classHolder import multiPDB
from scipy import stats
import numpy as np
import string

class processedAtom(multiPDB):
	# A subclass extension for a collection of multiple different dose pdb file structures as defined by the multiPDB class. This class adds additonal 
	# attributes and methods
	def __init__(self,atomnum=0,residuenum=0,atomtype="",basetype="",chaintype="",
				 X_coord=0,Y_coord=0,Z_coord=0,Bfactor=[],Occupancy=[],meandensity=[],
				 maxdensity=[],mindensity=[],mediandensity=[],atomidentifier="",
				 numsurroundatoms=0,numsurroundprotons=0,bdam=[],bdamchange=[],Bfactorchange=[],
				 meandensity_norm=[],maxdensity_norm=[],mindensity_norm=[],
				 mediandensity_norm=[],numvoxels=[],stddensity=[],min90tile=[],max90tile=[],
				 min95tile=[],max95tile=[],modedensity=[],rsddensity=[],dipstat=[],rangedensity=[]):

		super(processedAtom, self).__init__(atomnum,residuenum,atomtype,basetype,chaintype,
											X_coord,Y_coord,Z_coord,atomidentifier,numsurroundatoms,
											numsurroundprotons,bdam,bdamchange,Bfactorchange,
											meandensity_norm,maxdensity_norm,mindensity_norm,
											mediandensity_norm,numvoxels,stddensity,min90tile,
											max90tile,min95tile,max95tile,modedensity,rsddensity,
											dipstat,rangedensity)   

	def cloneInfo(self,atom):
		# clone the information from a specific PDBmulti atom to a 
		# new processedAtom
		self.__dict__ = atom.__dict__.copy()

		# these attributes are dictionaries and will contain values for multiple
		# variations of the density change metrics
		self.densMetric	= {}
		for metricType in ('loss','gain','mean','bfactor','bdamage','|loss|','median','median-simple','max-simple'):
			self.densMetric[metricType] = {}
			self.densMetric[metricType]['Standard'] = {}
		self.densMetric['loss']['Standard']['values'] = self.mindensity
		self.densMetric['gain']['Standard']['values'] = self.maxdensity
		self.densMetric['mean']['Standard']['values'] = self.meandensity
		self.densMetric['median']['Standard']['values'] = self.mediandensity
		self.densMetric['bfactor']['Standard']['values'] = self.Bfactorchange
		self.densMetric['bdamage']['Standard']['values'] = self.bdamchange
		self.densMetric['|loss|']['Standard']['values'] = [-val for val in self.mindensity]

		self.densMetric['median-simple']['Standard']['values'] = np.array(self.mediandensity)/self.mediandensity[0]
		self.densMetric['max-simple']['Standard']['values'] = np.array(self.maxdensity)/self.maxdensity[0]

	def boundOrUnbound(self):
		# for the TRAP-RNA study, determine whether an atom is a member of an RNA- bound or
		# unbound protein chain.
		unboundChains = list(string.ascii_lowercase.upper())[:11]
		boundChains = list(string.ascii_lowercase.upper())[11:22]
		rnaChains = list(string.ascii_lowercase.upper())[22:23]
		if self.chaintype in unboundChains:
			return 'unbound protein'
		elif self.chaintype in boundChains:
			return 'bound protein'
		elif self.chaintype in rnaChains:
			return 'rna'
		else:
			return 'something else'

	def calculateLinReg(self,numDatasets,type):
		# Calculates linear regression for the density change metrics and
		# determines the linear slope of density change
		# 'numDatasets' in the number of difference map datasets across which 
		# linear regression will be preformed. 'type' specifies whether 'Standard'
		# or 'Calpha normalised' metric values should be used
		x = np.array(range(2,numDatasets+1))

		# calculate lin reg for max dens loss metric
		for metricType in ('loss','gain','mean','net','bfactor','bdamage'):
			# check each metric key is defined and skip if not
			try:
				self.densMetric[metricType]
			except KeyError:
				continue
			# this ensures metrics not included if empty list
			if len(self.densMetric[metricType]['Standard']['values']) == 0:
				continue

			self.densMetric[metricType][type]['lin reg'] = {}
			y = np.array(self.densMetric[metricType][type]['values'][0:len(x)])
			slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

			self.densMetric[metricType][type]['lin reg']['slope'] 		= slope
			self.densMetric[metricType][type]['lin reg']['intercept'] 	= intercept
			self.densMetric[metricType][type]['lin reg']['r_squared'] 	= r_value**2
			self.densMetric[metricType][type]['lin reg']['p_value'] 	= p_value
			self.densMetric[metricType][type]['lin reg']['std_err'] 	= std_err

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
		# B-factor, BDamage etc)

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




