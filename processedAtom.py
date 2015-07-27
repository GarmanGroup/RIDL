# processedAtom class to refine the PDBmulti list of objects with additional 
# attributes and methods

from classHolder import multiPDB
from scipy import stats
import numpy as np

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
		self.maxDensLoss	= {}     
		self.maxDensGain	= {}
		self.meanDensChange	= {}

		self.maxDensLoss['Standard'] 	= {}
		self.maxDensGain['Standard'] 	= {}
		self.meanDensChange['Standard'] = {}

		self.maxDensLoss['Standard']['values']		= self.mindensity     
		self.maxDensGain['Standard']['values']		= self.maxdensity
		self.meanDensChange['Standard']['values']	= self.meandensity

	def calculateLinReg(self,numDatasets,type):
		# Calculates linear regression for the density change metrics and
		# determines the linear slope of density change
		# 'numDatasets' in the number of difference map datasets across which 
		# linear regression will be preformed. 'type' specifies whether 'Standard'
		# or 'Calpha normalised' metric values should be used
		x = np.array(range(2,numDatasets+1))

		# calculate lin reg for max dens loss metric
		self.maxDensLoss[type]['lin reg'] = {}
		y = np.array(self.maxDensLoss[type]['values'][0:numDatasets-1])
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

		self.maxDensLoss[type]['lin reg']['slope'] 		= slope
		self.maxDensLoss[type]['lin reg']['intercept'] 	= intercept
		self.maxDensLoss[type]['lin reg']['r_squared'] 	= r_value**2
		self.maxDensLoss[type]['lin reg']['p_value'] 	= p_value
		self.maxDensLoss[type]['lin reg']['std_err'] 	= std_err

		# calculate lin reg for max dens gain metric
		self.maxDensGain[type]['lin reg'] = {}
		y = np.array(self.maxDensGain[type]['values'][0:numDatasets-1])
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

		self.maxDensGain[type]['lin reg']['slope'] 		= slope
		self.maxDensGain[type]['lin reg']['intercept'] 	= intercept
		self.maxDensGain[type]['lin reg']['r_squared'] 	= r_value**2
		self.maxDensGain[type]['lin reg']['p_value'] 	= p_value
		self.maxDensGain[type]['lin reg']['std_err'] 	= std_err

		# calculate lin reg for mean dens change metric
		self.meanDensChange[type]['lin reg'] = {}
		y = np.array(self.meanDensChange[type]['values'][0:numDatasets-1])
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

		self.meanDensChange[type]['lin reg']['slope'] 		= slope
		self.meanDensChange[type]['lin reg']['intercept'] 	= intercept
		self.meanDensChange[type]['lin reg']['r_squared'] 	= r_value**2
		self.meanDensChange[type]['lin reg']['p_value']		= p_value
		self.meanDensChange[type]['lin reg']['std_err'] 	= std_err

	def CalphaWeightedDensChange(self,CalphaWeights):
		# calculates the Calpha weighted density metrics, compensating for 
		# the overall background degregation of the electron density maps 
		# with increasing dose/dataset number. Since the same CalphaWeights
		# (per dataset) are applied to all atoms within the structure,
		# these are an input (and thus not recalculated for every single atom)
		# - see the CalphaWeights class for details on how to calculate these

		try:
			CalphaWeights.weight_MaxDensLoss
		except AttributeError:
			print 'Calpha weights not yet calculated.. need to calculate first, see CalphaWeight class'
			return

		self.maxDensLoss['Calpha normalised'] 		= {}
		self.maxDensGain['Calpha normalised'] 		= {}
		self.meanDensChange['Calpha normalised'] 	= {}

		self.maxDensLoss['Calpha normalised']['values'] 	= list(np.divide(self.maxDensLoss['Standard']['values'],CalphaWeights.weight_MaxDensLoss))
		self.maxDensGain['Calpha normalised']['values'] 	= list(np.divide(self.maxDensGain['Standard']['values'],CalphaWeights.weight_MaxDensGain))
		self.meanDensChange['Calpha normalised']['values'] 	= list(np.divide(self.meanDensChange['Standard']['values'],CalphaWeights.weight_MeanDensChange))



		
