# processedAtom class to refine the PDBmulti list of objects with additional 
# attributes and methods

from classHolder import multiPDB
from scipy import stats

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

		# these attributes are dictionaries and will contain values for multiple
		# variations of the density change metrics
		self.maxDensLoss	= {}     
		self.maxDensGain	= {}
		self.meanDensChange	= {}

	def cloneInfo(self,atom):
		# clone the information from a specific PDBmulti atom to a 
		# new processedAtom
		self.__dict__ = atom.__dict__.copy()

		self.maxDensLoss['Standard']	= self.mindensity     
		self.maxDensGain['Standard']	= self.maxdensity
		self.meanDensChange['Standard']	= self.meandensity

	def calculateLinReg(self,numDatasets):
		# Calculates linear regression for the density change metrics and
		# determines the linear slope of density change
		# 'numDatasets' in the number of difference map datasets across which 
		# linear regression will be preformed
		x = np.array(range(2,numDatasets+1))

		# calculate lin reg for max dens loss metric
		self.maxDensLoss['lin reg'] = {}
		y = np.array(atom.mindensity[0:numDatasets-1])
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

		self.maxDensLoss['lin reg']['slope'] 		= slope
		self.maxDensLoss['lin reg']['intercept'] 	= intercept
		self.maxDensLoss['lin reg']['r_squared'] 	= r_value**2
		self.maxDensLoss['lin reg']['p_value'] 		= p_value
		self.maxDensLoss['lin reg']['std_err'] 		= std_err

		# calculate lin reg for max dens gain metric
		self.maxDensGain['lin reg'] = {}
		y = np.array(atom.maxdensity[0:numDatasets-1])
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

		self.maxDensGain['lin reg']['slope'] 		= slope
		self.maxDensGain['lin reg']['intercept'] 	= intercept
		self.maxDensGain['lin reg']['r_squared'] 	= r_value**2
		self.maxDensGain['lin reg']['p_value'] 		= p_value
		self.maxDensGain['lin reg']['std_err'] 		= std_err

		# calculate lin reg for mean dens change metric
		self.meanDensChange['lin reg'] = {}
		y = np.array(atom.meandensity[0:numDatasets-1])
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

		self.meanDensChange['lin reg']['slope'] 	= slope
		self.meanDensChange['lin reg']['intercept'] = intercept
		self.meanDensChange['lin reg']['r_squared'] = r_value**2
		self.meanDensChange['lin reg']['p_value'] 	= p_value
		self.meanDensChange['lin reg']['std_err'] 	= std_err

	def CalphaWeightedDensChange(self,CalphaWeights):
		# calculates the Calpha weighted density metrics, compensating for 
		# the overall background degregation of the electron density maps 
		# with increasing dose/dataset number. Since the same CalphaWeights
		# (per dataset) are applied to all atoms within the structure,
		# these are an input (and thus not recalculated for every single atom)
		# - see the CalphaWeights class for details on how to calculate these

		self.maxDensLoss['Calpha normalised'] 		= list(np.divide(self.self.maxDensLoss['Standard'],CalphaWeights))
		self.maxDensGain['Calpha normalised'] 		= list(np.divide(self.self.maxDensGain['Standard'],CalphaWeights))
		self.meanDensChange['Calpha normalised'] 	= list(np.divide(self.self.meanDensChange['Standard'],CalphaWeights))



		
