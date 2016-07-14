from classHolder import StructurePDB
from scipy import stats
import numpy as np
import string

class combinedAtom(StructurePDB):

	# A subclass extension for a collection of multiple different
	#  dose pdb file structures as defined by the StructurePDB 
	# class. This class adds additonal attributes and methods

	def __init__(self,
				 atomnum            = 0,
				 residuenum         = 0,
				 atomtype           = "",
				 basetype           = "",
				 chaintype          = "",
				 X_coord            = 0,
				 Y_coord            = 0,
				 Z_coord            = 0,
				 atomID             = "",
				 numsurroundatoms   = 0,
				 numsurroundprotons = 0,
				 densMetric         = {},
				 partialInfo        = False):

		super(combinedAtom, self).__init__(atomnum,
										   residuenum,
										   atomtype,
										   basetype,
										   chaintype,
										   X_coord,
										   Y_coord,
										   Z_coord,
										   atomID,
										   numsurroundatoms,
										   numsurroundprotons)

		self.densMetric = {} # dictionary of density metrics to be filled

	def getPresentDatasets(self):

		# for atoms only present within subset of datasets, 
		# get the datasets for which present
		presentList = []
		i = -1
		for val in self.densMetric['loss']['Standard']['values']:
			i += 1
			if not np.isnan(val):
				presentList.append(i)
		return presentList

	def getDensMetricInfo(self,
						  metric   = 'loss',
						  normType = 'Standard',
						  values   = []):

		# these attributes are dictionaries and will contain 
		# values for multiple variations of the density 
		# change metrics 'normType' is 'standard' or 
		# 'Calpha weighted'

		# special case for Dloss metric (change sign of values)
		if metric == 'loss' and normType == 'Standard':
			values = list(-np.array(values))

		try:
			self.densMetric[metric]
		except KeyError: 
			self.densMetric[metric] = {}

		self.densMetric[metric][normType] = {}
		self.densMetric[metric][normType]['values'] = values

	def calcAvMetric(self,type,densMetric):

		# calculate the average of a given
		#  metric over series of datasets

		densVals = self.densMetric[densMetric][type]['values']
		self.densMetric[densMetric][type]['average'] = np.nanmean(densVals)
		
	def calcLinReg(self,
				   numLinRegDsets = 1,
				   normType       = 'Standard',
				   metric         = 'loss'):

		# Calculates linear regression for the density metric 
		# 'metric' and determines the linear slope of density 
		# change 'numLinRegDsets' in the number of difference 
		# map datasets across which linear regression will be 
		# preformed. 'normType' specifies whether 'Standard'
		# or 'Calpha normalised' metric values should be used

		x = np.array(range(2,numLinRegDsets+1))

		# this ensures no calculation performed if metric values are empty list
		if len(self.densMetric[metric][normType]['values']) == 0: return

		self.densMetric[metric][normType]['lin reg'] = {}
		y = np.array(self.densMetric[metric][normType]['values'][0:len(x)])
		linRegRslts = stats.linregress(x,y)

		for v1,v2 in zip(['slope','intercept','r_squared','p_value','std_err'],linRegRslts):
			if v1 == 'r_value': v2 = v2**2
			self.densMetric[metric][normType]['lin reg'][v1] = v2

	def CalphaWeightedDensChange(self,
								 CalphaWeights = [],
								 metric        = 'loss'):

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

		weight = CalphaWeights.weight[metric]
		metricVals = self.densMetric[metric]['Standard']['values']
		self.getDensMetricInfo(metric   = metric,
							   normType = 'Calpha normalised',
							   values   = list(np.sign(weight)*np.divide(metricVals-weight,weight)))

	def calcFirstDatasetSubtractedMetric(self,
										 normType = 'Standard',
										 metric   = 'loss'):

		# for a specified metric calculate a new metric: 
		# metric(dn)-metric(d1) where dn is dataset n

		metricVals = self.densMetric[metric][normType]['values']
		self.getDensMetricInfo(metric   = metric,
							   normType = 'dataset 1 subtracted',
							   values   = metricVals-metricVals[0]*(np.ones(len(metricVals))))

	def getNumDatasets(self,*metric):

		# get number of datasets in current series

		if len(metric) == 0:
			return len(self.densMetric['loss']['Standard']['values'])
		else:
			return len(self.densMetric[metric[0]]['Standard']['values'])

	def calcDiffFromMeanMetric(self,
							   metric   = 'loss',
							   normType = 'Standard',
							   avMetric = []):

		# for a specified metric calculate difference between 
		# metric for atom and from mean of structure 'avMetric' 
		# is average metric over whole structure

		metricVals = self.densMetric[metric][normType]['values']
		self.getDensMetricInfo(metric   = metric,
							   normType = 'subtract mean',
							   values   = list(np.array(metricVals)-np.array(avMetric)))

	def calcRatioToMeanMetric(self,
							  metric   = 'loss',
							  normType = 'Standard',
							  avMetric = []):

		# for a specified metric calculate difference between
		#  metric for atom and from mean of structure 'avMetric' 
		# is average metric over whole structure

		metricVals = self.densMetric[metric][normType]['values']
		self.getDensMetricInfo(metric   = metric,
							   normType = 'divide mean',
							   values   = list(np.array(metricVals)/np.array(avMetric)))

	def calcNumStdFromMeanMetric(self,
								 metric    = 'loss',
								 normType  = 'Standard',
								 avMetric  = [],
								 stdMetric = []):

		# for a specified metric calculate number of structure-wide
		# standard deviations between metric for atom and from mean
		# of structure. 'avMetric' is average metric over whole structure

		self.calcDiffFromMeanMetric(metric   = metric,
									normType = normType,
									avMetric = avMetric)

		metricVals = self.densMetric[metric]['subtract mean']['values']
		self.getDensMetricInfo(metric   = metric,
							   normType = 'num stds',
							   values   = list(np.array(metricVals)/np.array(stdMetric)))

	def calcNetChangeMetric(self,
							normType = 'Standard'):

		# the following metric attempts to determine whether there 
		# is a net loss, gain or disordering of density associated 
		# with a specific atom

		self.getDensMetricInfo(metric   = 'net',
							   normType = normType,
							   values   = [])

		for dataset in range(0,self.getNumDatasets()):
			abs_maxDensLoss = np.abs(self.densMetric['loss'][normType]['values'][dataset]) 
			abs_maxDensGain = np.abs(self.densMetric['gain'][normType]['values'][dataset]) 
			metricVal = (abs_maxDensLoss - abs_maxDensGain)
			self.densMetric['net'][normType]['values'].append(metricVal)

	def calcVectorWeightedMetric(self,
								 metric   = 'loss',
								 normType = 'Standard',
								 vector   = []):

		# for a specified metric calculate new vector-weighted 
		# values, where the metric over a series of doses
		# is multiplied by a vector of per-dataset scalings

		metricVals = self.densMetric[metric][normType]['values']
		if len(vector) != len(metricVals):
			print 'Incompatible metric and per-dataset scale vector lengths'
			return
		self.getDensMetricInfo(metric   = metric,
							   normType = 'vector-weighted',
							   values   = metricVals/np.array(vector))

	def calcVectorSubtractedMetric(self,
								   metric   = 'loss',
								   normType = 'Standard',
								   vector   = []):

		# for a specified metric calculate new vector-weighted 
		# values, where the metric over a series of doses
		# is multiplied by a vector of per-dataset scalings

		metricVals = self.densMetric[metric][normType]['values']
		if len(vector) != len(metricVals):
			print 'Incompatible metric and per-dataset scale vector lengths'
			return
		self.getDensMetricInfo(metric   = metric,
							   normType = 'vector-subtracted',
							   values   = metricVals - np.array(vector))

	def findSolventAccessibility(self,
								 inputPDBfile = 'untitled.pdb',
								 printText    = False):

		# read in a pdb file output by ccp4 program 'areaimol' 
		# (run independently) which calculates solvent accessibility
		# for each atom within a pdb file and writes value in Bfactor column

		openFile =  open(str(inputPDBfile), "r")
		for line in openFile.readlines():
			if (self.atomtype == str(line[12:16].strip())
				and self.residuenum == str(line[22:27].strip())
				and self.chaintype == str(line[21])  
				and self.basetype == str(line[17:20].strip())):
				self.solventAccess = str(line[60:66].strip()) 
				break

		if printText is True:
			print 'Solvent Accessibility: {}'.format(self.solventAccess)
		openFile.close()
		return self.solventAccess
