# script to convert multiPDB list of atom objects to list of processedAtom class objects
# (which have additional attributes and methods over those of multiPDB)
from processedAtom import processedAtom
from CalphaWeight import CalphaWeight

class processedAtomList(object):
	# class for list of atom objects defined by processedAtom class

	def __init__(self,unprocessedAtomList=[],numDatasets=0):
		self.unprocessedAtomList = unprocessedAtomList
		self.numDatasets = numDatasets # number of datasets to perform linear regression over

		# if number of datasets to perform lin reg not stated then default
		# to number of datasets present within unprocessedAtomList
		if self.numDatasets == 0:
			self.numDatasets = len(self.unprocessedAtomList[0].mindensity)

	def processAtomList(self):
		# process the input multiPDB list of atom objects to create new
		# list of atom objects from processedAtom class
		processedList = []

		# calculate the Calpha weights for each dataset (see CalphaWeight class for details)
		print 'Calculating Calpha weights at each dataset...'
		CAweights = CalphaWeight(self.unprocessedAtomList)
		CAweights.calculateWeights()

		# loop over all atoms in list and determine new atom info (defined by processedAtom class)
		print 'Creating new list of atom objects within class processedAtom...'
		for oldAtom in self.unprocessedAtomList:
			newAtom = processedAtom()
			newAtom.cloneInfo(oldAtom)
			newAtom.CalphaWeightedDensChange(CAweights)
			newAtom.calculateAdditionalMetrics()
			newAtom.calculateLinReg(self.numDatasets,'Standard')
			newAtom.calculateLinReg(self.numDatasets,'Calpha normalised')

			processedList.append(newAtom)
		self.processedAtomList = processedList


	def graphMetric(self):
		# produce a graph of selected metric against dataset number for a specified atom
		import seaborn as sns
		import matplotlib.pyplot as plt

		# only run this function if atom list has been processed already
		try:
			self.processedAtomList
		except AttributeError:
			print 'Need to process atom list before this function can be run'
			return

		self.atomType 	= raw_input("Atom type: ")
		self.baseType 	= raw_input("Residue/nucleotide type: ")
		self.residueNum = int(raw_input("Residue number: "))
		self.densityMetric = raw_input("Density metric type: ")
		self.normalise  = raw_input("Normalisation type: ")

		# find the correct atom specified by above properties
		equivAtoms = []
		for atom in self.processedAtomList:
			if (atom.atomtype == self.atomType and atom.basetype == self.baseType 
				and atom.residuenum == self.residueNum):
				equivAtoms.append(atom)

		# if an atom with specified properties not found, don't proceed
		if len(equivAtoms) == 0:
			print 'Atom not found within structure.. make sure atom properties correctly specified'
			return

		# define x range here (damage set numbers)
		x = range(2,len(equivAtoms[0].mindensity)+2)

		sns.set(style="white", context="talk")
		f, axes = plt.subplots(1, 1, figsize=(12, 12), sharex=True)

		# determine y values here dependent on density metric type specified 
		for foundAtom in equivAtoms:
			y = foundAtom.densMetric[self.densityMetric][self.normalise]['values']
			try:
				y
			except NameError:
				print 'Unexpected density metric name.. please try again'
				return

			# plot the graph here
			plt.plot(x,y)

		plt.xlabel('Dataset', fontsize=18)
		plt.ylabel('{} D{}'.format(self.normalise,self.densityMetric), fontsize=18)
		f.suptitle('{} D{}: {} {} {}'.format(self.normalise,self.densityMetric,
												self.baseType,self.residueNum,
												self.atomType),fontsize=24)

		plt.show()




		





