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
		CAweights = CalphaWeight(self.unprocessedAtomList)
		CAweights.calculateWeights()

		# loop over all atoms in list and determine new atom info (defined by processedAtom class)
		for oldAtom in self.unprocessedAtomList:
			newAtom = processedAtom()
			newAtom.cloneInfo(oldAtom)
			newAtom.calculateLinReg(self.numDatasets)
			newAtom.CalphaWeightedDensChange(CAweights)
			processedList.append(newAtom)
		self.processedAtomList = processedList

		





