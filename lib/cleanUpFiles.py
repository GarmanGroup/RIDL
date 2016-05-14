import os

class cleanUpFinalFiles():

	def __init__(self,
				 autoRun        = True,
				 outputDir      = './',
				 keepFCmap      = False,
				 keepAtomTagMap = False,
				 keepDensityMap = True,
				 keepPdbs       = True,
				 printText      = True):

		# class to clean up output directories 
		# after a run has been performed

		self.outputDir      = outputDir
		self.keepFCmap      = keepFCmap
		self.keepAtomTagMap = keepAtomTagMap
		self.keepDensityMap = keepDensityMap
		self.keepPdbs       = keepPdbs

		self.mapProcessDir  = '{}maps/'.format(self.outputDir)
		self.calculationDir = '{}metrics/'.format(self.outputDir)

		if autoRun is True:
			if self.checkCorrectOutput() is True:
				self.cleanUpMapProcessingDir()
				if printText is True:
					self.statePurpose()

	def statePurpose(self):

		# print purpose to command line

		print 'Unnecessary files removed from output directories'

	def checkCorrectOutput(self):

		# check that the output directory 
		# format seems reasonable for a run

		subdirs = [self.mapProcessDir,
				   self.calculationDir]
		for subdir in subdirs:	  
			if os.path.isdir(subdir) is False:
				print 'Error! Subdirectory {} not found.'.format(subdir)
				return False

		return True

	def cleanUpMapProcessingDir(self):

		# clean up the map processing directory by 
		# removing unwanted end files

		for f in os.listdir(self.mapProcessDir):
			self.removeFCmap(fName = f)
			self.removeAtomtaggedMap(fName = f)
			self.removeDensityMap(fName = f)

	def removeFCmap(self,
					fName  = ''):

		# remove file if it is a FC map

		if self.keepFCmap is False:
			if fName.endswith('_FC.map') is True:
				os.remove(self.mapProcessDir+fName)

	def removeAtomtaggedMap(self,
							fName  = ''):

		# remove file if it is an 'atom-tagged' map

		if self.keepAtomTagMap is False:
			if fName.endswith('_atoms.map') is True:
				os.remove(self.mapProcessDir+fName)

	def removeDensityMap(self,
						 fName  = ''):

		# remove file if it is an 'atom-tagged' map

		if self.keepDensityMap is False:
			if fName.endswith('_density.map') is True:
				os.remove(self.mapProcessDir+fName)






