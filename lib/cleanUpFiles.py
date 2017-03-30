import os
import shutil

class cleanUpFinalFiles():

	def __init__(self,
				 cleanMapDir    = True,
				 outputDir      = './',
				 keepFCmap      = False,
				 keepAtomTagMap = False,
				 keepDensityMap = True,
				 keepMapDir     = True,
				 keepPdbs       = True,
				 printText      = True):

		# class to clean up output directories 
		# after a run has been performed

		self.outputDir      = outputDir
		self.keepFCmap      = keepFCmap
		self.keepAtomTagMap = keepAtomTagMap
		self.keepDensityMap = keepDensityMap
		self.keepPdbs       = keepPdbs

		self.mapDir    = '{}RIDL-maps/'.format(self.outputDir)
		self.metricDir = '{}RIDL-metrics/'.format(self.outputDir)

		if cleanMapDir:
			if self.checkSubdirExists(self.mapDir):
				self.cleanUpMapDir()
				if printText:
					self.statePurpose()

		if not keepMapDir:
			self.removeMapDir()

	def statePurpose(self):

		# print purpose to command line

		print 'Unnecessary files removed from output directories'

	def checkSubdirExists(self,
						  subdir   = './',
						  printTxt = True):

		# check that the output directory exists

		if not os.path.isdir(subdir):
			if printTxt:
				print 'Error! Subdirectory {} not found.'.format(subdir)
			return False

		return True

	def cleanUpMapDir(self):

		# clean up the map processing directory by 
		# removing unwanted end files

		for f in os.listdir(self.mapDir):
			self.removeFCmap(fName = f)
			self.removeAtomtaggedMap(fName = f)
			self.removeDensityMap(fName = f)

	def removeFCmap(self,
					fName  = ''):

		# remove file if it is a FC map

		if not self.keepFCmap:
			if fName.endswith('_FC.map'):
				os.remove(self.mapDir+fName)

	def removeAtomtaggedMap(self,
							fName  = ''):

		# remove file if it is an 'atom-tagged' map

		if not self.keepAtomTagMap:
			if fName.endswith('_atoms.map'):
				os.remove(self.mapDir+fName)

	def removeDensityMap(self,
						 fName  = ''):

		# remove file if it is an 'atom-tagged' map

		if not self.keepDensityMap:
			if fName.endswith('_density.map'):
				os.remove(self.mapDir+fName)

	def removeMapDir(self):

		# give option to remove the RIDL-maps/ directory
		# following a successful run of the pipeline

		if self.checkSubdirExists(self.mapDir,printTxt=False):
			print 'Removing RIDL-maps/ directory'
			shutil.rmtree(self.mapDir)






