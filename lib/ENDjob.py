# run shell script code to create END map for dataset
# NOTE: model must have been through one round of phenix refinement and have .eff file available
from mapTools import mapTools
from ccp4Job import checkInputsExist
import os

class ENDjob():

	def __init__(self,
				 pdbFile   = 'untitled.pdb',
				 mtzFile   = 'untitled.mtz',
				 effFile   = 'untitled.eff',
				 outputDir = './',
				 gridSamps = [0,0,0],
				 runLog    = ''):

		self.pdbFile		= pdbFile
		pdbName 			= (pdbFile.split('/')[-1]).split('.pdb')[0]
		self.mtzFile 		= mtzFile
		self.effFile		= effFile
		self.outputDir		= outputDir
		self.gridSamps		= gridSamps
		self.outputMapFile 	= '{}/{}_END.map'.format(outputDir,pdbName)
		self.runLog 		= runLog

		self.runLog.writeToLog('Running END map generation job')

	def run(self):
		inputFiles = [self.pdbFile,self.mtzFile,self.effFile]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False
		self.runEND()
		if self.jobSuccess is True:
			self.provideFeedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMapFile))
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

	def runEND(self):
		self.outputLogfile = 'ENDlogfile.txt'
		cmdLineStr = 'END_RAPID.com {} cycles=1 seeds=5 cpus=1 g1={} g2={} g3={} -nofofc -nosigf -norapid'.format(self.effFile,*self.gridSamps)
		os.system('{} > {}'.format(cmdLineStr,self.outputLogfile))
		os.system('mv 2FoFc_END.map {}'.format(self.outputMapFile))
		os.system('mv {} {}/{}'.format(self.outputLogfile,self.outputDir,self.outputLogfile))

		if 'ENDmapFiles' not in os.listdir('./'):
			os.makedirs('ENDmapFiles')

		for f in self.filesForDir():
			os.system('mv {} ENDmapFiles/{}'.format(f,f))
		os.system('mv ENDmapFiles {}/ENDmapFiles'.format(self.outputDir))

		self.jobSuccess = self.checkJobSuccess()

	def checkJobSuccess(self):
		# job success checked, based on whether output files exist
		if os.path.isfile(self.outputMapFile) is False:
			print 'END map generation did not proceed to completion'
			return False
		else:
			return True

	def provideFeedback(self):
		# provide some feedback
		print '--------------------------'
		print 'END map Summary:'
		print 'Input pdb file: {}'.format(self.pdbFile)
		print 'Input mtz file: {}'.format(self.mtzFile)
		print 'Output map file: {}'.format(self.outputMapFile)
		Map = mapTools(self.outputMapFile)
		Map.printMapInfo()
		print '--------------------------'

	def filesForDir(self):
		filesForDir = ['2FoFc_scaled.map',
					   '2FoFc.map',
					   'cadded.mtz','dens.txt',
					   'find_F000_001_f_model.mtz',
					   'find_F000_001.eff',
					   'find_F000_001.geo',
					   'find_F000_001.log',
					   'find_F000_001.mtz',
					   'find_F000_001.pdb',
					   'find_F000_data.mtz',
					   'find_F000_input.eff',
					   'find_F000.eff',
					   'fmodel_ss.map',
					   'fmodel.mtz',
					   'FoFc_scaled.map',
					   'FoFc.map',
					   'kickme.mtz',
					   'neg.map',
					   'nobulk.map',
					   'nobulk.mtz',
					   'phenix_fmodel.log',
					   'refined.eff',
					   'refined.pdb',
					   'scaled_map_coeffs.mtz',
					   'scaled.mtz',
					   'scaleit_fobs.log',
					   'scaleit_map.log',
					   'sfall.log',
					   'sfall.map',
					   'sfallme.pdb',
					   'sharpsolvent_vac.log',
					   'sharpsolvent.map',
					   'ss.mtz',
					   'vacuum_zero.map']
					   
		return filesForDir



