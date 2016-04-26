import sys
sys.path.insert(0,'./lib')
from eTrack_RUN import eTrack
import os

class run():
	# run the main ETRACK scripts to calculate per-atom density 
	# metrics for a series of increasing doses
	def __init__(self,
				 runAll       = True,
				 inputFileLoc = ''):

		self.inputFileName = inputFileLoc+'e_Track_inputfile.txt'
		
		if runAll is True:
			self.runETRACK()

	def runETRACK(self,
				  mapProcess  = True,
				  postProcess = True,
				  retrieve    = False):
		# run the ETRACK processing for the currently defined input file

		print 'Running main ETRACK script...'
		exists = self.checkInputFileExists()
		if exists is False:
			return
		eT = eTrack()
		eT.runPipeline(map_process   = mapProcess,
					   post_process  = postProcess,
					   retrieve      = retrieve,
					   inputFileName = self.inputFileName)
		self.et = eT

	def defineDoseList(self,doses,names,version):
		# do not include first dataset dose if difference maps chosen
		if version == 'DIFF':
			if len(doses.split(',')[1:])>1:
				dosesOut = ','.join(doses.split(',')[1:])
				namesOut = ','.join(names.split(',')[1:])
			else:
				dosesOut = doses.split(',')[1]
				namesOut = names.split(',')[1]
			return dosesOut,namesOut
		else: return doses,names

	def writeBlankInputFile(self):
		# Need to create input file for a generic damage series to be completed by the user

		args = ['inDir <input file location>',
				'outDir <output file location>',
				'damageset_name <consistent name of series, e.g. for TRAP1.pdb, TRAP2.pdb, this is TRAP>',
				'damageset_num <dataset numbers, e.g. 1,2 for TRAP. Can be letters corresponding to pdb series>',
				'initialPDB <initial dataset pdb file e.g. TRAP1.pdb>',
				'doses <list of doses, length must match length of damageset_num above>',
				'PKLMULTIFILE <if already processed in ETRACK, can specify single output .pkl for series>']
		inputString = self.writeInputFile(*args)	

	def writeInputFile(self,
					   inDir         = '',
					   outDir        = '',
					   damSetName    = '',
					   laterDatasets = '',
					   initialPDB    = '',
					   doses         = '',
					   PKLMULTIFILE  = ''):
		# write a generic input file for a damage series here

		inputString = 'inDir {}\n'.format(inDir)+\
					  'outDir {}\n'.format(outDir)+\
					  'damageset_name {}\n'.format(damSetName)+\
					  'laterDatasets {}\n'.format(laterDatasets)+\
					  'initialPDB {}\n'.format(initialPDB)+\
					  'doses {}\n'.format(doses)
		if PKLMULTIFILE != '':
			inputString += 'PKLMULTIFILE {}'.format(PKLMULTIFILE)
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()	

	def checkInputFileExists(self):
		if os.path.isfile(self.inputFileName) and os.access(self.inputFileName,os.R_OK):
			print "Input file exists and is readable"
			return True
		else:
			print "Either input file is missing or is not readable"
			print 'Check whether "{}" is a suitable input file'.format(self.inputFileName)
			return False


