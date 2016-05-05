import sys,os
sys.path.insert(0,'./lib')
from processFiles import processFiles

class process():

	# a class to overview the generation of atom-tagged and
	# density maps from pdb & mtz files for a specified damage
	# series. This class contains some methods to write input 
	# files, and can call the class processFiles, in which the 
	# actual map-generation (and subsequent calculation of 
	# per-atom damage metrics) can be performed

	def __init__(self,
				 inputFile         = 'fullInput.txt',
				 run               = True,
				 proceedToETRACK   = False,
				 skipToETRACK      = False,
				 outputGraphs      = True,
				 cleanUpFinalFiles = False):

		self.titleCaption('ETRACK file preparation')

		self.inputFile         = inputFile
		self.proceedToETRACK   = proceedToETRACK
		self.skipToETRACK      = skipToETRACK
		self.outputGraphs      = outputGraphs
		self.cleanUpFinalFiles = cleanUpFinalFiles

		if run is True:
			self.run()

	def setInputFile(self,name):

		# specify input file name as required for the run

		self.inputFile = name
		print 'Input file name set as "{}"'.format(self.inputFile)

	def checkInputFileExists(self):

		# check if input file exists and is readable 

		if os.path.isfile(self.inputFile) and os.access(self.inputFile,os.R_OK):
			print "Input file exists and is readable"
			return True
		else:
			print "Either input file is missing or is not readable"
			return False

	def run(self):

		# run to process files with the specified input file

		success = self.checkInputFileExists()
		if success is False: 
			return success
		pro = processFiles(inputFile       = self.inputFile,
						   proceedToETRACK = self.proceedToETRACK,
						   skipToETRACK    = self.skipToETRACK,
						   outputGraphs    = self.outputGraphs,
						   cleanFinalFiles = self.cleanUpFinalFiles)
		return pro.jobSuccess

	def printInputFile(self):

		# print the contents of the specified input file

		fileIn = open(self.inputFile,'r')
		for line in fileIn.readlines():
			if len(line.strip()) != 0:
				print line.strip()
		fileIn.close()

	def writeTemplateInputFile(self,numHigherDoseDatasets = 1):

		# write a template input file to current directory to be completed

		f = open(self.inputFile,'w')
		highDoseStrs = [','.join(['???']*numHigherDoseDatasets)]*4
		string   = 	'dir ???\n'+\
					'\nINITIALDATASET\nname1 ???\nmtz1 ???\nmtzlabels1 ???\npdb1 ???\nRfreeFlag1 ???\n'+\
					'\nLATERDATASET\nname2 {}\nmtz2 {}\nmtzlabels2 {}\npdb2 {}\n'.format(*highDoseStrs)+\
					'\nPHASEDATASET\nname3 ???\nmtz3 ???\nmtzlabels3 ???\n'+\
					'\nMAPINFO\nsfall_VDWR 1\ndensMapType DIFF\nFFTmapWeight False\ndeleteIntermediateFiles TRUE'
		f.write(string)
		f.close()
		print 'Template input file "{}" written'.format(self.inputFile)

	def writePDBredoInputFile(self,pdb1,pdb2,inputFileDir,outputDir):

		# write an input file that is suitable for a pdb_redo-downloaded damage series
		# inputFileDir = '/Users/charlie/DPhil/PDBredo_damageSeries/??'
		# outputDir = '/Users/charlie/DPhil/YEAR2/JAN/??'

		self.setInputFile('fullInput_{}-{}DIFF.txt'.format(pdb2,pdb1))
		f = open(self.inputFile,'w')
		string 	 = 	'dir {}\n'.format(outputDir)+\
					'\nINITIALDATASET\n'+\
					'name1 {}init\n'.format(pdb1)+\
					'mtz1 {}{}/{}.mtz\n'.format(inputFileDir,pdb1,pdb1)+\
					'mtzlabels1 P_{}\n'.format(pdb1)+\
					'pdb1 {}{}/{}.pdb\n'.format(inputFileDir,pdb1,pdb1)+\
					'RfreeFlag1 FreeR_flag\n'+\
					'\nLATERDATASET\n'+\
					'name2 {}\n'.format(pdb2)+\
					'mtz2 {}{}/{}.mtz\n'.format(inputFileDir,pdb2,pdb2)+\
					'mtzlabels2 P_{}\n'.format(pdb2)+\
					'pdb2 {}{}/{}.pdb\n'.format(inputFileDir,pdb2,pdb2)+\
					'\nPHASEDATASET\n'+\
					'name3 {}init\n'.format(pdb1)+\
					'mtz3 {}{}/{}.mtz\n'.format(inputFileDir,pdb1,pdb1)+\
					'mtzlabels3 C_{}\n'.format(pdb1)+\
					'\nMAPINFO\n'+\
					'sfall_VDWR 1\n'+\
					'densMapType DIFF\n'+\
					'FFTmapWeight True'
		f.write(string)
		f.close()

	def howToWriteInputFile(self,
							printStr = True):

		# some information on how write the input file

		infoString 	= """
*** Information on how to successfully write an input file for job ***

* Note:
* If multiple higher-dose datasets are to be processed successively in a single batch,
* then please input a comma-separated list for each argument within the LATERDATASET
* section below:

FILELOCATION
	dir 		: full path to output working directory

INITIALDATASET
	name1 		: assign a name to your low dose damage set
	mtz1  		: full path to the low dose mtz file
	mtzlabels1	: F & SIGF column labels (look in low dose mtz file), if "FP_X" then type "P_X"
	pdb1 		: full path to the low dose pdb file
	RfreeFlag1 	: the Rfree flag label within the low dose mtz file (e.g. "FreeR_flag")

LATERDATASET
	name2 		: assign a name to your high dose damage set (e.g. a pdb code "1qid")
	mtz2 		: full path to the high dose mtz file
	mtzlabels2 	: F & SIGF column labels (look in high dose mtz file), if "FP_X" then type "P_X"
	pdb2 		: full path to the high dose pdb file

PHASEDATASET
	name3 		: assign a name to your low dose damage set (same as INITIALDATASET above)
	mtz3 		: full path to the low dose mtz file (same as INITIALDATASET above)
	mtzlabels3 	: PHI phase column labels (look in low dose mtz file), if "PHIC_X" then type "C_X"

MAPINFO
	For difference map analysis, do not change these parameters
	"""
		if printStr is True:
			print infoString
		return infoString

	def writeTestInputFile(self,dataset):

		# write an input file for the test case on the github README.md

		inds = ['d','e','f']
		if dataset not in range(1,4):
			print 'Select dataset between 1-3'
			return
		ind = inds[dataset-1]

		self.setInputFile('testInput{}.txt'.format(dataset))
		f = open(self.inputFile,'w')
		string 	 =  'dir ./testOutput/ETRACK/\n'+\
					'INITIALDATASET\n'+\
					'name1 1qidinit\n'+\
					'mtz1 ./testOutput/1qid/1qid.mtz\n'+\
					'mtzlabels1 P_1qid\n'+\
					'pdb1 ./testOutput/1qid/1qid.pdb\n'+\
					'RfreeFlag1 FreeR_flag\n'+\
					'LATERDATASET\n'+\
					'name2 1qi{}\n'.format(ind)+\
					'mtz2 ./testOutput/1qi{}/1qi{}.mtz\n'.format(ind,ind)+\
					'mtzlabels2 P_1qi{}\n'.format(ind)+\
					'pdb2 ./testOutput/1qi{}/1qi{}.pdb\n'.format(ind,ind)+\
					'PHASEDATASET\n'+\
					'name3 1qidinit\n'+\
					'mtz3 ./testOutput/1qid/1qid.mtz\n'+\
					'mtzlabels3 C_1qid\n'+\
					'MAPINFO\n'+\
					'sfall_VDWR 1\n'+\
					'densMapType DIFF\n'+\
					'FFTmapWeight True'
		f.write(string)
		f.close()

	def fillerLine(self):

		# a simple filler line to print to command line

		print '---------------------------------------------------------------'	

	def titleCaption(self,title):

		# a template title line to print to command line

		print '\n\n||========================== {} ==========================||'.format(title)

