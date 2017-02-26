import sys,os,time
sys.path.insert(0,'./lib')
from errors import error
from processFiles import processFiles
from logFile import logFile

class process():

	# a class to overview the generation of atom-tagged and
	# density maps from pdb & mtz files for a specified damage
	# series. This class contains some methods to write input 
	# files, and can call the class processFiles, in which the 
	# actual map-generation (and subsequent calculation of 
	# per-atom damage metrics) can be performed

	def __init__(self,
				 inputFile           = 'fullInput.txt',
				 run                 = True,
				 proceedToMetricCalc = False,
				 skipToMetricCalc    = False,
				 outputGraphs        = 'yes',
				 cleanUpFinalFiles   = False,
				 printOutput         = True,
				 skipToSummaryFiles  = False,
				 writeSummaryFiles   = False):

		self.inputFile           = inputFile
		self.proceedToMetricCalc = proceedToMetricCalc
		self.skipToMetricCalc    = skipToMetricCalc
		self.outputGraphs        = outputGraphs # takes 'yes', 'no' or 'slim'
		self.cleanUpFinalFiles   = cleanUpFinalFiles
		self.printOutput         = printOutput
		self.skipToSummaryFiles  = skipToSummaryFiles
		self.writeSummaryFiles   = writeSummaryFiles

		if run:
			self.run()

	def setInputFile(self,name):

		# specify input file name as required for the run

		self.inputFile = name
		print 'Input file name set as "{}"'.format(self.inputFile)

	def checkInputFileExists(self):

		# check if input file exists and is readable 

		if os.path.isfile(self.inputFile) and os.access(self.inputFile,os.R_OK):
			# print "Input file exists and is readable"
			return True
		else:
			err = "Either input file is missing or is not readable"
			error(text = err)
			return False

	def run(self):

		# run to process files with the specified input file

		success = self.checkInputFileExists()
		if not success: 
			return success

		self.startLogFile()
		self.titleCaption()
		self.info()

		pro = processFiles(inputFile           = self.inputFile,
						   proceedToMetricCalc = self.proceedToMetricCalc,
						   skipToMetricCalc    = self.skipToMetricCalc,
						   outputGraphs        = self.outputGraphs,
						   cleanFinalFiles     = self.cleanUpFinalFiles,
						   logFileObj          = self.logFile,
						   skipToSummaryFiles  = self.skipToSummaryFiles,
						   writeSummaryFiles   = self.writeSummaryFiles)

		return pro.jobSuccess

	def printInputFile(self):

		# print the contents of the specified input file

		fileIn = open(self.inputFile,'r')
		for line in fileIn.readlines():
			if len(line.strip()) != 0:
				print line.strip()
		fileIn.close()

	def writeTemplateInputFile(self,
							   numHigherDoseDatasets = 1):

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
		if printStr:
			print infoString

		return infoString

	def fillerLine(self):

		# a simple filler line to print to command line

		ln = '---------------------------------------------------------------'	
		self.logFile.writeToLog(str = ln)


	def titleCaption(self,
					 whichString = 'new'):

		# a title line to print to command line

		if whichString == 'old':
			txt = '\n||========================== RIDL ==========================||'+\
			      '||======== Radiation-Induced Density Loss Analysis =========||\n'
		else:
			txt = '\n----------------------------------------------------\n'+\
				  ' '*24+'RIDL\n'+\
				  ' '*4+'(Radiation-Induced Density Loss Analysis)\n'+\
				  '----------------------------------------------------'

		self.logFile.writeToLog(str = txt)

	def info(self):

		# print a small ammount of information to the 
		# command line (date and email)

		txt = 'date: {}\n'.format(time.strftime("%c"))+\
			  'email queries to charles.bury@dtc.ox.ac.uk\n'
		self.logFile.writeToLog(str = txt)

	def quickParseInputFile(self):

		# parse input file to retrieve directory name for output

		f = open(self.inputFile,'r')
		for l in f.readlines():
			if len(l) == 0:
				continue
			if 'dir' == l.split()[0]:
				self.outputDir = l.split()[1]
				if self.outputDir[-1] not in ('/'):
					self.outputDir  += '/'
				return

	def startLogFile(self):

		# create log file for the current job

		self.quickParseInputFile()
		self.checkOutputDirExists(printToScreen = self.printOutput)

		logDir = '/'+self.outputDir.strip('/')+'/RIDL-log/'
		if not os.path.exists(logDir):
			os.makedirs(logDir)

		name = '{}RIDLjob.log'.format(logDir)

		logName = lambda x: name.replace('.log','_{}{}'.format(x,'.log'))
		i = 1
		while os.path.isfile(logName(i)): 
			i += 1 
		uniqLogName = logName(i)

		log = logFile(fileName      = uniqLogName,
					  fileDir       = self.outputDir,
					  printToScreen = self.printOutput)
		self.logFile = log

	def checkOutputDirExists(self,
							 printToScreen = True):

		# check whether output directory exists and make if not

		if not self.outputDir.endswith('/'):
			if printToScreen:
				print 'Working directory specified in input '+\
					  'file must end in "/" - appending.'
			self.outputDir += '/'

		if not os.path.exists(self.outputDir):
			if printToScreen:
				print 'Output directory "{}" not found, making directory'.format(self.outputDir)
			os.makedirs(self.outputDir)




