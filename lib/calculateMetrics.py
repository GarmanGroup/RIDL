print 'Importing packages...'
import time
t0 = time.time()
from combinedAtomList import combinedAtomList
t1 = time.time()
from checkSeabornPresent import checkSeabornPresent as checkSns
t1 = time.time()
from savevariables import retrieve_objectlist,save_objectlist,saveGenericObject,retrieveGenericObject
from PDBFileManipulation import PDBtoList,writePDBline
from readAtomMap import maps2DensMetrics
from time import gmtime, strftime
from numpy import isnan
from shutil import move
from os import path,makedirs,listdir,system
print 'Total import time: {}'.format(time.time()-t0)

class calculateMetrics(object):

	# a class for retrieving the RIDL input text file information and running 
	# the RIDL pipeline to calculate per-atom damage metrics for a specified 
	# damage series. This code requires previously generated atom-tagged and
	# density maps (typically Fobs(n) - Fobs(1) Fourier difference maps) to 
	# have been generated for the current damage series (as specified within 
	# the input .txt file parsed below). If run as part of the full RIDL 
	# pipeline (by running 'python runRIDL.py -i [inputfile.txt] -pc') then 
	# this will automatically run directly after the suitable map files have
	# been generated, with no need to explicitly write a new input file for
	# this class to work.

	def __init__(self,
				 inDir      = './',
				 outDir     = './',
				 pdbNames   = [],
				 pklFiles   = [],
				 initialPDB = "",
				 seriesName = "untitled-series",
				 pklSeries  = "",
				 doses      = [],
				 plot       = 'no',
				 output     = 'simple',
				 logFile    = ''):

		self.inDir 		= inDir 		# the input file directory
		self.outDir  	= outDir 		# the output file directory
		self.pdbNames 	= pdbNames 		# the list of pdb codes for series
		self.pklFiles 	= pklFiles 		# list of pkl files from map_processing
		self.initialPDB = initialPDB 	# the first dataset pdb code
		self.seriesName = seriesName 	# the general series name 
		self.pklSeries	= pklSeries 	# the combined series pkl file from post_processing
		self.doses 		= doses 		# list of increasing doses
		self.plot       = plot 			# (bool) decide to plot per-residue summary plots per dataset
		self.output     = output        # the amount of output to provide (either 'simple' for just Dloss
										# info or 'full' larger selection of output files)
		self.logFile    = logFile
		
	def runPipeline(self,
					map_process   = True,
					post_process  = True,
					retrieve      = True,
					inputFileName = ''):

		# the following function reads in the above functions one by one in a scripted 
		# pipeline. Takes inputs to specify which parts of above pipeline are to be
		# included. For each function input specify True if this part is to be performed
		# and False otherwise.

		self.inputFileName = inputFileName

		# check whether valid inputs to function
		valid = self.checkValidInputs(map_process,post_process,retrieve)
		if valid is False: 
			return

		# first need to run function above to read in input file containing info
		# on where input files are and where output files should be written
		ln = 'Reading input file: {}'.format(self.inputFileName)
		self.logFile.writeToLog(str = ln)

		success = self.readInputFile()
		if success is False: 
			return
		success = self.checkInOutDirExist()
		if success is False: 
			return
		self.setOutputDirs()

		# check whether seaborn exists
		self.checkSeaborn()

		if map_process is True:
			self.map_processing()
		else:
			ln = 'Map processing task not chosen...'
			self.logFile.writeToLog(str = ln)
		self.fillerLine()

		if post_process is True:
			self.post_processing()

			# save metric data to pkl file
			pklSeries = saveGenericObject(obj      = self.combinedAtoms,
										  fileName = self.seriesName)

			move(pklSeries,
				 '{}{}'.format(self.outputDir,pklSeries))
			self.pklSeries = pklSeries

			inputfile = open(self.inputFileName,'a')
			inputfile.write('\npklDataFile ' + pklSeries)
			inputfile.close()

			self.provideFeedback()

		else: 
			ln = 'Post processing job not chosen...'
			self.logFile.writeToLog(str = ln)

		if retrieve is True:
			self.PDBmulti_retrieve()
			self.provideFeedback()

		self.fillerLine(blank = True)

		if self.seabornFound is False:
			self.printNoSeabornWarning()

	def checkValidInputs(self,map_process,post_process,retrieve):

		# check the runPipeline inputs to make sure that they are valid

		for var in [map_process,post_process,retrieve]:
			if var not in (True,False):
				print 'Only Boolian parameters expected here'
				return False
		try:
			f = open(self.inputFileName,'r')
			f.close()
		except IOError:
			err = 'metric calculation input file "{}" not found'.format(self.inputFileName)
			self.logFile.writeToLog(str = err)
			return False
		return True

	def readInputFile(self,
					  printText = True):

		# read input file e_Track_input.txt to specify location 
		# of input files and where to write output files

		props = {'inDir'          : 'inDir',
				 'outDir'         : 'outDir',
				 'damageset_name' : 'seriesName',
				 'initialPDB'     : 'initialPDB',
				 'pklDataFile'    : 'pklSeries',
				 'laterDatasets'  : 'laterDatasets'}

		inputfile = open(self.inputFileName,'r')
		pklFiles = []
		for line in inputfile.readlines():
			l = line.split()
			if '#' == line[0]: 
				continue
			elif l[0] in props.keys():
				setattr(self,props[l[0]],l[1])
			elif 'damageset_num' in l[0]:
				datasetNums = l[1]
			elif 'PKLFILE' == l[0]:
				self.pklFiles.append(l[1])
			elif 'doses' == l[0]:
				self.doses 	= [float(d) for d in l[1].split(',')]
			elif 'plot' == l[0]:
				self.plot         = True
				self.plotHeatmaps = True
			elif 'slim-plot' == l[0]:
				self.plot         = True
				self.plotHeatmaps = False
		inputfile.close()

		if not self.initialPDB.endswith('.pdb'):
			self.initialPDB += '.pdb'

		# locate the correct format for the list of datasets within damage 
		# series. Currently two formats acceptable: (a) series-name + dataset-id
		# (per dataset), (b) input list of full dataset names (recommended).
		found = True
		try:
			datasetNums
		except UnboundLocalError:
			found = False
		if found is True:
			self.pdbNames = [self.seriesName+num for num in datasetNums.split(',')]
			return True
		found = True
		try:
			self.laterDatasets
		except AttributeError:
			found = False
		if found is True:
			self.pdbNames = self.laterDatasets.split(',')
			return True
		else:
			if printText is True:
				err = 'Error! Unable to extract list of dataset names from input file'
				self.logFile.writeToLog(str = err)
			return False

	def checkInOutDirExist(self):

		# check that an input/output directories have been 
		# found and make subdirectories if present

		for dir in ([[self.inDir,'Input'],[self.outDir,'Output']]):
			if path.isdir(dir[0]) == False:
				err = '{} file location: {} does not exist. '.format(dir[1],dir[0])+\
					  'Please select an appropriate directory'
				self.logFile.writeToLog(str = err)
				return False
		return True

	def makeOutputDir(self,
					  dirName   = './'):

		# if the above sub directory does not exist, make it

		if not path.exists(dirName):
			makedirs(dirName)
			ln = 'New sub directory "{}" created to contain output files'.format(dirName.replace(self.outputDir,''))
			self.logFile.writeToLog(str = ln)

	def makeNewPlotSubdir(self,
						  subdir = 'untitled/'):

		# make a new sub directory for plots within plots/ subdir
		outDir = self.outputPlotDir + subdir
		self.makeOutputDir(dirName = outDir)

		return outDir

	def setOutputDirs(self):

		# set the locations of the output directories

		self.outputDir 		= '{}RIDL-metrics/'.format(self.outDir)
		self.outputPlotDir 	= '{}plots/'.format(self.outputDir)

		# add pkl file names as attribute if specified in input file
		if len(self.pklFiles) != 0: 
			self.pklFiles = [self.outputDir+f for f in self.pklFiles]

	def map_processing(self):

		# combine the density map and atom-tagged map for a given dataset,
		# to calculate per-atom density metrics for each refined atom

		txt = 'Combining density maps and atom-tagged maps to calculate '+\
			  'per-atom density metrics for each refined atom in structure.\n'
		self.logFile.writeToLog(str = txt)

		txt = 'input directory:  {}\n'.format(self.inDir)+\
			  'output directory: {}'.format(self.outputDir)
		self.logFile.writeToLog(str   = txt,
								strip = False)

		# create additional subdirectories
		for oDir in (self.outputDir,self.outputPlotDir):
			self.makeOutputDir(dirName = oDir)

		# make pklFiles and dir to move all generated per-dataset pkl files to this
		pklFileDir = 'pklFiles-perDataset/'
		self.makeOutputDir(dirName = '{}{}'.format(self.outputDir,pklFileDir))

		txt = '{} higher dose datasets located within input file.'.format(len(self.pdbNames))+\
			  'Calculating per-atom density metrics for each dataset individually.\n'

		pklFileNames = []
		for i,dataset in enumerate(self.pdbNames):

			# derive per-atom density metrics from maps
			mapName1 = '{}_atoms.map'.format(dataset)
			mapName2 = '{}_density.map'.format(dataset)
			mapName3 = '{}_FC.map'.format(self.initialPDB.replace('.pdb',''))

			txt = '\n---------------------------------\n'+\
				  'Higher dose dataset {} starts here'.format(i+1)
			self.logFile.writeToLog(str = txt)

			maps2DensMets = maps2DensMetrics(filesIn     = self.inDir,
										     filesOut    = self.outputDir,
										     pdbName     = dataset,
										     atomTagMap  = mapName1,
										     densityMap  = mapName2,
										     FCmap 	     = mapName3,
										     plotScatter = False,
										     plotHist    = self.plot,
										     plotBar     = False,
										     logFile     = self.logFile)

   			maps2DensMets.maps2atmdensity()

			# move pkl file to working output directory
			pklFileName = maps2DensMets.pklFileName

			move(pklFileName,
				 '{}{}{}'.format(self.outputDir,
								 pklFileDir,
								 pklFileName))

			pklFileNames.append('{}{}{}'.format(self.outputDir,
												pklFileDir,
												pklFileName))

		self.pklFiles = pklFileNames

	def post_processing(self):

		# group the per-atom density metrics for each dataset together

		txt = 'Combining density metric information for each dataset '+\
			  'together within the damage series'
		self.logFile.writeToLog(str = txt)

		txt = 'Input pkl files for post processing chosen from input file:'
		for file in self.pklFiles: 
			txt += '\n\t{}'.format(file.replace(self.outDir,""))
		self.logFile.writeToLog(str = txt)

		# next read in the pdb structure file as list of atom objects
		initialPDBlist = PDBtoList(pdbFileName = self.get1stDsetPDB())

		# retrieve object lists of atoms for each damage set
		ln = 'Reading in damaged pkl files...'
		self.logFile.writeToLog(str = ln)

		dList = [] # dataset list
		for pkl_filename in self.pklFiles:
			ln = 'Damage file number: {}'.format(len(dList)+1)
			self.logFile.writeToLog(str = ln)
			PDB_ret = retrieve_objectlist(fileName = pkl_filename,
										  logFile  = self.logFile)

			# add new retrieved damage set list to dList
			dList.append(PDB_ret)

		# create a list of atom objects with attributes as lists varying over 
		# dose range, only including atoms present in ALL damage datasets
		ln = 'New list of atoms over full dose range calculated...'
		self.logFile.writeToLog(str = ln)
		combinedAtoms = combinedAtomList(datasetList    = dList,
										 numLigRegDsets = len(dList),
										 doseList       = self.doses,
										 initialPDBList = initialPDBlist,
										 outputDir      = self.outputDir,
										 seriesName     = self.seriesName)

		combinedAtoms.getMultiDoseAtomList()

		# calculate 'average' variant Dloss metrics
		combinedAtoms.calcAdditionalMetrics(newMetric = 'average')

		# calculate Calpha normalised metrics, if Calpha atoms exist
		if self.checkCalphasPresent(atomObjList = combinedAtoms) is True:
			for m in ('loss','mean','gain','Bfactor'):
				combinedAtoms.calcAdditionalMetrics(metric = m)
		
		self.combinedAtoms = combinedAtoms

	def provideFeedback(self,
						standardFeedback = False,
					    includeTests     = True,
					    writeCsvs        = False,
					    writeSumFile     = False,
					    writeTopDamSites = False,
					    plotHeatMaps     = False):

		# create series of output feedback files

		# define a 'standard' feedback for the program
		if standardFeedback is True:
			writeCsvs        = True
			writeSumFile     = True
			writeTopDamSites = True
			includeTests     = False

		for norm in ('Standard','Calpha normalised'):

			if norm == 'Calpha normalised':
				if self.checkCalphasPresent(atomObjList = self.combinedAtoms) is False:
					continue

			self.tests_FurtherAnalysis(perform  = includeTests,
									   normfype = norm)

		if writeCsvs is True:
			self.writeCsvFiles(output = self.output)

		# provide summary html file for Dloss metric per-dataset
		if writeSumFile is True:
			self.summaryFile(normType = 'Standard') 
			if self.checkCalphasPresent(atomObjList = self.combinedAtoms) is True:
				self.summaryFile(normType = 'Calpha normalised') 
		
		if writeTopDamSites is True:
			self.writeDamSitesToFile()

		# create heatmap plots (large files) if requested. Here
		# the plotting can be overridden within this method by 
		# setting 'plotHeatMaps' to True
		if self.plot is True:
			if self.plotHeatmaps is True and plotHeatMaps is True:

				subdir = 'metric_heatmap/'
				outDir = self.makeNewPlotSubdir(subdir = subdir)

				for norm in ('Standard','Calpha normalised'):
					if norm == 'Calpha normalised':
						if self.checkCalphasPresent(atomObjList = self.combinedAtoms) is False:
							continue
					self.combinedAtoms.densityMetricHeatMap(saveFig   = True,
													        metric    = 'loss',
													        normType  = norm,
													        outputDir = outDir)	

	def writeCsvFiles(self,
					  moveCsv     = True,
					  output      = 'simple',
					  inclGroupby = True):

		# write atom numbers and density metrics to simple 
		# csv files,one for each density metric separately.
		# inclGroupby = True also writes csv files for atoms
		# grouped by residue type and atom type

		CalphasPresent = self.checkCalphasPresent(atomObjList = self.combinedAtoms)
		self.fillerLine()
		ln = 'Writing .csv file for per-atom density metric:'
		self.logFile.writeToLog(str = ln)

		if output == 'simple':

			metrics = [['loss','Standard'],
					   ['loss','reliability']]

			if CalphasPresent is True:
				metrics += [['loss','Calpha normalised']]
		else:
			metrics = self.combinedAtoms.getDensMetrics()

		for densMet in metrics:
			ln = '\tmetric: {}, normalisation: {}'.format(*densMet)
			self.logFile.writeToLog(str = ln)
			self.combinedAtoms.writeMetric2File(where    = self.outputDir,
										        metric   = densMet[0],
										        normType = densMet[1])

		if inclGroupby is True:
			for groupBy in ('residue','atomtype'):
				self.combinedAtoms.writeMetric2File(where   = self.outputDir,
											   		groupBy = groupBy)
				if CalphasPresent is True:
					self.combinedAtoms.writeMetric2File(where    = self.outputDir,
											  	        groupBy  = groupBy,
										            	normType = 'Calpha normalised')

		# make csvFiles dir and move all generated csv files to this
		if moveCsv is True:
			self.makeOutputDir(dirName = '{}csvFiles/'.format(self.outputDir))
			self.makeOutputDir(dirName = '{}csvFiles/Calpha-normalised/'.format(self.outputDir))
			self.makeOutputDir(dirName = '{}csvFiles/Standard/'.format(self.outputDir))

			for fName in listdir(self.outputDir):
				if fName.endswith(".csv"):
					if '-Calphanormalised' in fName:
						loc = 'Calpha-normalised/'
					else:
						loc = 'Standard/'
					move('{}{}'.format(self.outputDir,fName),
						 '{}csvFiles/{}{}'.format(self.outputDir,loc,fName))

	def PDBmulti_retrieve(self):

		# retrieve list of atom objects from .pkl file

		self.fillerLine(blank = True)

		txt = 'Retrieving per-atom damage metric information from '+\
			  ' .pkl file for each dataset within damage series.'
		self.logFile.writeToLog(str = txt)

		txt = 'Input pkl file for data retrieval chosen from input file:\n'+\
			  '\t{}'.format(self.pklSeries)
		self.logFile.writeToLog(str = txt)

		# retrieve the combinedAtoms object from the pkl file
		self.combinedAtoms = retrieveGenericObject(fileName = self.outputDir+self.pklSeries)

	def summaryFile(self,
					fileType = 'html',
					metric   = 'loss',
					normType = 'Standard'):

		# write a summary output file (either html or plain text)

		ln = 'Writing {} summary output file for metric: {}, normalisation: {}'.format(fileType,metric,normType)
		self.logFile.writeToLog(str = ln)

		if fileType == 'html':
			self.summaryHTML(metric        = metric,
							 normType      = normType,
							 includeGraphs = self.plot)
		else:
			print 'Unknown file format. Only currently supported format is "html"'

	def summaryHTML(self, 
					metric        = 'loss', 
					normType      = 'Standard',
					includeGraphs = True,
					imWidth       = 750):

		# produce a selection of per-dataset summary statistics

		if normType == 'Calpha normalised':
			norm = 'C<sub>&#945</sub>-normalised'
		else:
			norm = normType

		numDsets = self.getNumDatasets()
		summaryFile = open('{}summaryFile-D{}-{}.html'.format(self.outputDir,metric,normType.replace(' ','-')),'w')
		summaryFile.write('<!DOCTYPE html>\n<html>\n')

		# create html head 
		headString = '<head>\n'+\
  					 '<meta name="viewport" content="width=device-width, initial-scale=1">\n'+\
  					 '<link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css">\n'+\
  					 '<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.2/jquery.min.js"></script>\n'+\
  					 '<script src="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js"></script>\n'+\
  					 '<title>RIDL summary file</title>\n'+\
  					 '<style>\ntable, th, td {\nborder: 1px solid black;\nborder-collapse: collapse;\n}\n'+\
					 'th, td {\npadding: 5px;\ntext-align: center;\n}\n</style>\n'+\
					 '</head>\n'

		summaryFile.write(headString)

		bodyString = '<body>\n'+\
					 '<div class="container">\n'+\
					 '<h1>D<sub>{}</sub> ({}) RIDL summary file</h1>\n'.format(metric,norm)+\
					 'Created: {}<br>\n'.format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))+\
					 'Summary information derived from {}<br>\n'.format(self.pklSeries)+\
					 'Email charles.bury@dtc.ox.ac.uk for queries<br>\n'+\
					 'Number of datasets reported in file: {}<br>\n'.format(numDsets)+\
					 'Providing analysis on each dataset individually below<br>\n'+\
					 '<h3>Key of reported quantities</h3>\n'+\
					 '<ul><li>mean: average of metric calculated over all atoms of a specified type</li>\n'+\
					 '<li>std: standard deviation of metric calculated over all atoms of a specified type</li>\n'+\
					 '<li>#atoms: total number of atoms of a specified type</li>\n'+\
					 '<li>outliers: assuming a symmetric distn around the mode, number of atoms that fall outside this domain</li>\n'+\
					 '<li>skew: skewness of metric distribution for atoms of specified type</li>\n'+\
					 '<li>normality: p-value for null hypothesis that distribution of metric values is normally distributed.\n'+\
					 'If not enough atoms are present to perform this test, "n/a" is given</li>\n</ul>'
		summaryFile.write(bodyString)

		# create plot of top n damage sites per dataset by residue type
		figName = self.combinedAtoms.getTopAtomsStackedBarplot(outputDir = self.outputPlotDir, n = 25)

		# provide some links to useful output files
		bodyString = '<h3>Links to useful output files</h3>\n'+\
					 '<ul><li><a href = "csvFiles/{}/{}-{}.csv">{} D<sub>{}</sub> csv file</a></li>\n'.format(normType.replace(' ','-'),metric,normType.replace(' ',''),normType,metric)+\
					 '<li><a href = "plots/{}">Top 25 damage sites per residue/nucleotide type</a></li>'.format(figName.split('/')[-1])
		

		# create heatmap plots (large files) if requested
		if self.plot is True:
			if self.plotHeatmaps is True:
				bodyString += '<li><a href = "plots/metric_heatmap/heatmap_metric-{}_normalisation-{}.svg">{} D<sub>{}</sub> per-atom heat map</a></li></ul>'.format(metric,normType.replace(' ',''),normType,metric)
			else: 
				bodyString += '</ul>'
		else:
			bodyString += '</ul>'

		summaryFile.write(bodyString)

		# find top overall damaged atoms and plot line plots 
		# versus dataset. Only plot if more than 1 high dose
		# dataset is present.
		if self.combinedAtoms.getNumDatasets() == 1:
			pass
		elif includeGraphs is True:
			numDamSites = 10

			figCalls = []
			for t in ['Top','Bottom']:

				topAtoms = self.combinedAtoms.getTopNAtoms(dataset  = 'all', 
														   n        = numDamSites,
														   topOrBot = t.lower())

				aInfoDic = {'resType'   : [],
							'atomType'  : [],
							'chainType' : [],
							'resNum'    : []}

				for atom in topAtoms:
					aInfo = atom.split('-')

					aInfoDic['resType'].append(aInfo[2])
					aInfoDic['atomType'].append(aInfo[3])
					aInfoDic['chainType'].append(aInfo[0])
					aInfoDic['resNum'].append(aInfo[1])

				figName = self.plotLinePlot(resType   = aInfoDic['resType'],
				 				  			atomType  = aInfoDic['atomType'],
				 				  			chainType = aInfoDic['chainType'],
				 				  			resNum    = aInfoDic['resNum'],
				 				  			metric    = metric,
				 				  			normType  = normType,
				 				  			figTitle  = '{} {} D{} damage sites'.format(t,numDamSites,metric),
				 				  			saveName  = 'Lineplot_Metric-D{}_Normalisation-{}_{}-atoms'.format(metric,normType,t))
				figCalls.append('<img class="img-responsive" src="plots/{}" width="{}">'.format(figName.split('/')[-1],400))

			info = '<div class = "row">\n'+\
			  	   '<div class = "col-sm-6">{}</div>\n'.format(figCalls[0])+\
			       '<div class = "col-sm-6">{}</div>\n'.format(figCalls[1])+\
			       '</div>\n'

			summaryFile.write(info)

		# get structure-wide metric average & std dev
		av,std = self.combinedAtoms.getAverageMetricVals(densMet  = metric,
														 normType = normType)
		for i in range(numDsets):
			summaryFile.write('<hr>')

			# start block of collapsible panels here
			summaryFile.write('<div class = "panel-group" id = "datasetinfo{}">'.format(i))
			summaryFile.write('<h2>Dataset {}</h2>\n'.format(i+1))

			t = 'Dataset info'
			c = 'Number in series : {}<br>\n'.format(i+1)+\
				'Dose (MGy)       : {}<br>\n'.format(self.doses[i])+\
				'Number of atoms  : {}<br>\n'.format(self.combinedAtoms.getNumAtoms())+\
				'Fourier diff map : <a href ="../RIDL-maps/{}_density.map">Download</a><br>\n'.format(self.pdbNames[i])
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)	
			summaryFile.write(panelStr)

			t = 'Structure-wide D<sub>{}</sub> summary'.format(metric)
			c = 'Average D<sub>{}</sub>   : {}<br>\n'.format(metric,round(av[i],3))+\
				'Std dev in D<sub>{}</sub>: {}<br>\n'.format(metric,round(std[i],3))
					 				
			if normType == 'Calpha normalised':
				CAweights = self.combinedAtoms.retrieveCalphaWeight(metric = metric)
				c += 'Calpha weight for current dataset: {}<br>\n'.format(round(CAweights.weight[metric][i],3))
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)	
			summaryFile.write(panelStr)

			if includeGraphs is True:

				subdir = 'metricDistn_allAtoms/'
				outDir = self.makeNewPlotSubdir(subdir = subdir)

				failedPlots = self.makeDistnPlots(densMet   = metric,
											      normType  = normType,
											      plotSet   = 4,
											      outputDir = outDir,
											      dataset   = i)

				t = 'Distribution of D<sub>{}</sub> for all refined atoms within structure'.format(metric)
				c = '<img class="img-responsive" src="plots/{}DistnPlot_Residues-all_Metric-D{}_Normalisation-{}_Dataset-{}.svg" width="{}">'.format(subdir,metric,normType.replace(' ',''),i,imWidth)
				panelStr = self.writeHtmlDropDownPanel(title   = t,
											           content = c,
										               dataset = i)
				summaryFile.write(panelStr)

			subdir = 'zScorePlots/'
			outDir = self.makeNewPlotSubdir(subdir = subdir)

			figName = self.combinedAtoms.plotNumAtomsWithMetricAboveStructureWideMean(metric   = metric,
																					  normType = normType,
																					  dataset  = i,
																					  outputLoc = outDir)

			t = '# atoms with {} D<sub>{}</sub> metric above N std dev of structure-wide mean'.format(norm,metric)
			c = '<img class="img-responsive" src="plots/{}{}" width="{}">'.format(subdir,figName.split('/')[-1],imWidth)
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)	
			summaryFile.write(panelStr)

			statsOut = self.combinedAtoms.getTopNAtomsString(metric   = metric,
															 normType = normType,
															 dataset  = i,
															 n        = 25)
			t = 'Top hits ranked by D<sub>{}</sub> metric'.format(metric)
			c = 'Top hits ranked by D<sub>{}</sub> metric.<br>D<sub>mean</sub> and D<sub>gain</sub> are mean '.format(metric)+\
			    'and maximum voxel difference density values, respectively, '+\
			    'assigned within a local region around each atom.<br>Proximity (from 0 '+\
			    'to 1) is a measure of the closeness of the voxel exhibiting the '+\
			    'maximum density loss D<sub>loss</sub> value from the specified atom (higher '+\
			    'values indicate smaller distances):<br><br>\n'
			c += self.convertPlainTxtTable2html(statsOut, width = '50%')
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)
			summaryFile.write(panelStr)

			statsOut = self.combinedAtoms.getPerAtmtypeStats(metric   = metric,
															 normType = normType,
															 dataset  = i)
			t = 'Per-atom-type Statistics'
			c = 'D<sub>{}</sub> metric ranked by mean value:<br><br>\n'.format(metric)
			c += self.convertPlainTxtTable2html(statsOut[0],width='60%')
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)
			summaryFile.write(panelStr)

			n = self.combinedAtoms.getNumAtoms()*0.05 # take top 5% of atoms here
			statsOut = self.combinedAtoms.breakdownTopNatomsBy(metric   = metric,
															   normType = normType,
															   dataset  = i,
															   n        = n)
			t = 'Most frequently damaged atom types'
			c = '{}<br><br>\n'.format(statsOut[0])
			c += self.convertPlainTxtTable2html(statsOut[1])
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)
			summaryFile.write(panelStr)

			statsOut = self.combinedAtoms.getPerResidueStats(metric   = metric,
															 normType = normType,
															 dataset  = i)
			t = 'Per-residue Statistics'
			c = 'D<sub>{}</sub> metric ranked by mean value:<br><br>\n'.format(metric)
			c += self.convertPlainTxtTable2html(statsOut[0],width='60%')
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)
			summaryFile.write(panelStr)

			statsOut = self.combinedAtoms.getPerChainStats(metric   = metric,
														   normType = normType,
														   dataset  = i,
														   n        = 'all')
			t = 'Per-chain Statistics'
			c = 'D<sub>{}</sub> metric ranked by mean value:<br><br>\n'.format(metric)
			c += self.convertPlainTxtTable2html(statsOut[0], width = '60%')
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)
			summaryFile.write(panelStr)

			if includeGraphs is True:

				params = [['Residues',1,'GLUASPCYSMETTYR','residue'],
						  ['nucleotides',3,'DADCDGDT','nucleotide']]

				for paramSet in params:
					subdir = 'metricDistn_key{}/'.format(paramSet[0])
					outDir = self.makeNewPlotSubdir(subdir = subdir)

					failedPlots = self.makeDistnPlots(densMet   = metric,
												      normType  = normType,
												      plotSet   = paramSet[1],
												      outputDir = outDir,
												      dataset   = i)

					if failedPlots[paramSet[2]] is False:

						args = [subdir,
								paramSet[2],
								metric,
								normType.replace(' ',''),
								i,
								imWidth]

						t = 'Distribution of D<sub>{}</sub> for known susceptible {} types\n'.format(metric,paramSet[3])
						c = '<img class="img-responsive" src="plots/{}DistnPlot_Residues-{}_Metric-D{}_Normalisation-{}_Dataset-{}.svg" width="{}"><br>'.format(*args)
						panelStr = self.writeHtmlDropDownPanel(title   = t,
													           content = c,
											                   dataset = i)
						summaryFile.write(panelStr)

			infoString = 'Atoms with unusually <font color="red">high</font> '+\
				  		 'or <font color="blue">low</font> D<sub>{}</sub> values '.format(metric)+\
				  		 'relative to the mean D<sub>{}</sub> value for that specific '.format(metric)+\
				  		 'atom type will be reported.<br><br>'

			lastInfoString = ''
			suspAtomLens   = []
			notDoneBefore  = True

			for t in [6,5,4,3,2]:
				suspAtoms,highOrLow = self.combinedAtoms.detectSuspiciousAtoms(dataset   = i,
							  										 		   metric    = metric,
							  										 		   normType  = normType,
							  										 		   threshold = t)

				tmpInfoString = '{} atoms found with unusually high/low D{} '.format(len(suspAtoms),metric)+\
							    'values compared to average of that atom type '+\
							    '(over {} standard deviations from average)'.format(t)

				if len(suspAtoms) > 0:

					if notDoneBefore is True:
						infoString += lastInfoString + tmpInfoString
						notDoneBefore = False
					else:
						infoString += tmpInfoString

					infoString += ':<br>'
					suspAtomLens.append(len(suspAtoms))
					for s,h in zip(suspAtoms,highOrLow):
						if h == 'high':
							c = 'red'
						else:
							c = 'blue'
						infoString += '<font color="{}">{}</font><br>'.format(c,s)
					infoString += '<br>'

				if len(suspAtomLens) > 1:
					if suspAtomLens[-1] > 0 and suspAtomLens[-2] > 0:
						break

				lastInfoString = tmpInfoString + '<br>'


			t = 'Suspicious Atoms'
			c = infoString
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)
			summaryFile.write(panelStr)

			# end block of collapsible panels here
			summaryFile.write('</div>')


		endString = '</div>\n'+\
					'</body>\n'+\
					'</html>'

		summaryFile.write(endString)
		summaryFile.close()

	def writeHtmlDropDownPanel(self,
							   title      = 'untitled',
							   header     = 'h3',
							   content    = 'no content',
							   panelColor = 'default',
							   dataset    = 0):

		# write the info for a html collapsible panel

		try:
			self.panelIndex
		except AttributeError:
			self.panelIndex = 1

		if self.panelIndex == 1:
			expanded = 'in'
		else:
			expanded = ''
		i   = dataset
		txt = '<div class="panel panel-{}">\n'.format(panelColor)+\
			  '<div class="panel-heading">\n'+\
			  '<{} class="panel-title">\n'.format(header)+\
			  '<a data-toggle="collapse" data-parent="#datasetinfo{}" href="#collapse{}">{}</a>\n'.format(i,self.panelIndex,title)+\
			  '</{}></div>\n'.format(header)+\
			  '<div id="collapse{}" class="panel-collapse collapse {}">'.format(self.panelIndex,expanded)+\
			  '<div class="panel-body">{}'.format(content)+\
			  '</div>\n</div>\n</div>\n'

		self.panelIndex += 1

		return txt

	def convertPlainTxtTable2html(self,
								  plainText      = '',
								  numHeaderLines = 1,
								  width          = '40%',
								  border         = '1'):

		# convert a plain text table of values into a html format table
 # border="{}" style="width:{}"
		htmlTable = '<table class="table table-striped">\n'.format(border,width)
		for i, l in enumerate(plainText.split('\n')):
			if i < numHeaderLines:
				newStr = '<tr><th>'+' '.join(l.split()).replace(' ','</th><th>')+'</th></tr>'
			else:
				newStr = '<tr><td>'+' '.join(l.split()).replace(' ','</td><td>')+'</td></tr>'
			htmlTable += newStr+'\n'
		htmlTable += '</table>\n'
		return htmlTable

	def getNumDatasets(self):

		# get number of datasets in damage series

		numDsets = self.combinedAtoms.atomList[0].getNumDatasets()
		return numDsets

	def writeDamSitesToFile(self, 
							metric      = 'loss', 
							normType    = 'Standard', 
							numDamSites = 'all'):

		# write top damage sites to .pdb file for each dataset

		self.damSitesPDB = []
		for i in range(self.getNumDatasets()): 
			damPDB = self.combinedAtoms.getTopNAtomsPDBfile(metric   = metric,
															normType = normType,
															dataset  = 'all',
															n        = numDamSites,
															pdbFile  = self.get1stDsetPDB())

			self.damSitesPDB.append(damPDB)

	def colorByMetric(self, 
					  metric    = 'loss', 
					  normType  = 'Standard', 
					  dataset   = 0, 
					  singleRes = ''):

		# for the initial pdb file in damage series, convert the Bfactor 
		# column to values for the specified metric. If 'singleRes' is 
		# specified (use 3-letter residue code) then the average density 
		# metric for each atom within this residue is calculated within 
		# the structure and only this is output to resulting PDB file, 
		# otherwise use ''

		if normType == 'Calpha normalised': 
			self.combinedAtoms.calcAdditionalMetrics(metric=metric,normType=normType)

		pdbIn = open(self.get1stDsetPDB(), 'r')
		fileOut = self.outputDir+self.initialPDB.strip('.pdb')+\
				  '_{}D{}_{}.pdb'.format(normType.replace(" ",""),metric,dataset)
		if singleRes != '':
			fileOut = fileOut.strip('.pdb')+'-{}.pdb'.format(singleRes)

		pdbOut = open(fileOut,'w')
		pdbOut.write('REMARK\tBfactor column replaced by {} D{} metric values\n'.format(normType,metric))
		for l in pdbIn.readlines():
			if l.split()[0] in ('CRYST1','SCALE1','SCALE2','SCALE3'):
				pdbOut.write(l)
			elif 'ATOM' in l.split()[0]:
				break
		pdbIn.close()

		if singleRes == '':
			for atm in self.combinedAtoms.atomList:
				dens = atm.densMetric[metric][normType]['values'][dataset]
				if not isnan(dens): # don't include atoms for which not calculated properly
					l = writePDBline(atm,dens)
					pdbOut.write(l+'\n')
		else:
			atmDic = self.combinedAtoms.getAvMetricPerAtmInRes(singleRes,metric,normType,dataset)
			for atm in self.combinedAtoms.atomList:
				if atm.basetype == singleRes:
					resNum = atm.residuenum # save the residue number now
					break
			for atm in self.combinedAtoms.atomList:
				if atm.residuenum == resNum:
					dens = atmDic[atm.atomtype]
					l = writePDBline(atm,dens)
					pdbOut.write(l+'\n')
		pdbOut.write('END')
		pdbOut.close()

	def visualiseDamSites(self, 
						  dataset  = 0, 
						  metric   = 'loss', 
						  software = 'pymol', 
						  size     = 1,
						  savePic  = True):

		# open coot/pymol to view top damage sites. size is the 
		# density metric scale factor used when visualising damage 
		# sites as spheres of vdw = size*{density metric} within pymol

		if software not in ('coot','pymol'):
			print 'Damage sites can be visualised in either "coot" or "pymol"'
			print 'Please make sure the paths to these programs are '+\
				  'correctly set up before proceeding'
			return
		try: self.damSitesPDB
		except AttributeError:
			return "Must run .writeDamSitesToFile() before damage sites can be read in {}".format(software)
		if software == 'coot':
			system('coot -pdb {} -pdb2 {}'.format(self.get1stDsetPDB(),
													 self.damSitesPDB[dataset]))
		else:
			# need to write script for pymol to run with
			damSitesTag = (self.damSitesPDB[dataset].split('/')[-1]).strip('.pdb')
			structureTag = self.initialPDB.strip('.pdb')
			scriptName = self.outputDir+'runPymolScript.pml'
			pymolScript = open(scriptName,'w')

			pymolStr = 'load {}\n'.format(self.get1stDsetPDB())+\
					   'load {}\n'.format(self.damSitesPDB[dataset])+\
					   'hide lines\nshow cartoon\n'+\
					   'set antialias, 1\n'+\
					   'set ray_trace_mode, 0\n'+\
					   'set cartoon_fancy_helices, 1\n'+\
					   'set cartoon_side_chain_helper, on\n'+\
					   'set ray_opaque_background, 0\n'+\
					   'show sphere, {}\n'.format(damSitesTag)+\
					   'alter {}, vdw=b*{}\n'.format(damSitesTag,size)+\
					   'rebuild\n'+\
					   'select nearDamage, {} w. 4 of {}\n'.format(structureTag,damSitesTag)+\
					   'show sticks, nearDamage\n'+\
					   'set stick_transparency, 0.6\n'+\
					   'color red, {}\n'.format(damSitesTag)+\
					   'color white, {}\n'.format(structureTag)

			if savePic is True:
				# DOESN'T WORK CURRENTLY
				pymolStr += 'ray\n png {}\nquit'.format(self.damSitesPDB[dataset].replace('.pdb','.png'))

			pymolScript.write(pymolStr)
			pymolScript.close()
			system('pymol -q {}'.format(scriptName))

	def tests_FurtherAnalysis(self,
							  endAfter     = False,
							  normType     ='Standard',
							  metric       = 'loss',
							  perform      = False,
							  metSigPlots  = False,
							  mainSideCmp  = True,
							  distToDiSulp = False,
							  tTests       = False,
							  TopOfType    = False,
							  rankings     = False,
							  sigKSstats   = False,
							  atmCorels    = False,
							  HbondCorr    = False):
	
		# plot line plots of statistics against dataset number and also
		# plot scatter plots to compare two statistics for each dataset
		# for each residue present within a structure. The currently
		# plotted statistics are mean, std dev, skewness and kurtosis.
		# This plot should NOT be output in a default run of program

		if perform is False:
			return

		if self.getNumDatasets() == 1:
			return

		if perform is False:
			return

		if metSigPlots is True:
			for t,u in zip(['','-sideOnly'],[False,True]):
				outDir = self.makeNewPlotSubdir(subdir = 'metricSignaturePlots{}/'.format(t))
				for i in range(self.getNumDatasets()):
					self.makeDistnPlots(densMet    = metric,
										normType   = normType,
									    plotSet    = 2,
									    calcKSstat = True,
									    calcADstat = True,
									    sideOnly   = u,
									    outputDir  = outDir,
									    dataset    = i,
									    requireAll = True)

		if mainSideCmp is True:
			outDir = self.makeNewPlotSubdir(subdir = 'side-main-compare-distn/')

			for d in range(self.getNumDatasets()):
				for s in [4,5]:
					self.makeDistnPlots(densMet    = metric,
								        normType   = normType,
								        plotSet    = s,
								        calcKSstat = True,
								   	 	calcADstat = True,
								        outputDir  = outDir,
								        dataset    = d)

		# COMMENTED CODE BELOW ONLY WORKS FOR MYROSINASE (BURMEISTER, 2000)
		# loc = '/Users/charlie/DPhil/YEAR2/JUN/setOccupanciesTo1/'
		# for bmet in ('Bdamage','Bfactor'):
		# 	for d,i in zip(range(self.getNumDatasets()),['f','g','h','i']):
		# 		self.combinedAtoms.compareDensToBdamChange(BDamageFile1  = loc + '1dwa_editted_refmac1Bdamage.txt',
		# 												   BDamageFile2  = loc + '1dw{}_editted_refmac1Bdamage.txt'.format(i),
		# 												   BdamChange    = True,
		# 												   resType       = 'TYR',
		# 												   atomType      = 'OH',
		# 												   Bmetric       = bmet,
		# 												   percentChange = True,
		# 												   densMet       = metric,
		# 												   normType      = normType,
		# 												   outputDir     = loc,
		# 												   dataset       = d)

		if distToDiSulp is True:
			print 'Finding distances to nearest disulphide bridges'
			self.combinedAtoms.scatterAtmsByDistToOtherAtms(outputDir = self.outputPlotDir)

		# COMMENTED CODE BELOW ONLY WORKS FOR TRAP (BURY ET AL., 2016)
		# for resNum in ['62']:
		# 	self.combinedAtoms.metricValsAtDistFromAtom(atomType  = 'OH',
		# 												resType   = 'TYR',
		# 												resNum    = resNum,
		# 												chainType = 'A',
		# 												metric    = metric,
		# 												normType  = normType,
		# 												outputDir = self.outputPlotDir)

		if TopOfType is True:
			t = ['TYR','OH']
			print 'Finding top atom of type: {}'.format('-'.join(t))
			self.combinedAtoms.findTopAtomOfType(normType = normType,
								  				 atomType  = t[1]
						  						 resType   = t[0])

		if tTests is True:
			print 'Performing t-test to compare metric correlation between atom types'
			for d in range(self.getNumDatasets()):
				stats = self.combinedAtoms.twoAtomTypeTtest(dataset   = d,
															atomType2 = ['OE1','OE2','OD1','OD2'],
															resType2  = ['GLU','GLU','ASP','ASP'])
				print 'Glu/Asp) Dataset: {}, stats: {}'.format(d,stats)

				stats = self.combinedAtoms.twoAtomTypeTtest(dataset   = d,
															atomType2 = ['OE1','OD1'],
															resType2  = ['GLN','ASN'])
				print 'Gln/Asn) Dataset: {}, stats: {}'.format(d,stats)

		res = ['CYS','TYR','TYR','GLU','GLU','GLU','GLU','ASP','ASP','ASP','ASP','GLN','GLN','ASN','ASN','SER','THR']
		atm = ['SG','OH','CZ','CD','OE1','OE2','CG','CG','CB','OD1','OD2','OE1','NE2','OD1','ND2','OG','OG1']
		if rankings is True:
			print 'Determining per atom-type rankings for selected atoms in structure'
			rank = self.combinedAtoms.getAtomtypeRanking(metric    = metric,
												         normType  = normType,
												         dataset   = 'all',
												         residue   = res,
												         atomtype  = atm,
												         printText = True)
			for key in rank.keys():
				print '{}: {}'.format(key,rank[key])

			for r,a in zip(res,atm):
				num = self.combinedAtoms.numAtmsWithMetricAboveLevel(dataset   = 'all',
																     metric    = metric,
																     normType  = normType,
																     threshold = 1,
																     firstTime = True,
																     atomType  = a,
																     resType   = r)
				print '{}-{}: {}'.format(r,a,num)

			if self.getNumDatasets() == 1:
				return

			outDir = self.makeNewPlotSubdir(subdir = 'TyrOH-ranking/')
			self.combinedAtoms.plotAtomtypeRankingWithDataset(outputDir = outDir,
															  percent   = False,
															  normType  = normType)

			for r,a in zip(res,atm):
				self.combinedAtoms.plotNumAtomsWithMetricAboveStructureWideMean(metric    = metric,
																				normType  = normType,
																				dataset   = 'all',
																				outputLoc = outDir,
																				atomType  = a,
																				resType   = r)

		if sigKSstats is True:
			print 'Calculating Kolmogorov-Smirnov statistics to compare damage '+\
				  'signatures between atom types in current structure'
			outDir = self.makeNewPlotSubdir(subdir = 'MetricSignature-statistics/')
			self.combinedAtoms.plotStatVsDataset(normType  = normType,
											     outputDir = outDir)

			self.combinedAtoms.plotStatVsStat(normType  = normType,
											  outputDir = outDir)

			if self.checkCalphasPresent(atomObjList = self.combinedAtoms) is False:
				return

			plotType = 'box'
			outDir = self.makeNewPlotSubdir(subdir = 'KS-stat-values/')

			# for ref in ('GLY','ALA'):
			# 	saveName = '{}plot_Metric-D{}_Normalisation-{}-KSstat_allResidues-rel-to-{}'.format(plotType,metric,normType,ref)
			# 	self.combinedAtoms.plotKSstatVsDataset(normType  = normType,
			# 									  	   outputDir = outDir,
			# 									  	   reference = ref,
			# 									  	   plotType  = plotType,
			# 									  	   saveName  = saveName)

			res = ['GLU','ASP','TYR','LYS','ILE','THR']
			ref = ['GLN','ASN','PHE','ARG','LEU','SER']
			saveName = '{}plot_Metric-D{}_Normalisation-{}-KSstat_keyResidues'.format(plotType,metric,normType)

			sideOnly = True
			if sideOnly is True:
				res += ['GLU']
				ref += ['ASP']
				saveName += '_sideChainsOnly'

			self.combinedAtoms.plotKSstatVsDataset(normType  = normType,
											  	   outputDir = outDir,
											  	   reference = ref,
											  	   residues  = res,
											  	   plotType  = plotType,
											  	   inclLine  = False,
											  	   sideOnly  = sideOnly,
											  	   saveName  = saveName)

		if atmCorels is True:
			print 'Determining correlation between atom types in same residues'
			outDir = self.makeNewPlotSubdir(subdir = 'atomtype-correlationPlots/')
			self.combinedAtoms.compareSensAtoms_RsquaredLinePlotWithDataset(normType  = normType,
							 											    outputDir = outDir,
							 											    plotType  = 'box')

		if HbondCorr is True:
			t = ['TYR','OH']
			print 'Determining whether a correlation exists between {} and '.format('-'.join(t))+\
				  'the presence of carboxylate salt bridge interactions'
			outDir = self.makeNewPlotSubdir(subdir = 'carboxylate-TyrOH-correlation/')
			self.combinedAtoms.plot_densMetSurroundAtmsCorrel(pdbName       = self.get1stDsetPDB(),
															  symmetrygroup = self.getSpaceGroup(),
															  restype       = t[0],
								  	   						  atomtype      = t[1],
															  normType      = normType,
															  outputDir     = outDir,
															  errorbars     = False,
															  printText     = False)

		if endAfter is True:
			self.endAfter()

	def plotLinePlot(self,
					 resType   = 'GLU',
					 errorBars = 'NONE',
					 atomType  = '',
					 chainType = 'A',
					 resNum    = '',
					 metric    = 'loss',
					 normType  = 'Standard',
					 figTitle  = '',
					 saveName  = ''):

		# create a line plot for Dloss density metric vs dataset number

		figName = self.combinedAtoms.graphMetric(atomType  = atomType,
											     resType   = resType,
											     chainType = chainType,
											     resiNum   = resNum,
											     densMet   = metric,
											     normType  = normType,
											     errorBars = errorBars,
											     outputDir = self.outputPlotDir,
											     saveFig   = True,
											     figTitle  = figTitle,
											     saveName  = saveName)
		return figName

	def metricBarplots(self):

		# plot barplots of damage metric for each susceptible residue type

		numDsets = self.getNumDatasets()
		for i in range(numDsets):
			for set in [1,2]:
				for b in ('Box','Bar'):
					self.combinedAtoms.susceptAtmComparisonBarplot(metric,normType,i,set,b)

	def makeDistnPlots(self,
					   densMet    = 'loss',
					   normType   = 'Standard',
					   plotType   = 'both',
					   plotSet    = 1,
					   calcKSstat = False,
					   calcADstat = False,
					   sideOnly   = False,
					   resSplit   = False,
					   outputDir  = '',
					   dataset    = 0,
					   requireAll = False):

		# retrieve the distribution for damage metric values for all 
		# atoms of specified residues (see 'resList' below), and then 
		# plot a kde plot for each.
		# If 'calcKSstat' is True, perform the 2-sample Kolmogorov-Smirnov 
		# test (will only calculate when two histograms plotted on 1 axis)

		if outputDir == '':
			outputDir = self.outputPlotDir

		title      = ''
		if plotSet == 1:
			resList = [['GLU','ASP','CYS','MET','TYR']]
			title   = 'Predicted susceptible residue types'
		elif plotSet == 2:
			resList = [['GLU','GLN'],
					   ['ASP','ASN'],
					   ['ILE','LEU'],
					   ['TYR','PHE'],
					   ['GLU','ASP'],
					   ['GLU','GLY'],
					   ['ASP','GLY'],
					   ['ALA','GLY'],
					   ['TYR','GLY'],
					   ['CYS','GLY'],
					   ['MET','GLY'],
					   ['ARG','LYS'],
					   ['A','GLY'],
					   ['G','GLY'],
					   ['U','GLY']]
		elif plotSet == 3:
			resList = [['DA','DC','DG','DT']]

		elif plotSet == 4:
			resList = ['all']
			resSplit  = True

		elif plotSet == 5:
			resList  = [['CYS'],['GLU'],['GLN'],['TYR'],
						['PHE'],['ASP'],['ASN'],['LYS'],
						['ARG'],['SER'],['GLY'],
						['DA'],['DC'],['DG'],['DT'],
						['A'],['C'],['G'],['U']]
			resSplit = True

		for resGroup in resList:
			k = ''.join(resGroup)
			failedPlots = {k : []}
			if len(resGroup) > 6 or plotSet == 4:
				plotType = 'kde'
			else:
				plotType = 'both'
			data = self.combinedAtoms.graphMetricDistn(metric     = densMet,
													   normType   = normType,
													   valType    = dataset,
													   plotType   = plotType,
													   resiType   = resGroup,
													   resSplit   = resSplit,
													   sideOnly   = sideOnly,
													   outputDir  = outputDir,
													   printText  = False,
													   plotTitle  = title,
													   calcKSstat = calcKSstat,
													   calcADstat = calcADstat,
													   requireAll = requireAll)	
			if data == {}:
				failedPlots[k] = True
			else:
				failedPlots[k] = False

		return failedPlots

	def sensAtomPlots(self,
					  metric   = 'loss',
					  normType = 'Standard'):

		# set of plots investigating damage progression 
		# for sensitive atoms within structure.

		pairs = {'susceptRes': [['GLU','CD','CG'],
							    ['GLU','CD','OE1'],
							    ['ASP','CG','CB'],
							    ['ASP','CG','OD1'],
							    ['TYR','OH','CZ'],
							    ['TYR','CZ','CE2'],
							    ['MET','SD','CG'],
							    ['MET','SD','CE'],
							    ['CYS','SG','CB'],
							    ['CYS','CA','CB']],
				'PHEvTYR'    : [['TYR','CZ','OH'],
					            ['TYR','CZ','CE2'],
					            ['TYR','CZ','CE1'],
							    ['TYR','CZ','CG'],
							    ['PHE','CZ','CE2'],
							    ['PHE','CZ','CE1'],
							    ['PHE','CZ','CG']],
				'TYR'        : [['TYR','CB','CG'],
					            ['TYR','CB','CD1'],
							    ['TYR','CB','CD2'],
							    ['TYR','CB','CE1'],
							    ['TYR','CB','CE2'],
							    ['TYR','CB','CZ'],
							    ['TYR','CB','OH']],
				'PHE'        : [['PHE','CB','CG'],
							    ['PHE','CB','CD1'],
							    ['PHE','CB','CD2'],
							    ['PHE','CB','CE1'],
							    ['PHE','CB','CE2'],
							    ['PHE','CB','CZ']]}

		# determine ratios of atoms within susceptible residue 
		# types and write output txt file and plot
		for met in ('distance','ratio'):

			for key in pairs.keys():
				self.combinedAtoms.findMetricRatioKeyResidues(metric   = 'loss',
															  normType = 'Standard',
															  rType    = met,
															  pairs    = pairs[key],
															  title    = key)
		atomtypes = [[['GLU','CD'],
					  	['ASP','CG'],
					  		['TYR','OH'],
					  			['CYS','SG'],
					  				['MET','SD']],
					 [['GLU','CD'],
					  	['ASP','CG'],
					  		['TYR','OH'],
					  			['CYS','SG'],
					  				['MET','SD'],
					 					['PHE','CZ'],
					  						['TYR','CZ']],
					 [['GLU','CG'],
					  	['ASP','CB'],
					  		['TYR','OH'],
					  			['CYS','CB'],
					  				['MET','CG'],
					  					['MET','CE'],
					  						['PHE','CZ'],
					  							['TYR','CZ']],
					 [['TYR','OH'],
					  	['PHE','CZ'],
					  		['TYR','CZ'],
					  			['PHE','CE1'],
					  				['TYR','CE1'],
					  					['PHE','CE2'],
					  						['TYR','CE2']]]

		for aTypes in atomtypes:
			self.combinedAtoms.plotSusceptibleAtoms(densMet  = metric,
													sus      = aTypes,
													normType = normType)

			self.combinedAtoms.plotSusceptibleAtoms(densMet   = metric,
													normType  = normType,
													errorbars = False,
													susAtms   = aTypes)

	def checkCalphasPresent(self,
							atomObjList = []):

		# check whether structure contains any Calpha 
		# protein backbone atoms within it

		return atomObjList.checkCalphaAtomsExist()

	def get1stDsetPDB(self):

		# retrieve name of first dataset pdb coordinate file

		pdbFile = self.inDir + self.initialPDB

		return pdbFile

	def getSpaceGroup(self):

		# parse the first dataset pdb file and retrieve the space group

		pdbFile = self.get1stDsetPDB()
		pdbin = open(pdbFile,'r')

		for line in pdbin.readlines():
			if line.split()[0] == 'CRYST1':
				self.spaceGroup = line[55:66].replace(' ','')
		pdbin.close()

		try: 
			self.spaceGroup
		except attributeError:
			err = 'Unable to find space group from file: {}'.format(pdbFile)
			self.logFile.writeToLog(str = err)
			return False
		return self.spaceGroup

	def fillerLine(self,
				   blank = False):

		# print a filler line to command line

		if blank is False:
			ln = '\n---------------------------------------------------------------'	
		else:
			ln = '\n'
		self.logFile.writeToLog(str = ln)

	def checkSeaborn(self):

		# force no plotting if seaborn library not found

		self.seabornFound = checkSns(printText = False,
									 logFile   = self.logFile)
		if self.seabornFound is False:
			self.plot = False

	def printNoSeabornWarning(self):

		# if seaborn library not found, print a warning to command line

		txt = '*** WARNING ***\n'+\
			  'Seaborn plotting library not found.\n'+\
			  'Some output plots could not be produced.\n'+\
			  'Use "pip install seaborn" to install package.'
		self.logFile.writeToLog(str = txt)


	def endAfter(self):

		# force quit of the class

		import sys
		sys.exit()
