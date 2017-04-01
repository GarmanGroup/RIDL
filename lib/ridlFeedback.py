from checkDependencies import checkDependencies
from os import path,makedirs,listdir,system
from time import gmtime, strftime
from shutil import move
import sys
import os

class provideFeedback(object):

	def __init__(self,
				 standardFeedback = True,
				 csvOnly          = False,
			     writeCsvs        = False,
			     writeSumFile     = False,
			     writeTopSites    = False,
			     plotHeatMaps     = False,
			     atmsObjs         = [],
			     outputDir        = './',
			     outputPlotDir    = './',
			     csvExtent        = 'simple',
			     plotGraphs       = True,
			     logFile          = [],
			     pklSeries        = '',
			     doses            = [],
			     pdbNames         = [],
			     inputDir         = './',
			     initialPDB       = 'untitled.pdb',
			     inclFCmetrics    = False,
			     autoRun          = True):

		# create series of output feedback files and graphs

		self.atmsObjs      = atmsObjs
		self.outputDir     = outputDir
		self.outputPlotDir = outputPlotDir
		self.plot          = plotGraphs
		self.logFile       = logFile
		self.pklSeries     = pklSeries
		self.plotHeatMaps  = False
		self.doses         = doses
		self.pdbNames      = pdbNames
		self.inDir         = inputDir
		self.initialPDB    = initialPDB
		self.inclFCmetrics = inclFCmetrics
		self.writeCsvs     = writeCsvs
		self.writeSumFile  = writeSumFile
		self.writeTopSites = writeTopSites
		self.csvExtent     = csvExtent

		if plotGraphs:
			if not os.path.exists(self.outputPlotDir):
				os.mkdir(self.outputPlotDir)

		# define a 'standard' feedback for the program
		if standardFeedback:
			self.writeCsvs     = True
			self.writeSumFile  = True
			self.writeTopSites = True

		if csvOnly:
			self.writeCsvs     = True
			self.writeSumFile  = False
			self.writeTopSites = False
			self.plotHeatMaps  = False

		if autoRun:

			# LEAVE THIS SECTION COMMENTED FOR A STANDARD RUN OF THE CODE (WORKING PROGRESS STILL)
			# atmsObjs.findProbAboveAvDam()

			# for n in [1,2,3]:
			# 	csvOut = open(self.outputDir+'CalphaDloss-probNeighbourHigh-{}.csv'.format(n),'w')
			# 	print '\nFor Calpha Dloss:'
			# 	for dist in [4,5,6,7,8,9,10]:
			# 		print '>>> Distance: {} Angstroms'.format(dist)
			# 		prbList = []
			# 		for d in range(self.getNumDatasets()):
			# 			prb = atmsObjs.findProbHighNeighbourGivenHighAtom(densMet = 'loss',
			# 															  normType = 'Calpha normalised',
			# 															  dataset  = d,
			# 															  distance = dist,
			# 															  criteria = '{}std'.format(n))
			# 			prbList.append(str(prb))
			# 		ln = '{},{}'.format(dist,','.join(prbList))
			# 		csvOut.write(ln+'\n')
			# 	csvOut.close()
			# for n in [1,2,3]:
			# 	csvOut = open(self.outputDir+'Bfactor-probNeighbourHigh-{}.csv'.format(n),'w')
			# 	print '\nFor Bfactor:'
			# 	for dist in [4,5,6,7,8,9,10]:
			# 		print '>>> Distance: {} Angstroms'.format(dist)
			# 		prbList = []
			# 		for d in range(self.getNumDatasets()):
			# 			prb = atmsObjs.findProbHighNeighbourGivenHighAtom(densMet = 'Bfactor',
			# 															  normType = 'Standard',
			# 															  dataset  = d,
			# 															  distance = dist,
			# 															  criteria = '{}std'.format(n))
			# 			prbList.append(str(prb))
			# 		ln = '{},{}'.format(dist,','.join(prbList))
			# 		csvOut.write(ln+'\n')
			# 	csvOut.close()

			self.run()

	def run(self):

		# main procedure for writing output files

		# no plotting if seaborn not found
		self.checkSeaborn()

		for norm in ('Standard','Calpha normalised'):

			if norm == 'Calpha normalised':
				if self.checkCalphasPresent(atomObjList = self.atmsObjs) is False:
					continue

		if self.writeCsvs:
			self.writeCsvFiles()

		# provide summary html file for Dloss metric per-dataset
		if self.writeSumFile:
			self.summaryFile() 
			# self.summaryFile(metric = 'density_weighted_mean_negOnly') 
		
		if self.writeTopSites:
			self.writeDamSitesToFile()

		# create heatmap plots (large files) if requested
		if self.plot:
			if self.plotHeatMaps:

				subdir = 'metric_heatmap/'
				outDir = self.makeNewPlotSubdir(subdir = subdir)

				for norm in ('Standard','Calpha normalised'):
					if norm == 'Calpha normalised':
						if self.checkCalphasPresent(atomObjList = self.atmsObjs) is False:
							continue
					self.atmsObjs.densityMetricHeatMap(saveFig   = True,
													   metric    = 'loss',
													   normType  = norm,
													   outputDir = outDir)	

	def writeCsvFiles(self,
					  moveCsv     = True,
					  inclGroupby = True,
					  inclGainMet = False,
					  inclMeanMet = False,
					  numDP       = 2):

		# write atom numbers and density metrics to simple 
		# csv files,one for each density metric separately.
		# inclGroupby = True also writes csv files for atoms
		# grouped by residue type and atom type

		CalphasPresent = self.checkCalphasPresent(atomObjList = self.atmsObjs)
		self.fillerLine()
		ln = 'Writing .csv file for per-atom density metric:'
		self.logFile.writeToLog(str = ln)

		if self.csvExtent == 'simple':
			n = 'Standard'
			metrics = [['loss',n],
					   ['mean',n]]

			if inclMeanMet:
				metrics += [['mean',n]]

			if inclGainMet:
				metrics += [['gain',n]]

			if self.inclFCmetrics:

				metrics += [['density_weighted_loss',n],
						    ['density_weighted_mean_negOnly',n]]	

				if inclGainMet:
					metrics += [['density_weighted_gain',n]]

				if inclMeanMet:
					metrics += [['density_weighted_mean',n]]

			if CalphasPresent:
				n = 'Calpha normalised'
				metrics += [['loss',n]]

				if self.inclFCmetrics:
					metrics += [['density_weighted_loss',n],
							    ['density_weighted_mean_negOnly',n]]
		else:
			metrics = self.atmsObjs.getDensMetrics()

		for densMet in metrics:
			ln = '\tmetric: {}, normalisation: {}'.format(*densMet)
			self.logFile.writeToLog(str = ln)
			self.atmsObjs.writeMetric2File(where    = self.outputDir,
										   metric   = densMet[0],
										   normType = densMet[1],
										   numDP    = numDP)

		if inclGroupby:
			for groupBy in ('residue','atomtype'):
				self.atmsObjs.writeMetric2File(where   = self.outputDir,
											   groupBy = groupBy)
				if CalphasPresent:
					self.atmsObjs.writeMetric2File(where    = self.outputDir,
											  	   groupBy  = groupBy,
										           normType = 'Calpha normalised')

		# make csvFiles dir and move all generated csv files to this
		if moveCsv:
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

	def summaryFile(self,
					fileType = 'html',
					metric   = 'loss',
					normType = 'Standard'):

		# write a summary output file (only html currently available)

		ln = 'Writing {} summary output file for metric: {}, normalisation: {}'.format(fileType,metric,normType)
		self.logFile.writeToLog(str = ln)

		if fileType == 'html':

			self.summaryHTML(primaryMetric = metric)
		else:
			print 'Unknown file format. Only currently supported format is "html"'

	def summaryHTML(self, 
					primaryMetric = 'loss',
					primaryNorm   = 'Calpha normalised'):

		# produce a selection of per-dataset summary statistics

		calphaNorm = 'C<sub>&#945</sub>-normalised'
		norms      = ['Standard','Calpha normalised']
		plotNorm   = primaryNorm

		numDsets = self.getNumDatasets()
		summaryFile = open('{}summaryFile.html'.format(self.outputDir),'w')
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
					 '<h1>RIDL summary file</h1>\n'+\
					 'Created: {}<br>\n'.format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))+\
					 'Summary information derived from {}<br>\n'.format(self.pklSeries)+\
					 'Email charles.bury@dtc.ox.ac.uk for queries<br>\n'+\
					 'Number of datasets reported in file: {}<br>\n'.format(numDsets)

		summaryFile.write(bodyString)

		###########################################################################
		# Plot of top n damage sites per dataset by residue type created here

		subdir  = 'topDamSites/'
		plotDir = self.makeNewPlotSubdir(subdir = subdir)
		figName = self.atmsObjs.getTopAtomsStackedBarplot(outputDir = plotDir,
														  metric    = primaryMetric)

		# provide some links to useful output files
		bodyString = '<h3>Links to useful output files</h3><ul>\n'

		for m,norm in zip([primaryMetric]*2,norms):
			bodyString += '<li><a href = "csvFiles/{}/{}-'.format(norm.replace(' ','-'),m)+\
					 	  '{}.csv">{} D<sub>{}</sub> csv file</a></li>\n'.format(norm.replace(' ',''),norm,m)
		
		bodyString += '<li><a href = "plots/{}{}">Top 25 damage sites per residue/nucleotide type</a></li>'.format(subdir,figName.split('/')[-1])

		# create heatmap plots (large files) if requested
		if self.plotHeatMaps:
			for m,norm in zip([primaryMetric]*2,norms):
				bodyString += '<li><a href = "plots/metric_heatmap/heatmap_metric-{}_normalisation-'.format(m)+\
						  	  '{}.svg">{} D<sub>{}</sub> per-atom heat map</a></li></ul>'.format(norm.replace(' ',''),norm,m)
		else: 
			bodyString += '</ul>'

		summaryFile.write(bodyString)

		# find top overall damaged atoms and plot line plots 
		# versus dataset. Only plot if more than 1 high dose
		# dataset is present.
		if numDsets == 1:
			pass

		else:

			###########################################################################
			# Line plots of metric values versus dose provided for top and bottom 10
			# atoms located in structure (top and bottom determined over whole dataset range)

			subdir  = 'topDamSites/'
			plotDir = self.makeNewPlotSubdir(subdir = subdir)

			numDamSites = 10
			figCalls    = []

			for metric, norm in zip([primaryMetric]*2,norms):

				topAtoms = self.atmsObjs.getTopNAtoms(dataset  = 'all', 
													  n        = numDamSites)

				keys = ['chain','num','res','atm']
				info = {k : [] for k in keys}

				for atm in topAtoms:
					aInfo = atm.split('-')
					for i,j in zip(keys,aInfo):
						info[i].append(j)

				figName = self.atmsObjs.graphMetric(atomType  = info['atm'],
													resType   = info['res'],
													chainType = info['chain'],
													resiNum   = info['num'],
													densMet   = metric,
													normType  = norm,
													outputDir = self.outputPlotDir+subdir,
													saveFig   = True,
													figTitle  = 'Top {} D{} damage sites'.format(numDamSites,metric),
													saveName  = 'Lineplot_Metric-D{}_Normalisation-{}_topDamSites'.format(metric,norm))

				figCalls.append('<img class="img-responsive" src="plots/{}{}">'.format(subdir,figName.split('/')[-1]))

			# make distn plots over entire structure to indicate
			# the role of Calpha normalisation
			subdir = 'metricDistn_allAtoms/'
			outDir = self.makeNewPlotSubdir(subdir = subdir)

			for n in ['Standard','Calpha normalised']:

				data,saveName = self.makeDistnPlots(densMet   = primaryMetric,
													normType  = n,
											     	plotSet   = 4,
											        outputDir = outDir,
											        dataset   = 'all')

				figCalls.append('<img class="img-responsive" src="plots/{}{}">'.format(subdir,saveName.split('/')[-1]))

			info = '<h3>Metric distribution</h3>\n'+\
				   '<div class = "row">\n'+\
			  	   '<div class = "col-sm-6">{}</div>\n'.format(figCalls[2])+\
			       '<div class = "col-sm-6">{}</div>\n'.format(figCalls[3])+\
			       '</div>\n'
			summaryFile.write(info)

			info = '<h3>Top 10 density loss sites with dose:</h3>\n'+\
				   '<div class = "row">\n'+\
			  	   '<div class = "col-sm-6">{}</div>\n'.format(figCalls[0])+\
			       '<div class = "col-sm-6">{}</div>\n'.format(figCalls[1])+\
			       '</div>\n'
			summaryFile.write(info)

		# make a set of tabs for easy navigation to each dataset
		summaryFile.write('<h2>Per-dataset analysis</h2>\n')
		tabs = '<ul class="nav nav-tabs">'+\
			   '<li class="active"><a data-toggle="tab" href="#dset-tab0">Dataset 1</a></li>'

		for i in range(1,numDsets):
			tabs += '<li><a data-toggle="tab" href="#dset-tab{}">Dataset {}</a></li>'.format(i,i+1)
		tabs += '</ul>'

		summaryFile.write(tabs)
		summaryFile.write('<div class="tab-content">')

		###########################################################################
		# Per-dataset breakdown analysis begins at this point

		for i in range(numDsets):
			if i != 0:
				str = '<div id="dset-tab{}" class="tab-pane fade">'.format(i)
			else:
				str = '<div id="dset-tab{}" class="tab-pane fade in active">'.format(i)
			summaryFile.write(str + '<hr>')

			# start block of collapsible panels here
			summaryFile.write('<div class = "panel-group" id = "datasetinfo{}">'.format(i))

			###########################################################################
			# General information regarding current dataset provided here

			t = 'Dataset info'
			c = 'Number in series : {}<br>\n'.format(i+1)+\
				'Dose (MGy)       : {}<br>\n'.format(self.doses[i])+\
				'Number of atoms  : {}<br>\n'.format(self.atmsObjs.getNumAtoms())+\
				'Fourier diff map : <a href ="../RIDL-maps/{}_density.map">Download</a><br>\n'.format(self.pdbNames[i])
			
			CAweights = self.atmsObjs.retrieveCalphaWeight(metric = primaryMetric)
			c += 'Calpha weight for current dataset: {}<br>\n'.format(round(CAweights.weight[metric][i],3))
			self.writeHtmlDropDownPanel(title   = t,
										content = c,
								        dataset = i,
								        sumFile = summaryFile)

			###########################################################################
			# Create distribution plots for metric values over whole structure
			# and key residue types that are typically damaged

			t = 'D<sub>{}</sub> Distribution Plots'.format(primaryMetric)

			subdir = 'metricDistn_allAtoms/'
			outDir = self.makeNewPlotSubdir(subdir = subdir)

			failedPlots = self.makeDistnPlots(densMet   = primaryMetric,
										      normType  = plotNorm,
										      plotSet   = 4,
										      outputDir = outDir,
										      dataset   = i)

			figCall = '<img class="img-responsive" src="plots/{}DistnPlot_Residues-all_Metric-D{}_Normalisation-{}_Dataset-{}.svg">'.format(subdir,primaryMetric,plotNorm.replace(' ',''),i)
			figInfo = '<div class = "row">\n'+\
		  	          '<div class = "col-sm-6">{}</div>\n'.format(figCall)

			params = [['Residues',1,'GLUASPCYSMETTYR','residue'],
					  ['nucleotides',3,'DADCDGDT','nucleotide']]

			for paramSet in params:
				subdir = 'metricDistn_key{}/'.format(paramSet[0])
				outDir = self.makeNewPlotSubdir(subdir = subdir)

				failedPlots,saveName = self.makeDistnPlots(densMet   = primaryMetric,
													       normType  = plotNorm,
													       plotSet   = paramSet[1],
													       outputDir = outDir,
													       dataset   = i)

				if not failedPlots[paramSet[2]]:
					figCall = '<img class="img-responsive" src="plots/{}{}"><br>'.format(subdir,saveName.split('/')[-1])
					figInfo += '<div class = "col-sm-6">{}</div>\n'.format(figCall)

			figInfo += '</div>\n'

			self.writeHtmlDropDownPanel(title   = t,
										content = figInfo,
								        dataset = i,
								        sumFile = summaryFile)	

			###########################################################################
			# Z score plot for current dataset of number of atoms with metric above 
			# number of standard deviations of mean value given here

			# subdir = 'zScorePlots/'
			# outDir = self.makeNewPlotSubdir(subdir = subdir)

			# figName = self.atmsObjs.plotNumAtomsWithMetricAboveStructureWideMean(metric   = primaryMetric,
			# 																	   normType = plotNorm,
			# 																       dataset  = i,
			# 																	   outputLoc = outDir)

			# t = '# atoms with D<sub>{}</sub> metric above N std dev of structure-wide mean'.format(primaryMetric)
			# c = '<img class="img-responsive" src="plots/{}{}">'.format(subdir,figName.split('/')[-1])
			# self.writeHtmlDropDownPanel(title   = t,
			# 							content = c,
			# 					        dataset = i,
			# 					        sumFile = summaryFile)

			# retrieve top N atoms and ranked by specified metric
			# statsOut = self.atmsObjs.getTopNAtomsString(dataset  = i,
			# 											metric   = primaryMetric,
			# 										    n        = 25)
			# t = 'Top hits ranked by D<sub>{}</sub> metric'.format(primaryMetric)
			# c = 'Top hits ranked by D<sub>{}</sub> metric.<br>D<sub>mean</sub> and D<sub>gain</sub> are mean '.format(primaryMetric)+\
			#     'and maximum voxel difference density values, respectively, '+\
			#     'assigned within a local region around each atom.<br>Proximity (from 0 '+\
			#     'to 1) is a measure of the closeness of the voxel exhibiting the '+\
			#     'maximum density loss D<sub>loss</sub> value from the specified atom (higher '+\
			#     'values indicate smaller distances):<br><br>\n'
			# c += self.convertPlainTxtTable2html(statsOut, width = '50%', splitBy = ',')
			# self.writeHtmlDropDownPanel(title   = t,
			# 							content = c,
			# 					        dataset = i,
			# 					        sumFile = summaryFile)


			###########################################################################
			# Determination of top N=25 damage sites for current dataset given here

			subdir = 'topNatoms/'
			outDir = self.makeNewPlotSubdir(subdir = subdir)

			# choose the order to plot the metrics for the dot plot
			# (including normalisation types). atoms will be ordered
			# by the first metric in the list

			if self.inclFCmetrics:
				metsToPlot = [['density_weighted_mean_negOnly','Calpha normalised'],
							  ['loss','Calpha normalised'],
							  ['density_weighted_loss','Calpha normalised'],
							  ['Bfactor','Standard']]
			else:
				metsToPlot = [['loss','Calpha normalised'],
							  ['Bfactor','Standard']]

			firstToPlot = [primaryMetric,plotNorm]
			if firstToPlot in metsToPlot:
				metsToPlot.remove(firstToPlot)
			metsToPlot = [firstToPlot] + metsToPlot
			metsToPlot = [[m[j] for m in metsToPlot] for j in range(2)]

			saveName = self.atmsObjs.getTopNAtomsDotPlot(dataset   = i,
														 metrics   = metsToPlot[0],
														 normTypes = metsToPlot[1],
														 numHits   = 25)

			t = 'Top 25 damage sites'
			c = '<img class="img-responsive" src="plots/{}{}">'.format(subdir,saveName.split('/')[-1])	

			self.writeHtmlDropDownPanel(title   = t,
										content = c,
								        dataset = i,
								        sumFile = summaryFile)	

			###########################################################################
			# statistics breakdown (per-atom, per-chain, per-residue, full structure) starts here

			c = 'Key:\n'+\
				'<ul><li>mean: average of metric calculated over all atoms of a specified type</li>\n'+\
				'<li>std: standard deviation of metric calculated over all atoms of a specified type</li>\n'+\
				'<li>#atoms: total number of atoms of a specified type</li>\n'+\
				'<li>outliers: assuming a symmetric distn around the mode, number of atoms that fall outside this domain</li>\n'+\
				'<li>skew: skewness of metric distribution for atoms of specified type</li>\n'+\
				'<li>kurtosis: kurtosis of metric distribution for atoms of specified type</li>\n</ul>'

			c += '<ul class="nav nav-tabs">'+\
				 '<li class="active"><a data-toggle="tab" href="#byatom{}">By Atom</a></li>'.format(i)+\
				 '<li><a data-toggle="tab" href="#byresidue{}">By Residue</a></li>'.format(i)+\
			     '<li><a data-toggle="tab" href="#bychain{}">By Chain</a></li>'.format(i)+\
			     '<li><a data-toggle="tab" href="#bystructure{}">Full Structure</a></li>'.format(i)+\
				 '</ul>'

			c += '<div class="tab-content">'

			# per-atom statistics
			c += '<div id="byatom{}" class="tab-pane fade in active">'.format(i)
			statsOut = self.atmsObjs.getPerAtmtypeStats(metric   = primaryMetric,
														normType = 'Calpha normalised',
														dataset  = i)
			c += '{} D<sub>{}</sub> metric ranked by mean value:<br><br>\n'.format(calphaNorm,primaryMetric)
			c += self.convertPlainTxtTable2html(statsOut[0], width = '60%')
			c += '</div>'

			# per-residue statistics
			c += '<div id="byresidue{}" class="tab-pane fade">'.format(i)
			statsOut = self.atmsObjs.getPerResidueStats(metric   = primaryMetric,
													    normType = 'Calpha normalised',
													    dataset  = i)
			c += '{} D<sub>{}</sub> metric ranked by mean value:<br><br>\n'.format(calphaNorm,primaryMetric)
			c += self.convertPlainTxtTable2html(statsOut[0], width = '60%')
			c += '</div>'

			# per-chain statistics
			c += '<div id="bychain{}" class="tab-pane fade">'.format(i)
			statsOut = self.atmsObjs.getPerChainStats(metric   = metric,
													  normType = 'Calpha normalised',
													  dataset  = i,
													  n        = 'all')
			c += '{} D<sub>{}</sub> metric ranked by mean value:<br><br>\n'.format(calphaNorm,primaryMetric)
			c += self.convertPlainTxtTable2html(statsOut[0], width = '60%')
			c += '</div>'
			
			# full-structure statistics
			c += '<div id="bystructure{}" class="tab-pane fade">'.format(i)
			statsOut = self.atmsObjs.getStructureStats(metric   = primaryMetric,
													  normType = 'Calpha normalised',
													  dataset  = i)
			c += '{} D<sub>{}</sub> metric ranked by mean value:<br><br>\n'.format(calphaNorm,primaryMetric)
			c += self.convertPlainTxtTable2html(statsOut[0], width = '60%')
			c += '</div>'
			c += '</div>'

			self.writeHtmlDropDownPanel(title   = 'Statistics Breakdown',
										content = c,
								        dataset = i,
								        sumFile = summaryFile)


			###########################################################################
			# Detection of any atoms that behaviour unusually (unlike rest of that atom
			# type) is presented here

			infoString = 'Atoms with unusually <font color="red">high</font> '+\
				  		 'or <font color="blue">low</font> D<sub>{}</sub> values '.format(primaryMetric)+\
				  		 'relative to the mean D<sub>{}</sub> value for that specific '.format(primaryMetric)+\
				  		 'atom type will be reported.<br><br>'

			lastInfoString = ''
			suspAtomLens   = []
			notDoneBefore  = True

			for t in [6,5,4,3]:
				suspAtoms,highOrLow = self.atmsObjs.detectSuspiciousAtoms(dataset   = i,
						  										 		  metric    = primaryMetric,
						  										 		  threshold = t)

				tmpInfoString = '{} atoms found with unusually high/low D{} '.format(len(suspAtoms),primaryMetric)+\
							    'values compared to average of that atom type '+\
							    '(over {} standard deviations from average)'.format(t)

				if len(suspAtoms) > 0:

					if notDoneBefore:
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
			self.writeHtmlDropDownPanel(title   = t,
										content = c,
										dataset = i,
										sumFile = summaryFile)

			# end block of collapsible panels here
			summaryFile.write('</div>')
			summaryFile.write('</div>')

		summaryFile.write('<div class="tab-content">')


		endString = '</div>\n'+\
					'</body>\n'+\
					'</html>'

		summaryFile.write(endString)
		summaryFile.close()

	def writeHtmlDropDownPanel(self,
							   title      = 'untitled',
							   header     = 'h3',
							   content    = 'no content',
							   panelColor = 'alternate',
							   dataset    = 0,
							   sumFile    = 'untitled.txt'):

		# write the info for a html collapsible panel

		try:
			self.panelIndex
		except AttributeError:
			self.panelIndex = 1

		if panelColor == 'alternate':
			if self.panelIndex%2 == 0:
				panelColor = 'link'
			else:
				panelColor = 'info'

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

		if sumFile != '':
			sumFile.write(txt)
		else:
			return txt

	def convertPlainTxtTable2html(self,
								  plainText      = '',
								  numHeaderLines = 1,
								  width          = '40%',
								  border         = '1',
								  splitBy        = '\t\t'):

		# convert a plain text table of values into a html format table

		htmlTable = '<table class="table table-striped">\n'.format(border,width)
		for i, l in enumerate(plainText.split('\n')):
			if i < numHeaderLines:
				newStr = '<tr><th>'+'</th><th>'.join(l.split(splitBy))+'</th></tr>'
			else:
				newStr = '<tr><td>'+'</td><td>'.join(l.split(splitBy))+'</td></tr>'
			htmlTable += newStr+'\n'
		htmlTable += '</table>\n'
		return htmlTable

	def getNumDatasets(self):

		# get number of datasets in damage series

		numDsets = self.atmsObjs.atomList[0].getNumDatasets()
		return numDsets

	def writeDamSitesToFile(self, 
							metric      = 'loss', 
							normType    = 'Standard', 
							numDamSites = 25):

		# write top damage sites to .pdb file for each dataset
		# 'numDamSites' takes 'all' or integer

		pdbTemplate = self.get1stDsetPDB()
		if not os.path.exists(pdbTemplate):
			print 'Warning: Could not coordinate file "{}"'.format(pdbTemplate)
			print '--> Not write top damage sites to .pdb file..'
			return 

		self.damSitesPDB = []
		for i in range(self.getNumDatasets()): 
			damPDB = self.atmsObjs.getTopNAtomsPDBfile(metric   = metric,
													   normType = normType,
													   dataset  = 'all',
													   n        = numDamSites,
													   pdbFile  = pdbTemplate)

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
			self.atmsObjs.calcAdditionalMetrics(metric   = metric,
												normType = normType)

		pdbIn = open(self.get1stDsetPDB(), 'r')
		fileOut = self.outputDir+self.initialPDB.replace('.pdb','')+\
				  '_{}D{}_{}.pdb'.format(normType.replace(" ",""),metric,dataset)
		if singleRes != '':
			fileOut = fileOut.replace('.pdb','-{}.pdb'.format(singleRes))

		pdbOut = open(fileOut,'w')
		pdbOut.write('REMARK\tBfactor column replaced by {} D{} metric values\n'.format(normType,metric))
		for l in pdbIn.readlines():
			if l.split()[0] in ('CRYST1','SCALE1','SCALE2','SCALE3'):
				pdbOut.write(l)
			elif 'ATOM' in l.split()[0]:
				break
		pdbIn.close()

		if singleRes == '':
			for atm in self.atmsObjs.atomList:
				dens = atm.densMetric[metric][normType]['values'][dataset]
				if not isnan(dens): # don't include atoms for which not calculated properly
					l = writePDBline(atm,dens)
					pdbOut.write(l+'\n')
		else:
			atmDic = self.atmsObjs.getAvMetricPerAtmInRes(singleRes,metric,normType,dataset)
			for atm in self.atmsObjs.atomList:
				if atm.basetype == singleRes:
					resNum = atm.residuenum # save the residue number now
					break
			for atm in self.atmsObjs.atomList:
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
			damSitesTag = (self.damSitesPDB[dataset].split('/')[-1]).replace('.pdb','')
			structureTag = self.initialPDB.replace('.pdb','')
			scriptName = self.outputDir + 'runPymolScript.pml'
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

			if savePic:
				pass
				# DOESN'T WORK CURRENTLY
				pymolStr += 'ray\n png {}\nquit'.format(self.damSitesPDB[dataset].replace('.pdb','.png'))

			pymolScript.write(pymolStr)
			pymolScript.close()
			system('pymol -q {}'.format(scriptName))

	def metricBarplots(self):

		# plot barplots of damage metric for each susceptible residue type

		numDsets = self.getNumDatasets()
		for i in range(numDsets):
			for set in [1,2]:
				for b in ('Box','Bar'):
					self.atmsObjs.susceptAtmComparisonBarplot(metric,normType,i,set,b)

	def makeDistnPlots(self,
					   densMet    = 'loss',
					   normType   = 'Standard',
					   plotType   = 'both',
					   inclTitle  = False,
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
			data,saveName = self.atmsObjs.graphMetricDistn(metric     = densMet,
														   normType   = normType,
														   valType    = dataset,
														   plotType   = plotType,
														   resiType   = resGroup,
														   resSplit   = resSplit,
														   sideOnly   = sideOnly,
														   outputDir  = outputDir,
														   printText  = False,
														   plotTitle  = title,
														   inclTitle  = inclTitle,
														   calcKSstat = calcKSstat,
														   calcADstat = calcADstat,
														   requireAll = requireAll)	
			if data == {}:
				failedPlots[k] = True
			else:
				failedPlots[k] = False

		return failedPlots,saveName

	def checkCalphasPresent(self,
							atomObjList = []):

		# check whether structure contains any Calpha 
		# protein backbone atoms within it

		return atomObjList.checkCalphaAtomsExist()

	def get1stDsetPDB(self):

		# retrieve name of first dataset pdb coordinate file

		pdbFile = self.inDir + str(self.initialPDB[0])

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

		if not blank:
			ln = '\n---------------------------------------------------------------'	
		else:
			ln = '\n'
		self.logFile.writeToLog(str = ln)

	def checkSeaborn(self):

		# force no plotting if seaborn library not found

		c = checkDependencies()
		if not c.checkSeaborn(logFile = self.logFile):
			self.plot = False
			self.printNoSeabornWarning()

	def printNoSeabornWarning(self):

		# if seaborn library not found, print a warning to command line

		txt = '*** WARNING ***\n'+\
			  'Seaborn plotting library not found.\n'+\
			  'Some output plots could not be produced.\n'+\
			  'Use "pip install seaborn" to install package.'
		self.logFile.writeToLog(str = txt)

	def makeOutputDir(self,
					  dirName   = './'):

		# if the above sub directory does not exist, make it

		if not path.exists(dirName):
			makedirs(dirName)
			ln = 'New sub directory "{}" '.format(dirName.replace(self.outputDir,''))+\
				 'created to contain output files'
			self.logFile.writeToLog(str = ln)

	def makeNewPlotSubdir(self,
						  subdir = 'untitled/'):

		# make a new sub directory for plots within plots/ subdir
		
		outDir = self.outputPlotDir + subdir
		self.makeOutputDir(dirName = outDir)

		return outDir

	def endAfter(self):

		# force quit of the class

		sys.exit()



