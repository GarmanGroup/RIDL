from checkSeabornPresent import checkSeabornPresent as checkSns
from os import path,makedirs,listdir,system
from time import gmtime, strftime
from shutil import move
import sys

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
			     autoRun          = True):

		# create series of output feedback files and graphs

		self.atmsObjs      = atmsObjs
		self.outputDir     = outputDir
		self.outputPlotDir = outputPlotDir
		self.plot          = plotGraphs
		self.logFile       = logFile
		self.pklSeries     = pklSeries
		self.plotHeatMaps  = plotHeatMaps
		self.doses         = doses
		self.pdbNames      = pdbNames
		self.inDir         = inputDir
		self.initialPDB    = initialPDB
		self.writeCsvs     = writeCsvs
		self.writeSumFile  = writeSumFile
		self.writeTopSites = writeTopSites
		self.csvExtent     = csvExtent

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
			self.run()

	def run(self):

		# main procedure for writing output files

		# no plotting if seaborn not found
		self.checkSeaborn()
		if not self.seabornFound:
			self.printNoSeabornWarning()

		for norm in ('Standard','Calpha normalised'):

			if norm == 'Calpha normalised':
				if self.checkCalphasPresent(atomObjList = self.atmsObjs) is False:
					continue

		if self.writeCsvs:
			self.writeCsvFiles()

		# provide summary html file for Dloss metric per-dataset
		if self.writeSumFile:
			self.summaryFile(normType = 'Standard') 
			if self.checkCalphasPresent(atomObjList = self.atmsObjs) is True:
				self.summaryFile(normType = 'Calpha normalised') 
		
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
					  inclGroupby = True):

		# write atom numbers and density metrics to simple 
		# csv files,one for each density metric separately.
		# inclGroupby = True also writes csv files for atoms
		# grouped by residue type and atom type

		CalphasPresent = self.checkCalphasPresent(atomObjList = self.atmsObjs)
		self.fillerLine()
		ln = 'Writing .csv file for per-atom density metric:'
		self.logFile.writeToLog(str = ln)

		if self.csvExtent == 'simple':

			metrics = [['loss','Standard'],
					   ['loss','reliability']]

			if CalphasPresent:
				metrics += [['loss','Calpha normalised']]
		else:
			metrics = self.atmsObjs.getDensMetrics()

		for densMet in metrics:
			ln = '\tmetric: {}, normalisation: {}'.format(*densMet)
			self.logFile.writeToLog(str = ln)
			self.atmsObjs.writeMetric2File(where    = self.outputDir,
										   metric   = densMet[0],
										   normType = densMet[1])

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
		figName = self.atmsObjs.getTopAtomsStackedBarplot(outputDir = self.outputPlotDir)

		# provide some links to useful output files
		bodyString = '<h3>Links to useful output files</h3>\n'+\
					 '<ul><li><a href = "csvFiles/{}/{}-'.format(normType.replace(' ','-'),metric)+\
					 '{}.csv">{} D<sub>{}</sub> csv file</a></li>\n'.format(normType.replace(' ',''),normType,metric)+\
					 '<li><a href = "plots/{}">Top 25 damage sites per residue/nucleotide type</a></li>'.format(figName.split('/')[-1])

		# create heatmap plots (large files) if requested
		if self.plotHeatMaps is True:
			bodyString += '<li><a href = "plots/metric_heatmap/heatmap_metric-{}_normalisation-'.format(metric)+\
						  '{}.svg">{} D<sub>{}</sub> per-atom heat map</a></li></ul>'.format(normType.replace(' ',''),normType,metric)
		else: 
			bodyString += '</ul>'

		summaryFile.write(bodyString)

		# find top overall damaged atoms and plot line plots 
		# versus dataset. Only plot if more than 1 high dose
		# dataset is present.
		if self.atmsObjs.getNumDatasets() == 1:
			pass
		elif includeGraphs is True:
			numDamSites = 10

			figCalls = []
			for t in ['Top','Bottom']:

				topAtoms = self.atmsObjs.getTopNAtoms(dataset  = 'all', 
													  n        = numDamSites,
												      topOrBot = t.lower())

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
													normType  = normType,
													outputDir = self.outputPlotDir,
													saveFig   = True,
													figTitle  = '{} {} D{} damage sites'.format(t,numDamSites,metric),
													saveName  = 'Lineplot_Metric-D{}_Normalisation-{}_{}-atoms'.format(metric,normType,t))

				figCalls.append('<img class="img-responsive" src="plots/{}" width="{}">'.format(figName.split('/')[-1],400))

			info = '<div class = "row">\n'+\
			  	   '<div class = "col-sm-6">{}</div>\n'.format(figCalls[0])+\
			       '<div class = "col-sm-6">{}</div>\n'.format(figCalls[1])+\
			       '</div>\n'

			summaryFile.write(info)

		# get structure-wide metric average & std dev
		av,std = self.atmsObjs.getAverageMetricVals(densMet  = metric,
												    normType = normType)
		for i in range(numDsets):
			summaryFile.write('<hr>')

			# start block of collapsible panels here
			summaryFile.write('<div class = "panel-group" id = "datasetinfo{}">'.format(i))
			summaryFile.write('<h2>Dataset {}</h2>\n'.format(i+1))

			t = 'Dataset info'
			c = 'Number in series : {}<br>\n'.format(i+1)+\
				'Dose (MGy)       : {}<br>\n'.format(self.doses[i])+\
				'Number of atoms  : {}<br>\n'.format(self.atmsObjs.getNumAtoms())+\
				'Fourier diff map : <a href ="../RIDL-maps/{}_density.map">Download</a><br>\n'.format(self.pdbNames[i])
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)	
			summaryFile.write(panelStr)

			t = 'Structure-wide D<sub>{}</sub> summary'.format(metric)
			c = 'Average D<sub>{}</sub>   : {}<br>\n'.format(metric,round(av[i],3))+\
				'Std dev in D<sub>{}</sub>: {}<br>\n'.format(metric,round(std[i],3))
					 				
			if normType == 'Calpha normalised':
				CAweights = self.atmsObjs.retrieveCalphaWeight(metric = metric)
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

			figName = self.atmsObjs.plotNumAtomsWithMetricAboveStructureWideMean(metric   = metric,
																				 normType = normType,
																			     dataset  = i,
																				 outputLoc = outDir)

			t = '# atoms with {} D<sub>{}</sub> metric above N std dev of structure-wide mean'.format(norm,metric)
			c = '<img class="img-responsive" src="plots/{}{}" width="{}">'.format(subdir,figName.split('/')[-1],imWidth)
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)	
			summaryFile.write(panelStr)

			statsOut = self.atmsObjs.getTopNAtomsString(metric   = metric,
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

			statsOut = self.atmsObjs.getPerAtmtypeStats(metric   = metric,
														normType = normType,
														dataset  = i)
			t = 'Per-atom-type Statistics'
			c = 'D<sub>{}</sub> metric ranked by mean value:<br><br>\n'.format(metric)
			c += self.convertPlainTxtTable2html(statsOut[0],width='60%')
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)
			summaryFile.write(panelStr)

			n = self.atmsObjs.getNumAtoms()*0.05 # take top 5% of atoms here
			statsOut = self.atmsObjs.breakdownTopNatomsBy(metric   = metric,
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

			statsOut = self.atmsObjs.getPerResidueStats(metric   = metric,
													    normType = normType,
													    dataset  = i)
			t = 'Per-residue Statistics'
			c = 'D<sub>{}</sub> metric ranked by mean value:<br><br>\n'.format(metric)
			c += self.convertPlainTxtTable2html(statsOut[0],width='60%')
			panelStr = self.writeHtmlDropDownPanel(title   = t,
										           content = c,
										           dataset = i)
			summaryFile.write(panelStr)

			statsOut = self.atmsObjs.getPerChainStats(metric   = metric,
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
				suspAtoms,highOrLow = self.atmsObjs.detectSuspiciousAtoms(dataset   = i,
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

		numDsets = self.atmsObjs.atomList[0].getNumDatasets()
		return numDsets

	def writeDamSitesToFile(self, 
							metric      = 'loss', 
							normType    = 'Standard', 
							numDamSites = 'all'):

		# write top damage sites to .pdb file for each dataset

		self.damSitesPDB = []
		for i in range(self.getNumDatasets()): 
			damPDB = self.atmsObjs.getTopNAtomsPDBfile(metric   = metric,
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
			self.atmsObjs.calcAdditionalMetrics(metric   = metric,
												normType = normType)

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
			damSitesTag = (self.damSitesPDB[dataset].split('/')[-1]).strip('.pdb')
			structureTag = self.initialPDB.strip('.pdb')
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

			if savePic is True:
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
			data = self.atmsObjs.graphMetricDistn(metric     = densMet,
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

	def endAfter(self):

		# force quit of the class

		sys.exit()



