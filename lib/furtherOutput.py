from ridlFeedback import provideFeedback

class furtherAnalysis(provideFeedback):

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
				 metSigPlots      = False,
				 mainSideCmp      = False,
				 distToDiSulp     = False,
				 tTests           = False,
				 topOfType        = False,
				 rankings         = False,
				 sigKSstats       = False,
				 atmCorels        = False,
				 HbondCorr        = False,
				 sensAtmPlts      = False,
				 doAll            = True,
				 autoRun          = False):

		# collection of statistics and further analysis.
		# This plot should NOT be output in a default run of program

		super(furtherAnalysis, self).__init__( standardFeedback,
											   writeCsvs,
											   csvOnly,
											   writeSumFile,
											   writeTopSites,
											   plotHeatMaps,
											   atmsObjs,
											   outputDir,
											   outputPlotDir,
											   csvExtent,
											   plotGraphs,
											   logFile,
											   pklSeries,
											   doses,
											   pdbNames,
											   inputDir,
											   initialPDB,
											   autoRun )

		# first perform standard output
		self.run()

		# now include extra analysis
		if doAll:
			metSigPlots  = True,
			mainSideCmp  = True,
			distToDiSulp = True,
			tTests       = True,
			topOfType    = True,
			rankings     = True,
			sigKSstats   = True,
			atmCorels    = True,
			HbondCorr    = True,
			sensAtmPlts  = True,

		self.metric = 'loss'
		for norm in ('Standard','Calpha normalised'):
			self.normType = norm
			if metSigPlots:
				self.metSigPlots()
			if mainSideCmp:
				self.mainSideCmp()
			if distToDiSulp:
				self.distToDiSulp()
			if tTests:
				self.tTests()
			if topOfType:
				self.topOfType()
			if rankings:
				self.rankings()
			if sigKSstats:
				self.sigKSstats()
			if atmCorels:
				self.atmCorels()
			if HbondCorr:
				self.HbondCorr()
			if sensAtmPlts:
				self.sensAtmPlts()

	def metSigPlots(self):

		for t,u in zip(['','-sideOnly'],[False,True]):
			outDir = self.makeNewPlotSubdir(subdir = 'metricSignaturePlots{}/'.format(t))
			for i in range(self.getNumDatasets()):
				self.makeDistnPlots(densMet    = self.metric,
							   		normType   = self.normType,
								    plotSet    = 2,
								    calcKSstat = True,
								    calcADstat = True,
								    sideOnly   = u,
								    outputDir  = outDir,
								    dataset    = i,
								    requireAll = True)

	def mainSideCmp(self):

		outDir = self.makeNewPlotSubdir(subdir = 'side-main-compare-distn/')

		for d in range(self.getNumDatasets()):
			for s in [4,5]:
				self.makeDistnPlots(densMet    = self.metric,
							        normType   = self.normType,
							        plotSet    = s,
							        calcKSstat = True,
							   	    calcADstat = True,
							        outputDir  = outDir,
							        dataset    = d)

	# def BdamageAnalysis(self):

	# 	COMMENTED CODE BELOW ONLY WORKS FOR MYROSINASE (BURMEISTER, 2000)
	# 	loc = '/Users/charlie/DPhil/YEAR2/JUN/setOccupanciesTo1/'
	# 	for bmet in ('Bdamage','Bfactor'):
	# 		for d,i in zip(range(self.getNumDatasets()),['f','g','h','i']):
	# 			self.atmsObjs.compareDensToBdamChange(BDamageFile1  = loc + '1dwa_editted_refmac1Bdamage.txt',
	# 												  BDamageFile2  = loc + '1dw{}_editted_refmac1Bdamage.txt'.format(i),
	# 												  BdamChange    = True,
	# 												  resType       = 'TYR',
	# 												  atomType      = 'OH',
	# 												  Bmetric       = bmet,
	# 												  percentChange = True,
	# 												  densMet       = metric,
	# 												  normType      = normType,
	# 												  outputDir     = loc,
	# 												  dataset       = d)

	def distToDiSulp(self):

		print 'Finding distances to nearest disulphide bridges'

		self.atmsObjs.scatterAtmsByDistToOtherAtms(outputDir = self.outputPlotDir)

	# def metValsAtDistFromAtm(self):

	# 	COMMENTED CODE BELOW ONLY WORKS FOR TRAP (BURY ET AL., 2016)
	# 	for resNum in ['62']:
	# 		self.atmsObjs.metricValsAtDistFromAtom(atomType  = 'OH',
	# 											   resType   = 'TYR',
	# 											   resNum    = resNum,
	# 											   chainType = 'A',
	# 											   metric    = metric,
	# 											   normType  = normType,
	# 											   outputDir = self.outputPlotDir)

	def topOfType(self,
				  res = 'TYR',
				  atm = 'OH'):

		print 'Finding top atom of type: {}'.format('-'.join([res,atm]))

		self.atmsObjs.findTopAtomOfType(normType  = self.normType,
								  		atomType  = res,
						  				resType   = atm)

	def tTests(self):

		print 'Performing t-test to compare metric correlation between atom types'

		for d in range(self.getNumDatasets()):
			stats = self.atmsObjs.twoAtomTypeTtest(dataset   = d,
												   atomType2 = ['OE1','OE2','OD1','OD2'],
												   resType2  = ['GLU','GLU','ASP','ASP'])
			print 'Glu/Asp) Dataset: {}, stats: {}'.format(d,stats)

			stats = self.atmsObjs.twoAtomTypeTtest(dataset   = d,
												   atomType2 = ['OE1','OD1'],
												   resType2  = ['GLN','ASN'])
			print 'Gln/Asn) Dataset: {}, stats: {}'.format(d,stats)

	def rankings(self):

		res = ['CYS','TYR','TYR','GLU','GLU','GLU','GLU','ASP','ASP','ASP','ASP','GLN','GLN','ASN','ASN','SER','THR']
		atm = ['SG','OH','CZ','CD','OE1','OE2','CG','CG','CB','OD1','OD2','OE1','NE2','OD1','ND2','OG','OG1']

		print 'Determining per atom-type rankings for selected atoms in structure'

		rank = self.atmsObjs.getAtomtypeRanking(metric    = self.metric,
											    normType  = self.normType,
											    dataset   = 'all',
											    residue   = res,
											    atomtype  = atm,
											    printText = True)
		for key in rank.keys():
			print '{}: {}'.format(key,rank[key])

		for r,a in zip(res,atm):
			num = self.atmsObjs.numAtmsWithMetricAboveLevel(dataset   = 'all',
															metric    = self.metric,
															normType  = self.normType,
															threshold = 1,
															firstTime = True,
															atomType  = a,
															resType   = r)
			print '{}-{}: {}'.format(r,a,num)

		if self.getNumDatasets() == 1:
			return

		outDir = self.makeNewPlotSubdir(subdir = 'TyrOH-ranking/')
		self.atmsObjs.plotAtomtypeRankingWithDataset(outputDir = outDir,
												     percent   = False,
												     normType  = self.normType)

		for r,a in zip(res,atm):
			self.atmsObjs.plotNumAtomsWithMetricAboveStructureWideMean(metric    = self.metric,
																       normType  = self.normType,
																	   dataset   = 'all',
																	   outputLoc = outDir,
																	   atomType  = a,
																	   resType   = r)

	def sigKSstats(self):

		print 'Calculating Kolmogorov-Smirnov statistics to compare damage '+\
			  'signatures between atom types in current structure'

		outDir = self.makeNewPlotSubdir(subdir = 'MetricSignature-statistics/')
		self.atmsObjs.plotStatVsDataset(normType  = self.normType,
										outputDir = outDir)

		self.atmsObjs.plotStatVsStat(normType  = self.normType,
								     outputDir = outDir)

		plotType = 'box'
		outDir = self.makeNewPlotSubdir(subdir = 'KS-stat-values/')

		# for ref in ('GLY','ALA'):
		# 	saveName = '{}plot_Metric-D{}_Normalisation-{}-KSstat_allResidues-rel-to-{}'.format(plotType,metric,normType,ref)
		# 	self.atmsObjs.plotKSstatVsDataset(normType  = normType,
		# 									  outputDir = outDir,
		# 									  reference = ref,
		# 									  plotType  = plotType,
		# 									  saveName  = saveName)

		res = ['GLU','ASP','TYR','LYS','ILE','THR']
		ref = ['GLN','ASN','PHE','ARG','LEU','SER']
		saveName = '{}plot_Metric-D{}_Normalisation-{}-KSstat_keyResidues'.format(plotType,self.metric,self.normType)

		sideOnly = True
		if sideOnly is True:
			res += ['GLU']
			ref += ['ASP']
			saveName += '_sideChainsOnly'

		self.atmsObjs.plotKSstatVsDataset(normType  = self.normType,
										  outputDir = outDir,
										  reference = ref,
										  residues  = res,
										  plotType  = plotType,
										  inclLine  = False,
										  sideOnly  = sideOnly,
										  saveName  = saveName)

	def atmCorels(self):

		print 'Determining correlation between atom types in same residues'

		outDir = self.makeNewPlotSubdir(subdir = 'atomtype-correlationPlots/')
		self.atmsObjs.compareSensAtoms_RsquaredLinePlotWithDataset(normType  = self.normType,
						 										   outputDir = outDir,
						 										   plotType  = 'box')

	def HbondCorr(self,
				  res = 'TYR',
				  atm = 'OH'):

		print 'Determining whether a correlation exists between {} and '.format('-'.join([res,atm]))+\
			  'the presence of carboxylate salt bridge interactions'

		outDir = self.makeNewPlotSubdir(subdir = 'carboxylate-TyrOH-correlation/')
		self.atmsObjs.plot_densMetSurroundAtmsCorrel(pdbName       = self.get1stDsetPDB(),
												     symmetrygroup = self.getSpaceGroup(),
													 restype       = res,
							  	   					 atomtype      = atm,
													 normType      = self.normType,
													 outputDir     = outDir,
													 errorbars     = False,
													 printText     = False)

	def sensAtmPlts(self):

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
				self.atmsObjs.findMetricRatioKeyResidues(metric   = 'loss',
														 normType = self.normType,
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
			for eBars in (True,False):
				self.atmsObjs.plotSusceptibleAtoms(densMet   = self.metric,
												   susAtms   = aTypes,
												   errorbars = eBars,
												   normType  = self.normType)
