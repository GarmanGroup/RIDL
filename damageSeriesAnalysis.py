from runETRACK import run as ETRACK
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt

class damageSeriesAnalysis():
	# For damage series from PDB run the ETRACK scripts to calculate per-atom density metrics.
	# Below are a series of methods for processing specific damage series within the pdb.

	def __init__(self):
		self.ETRACK = ETRACK()

	def runDataseries(self,name,version,process,postprocess,retrieve):
		if name == 'TRAP':
			self.writeETRACKinputfile_TRAP(version)
		elif name == 'BURM':
			self.writeETRACKinputfile_Burm2000(version)
		elif name == 'FIOR':
			self.writeETRACKinputfile_Fior2007()
		elif name == 'DELAMORA':
			self.writeETRACKinputfile_DelaMora2011(version)
		elif name == 'FRANK':
			self.writeETRACKinputfile_Frankaer2014()
		elif name == 'JUERS100':
			self.writeETRACKinputfile_Juers2011_100K()
		elif name == 'JUERS160':
			self.writeETRACKinputfile_Juers2011_160K()		
		elif name == 'BURY':
			self.writeETRACKinputfile_Bury2015(version)
		elif name == 'WEIK':
			self.writeETRACKinputfile_Weik2000(version)
		elif name == 'DIXON':
			self.writeETRACKinputfile_TDIXONinsulin(version)
		elif name == 'SUTTON':
			self.writeETRACKinputfile_Sutton2013()
		elif name == 'PETROVA':
			self.writeETRACKinputfile_Petrova2010(version)
		elif name == 'NANAO':
			self.writeETRACKinputfile_Nanao2005(version)
		else:
			'Data series not recognised'
			return

		self.ETRACK.runETRACK(process,postprocess,retrieve)

	def getDensSeries(self):
		dSeries = [['DELAMORA','DIFF'],['JUERS100',''],['JUERS160',''],
 					['BURY','initial'],['WEIK','initial'],['DIXON','DIFF'],
 					['FIOR',''],['BURM','final'],['TRAP','DIFF'],['SUTTON',''],
 					['PETROVA','DIFF'],['NANAO','Elastase'],['NANAO','Thaumatin'],
 					['NANAO','Trypsin'],['NANAO','Lysozyme'],['NANAO','Insulin'],
 					['NANAO','RibonucleaseA']]
 		return dSeries

 	def writeETRACKinputfile_TRAP(self,version):
		# Need to create input file for this test TRAP damage series
		doses = '1.31,3.88,6.45,9.02,11.58,14.15,16.72,19.29,21.86,24.98'
		dNames = '1,2,3,4,5,6,7,8,9,10'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,version)
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/TRAP_ETRACK/{}/'.format(version),
									      		 'TRAP',dNames,'TRAP1.pdb',doses,'TRAP_data.pkl')					

	def writeETRACKinputfile_Burm2000(self,version):
		# Need to create input file for the Burmeister 2000 damage series
		# version in ('final','initial') depending on which PDB_REDO version required
		doses = '1,2,3,4'
		dNames = 'a,f,g,h'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,'DIFF')
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Burm2000_ETRACK/{}/'.format(version),
									      '1dw',dNames,'1dwa.pdb',doses,'1dw_data.pkl')							   

	def writeETRACKinputfile_Sutton2013(self):
		# Need to create input file for the Sutton 2013 damage series
		doses = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15'
		dNames = '8x,8y,8z,9a,9b,9c,9e,9f,9h,9i,90,91,92,93,94'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,'DIFF')
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Sutton2013_ETRACK/DIFF/',
									      '4h',dNames,'4h8x.pdb',doses,'4h_data.pkl')

	def writeETRACKinputfile_Petrova2010(self,version):
		# Need to create input file for the Petrova 2010 damage series
		dNames = '1,2,3,4,5,6,7,8'
		doses = '1.2,14.2,15.4,28.4,29.6,42.6,43.8,56.8'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,version)
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Petrova2010_ETRACK/M_100K/{}/'.format(version),
									      'PET',dNames,'PET1.pdb',doses,'PET_data.pkl')

	def writeETRACKinputfile_Nanao2005(self,protein):
		# Need to create input file for the Nanao 2005 damage sets
		# 'protein' in ('Elastase','Thaumatin','Trypsin','Lysozyme','Insulin','RibonucleaseA')
		pInfo = {'Trypsin':'lv,lw','Thaumatin':'lr,lu','Elastase':'lo,lq',
				 'Lysozyme':'lx,ly','RibonucleaseA':'lp,lz','Insulin':'n3,n1'}
		if protein not in pInfo.keys(): 
			return

		doses = '1,2'
		dNames = pInfo[protein]
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,'DIFF')
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Nanao2005_ETRACK/{}/DIFF/'.format(protein),
									      '2b',dNames,'2b{}.pdb'.format(pInfo[protein].split(',')[0]),doses,'2b_data.pkl')

	def writeETRACKinputfile_DelaMora2011(self,version):
		# Need to create input file for the DelaMora 2011 damage series
		# type in ('DIFF','SIMPLE','END')
		doses = '2.31,6.62,12.31,17.9,23.3,28.6'
		dNames = 'h,i,j,l,m,n'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,version)
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/DelaMora2011_ETRACK/{}/'.format(version),
									      '2yb',dNames,'2ybh.pdb',doses,'2yb_data.pkl')

	def writeETRACKinputfile_Fior2007(self):
		# Need to create input file for the Fioravanti 2007 series
		doses = '1,2,3'
		dNames = 'k,q,r'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,'DIFF')
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Fioravanti2007_ETRACK/',
									      '2j5',dNames,'2j5k.pdb',doses,'2j5_data.pkl')

	def writeETRACKinputfile_Frankaer2014(self):
		# Need to create input file for the Frankaer 2014 series
		doses = '2,3,4'
		dNames = 'h,i,j'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,'DIFF')
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Frankaer2014_ETRACK/',
									      '4m4',dNames,'4m4f.pdb',doses,'4m4_data.pkl')

	def writeETRACKinputfile_Juers2011_100K(self):
		# Need to create input file for the Frankaer 2014 100K series
		doses = '1,2,3,4'
		dNames = 'p,q,r,s'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,'DIFF')
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Juers2011-100K_ETRACK/',
									      '3p7',dNames,'3p7p.pdb',doses,'3p7_data.pkl')

	def writeETRACKinputfile_Juers2011_160K(self):
		# Need to create input file for the Frankaer 2014 160K series
		doses = '1,2,3,4'
		dNames = 't,u,v,w'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,'DIFF')
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Juers2011-160K_ETRACK/',
									      '3p7',dNames,'3p7t.pdb',doses,'3p7_data.pkl')

	def writeETRACKinputfile_Bury2015(self,version):
		# Need to create input file for the Bury2015 series
		# version in ('final','initial') depending on which PDB_REDO version required
		doses = '2.1,6.2,10.3,14.4,20.6,26.8,35.7,44.6'
		dNames = 'b,c,d,e,f,g,h,i'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,'DIFF')
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Bury2015_ETRACK/{}/'.format(version),
									      '4x4',dNames,'4x4b.pdb',doses,'4x4_data.pkl')

	def writeETRACKinputfile_Weik2000(self,version):
		# Need to create input file for the Weik2000 series
		# version in ('final','initial') depending on which PDB_REDO version required
		doses = '1,2,3,4,5,6,7,8,9'
		dNames = 'd,e,f,g,h,i,j,k,m'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,'DIFF')
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Weik2000_ETRACK/{}/'.format(version),
									      '1qi',dNames,'1qid.pdb',doses,'1qi_data.pkl')

	def writeETRACKinputfile_TDIXONinsulin(self,version):
		# Need to create input file for the TDIXON-insulin series
		doses = '0.89,2.74,4.59,6.44,8.29,10.14,11.98,13.84,15.68,17.54'
		dNames = '1,2,3,4,5,6,7,8,9,10'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,version)
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/TDIXON_InsulinSeries_ETRACK/{}/'.format(version),
									      'insu',dNames,'insu1.pdb',doses,'insu_data.pkl')

	def writeETRACKinputfile_TDIXONinsulin_FEB(self,version):
		# Need to create input file for the TDIXON-insulin series
		doses = '0.89,2.74,4.59,6.44,8.29,10.14,11.98,13.84,15.68,17.54'
		dNames = '1,2,3,4,5,6,7,8,9,10'
		doses,dNames = self.ETRACK.defineDoseList(doses,dNames,version)
		inputString = self.ETRACK.writeInputFile('/Users/charlie/DPhil/YEAR2/FEB/ETRACK-testing/2FOFC-testing/TDinsulin/',
									      'insu',dNames,'insu1.pdb',doses,'insu_data.pkl')


	def residueMetricDistributionPlots(self,mapPro,postPro,retr,densMet,normType):
		# for each damage series retrieve the distribution for damage metric values for all atoms of 
		# specified residues (see 'resList' below), and then plot a kde plot for each
		resList = [['GLU','GLN'],['ASP','ASN'],['ILE','LEU'],['TYR','PHE'],
				   ['TYR','ASP','GLU'],['TYR','PHE','GLY'],
				   ['GLU','GLY'],['ASP','GLY'],['TYR','GLY'],['PHE','GLY'],
				   ['CYS','GLY'],['MET','GLY'],['GLU','ASP','CYS','MET','TYR']]
		dSeries = self.getDensSeries()
		plotData = {'-'.join(res):{} for res in resList}
		for key in plotData.keys():
			for k in key.split('-'):
				plotData[key][k] = []
		for dSer in dSeries:
			self.runDataseries(dSer[0],dSer[1],mapPro,postPro,retr)
			for resGroup in resList:
				if len(resGroup) > 4:
					plotType = 'kde'
				else:
					plotType = 'both'
				data = self.ETRACK.et.combinedAtoms.graphMetricDistn(densMet,normType,True,plotType,resGroup,True)

				for r in resGroup:
					plotData['-'.join(resGroup)][r] += data[r]

		for resGroup in resList:
			# sns.set_palette("deep", desat=.6)
			sns.set_style("whitegrid")
			sns.set_context(rc={"figure.figsize": (10, 6)})
			fig = plt.figure()
			ax = plt.subplot(111)
			if len(resGroup) > 4:
				plotType = 'kde'
			else:
				plotType = 'both'
			for j,(res, color) in enumerate(zip(resGroup, sns.color_palette('hls', len(resGroup)))):
				datax = plotData['-'.join(resGroup)][res]
				print 'number of {} atoms = {}'.format(res,len(datax))
				self.ETRACK.et.combinedAtoms.plotHist(plotType,300,datax,'average,{}'.format(res),color)
			plt.legend()
			plt.xlabel('{} D{} per atom'.format(normType,densMet),fontsize=18)
			plt.ylabel('Norm-frequency',fontsize=18)
			plt.title('{} D{} per atom, residues: {}'.format(normType,densMet,','.format(resGroup)))
			saveName = '{}_{}D{}-{}-COMBINED.png'.format(''.join(resGroup),normType.replace(" ",""),densMet,plotType)
			fig.savefig(saveName)

	def findTyrDlossVNeighbourhood(self,mapPro,postPro,retr,densMet,normType,resType,atomType,distance,weighted):
		# for each damage series determine the speicified damage metric for all resType-atomType atoms and compare
		# to the average density metric within a local distance of each atom.
		# 'distance' is distance from an atom in Angstrom
		# 'weighted' is True if distance-weighted average metric taken over local environment, otherwise, 
		# if False then standard average of metric taken.
		combPlotData = {'x':[],'y':[]}
		dSeries = self.getDensSeries()
		for dSer in dSeries:
			self.runDataseries(dSer[0],dSer[1],mapPro,postPro,retr)
			if normType == 'Calpha normalised': 
				self.et.combinedAtoms.calcAdditionalMetrics(densMet,normType,'Calpha')
			self.ETRACK.et.combinedAtoms.calcAdditionalMetrics(densMet,normType,'average')	
			plotData = 	self.ETRACK.et.combinedAtoms.calculateLocalDloss(resType,atomType,distance,densMet,normType,weighted)
			for k in plotData.keys():
				combPlotData[k] += plotData[k]

		# plot the relationship between atom Dloss and local environment average Dloss
		self.ETRACK.et.combinedAtoms.plotScatterPlot(combPlotData['x'],combPlotData['y'],
							 '{}-{} D{}'.format(resType,atomType,densMet),
							 'local environment D{}'.format(densMet),
							 'Average D{} within {} Angstrom of {}-{}'.format(densMet,distance,resType,atomType),
							 'D{}Scatter_{}-{}_localEnvironmentComparison-COMBINED.png'.format(densMet,resType,atomType),
							  True,False,'')

	def rankTyrOHdamage(self,mapPro,postPro,retr,densMet,normType):
		# for each damage series determine the specified damage metric for all TYR-OH atoms and rank them
		dSeries = self.getDensSeries()
		tyrOHatms = {'atmID':[],'metric':[]}
		for dSer in dSeries:
			self.runDataseries(dSer[0],dSer[1],mapPro,postPro,retr)
			if normType == 'Calpha normalised': 
				self.et.combinedAtoms.calcAdditionalMetrics(densMet,normType,'Calpha')
			self.ETRACK.et.combinedAtoms.calcAdditionalMetrics(densMet,normType,'average')

			atms = self.ETRACK.et.combinedAtoms.getAtom('','TYR','','OH')
			for a in atms:
				tyrOHatms['atmID'].append(dSer[0]+'-'+dSer[1]+'-'+a.getAtomID())
				tyrOHatms['metric'].append(a.densMetric[densMet][normType]['average'])

		# sort the Tyr-OH atms in order of metric value:
		list1, list2 = (list(t) for t in zip(*sorted(zip(tyrOHatms['metric'], tyrOHatms['atmID']))))
		for i in range(len(list1)):
			print '{} ---> {}'.format(list2[i],list1[i])


	def runBatchSeries(self,mapPro,postPro,retr,densMet,normType):
		# plot scatter plots of damage to GLU/ASP/TYR side-chain atoms to compare the relative damage to adjacent atoms
		# e.g. TYR-OH vs TYR-CZ or GLU-CD vs GLU-CG etc.
 		dSeries = self.getDensSeries()

		rSquaredDic,numPairsDic,dataDic = {},{},{}
		for dSer in dSeries:
			self.runDataseries(dSer[0],dSer[1],mapPro,postPro,retr)

			# calculate additional metrics as required
			self.ETRACK.et.combinedAtoms.calcMetricDiffFromStructureMean(densMet,'Standard','std-devs')
			if normType == 'Calpha normalised': 
				self.ETRACK.et.combinedAtoms.calcAdditionalMetrics(densMet,normType,'Calpha')

			rSquared,numPairs,data = self.ETRACK.et.combinedAtoms.compareSensAtoms(densMet,normType)
			rSquaredDic[''.join(dSer)] = rSquared
			numPairsDic[''.join(dSer)] = numPairs
			for key in data.keys():
				if key not in dataDic.keys(): dataDic[key] = {'x':[],'y':[]}
				for i in ('x','y'): 
					dataDic[key][i] += data[key][i]#+= list((np.array(data[key][i])-avMetric*np.ones(len(data[key][i])))/stdMetric)

		# plot this data as scatter plot
		for key in dataDic.keys():
			xData  = dataDic[key]['x']
			yData  = dataDic[key]['y']
			k 	   = key.split('_')
			xLabel = '{}-{} D{}'.format(k[0],k[1],densMet)
			yLabel = '{}-{} D{}'.format(k[0],k[2],densMet)
			rSquared = self.ETRACK.et.combinedAtoms.plotScatterPlot(xData,yData,xLabel,
							 yLabel,'D{} for {} atoms'.format(densMet,key),
							 'D{}_{}_scatterplotCOMBINED.svg'.format(densMet,key),True,True,'')

		return rSquaredDic,numPairsDic,dataDic

	def runBatchSeries2(self,mapPro,postPro,retr,densMet,normType,rType):
		pairs = [['GLU','CD','CG'],['GLU','CD','OE1'],['ASP','CG','CB'],
				 ['ASP','CG','OD1'],['TYR','OH','CZ']]

		dSeries = self.getDensSeries()
		seriesData = {}
		for dSer in dSeries:
			self.runDataseries(dSer[0],dSer[1],mapPro,postPro,retr)
			self.ETRACK.et.combinedAtoms.calcMetricDiffFromStructureMean(densMet,'Standard','std-devs')
			data = self.ETRACK.et.combinedAtoms.findMetricRatioKeyResidues_scatterplot(densMet,normType,rType,pairs,dSer)
			seriesData['-'.join(dSer)] = data

		# plot a combined dataseries plot here	
		sns.set_palette("deep", desat=.6)
		sns.set_context(rc={"figure.figsize":(40, 10)})
		fig = plt.figure()
		ax = plt.subplot(111)
		colors = ['#737474','#409cd6','#58bb6b','#faa71a','#ff6a6a']
		i = -2
		for key in seriesData.keys():
			i += 2
			for pair in seriesData[key].keys():
				j = (['-'.join(p) for p in pairs]).index(pair) 
				x = seriesData[key][pair]['x']
				y = seriesData[key][pair]['y']
				xNew = np.array(x)+i
				if i == 0:
					plt.scatter(xNew,y,marker='o',s=100,c=colors[j],edgecolors='#FFFFFF',label=pair)
				else:
					plt.scatter(xNew,y,marker='o',s=100,c=colors[j],edgecolors='#FFFFFF')

		plt.plot([-1, i+2],[0,0],':',color='#088DA5') # a horizontal line at y=0
		tickPts = list(2*np.array(range(len(seriesData.keys())))+0.5)
		ax.set_xlim([-1, i+2])
		plt.xticks(tickPts,seriesData.keys())
		plt.xlabel('Damage Series', fontsize=24)
		plt.ylabel('D{} {}'.format(densMet,rType), fontsize=24)
		figtitle = '{} D{} {}'.format(normType,densMet,rType)

		# place legend outside to right of plot
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.95, box.height])
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=24)

		fig.suptitle(figtitle,fontsize=28)
		saveTitle = figtitle.replace(' ','_')
		fname = lambda x: saveTitle.strip('.png')+'_{}.png'.format(x)
		i = 0
		while os.path.isfile(fname(i)): i += 1 
		fig.savefig(fname(i))

	def runBatchSeries3(self,mapPro,postPro,retr,densMet,normType,distLim):
		seriesData = {}
		dSeries = self.getDensSeries()

		plotData1 = {'x':[],'y':[],'colors':[]}
		plotData2 = {'x':[],'y':[]}
		for dSer in dSeries:
			self.runDataseries(dSer[0],dSer[1],mapPro,postPro,retr)
			for diff in ('std-devs','ratio'):
				self.ETRACK.et.combinedAtoms.calcMetricDiffFromStructureMean(densMet,'Standard',diff)

			if normType == 'Calpha normalised': 
				self.ETRACK.et.combinedAtoms.calcAdditionalMetrics(densMet,normType,'Calpha')

			self.ETRACK.et.combinedAtoms.calcAdditionalMetrics('loss',normType,'average')
			self.ETRACK.et.combinedAtoms.seriesName = '-'.join(dSer)

			# get pdb file information and space group
			pdbName = self.ETRACK.et.where+self.et.initialPDB 
			sGroup = self.ETRACK.et.getSpaceGroup()
			if sGroup is False: return

			# determine correlation between TYR-OH damage and carboxyl contacts
			seriesData['-'.join(dSer)],TyrDam,CO2Dam,scatterColor = self.ETRACK.et.combinedAtoms.densMetSurroundAtmsCorrel('TYR','OH',distLim,normType,densMet,
																						 			   True,pdbName,sGroup)
			# data for Glu-CD/Asp-CG Dloss vs Tyr-OH Dloss
			plotData1['x'] += TyrDam
			plotData1['y'] += CO2Dam
			plotData1['colors'] += scatterColor

			# determine the per-dataset correlation between solvent accessibility and Tyr-OH Dloss
			solvAccDic,plotData = self.ETRACK.et.combinedAtoms.compareSolvAccessWithAndWithoutGluAspGroups(pdbName,sGroup,'TYR','OH',[],densMet,normType)
			plotData2['x'] += plotData['x']
			plotData2['y'] += plotData['y']

		for key in seriesData.keys():
			print 'For dataseries: {}'.format(key)
			print seriesData[key]

		# plot the relationship between TYR-OH density and nearby carboxyl atom density
		self.ETRACK.et.combinedAtoms.plotScatterPlot(plotData1['x'],
													 plotData1['y'],
											  		 'TYR-OH D{}'.format(densMet),
											  		 'Carboxyl-contact D{}'.format(densMet),
											  		 'TYR-OH vs carboxyl-contact D{}'.format(densMet),
											  		 'D{}Scatter_TYR-OH_carboxylContacts-COMBINED.png'.format(densMet),
											  		 True,
											  		 False,
											  		 plotData1['colors'])

		# plot the relationship between TYR-OH density and solvent accessibility
		self.ETRACK.et.combinedAtoms.plotScatterPlot(plotData2['x'],
													 plotData2['y'],
											  		 'Solvent Accessibility',
											  		 'TYR-OH D{}'.format(densMet),
											  		 'TYR-OH D{} vs solvent accessibility'.format(densMet),
											  		 '{}D{}Scatter_TYR-OH_solvAccess-COMBINED.png'.format(normType,densMet),
											  		 True,False,
											  		 '')

