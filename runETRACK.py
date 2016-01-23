import sys
sys.path.insert(0,'./lib')
from eTrack_RUN import eTrack
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt

class run():
	# run the main ETRACK scripts to calculate per-atom density metrics for a series of increasing doses
	def __init__(self):
		self.inputFileName = 'e_Track_inputfile.txt'

	def runETRACK(self,mapProcess,postProcess,retrieve):
		# run the ETRACK processing for the currently defined input file
		eT = eTrack()
		eT.runPipeline(mapProcess,postProcess,retrieve,self.inputFileName)
		self.et = eT

	def defineDoseList(self,doses,names,version):
		# do not include first dataset dose if difference maps chosen
		if version == 'DIFF':
			dosesOut = ','.join(doses.split(',')[1:])
			namesOut = ','.join(names.split(',')[1:])
			return dosesOut,namesOut
		else: return doses,names

	def writeBlankInputFile(self):
		# Need to create input file for a generic damage series to be completed by the user
		inputString = self.writeInputFile('<input file location>',
										  '<consistent name of series, e.g. for TRAP1.pdb, TRAP2.pdb, this is TRAP>',
										  '<dataset numbers, e.g. 1,2 for TRAP. Can be letters corresponding to pdb series>',
										  '<initial dataset pdb file e.g. TRAP1.pdb>',
										  '<list of doses, length must match length of damageset_num above>',
										  '<if already processed in ETRACK, can specify single output .pkl for series>')	

	def writeInputFile(self,where,damageset_name,damageset_num,initialPDB,doses,PKLMULTIFILE):
		# write a generic input file for a damage series here
		inputString = 'where {}\n'.format(where)+\
					  'damageset_name {}\n'.format(damageset_name)+\
					  'damageset_num {}\n'.format(damageset_num)+\
					  'initialPDB {}\n'.format(initialPDB)+\
					  'doses {}\n'.format(doses)+\
					  'PKLMULTIFILE {}'.format(PKLMULTIFILE)
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()	

	#####################################################################################
	# series of additional methods for processing specific damage series within the pdb #
	#####################################################################################

	def runDataseries(self,name,version,process,postprocess,retrieve):
		if name == 'TRAP':
			self.writeETRACKinputfile_TRAP()
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
			'Dataseries name not recognised'
			return

		self.runETRACK(process,postprocess,retrieve)

	def getDensSeries(self):
		dSeries = [['DELAMORA','DIFF'],['JUERS100',''],['JUERS160',''],
 					['BURY','initial'],['WEIK','initial'],['DIXON',''],
 					['FIOR',''],['BURM','final'],['TRAP',''],['SUTTON',''],
 					['PETROVA',''],['NANAO','Elastase'],['NANAO','Thaumatin'],
 					['NANAO','Trypsin'],['NANAO','Lysozyme'],['NANAO','Insulin'],
 					['NANAO','RibonucleaseA']]
 		return dSeries

 	def writeETRACKinputfile_TRAP(self):
		# Need to create input file for this test TRAP damage series
		doses = '3.88,6.45,9.02,11.58,14.15,16.72,19.29,21.86,24.98'
		dNames = '2,3,4,5,6,7,8,9,10'
		doses,dNames = self.defineDoseList(doses,dNames,'DIFF')
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/TRAP_ETRACK/DIFF/',
									      'TRAP',dNames,'TRAP1.pdb',doses,'TRAP_data.pkl')					

	def writeETRACKinputfile_Burm2000(self,version):
		# Need to create input file for the Burmeister 2000 damage series
		# version in ('final','initial') depending on which PDB_REDO version required
		doses = '2,3,4'
		dNames = 'f,g,h'
		doses,dNames = self.defineDoseList(doses,dNames,'DIFF')
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Burm2000_ETRACK/{}/'.format(version),
									      '1dw',dNames,'1dwa.pdb',doses,'1dw_data.pkl')							   

	def writeETRACKinputfile_Sutton2013(self):
		# Need to create input file for the Sutton 2013 damage series
		doses = '2,3,4,5,6,7,8,9,10,11,12,13,14,15'
		dNames = '8y,8z,9a,9b,9c,9e,9f,9h,9i,90,91,92,93,94'
		doses,dNames = self.defineDoseList(doses,dNames,'DIFF')
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Sutton2013_ETRACK/DIFF/',
									      '4h',dNames,'4h8x.pdb',doses,'4h_data.pkl')

	def writeETRACKinputfile_Petrova2010(self,version):
		# Need to create input file for the Petrova 2010 damage series
		dNames = '1,2,3,4,5,6,7,8'
		doses = '1.2,14.2,15.4,28.4,29.6,42.6,43.8,56.8'
		doses,dNames = self.defineDoseList(doses,dNames,version)
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Petrova2010_ETRACK/M_100K/{}/'.format(version),
									      'PET',dNames,'PET1.pdb',doses,'PET_data.pkl')

	def writeETRACKinputfile_Nanao2005(self,protein):
		# Need to create input file for the Nanao 2005 damage sets
		# 'protein' in ('Elastase','Thaumatin','Trypsin','Lysozyme','Insulin','RibonucleaseA')
		pInfo = {'Trypsin':['lv','lw'],'Thaumatin':['lr','lu'],'Elastase':['lo','lq'],
				 'Lysozyme':['lx','ly'],'RibonucleaseA':['lp','lz'],'Insulin':['n3','n1']}
		if protein not in pInfo.keys(): return

		doses = '2'
		dNames = pInfo[protein][1]
		doses,dNames = self.defineDoseList(doses,dNames,'DIFF')
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Nanao2005_ETRACK/{}/DIFF/'.format(protein),
									      '2b',dNames,'2b{}.pdb'.format(pInfo[protein][0]),doses,'2b_data.pkl')

	def writeETRACKinputfile_DelaMora2011(self,version):
		# Need to create input file for the DelaMora 2011 damage series
		# type in ('DIFF','SIMPLE','END')
		doses = '2.31,6.62,12.31,17.9,23.3,28.6'
		dNames = 'h,i,j,l,m,n'
		doses,dNames = self.defineDoseList(doses,dNames,version)
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/DelaMora2011_ETRACK/{}/'.format(version),
									      '2yb',dNames,'2ybh.pdb',doses,'2yb_data.pkl')

	def writeETRACKinputfile_Fior2007(self):
		# Need to create input file for the Fioravanti 2007 series
		doses = '2,3'
		dNames = 'q,r'
		doses,dNames = self.defineDoseList(doses,dNames,'DIFF')
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Fioravanti2007_ETRACK/',
									      '2j5',dNames,'2j5k.pdb',doses,'2j5_data.pkl')

	def writeETRACKinputfile_Frankaer2014(self):
		# Need to create input file for the Frankaer 2014 series
		doses = '2,3,4'
		dNames = 'h,i,j'
		doses,dNames = self.defineDoseList(doses,dNames,'DIFF')
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Frankaer2014_ETRACK/',
									      '4m4',dNames,'4m4f.pdb',doses,'4m4_data.pkl')

	def writeETRACKinputfile_Juers2011_100K(self):
		# Need to create input file for the Frankaer 2014 100K series
		doses = '2,3,4'
		dNames = 'q,r,s'
		doses,dNames = self.defineDoseList(doses,dNames,'DIFF')
		inputString = self.writeInputFile('/insu/charlie/DPhil/YEAR2/JAN/Juers2011-100K_ETRACK/',
									      '3p7',dNames,'3p7p.pdb',doses,'3p7_data.pkl')

	def writeETRACKinputfile_Juers2011_160K(self):
		# Need to create input file for the Frankaer 2014 160K series
		doses = '2,3,4'
		dNames = 'u,v,w'
		doses,dNames = self.defineDoseList(doses,dNames,'DIFF')
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Juers2011-160K_ETRACK/',
									      '3p7',dNames,'3p7t.pdb',doses,'3p7_data.pkl')

	def writeETRACKinputfile_Bury2015(self,version):
		# Need to create input file for the Bury2015 series
		# version in ('final','initial') depending on which PDB_REDO version required
		doses = '6.2,10.3,14.4,20.6,26.8,35.7,44.6'
		dNames = 'c,d,e,f,g,h,i'
		doses,dNames = self.defineDoseList(doses,dNames,'DIFF')
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Bury2015_ETRACK/{}/'.format(version),
									      '4x4',dNames,'4x4b.pdb',doses,'4x4_data.pkl')

	def writeETRACKinputfile_Weik2000(self,version):
		# Need to create input file for the Weik2000 series
		# version in ('final','initial') depending on which PDB_REDO version required
		doses = '1,2,3,4,5,6,7,8,9'
		dNames = 'd,e,f,g,h,i,j,k,m'
		doses,dNames = self.defineDoseList(doses,dNames,'DIFF')
		inputString = self.writeInputFile('/Users/charlie/DPhil/YEAR2/JAN/Weik2000_ETRACK/{}/'.format(version),
									      '1qi',dNames,'1qid.pdb',doses,'1qi_data.pkl')

	def writeETRACKinputfile_TDIXONinsulin(self,version):
		# Need to create input file for the TDIXON-insulin series
		doses = '0.89,2.74,4.59,6.44,8.29,10.14,11.98,13.84,15.68,17.54'
		dNames = '1,2,3,4,5,6,7,8,9,10'
		doses,dNames = self.defineDoseList(doses,dNames,version)
		inputString = self.writeInputFile('/insu/charlie/DPhil/YEAR2/JAN/TDIXON_InsulinSeries_ETRACK/{}/'.format(version),
									      'insu',dNames,'insu1.pdb',doses,'insu_data.pkl')

	def runBatchSeries(self,mapPro,postPro,retr,densMet,normType):
 		dSeries = self.getDensSeries()

		rSquaredDic,numPairsDic,dataDic = {},{},{}
		for dSer in dSeries:
			self.runDataseries(dSer[0],dSer[1],mapPro,postPro,retr)
			self.et.combinedAtoms.calcMetricNumStdsFromStructureMean(densMet,'Standard')
			rSquared,numPairs,data = self.et.combinedAtoms.compareSensAtoms(densMet,normType)
			rSquaredDic[''.join(dSer)] = rSquared
			numPairsDic[''.join(dSer)] = numPairs
			# avMetric = np.mean(self.et.combinedAtoms.getAverageMetricVals(densMet,normType)[0]) # average metric over structure
			# stdMetric = np.mean(self.et.combinedAtoms.getAverageMetricVals(densMet,normType)[1])
			for key in data.keys():
				if key not in dataDic.keys(): dataDic[key] = {'x':[],'y':[]}
				for i in ('x','y'): dataDic[key][i] += data[key][i]#+= list((np.array(data[key][i])-avMetric*np.ones(len(data[key][i])))/stdMetric)

		# plot this data as scatter plot
		for key in dataDic.keys():
			xData = dataDic[key]['x']
			yData = dataDic[key]['y']
			rSquared = self.et.combinedAtoms.plotScatterPlot(xData,yData,'D{}'.format(densMet),
							 'D{}'.format(densMet),'D{} for {} atoms'.format(densMet,key),
							 'D{}_{}_scatterplotCOMBINED.png'.format(densMet,key),True)

		return rSquaredDic,numPairsDic,dataDic

	def runBatchSeries2(self,mapPro,postPro,retr,densMet,normType,rType):
		pairs = [['GLU','CD','CG'],['GLU','CD','OE1'],['ASP','CG','CB'],
				 ['ASP','CG','OD1'],['TYR','OH','CZ']]

		dSeries = self.getDensSeries()
		seriesData = {}
		for dSer in dSeries:
			self.runDataseries(dSer[0],dSer[1],mapPro,postPro,retr)
			self.et.combinedAtoms.calcMetricNumStdsFromStructureMean(densMet,'Standard')
			data = self.et.combinedAtoms.findMetricRatioKeyResidues_scatterplot(densMet,normType,rType,pairs,dSer)
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
		for dSer in dSeries:
			self.runDataseries(dSer[0],dSer[1],mapPro,postPro,retr)
			self.et.combinedAtoms.calcAdditionalMetrics('loss','Standard','average')
			seriesData['-'.join(dSer)] = self.et.combinedAtoms.densMetSurroundAtmsCorrel('TYR','OH',distLim,normType,densMet)

		for key in seriesData.keys():
			print 'For dataseries: {}'.format(key)
			print seriesData[key]
