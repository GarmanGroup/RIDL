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

	def writeETRACKinputfile_generic(self):
		# Need to create input file for a generic damage series to be completed by the user
		inputString =	'where <input file location>\n'+\
						'damageset_name <consistent name of series, e.g. for TRAP1.pdb, TRAP2.pdb, this is TRAP>\n'+\
						'damageset_num <dataset numbers, e.g. 1,2 for TRAP. Can be letters corresponding to pdb series>\n'+\
						'initialPDB <initial dataset pdb file e.g. TRAP1.pdb>\n'+\
						'doses <list of doses, length must match length of damageset_num above>\n'+\
						'PKLMULTIFILE <if already processed in ETRACK, can specify single output .pkl for series>'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def runETRACK(self,mapProcess,postProcess,retrieve):
		# run the ETRACK processing for the currently defined input file
		eT = eTrack()
		eT.runPipeline(mapProcess,postProcess,retrieve,False,self.inputFileName)
		self.et = eT


	####################################################################################

	####################################################################################

	####################################################################################
	# series of additional methods for processing specific damage series within the pdb

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
			self.writeETRACKinputfile_TDIXONinsulin()
		elif name == 'SUTTON':
			self.writeETRACKinputfile_Sutton2013()
		elif name == 'PETROVA':
			self.writeETRACKinputfile_Petrova2010()
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
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/TRAP_ETRACK/DIFF/\n'+\
						'damageset_name TRAP\n'+\
						'damageset_num 2,3,4,5,6,7,8,9,10\n'+\
						'initialPDB TRAP1.pdb\n'+\
						'doses 3.88,6.45,9.02,11.58,14.15,16.72,19.29,21.86,24.98\n'+\
						'PKLFILE 13830_TRAP2_data.pkl\n'+\
						'PKLFILE 13830_TRAP3_data.pkl\n'+\
						'PKLFILE 13830_TRAP4_data.pkl\n'+\
						'PKLFILE 13830_TRAP5_data.pkl\n'+\
						'PKLFILE 13830_TRAP6_data.pkl\n'+\
						'PKLFILE 13830_TRAP7_data.pkl\n'+\
						'PKLFILE 13830_TRAP8_data.pkl\n'+\
						'PKLFILE 13830_TRAP9_data.pkl\n'+\
						'PKLFILE 13830_TRAP10_data.pkl\n'+\
						'PKLMULTIFILE TRAP_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()	

	def writeETRACKinputfile_Burm2000(self,version):
		# Need to create input file for the Burmeister 2000 damage series
		# version in ('final','initial') depending on which PDB_REDO version required
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/Burm2000_ETRACK/{}/\n'.format(version)+\
						'damageset_name 1dw\n'+\
						'damageset_num f,g,h\n'+\
						'initialPDB 1dwa.pdb\n'+\
						'doses 2,3,4\n'+\
						'PKLMULTIFILE 1dw_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_Sutton2013(self):
		# Need to create input file for the Sutton 2013 damage series
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/Sutton2013_ETRACK/DIFF/\n'+\
						'damageset_name 4h\n'+\
						'damageset_num 8y,8z,9a,9b,9c,9e,9f,9h,9i,90,91,92,93,94\n'+\
						'initialPDB 4h8x.pdb\n'+\
						'doses 2,3,4,5,6,7,8,9,10,11,12,13,14,15\n'+\
						'PKLMULTIFILE 4h_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_Petrova2010(self):
		# Need to create input file for the Petrova 2010 damage series
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/Petrova2010_ETRACK/M_100K/DIFF/\n'+\
						'damageset_name 3m\n'+\
						'damageset_num nc,ns,nx,o3,o6,o9,oc\n'+\
						'initialPDB 3mnb.pdb\n'+\
						'doses 2,3,4,5,6,7,8\n'+\
						'PKLMULTIFILE 3m_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_Nanao2005(self,protein):
		# Need to create input file for the Nanao 2005 damage sets
		# 'protein' in ('Elastase','Thaumatin','Trypsin','Lysozyme','Insulin','RibonucleaseA')
		if protein == 'Trypsin': 
			initPDB,endPDB = 'lv','lw'
		elif protein == 'Thaumatin':
			initPDB,endPDB = 'lr','lu'	
		elif protein == 'Elastase':
			initPDB,endPDB = 'lo','lq'
		elif protein == 'Insulin':
			initPDB,endPDB = 'n3','n1'
		elif protein == 'Lysozyme':
			initPDB,endPDB = 'lx','ly'
		elif protein == 'RibonucleaseA':
			initPDB,endPDB = 'lp','lz'
		else: return

		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/Nanao2005_ETRACK/{}/DIFF/\n'.format(protein)+\
						'damageset_name 2b\n'+\
						'damageset_num {}\n'.format(endPDB)+\
						'initialPDB 2b{}.pdb\n'.format(initPDB)+\
						'doses 2\n'+\
						'PKLMULTIFILE 2b_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_DelaMora2011(self,type):
		# Need to create input file for the DelaMora 2011 damage series
		# type in ('DIFF','SIMPLE')
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/DelaMora2011_ETRACK/{}/\n'.format(type)+\
						'damageset_name 2yb\n'+\
						'damageset_num i,j,l,m,n\n'+\
						'initialPDB 2ybh.pdb\n'+\
						'doses 2,3,4,5,6\n'+\
						'PKLMULTIFILE 2yb_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_Fior2007(self):
		# Need to create input file for the Fioravanti 2007 series
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/Fioravanti2007_ETRACK/\n'+\
						'damageset_name 2j5\n'+\
						'damageset_num q,r\n'+\
						'initialPDB 2j5k.pdb\n'+\
						'doses 2,3\n'+\
						'PKLMULTIFILE 2j5_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_Frankaer2014(self):
		# Need to create input file for the Frankaer 2014 series
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/Frankaer2014_ETRACK/\n'+\
						'damageset_name 4m4\n'+\
						'damageset_num h,i,j\n'+\
						'initialPDB 4m4f.pdb\n'+\
						'doses 2,3,4\n'+\
						'PKLMULTIFILE 4m4_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_Juers2011_100K(self):
		# Need to create input file for the Frankaer 2014 100K series
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/Juers2011-100K_ETRACK/\n'+\
						'damageset_name 3p7\n'+\
						'damageset_num q,r,s\n'+\
						'initialPDB 3p7p.pdb\n'+\
						'doses 2,3,4\n'+\
						'PKLMULTIFILE 3p7_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_Juers2011_160K(self):
		# Need to create input file for the Frankaer 2014 160K series
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/Juers2011-160K_ETRACK/\n'+\
						'damageset_name 3p7\n'+\
						'damageset_num u,v,w\n'+\
						'initialPDB 3p7t.pdb\n'+\
						'doses 2,3,4\n'+\
						'PKLMULTIFILE 3p7_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_Bury2015(self,version):
		# Need to create input file for the Bury2015 series
		# version in ('final','initial') depending on which PDB_REDO version required
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/Bury2015_ETRACK/{}/\n'.format(version)+\
						'damageset_name 4x4\n'+\
						'damageset_num c,d,e,f,g,h,i\n'+\
						'initialPDB 4x4b.pdb\n'+\
						'doses 6.2,10.3,14.4,20.6,26.8,35.7,44.6\n'+\
						'PKLFILE 3925_4x4c_data.pkl\n'+\
						'PKLFILE 3925_4x4d_data.pkl\n'+\
						'PKLFILE 3925_4x4e_data.pkl\n'+\
						'PKLFILE 3925_4x4f_data.pkl\n'+\
						'PKLFILE 3925_4x4g_data.pkl\n'+\
						'PKLFILE 3925_4x4h_data.pkl\n'+\
						'PKLFILE 3925_4x4i_data.pkl\n'+\
						'PKLMULTIFILE 4x4_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_Weik2000(self,version):
		# Need to create input file for the Weik2000 series
		# version in ('final','initial') depending on which PDB_REDO version required
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/Weik2000_ETRACK/{}/\n'.format(version)+\
						'damageset_name 1qi\n'+\
						'damageset_num e,f,g,h,i,j,k,m\n'+\
						'initialPDB 1qid.pdb\n'+\
						'doses 2,3,4,5,6,7,8,9\n'+\
						'PKLMULTIFILE 1qi_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

	def writeETRACKinputfile_TDIXONinsulin(self):
		# Need to create input file for the TDIXON-insulin series
		inputString =	'where /Users/charlie/DPhil/YEAR2/JAN/TDIXON_InsulinSeries_ETRACK/\n'+\
						'damageset_name insu\n'+\
						'damageset_num 2,3,4,5,6,7,8,9,10\n'+\
						'initialPDB insu1.pdb\n'+\
						'doses 2,3,4,5,6,7,8,9,10\n'+\
						'PKLMULTIFILE insu_data.pkl'
		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()

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
