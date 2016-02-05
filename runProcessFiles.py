import sys
sys.path.insert(0,'./lib')
from processFiles import processFiles
from dependencies import addPaths

class process():
	# process pdb & mtz files to generate suitable density map files for subsequent analysis
	def __init__(self):
		self.inputFile = 'fullInput.txt'
		addPaths()

	def setInputFile(self,name):
		# specify input file name as required
		self.inputFile = name
		print 'Input file name set as "{}"'.format(self.inputFile)

	def run(self):
		# run to process files with the specified input file
		pro = processFiles(self.inputFile)
		success = pro.runProcessing()
		return success

	def printInputFile(self):
		# print the contents of the specified input file
		fileIn = open(self.inputFile,'r')
		for line in fileIn.readlines():
			if len(line.strip()) != 0:
				print line.strip()
		fileIn.close()

	def writeTemplateInputFile(self):
		# write a template input file to current directory to be completed
		f = open(self.inputFile,'w')
		string   = 	'dir ???\n'+\
					'INITIALDATASET\nname1 ???\nmtz1 ???\nmtzlabels1 ???\npdb1 ???\nRfreeFlag1 ???\n'+\
					'LATERDATASET\nname2 ???\nmtz2 ???\nmtzlabels2 ???\npdb2 ???\n'+\
					'PHASEDATASET\nname3 ???\nmtz3 ???\nmtzlabels3 ???\n'+\
					'MAPINFO\nsfall_VDWR 1\ndensMapType DIFF\nFFTmapWeight True'
		f.write(string)
		f.close()

	def writePDBredoInputFile(self,pdb1,pdb2,inputFileDir,outputDir):
		# write an input file that is suitable for a pdb_redo-downloaded damage series
		# inputFileDir = '/Users/charlie/DPhil/PDBredo_damageSeries/??'
		# outputDir = '/Users/charlie/DPhil/YEAR2/JAN/??'
		self.setInputFile('fullInput_{}-{}DIFF.txt'.format(pdb2,pdb1))
		f = open(self.inputFile,'w')
		string 	 = 	'dir {}\n'.format(outputDir)+\
					'\nINITIALDATASET\n'+\
					'name1 {}init\n'.format(pdb1)+\
					'mtz1 {}/{}/{}.mtz\n'.format(inputFileDir,pdb1,pdb1)+\
					'mtzlabels1 P_{}\n'.format(pdb1)+\
					'pdb1 {}/{}/{}.pdb\n'.format(inputFileDir,pdb1,pdb1)+\
					'RfreeFlag1 FreeR_flag\n'+\
					'\nLATERDATASET\n'+\
					'name2 {}\n'.format(pdb2)+\
					'mtz2 {}/{}/{}.mtz\n'.format(inputFileDir,pdb2,pdb2)+\
					'mtzlabels2 P_{}\n'.format(pdb2)+\
					'pdb2 {}/{}/{}.pdb\n'.format(inputFileDir,pdb2,pdb2)+\
					'\nPHASEDATASET\n'+\
					'name3 {}init\n'.format(pdb1)+\
					'mtz3 {}/{}/{}.mtz\n'.format(inputFileDir,pdb1,pdb1)+\
					'mtzlabels3 C_{}\n'.format(pdb1)+\
					'\nMAPINFO\n'+\
					'sfall_VDWR 1\n'+\
					'densMapType DIFF\n'+\
					'FFTmapWeight True'
		f.write(string)
		f.close()

	def howToWriteInputFile(self):
		# information on how write the input file
		infoString 	= 	'dir --> full path to output working directory\n'+\
						'INITIALDATASET\n'+\
						'name1 --> assign a name to your low dose damage set\n'+\
						'mtz1 --> full path to the low dose mtz file\n'+\
						'mtzlabels1 --> F & SIGF column labels (look in low dose mtz file), if "FP_X" then type "P_X" here for example\n'+\
						'pdb1 --> full path to the low dose pdb file\n'+\
						'RfreeFlag1 --> the Rfree flag label within the low dose mtz file (e.g. "FreeR_flag")\n\n'+\
						'LATERDATASET\n'+\
						'name2 --> assign a name to your high dose damage set (e.g. pdb code)\n'+\
						'mtz2 --> full path to the high dose mtz file\n'+\
						'mtzlabels2 --> F & SIGF column labels (look in high dose mtz file), if "FP_X" then type "P_X" here for example\n'+\
						'pdb2 --> full path to the high dose pdb file\n\n'+\
						'PHASEDATASET\n'+\
						'name3 --> assign a name to your low dose damage set (same as INITIALDATASET above)\n'+\
						'mtz3 --> full path to the low dose mtz file (same as INITIALDATASET above)\n'+\
						'mtzlabels3 --> PHI phase column labels (look in low dose mtz file), if "PHIC_X" then type "C_X" here for example\n\n'+\
						'MAPINFO\n'+\
						'For difference map analysis, do not change this parameters'
		print infoString

	def writeTestInputFile(self,dataset):
		# write an input file for the test case on the github README.md
		inds = ['d','e','f']
		if dataset not in range(1,4):
			print 'Select dataset between 1-3'
			return
		ind = inds[dataset-1]

		self.setInputFile('testInput{}.txt'.format(dataset))
		f = open(self.inputFile,'w')
		string 	 =  'dir ./testOutput/ETRACK\n'+\
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