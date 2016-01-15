import sys
sys.path.insert(0,'./lib')
from processFiles import processFiles

class process():
	# process pdb & mtz files to generate suitable density map files for subsequent analysis
	def __init__(self):
		self.inputFile = 'fullInput.txt'

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