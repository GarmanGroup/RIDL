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

	def runFull(self):
		# run ETRACK from start to finish (map processing & post processing)
		self.runETRACK(True,True,True)

	def runETRACK(self,mapProcess,postProcess,retrieve):
		# run the ETRACK processing for the currently defined input file
		eT = eTrack()
		eT.runPipeline(mapProcess,postProcess,retrieve,self.inputFileName)
		self.et = eT

	def defineDoseList(self,doses,names,version):
		# do not include first dataset dose if difference maps chosen
		if version == 'DIFF':
			if len(doses.split(',')[1:])>1:
				dosesOut = ','.join(doses.split(',')[1:])
				namesOut = ','.join(names.split(',')[1:])
			else:
				dosesOut = doses.split(',')[1]
				namesOut = names.split(',')[1]
			return dosesOut,namesOut
		else: return doses,names

	def writeBlankInputFile(self):
		# Need to create input file for a generic damage series to be completed by the user
		inputString = self.writeInputFile('where <input file location>',
										  'damageset_name <consistent name of series, e.g. for TRAP1.pdb, TRAP2.pdb, this is TRAP>',
										  'damageset_num <dataset numbers, e.g. 1,2 for TRAP. Can be letters corresponding to pdb series>',
										  'initialPDB <initial dataset pdb file e.g. TRAP1.pdb>',
										  'doses <list of doses, length must match length of damageset_num above>',
										  'PKLMULTIFILE <if already processed in ETRACK, can specify single output .pkl for series>')	

	def writeInputFile(self,where,damageset_name,damageset_num,initialPDB,doses,PKLMULTIFILE):
		# write a generic input file for a damage series here
		inputString = 'where {}\n'.format(where)+\
					  'damageset_name {}\n'.format(damageset_name)+\
					  'damageset_num {}\n'.format(damageset_num)+\
					  'initialPDB {}\n'.format(initialPDB)+\
					  'doses {}\n'.format(doses)
		if PKLMULTIFILE != '':
			inputString += 'PKLMULTIFILE {}'.format(PKLMULTIFILE)

		inputFile  = open(self.inputFileName,'w')
		inputFile.write(inputString)
		inputFile.close()	



