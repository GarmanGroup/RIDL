from ccp4Job import ccp4Job,checkInputsExist
from mapTools import mapTools
import os

class NCONTjob():

	def __init__(self,
				 inputPDBfile  = 'untitiled.pdb',
				 outputDir     = './',
				 symGroup      = 'P1',
				 seriesName    = 'untitled',
				 printText     = True):

		self.inputPDBfile 	= inputPDBfile
		self.outputDir 		= outputDir
		self.symGroup 		= symGroup
		self.seriesName		= seriesName
		self.printText      = printText

	def run(self):

		# run the ccp4 program NCONT

		if os.path.isfile(self.inputPDBfile) is False:
			return False
		self.runNCONT()
		if self.jobSuccess is True:
			self.provideFeedback()
			return True
		else: return False

	def runNCONT(self):

		# the actual running part here

		title = 'run of NCONT'
		self.commandInput1 = 'ncont '+\
				 			 'XYZIN {} '.format(self.inputPDBfile)
		self.commandInput2 ='source "*"\n'+\
							'target "*"\n'+\
							'mindist 0.0\n'+\
							'maxdist 4.0\n'+\
							'cells 1\n'+\
							'symm {}\n'.format(self.symGroup)+\
							'END'

		self.outputLogfile = '{}-NCONTlogfile.txt'.format(self.seriesName)

		# run NCONT job
		job = ccp4Job('NCONT',self.commandInput1,self.commandInput2,self.outputDir,
					   self.outputLogfile,self.outputDir+'/'+self.outputLogfile)
		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self):

		# provide some feedback

		if self.printText is False:
			return

		print '--------------------------'
		print 'NCONT Summary:'
		print 'Input pdb file: {}'.format(self.inputPDBfile)
		print 'Output log file: {}'.format(self.outputLogfile)
		print '--------------------------'
