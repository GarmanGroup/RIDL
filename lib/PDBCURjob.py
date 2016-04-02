from ccp4Job import ccp4Job,checkInputsExist
import os

class PDBCURjob():

	def __init__(self,inputPDBfile,outputDir,runLog):
		self.inputPDBfile 	= inputPDBfile
		self.outputPDBfile 	= '{}/{}_pdbcur.pdb'.format(outputDir,inputPDBfile.split('/')[-1].split('.pdb')[0])
		self.outputDir 		= outputDir
		self.runLog 			= runLog
		self.runLog.writeToLog('Running PDBCUR job')

	def run(self):
		inputFiles = [self.inputPDBfile]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False
		self.runPDBCUR()
		if self.jobSuccess is True:
			self.provideFeedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputPDBfile))
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

	def runPDBCUR(self):
		# PDBcur job to specify to remove hydrogen atoms and pick on the most probably conformation
		# (for two conformations with 0.5 occupancy, the first - A is chosen and occupancy
		# set to 1.00). Also remove all anisou info from file - since it is not needed for 
		# current analysis
		self.jobName = 'PDBCUR'

		# run PDBCUR from command line

		self.commandInput1 = 'pdbcur '+\
				 			 'XYZIN {} '.format(self.inputPDBfile)+\
				 			 'XYZOUT {} '.format(self.outputPDBfile)

		self.commandInput2 = 'delhydrogen\n'+\
				 			 'mostprob\n'+\
				 			 'noanisou\n'+\
				 			 'END'

		self.outputLogfile = 'PDBCURlogfile.txt'

		# run PDBCUR job
		job = ccp4Job('PDBCUR',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.outputPDBfile)
		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self,includeDir=False):
		# provide some feedback

		print '--------------------------'
		print 'PDBCUR Summary:'
		files = {'Input':self.inputPDBfile,
				 'Output':self.outputPDBfile}
		for k in files.keys():
			if includeDir is False:
				f = files[k].split('/')[-1]
			else:
				f = files[k]
			print '{} pdb file: {}'.format(k,f)

			# determine initial number of atoms in pdb file
			pdbin = open(files[k],'r')
			counter = 0
			for line in pdbin.readlines():
				if 'ATOM' in line[0:5]:
					counter += 1
			pdbin.close()
			print '\t# atoms in file: {}'.format(counter)
		print '--------------------------'
