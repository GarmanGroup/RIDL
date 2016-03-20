from ccp4Job import ccp4Job,checkInputsExist
import os

class SIGMAAjob():
	# run SIGMAA job to combine generate FOM weights if not present in input mtz file
	def __init__(self,inputMtz,MtzLabelNameIn,MtzLabelNameOut,RfreeFlag,inputPDB,outputDir,runLog):

		self.inputMtz 			= inputMtz # mtz containing initial dataset Fs
		self.LabelName 			= MtzLabelNameIn
		self.LabelRename 		= MtzLabelNameOut
		self.RfreeFlag			= RfreeFlag
		self.tmpMtz				= outputDir+'/'+(inputMtz.split('/')[-1]).split('.mtz')[0]+'.tmp'
		self.outputMtz			= self.tmpMtz.split('.tmp')[0]+'_sigmaa.mtz'
		self.inputPDB			= inputPDB
		self.outputDir			= outputDir
		self.runLog 			= runLog
		self.runLog.writeToLog('Running SIGMAA job')

	def run(self):
		inputFiles = [self.inputMtz,self.inputPDB]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False

		success = self.getSpaceGroup() # get space group from pdb file
		if success is False:
			return False

		self.runSIGMAA()
		if self.jobSuccess is True:
			self.provideFeedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMtz))
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

	def runSIGMAA(self):
		title = 'run of sigma-a'
										 
		self.commandInput1 	= 	'sfall '+\
								'HKLIN {} '.format(self.inputMtz)+\
								'HKLOUT {} '.format(self.tmpMtz)+\
								'XYZIN {}'.format(self.inputPDB)

		self.commandInput2 	= 'title {}\n'.format(title)+\
							  'LABIN  FP=F{} SIGFP=SIGF{} FREE={}\n'.format(self.LabelName,self.LabelName,self.RfreeFlag)+\
							  'labout -\n'+\
   							  'FC=FC PHIC=PHIC\n'+\
							  'MODE SFCALC -\n'+\
    						  'XYZIN -\n'+\
    						  'HKLIN\n'+\
							  'symmetry {}\n'.format(self.spaceGroup)+\
							  'NOSCALE\n'+\
							  'end'
					
		self.outputLogfile = 'SIGMAAlogfile1.txt'

		# run SFALL job within SIGMAA routine
		job = ccp4Job('SFALL',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.tmpMtz)
		self.jobSuccess = job.checkJobSuccess()
		if self.jobSuccess is False:
			return

		self.commandInput1 	= 	'sigmaa '+\
								'HKLIN {} '.format(self.tmpMtz )+\
								'HKLOUT {} '.format(self.outputMtz)

		self.commandInput2 	= 'title {}\n'.format(title)+\
							  'LABIN  FP=F{} SIGFP=SIGF{} FC=FC PHIC=PHIC\n'.format(self.LabelName,self.LabelName)+\
							  'labout -\n'+\
   							  'DELFWT=DELFWT_{} FWT=FWT_{} WCMB=FOM_{}\n'.format(self.LabelRename,self.LabelRename,self.LabelRename)+\
							  'ranges 20 -\n'+\
    						  '1000\n'+\
							  'symmetry "{}"\n'.format(self.spaceGroup)+\
							  'partial -\n'+\
    						  'damp 1.0\n'+\
							  'end'
					
		self.outputLogfile = 'SIGMAAlogfile2.txt'

		# run SIGMAA next
		job = ccp4Job('SIGMAA',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.outputMtz)
		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self):
		# provide some feedback
		print '--------------------------'
		print 'SIGMAA Summary:'
		print 'Input mtz file: {}'.format(self.inputMtz)
		print 'Output mtz file: {}'.format(self.outputMtz)
		print '--------------------------'

	def getSpaceGroup(self):
		pdbin = open(self.inputPDB,'r')
		for line in pdbin.readlines():
			if line.split()[0] == 'CRYST1':
				self.spaceGroup = line[55:66].replace(' ','')
				self.runLog.writeToLog('Retrieving space group from file: {}'.format(self.inputPDB))
				self.runLog.writeToLog('Space group determined to be {}'.format(self.spaceGroup))
		try:
			self.spaceGroup
		except attributeError:
			self.runLog.writeToLog('Unable to find space group from file: {}'.format(self.inputPDB))
			return False
		return True

