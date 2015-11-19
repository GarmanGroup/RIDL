from ccp4Job import ccp4Job
import os

class CADjob():
	# run CAD job to combine F and SIGF columns from two merged mtz files
	def __init__(self,inputMtz1,inputMtz2,inputMtz3,Mtz1LabelName,Mtz2LabelName,Mtz3LabelName,
				 Mtz1LabelRename,Mtz2LabelRename,Mtz3LabelRename,outputMtz,outputDir,runLog):

		self.inputMtz1 			= inputMtz1 # mtz containing initial dataset Fs
		self.inputMtz2 			= inputMtz2 # mtz containing later dataset Fs
		self.inputMtz3 			= inputMtz3 # mtz containing phases
		self.labels 			= [Mtz1LabelName,Mtz2LabelName,Mtz3LabelName]
		self.renameLabels 		= [Mtz1LabelRename,Mtz2LabelRename,Mtz3LabelRename]
		self.outputMtz			= outputMtz
		self.outputDir			= outputDir
		self.runLog 			= runLog
		self.runLog.writeToLog('Running CAD job')

	def run(self):
		if self.checkInputsExist() is False:
			return False
		self.runCAD()
		if self.jobSuccess is True:
			self.provideFeedback()
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

	def checkInputsExist(self):
		# check if input mtz files exist
		for fileName in (self.inputMtz1,self.inputMtz2,self.inputMtz3):
			if os.path.isfile(fileName) is False:
				self.runLog.writeToLog('Failed to find required input file "{}"'.format(fileName))
				return False
		else:
			self.runLog.writeToLog('Input files:')
			self.runLog.writeToLog('{}'.format(self.inputMtz1))	
			self.runLog.writeToLog('{}'.format(self.inputMtz2))
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMtz))
			return True

	def runCAD(self):
		title = 'run of cad'
										 
		self.commandInput1 	= 	'cad '+\
								'HKLIN1 {} '.format(self.inputMtz1)+\
								'HKLIN2 {} '.format(self.inputMtz2)+\
								'HKLIN3 {} '.format(self.inputMtz3)+\
								'HKLOUT {}'.format(self.outputMtz)

		self.commandInput2 	=	'title CAD JOB\n'+\
								'monitor BRIEF\n'+\
								'labin file 1 - \n'+\
								'E1 = F_{} - \n'.format(self.labels[0])+\
								'E2 = SIGF_{} \n'.format(self.labels[0])+\
								'labout file 1 - \n'+\
								'E1 = FP_{} - \n'.format(self.renameLabels[0])+\
								'E2 = SIGFP_{} \n'.format(self.renameLabels[0])+\
								'ctypin file 1 - \n'+\
								'E1 = F - \n'+\
								'E2 = Q \n'+\
								'labin file 2 - \n'+\
								'E1 = F_{} - \n'.format(self.labels[1])+\
								'E2 = SIGF_{} \n'.format(self.labels[1])+\
								'labout file 2 - \n'+\
								'E1 = FP_{} - \n'.format(self.renameLabels[1])+\
								'E2 = SIGFP_{} \n'.format(self.renameLabels[1])+\
								'ctypin file 2 - \n'+\
								'E1 = F - \n'+\
								'E2 = Q \n'+\
								'labin file 3 - \n'+\
								'E1 = PHI{} \n'.format(self.labels[2])+\
								'labout file 3 - \n'+\
								'E1 = PHIC_{} \n'.format(self.renameLabels[2])+\
								'ctypin file 3 - \n'+\
								'E1 = P \n'

		self.outputLogfile = 'CADlogfile.txt'

		# run CAD job
		job = ccp4Job('CAD',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.outputMtz)
		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self):
		# provide some feedback
		print '--------------------------'
		print 'CAD Summary:'
		print 'Input mtz files:\n{}\n{}\n{}'.format(self.inputMtz1,self.inputMtz2,self.inputMtz3)
		print 'Output mtz file:\n {}'.format(self.outputMtz)
		print '--------------------------'
