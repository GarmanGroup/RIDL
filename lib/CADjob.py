from ccp4Job import ccp4Job,checkInputsExist
import os

class CADjob():
	# run CAD job to combine F and SIGF columns from two merged mtz files
	def __init__(self,inputMtz1,inputMtz2,inputMtz3,Mtz1LabelName,Mtz2LabelName,Mtz3LabelName,
				 Mtz1LabelRename,Mtz2LabelRename,Mtz3LabelRename,outputMtz,outputDir,runLog,FOMWeight):

		self.inputMtz1 			= inputMtz1 # mtz containing initial dataset Fs
		self.inputMtz2 			= inputMtz2 # mtz containing later dataset Fs
		self.inputMtz3 			= inputMtz3 # mtz containing phases
		self.labels 			= [Mtz1LabelName,Mtz2LabelName,Mtz3LabelName]
		self.renameLabels 		= [Mtz1LabelRename,Mtz2LabelRename,Mtz3LabelRename]
		self.outputMtz			= outputMtz
		self.outputDir			= outputDir
		self.runLog 			= runLog
		self.runLog.writeToLog('Running CAD job')

		if FOMWeight == 'False':
			self.FOMtagIn = ''
			self.FOMtagOut = ''
		elif FOMWeight == 'recalculate':
			self.FOMtagIn = 'E3 = FOM{} \n'.format(self.labels[0])
			self.FOMtagOut = 'E3 = FOM_{} \n'.format(self.renameLabels[0])
		elif 'preset' in FOMWeight:
			lbl = FOMWeight.split(',')[-1]
			self.FOMtagIn = 'E3 = FOM{} \n'.format(lbl)
			self.FOMtagOut = 'E3 = FOM_{} \n'.format(self.renameLabels[0])

	def run(self):
		inputFiles = [self.inputMtz1,self.inputMtz2,self.inputMtz3]
		if checkInputsExist(inputFiles,self.runLog) is False:
			return False
		self.runCAD()
		if self.jobSuccess is True:
			self.provideFeedback()
			self.runLog.writeToLog('Output files:')	
			self.runLog.writeToLog('{}'.format(self.outputMtz))
			return True
		else:
			self.runLog.writeToLog('Job did not run successfully, see job log file "{}"'.format(self.outputLogfile))
			return False

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
								'E1 = F{} - \n'.format(self.labels[0])+\
								'E2 = SIGF{} - \n'.format(self.labels[0])+\
								'{}'.format(self.FOMtagIn)+\
								'labout file 1 - \n'+\
								'E1 = FP_{} - \n'.format(self.renameLabels[0])+\
								'E2 = SIGFP_{} - \n'.format(self.renameLabels[0])+\
								'{}'.format(self.FOMtagOut)+\
								'ctypin file 1 - \n'+\
								'E1 = F - \n'+\
								'E2 = Q - \n'+\
								'E3 = W \n'+\
								'labin file 2 - \n'+\
								'E1 = F{} - \n'.format(self.labels[1])+\
								'E2 = SIGF{} \n'.format(self.labels[1])+\
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

	def provideFeedback(self,includeDir=False):
		# provide some feedback
		if includeDir is False:
			fileIn1  = self.inputMtz1.split('/')[-1]
			fileIn2  = self.inputMtz2.split('/')[-1]
			fileIn3  = self.inputMtz3.split('/')[-1]
			fileOut = self.outputMtz.split('/')[-1]
		else:
			fileIn1  = self.inputMtz1
			fileIn2  = self.inputMtz2
			fileIn3  = self.inputMtz3
			fileOut = self.outputMtz
		print '--------------------------'
		print 'CAD Summary:'
		print 'Input mtz files: {}\n{}{}\n{}{}'.format(fileIn1,' '*17,fileIn2,' '*17,fileIn3)
		print 'Output mtz file: {}'.format(fileOut)
		print '--------------------------'
