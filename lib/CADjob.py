from ccp4Job import ccp4Job,checkInputsExist,fillerLine
import os

class CADjob():

	# run CAD job to combine F and SIGF columns from two merged mtz files
	
	def __init__(self,
				 inputMtz1       = '',
				 inputMtz2       = '',
				 inputMtz3       = '',
				 Mtz1LabelName   = '',
				 Mtz2LabelName   = '',
				 Mtz3phaseLabel  = '',
				 Mtz3FcalcLabel  = '',
				 Mtz1LabelRename = '',
				 Mtz2LabelRename = '',
				 Mtz3LabelRename = '',
				 outputMtz       = '',
				 outputDir       = './',
				 runLog          = '',
				 FOMWeight       = 'False'):

		self.inputMtz1 	  = inputMtz1 # mtz containing initial dataset Fs
		self.inputMtz2 	  = inputMtz2 # mtz containing later dataset Fs
		self.inputMtz3 	  = inputMtz3 # mtz containing phases
		self.labels 	  = [Mtz1LabelName,
							 Mtz2LabelName,
							 Mtz3phaseLabel,
							 Mtz3FcalcLabel]
		self.renameLabels = [Mtz1LabelRename,
							 Mtz2LabelRename,
							 Mtz3LabelRename]
		self.outputMtz	  = outputMtz
		self.outputDir	  = outputDir
		self.runLog 	  = runLog
		
		self.FOMtag = {'in'   : '',
					   'out'  : '',
					   'type' : ''}

		if FOMWeight != 'False':
			self.FOMtag['out'] = '- \nE3 = FOM_{}'.format(self.renameLabels[0])
			self.FOMtag['type'] = '- \nE3 = W'
		if FOMWeight == 'recalculate':
			self.FOMtag['in'] = '- \nE3 = FOM{}'.format(self.labels[0])
		elif 'preset' in FOMWeight:
			lbl = FOMWeight.split(',')[-1]
			self.FOMtag['in'] = '- \nE3 = FOM{}'.format(lbl)

	def run(self):
		
		inputFiles = [self.inputMtz1,
					  self.inputMtz2,
					  self.inputMtz3]

		if checkInputsExist(inputFiles,self.runLog) is False:
			return False

		self.runCAD()
		if self.jobSuccess is True:
			self.provideFeedback()
			return True

		else:
			err = 'Job did not run successfully, see job log file "{}"'.format(self.outputLogfile)
			self.runLog.writeToLog(err)
			return False

	def runCAD(self):

		# fillerLine()
		self.printPurpose()
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
								'E2 = SIGF{} '.format(self.labels[0])+\
								'{} \n'.format(self.FOMtag['in'])+\
								'labout file 1 - \n'+\
								'E1 = FP_{} - \n'.format(self.renameLabels[0])+\
								'E2 = SIGFP_{} '.format(self.renameLabels[0])+\
								'{} \n'.format(self.FOMtag['out'])+\
								'ctypin file 1 - \n'+\
								'E1 = F - \n'+\
								'E2 = Q '+\
								'{} \n'.format(self.FOMtag['type'])+\
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
								'E1 = {} -\n'.format(self.labels[2])+\
								'E2 = {} \n'.format(self.labels[3])+\
								'labout file 3 - \n'+\
								'E1 = PHIC_{} - \n'.format(self.renameLabels[2])+\
								'E2 = FC_{} \n'.format(self.renameLabels[2])+\
								'ctypin file 3 - \n'+\
								'E1 = P -\n'+\
								'E2 = F \n'

		self.outputLogfile = 'CADlogfile.txt'

		# run CAD job
		job = ccp4Job(jobName       = 'CAD',
					  commandInput1 = self.commandInput1,
					  commandInput2 = self.commandInput2,
					  outputDir     = self.outputDir,
					  outputLog     = self.outputLogfile,
					  outputFile    = self.outputMtz)

		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self,
					    includeDir = False):

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

		txt = 'CAD Summary:\n'+\
			  'Input mtz files: {}\n{}{}\n{}{}\n'.format(fileIn1,' '*17,fileIn2,' '*17,fileIn3)+\
			  'Output mtz file: {}'.format(fileOut)
		self.runLog.writeToLog(txt)


	def printPurpose(self,
					 include = True):

		# provide a summary of what this does 
		# (within ETRACK) to the command line

		txt = 'Combining relevant columns from multiple '+\
			  'input mtz files into 1 single mtz file'
		self.runLog.writeToLog(txt)
