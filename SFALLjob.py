from ccp4Job import ccp4Job
from mapTools import mapTools

class SFALLjob():

	def __init__(self,inputPDBfile,outputDir,VDWR,symmetrygroup,gridDimensions,mapoutType):
		self.inputPDBfile = inputPDBfile
		self.outputMapFile = '{}_sfall.map'.format(inputPDBfile.split('_reordered.pdb')[0])
		self.outputDir = outputDir
		self.symGroup = symmetrygroup
		self.VDWR = VDWR
		self.gridDims = gridDimensions
		self.mapoutType = mapoutType # atom-map ATMMOD or solvent-map SOLVMAP

	def run(self):
		self.runSFALL()
		if self.jobSuccess is True:
			self.provideFeedback()
		else:
			return

	def runSFALL(self):
		title = 'run of sfall'

		self.commandInput1 = '/Applications/ccp4-6.4.0/bin/sfall '+\
				 'XYZIN {} '.format(self.inputPDBfile)+\
				 'ATOMSF /Users/charlie/Desktop/atomsf_addNplus1.lib '+\
			 	 'MAPOUT {} '.format(self.outputMapFile)+\
				 'SYMINFO /Applications/ccp4-6.4.0/lib/data/syminfo.lib '

		if len(self.gridDims) == 0:
			self.commandInput2 = 'MODE ATMMAP {}\n'.format(self.mapoutType)+\
							 	 'SYMMETRY {}\n'.format(self.symGroup)+\
							 	 'VDWR {}\n'.format(self.VDWR)+\
							 	 'title {}\n'.format(title)+\
				   			 	 'END'
		else:
			self.commandInput2 = 'MODE ATMMAP {}\n'.format(self.mapoutType)+\
							 	 'SYMMETRY {}\n'.format(self.symGroup)+\
							 	 'VDWR {}\n'.format(self.VDWR)+\
							 	 'title {}\n'.format(title)+\
							 	 'GRID {} {} {}\n'.format(self.gridDims[0],self.gridDims[1],self.gridDims[2])+\
							 	 'END'

		self.outputLogfile = 'SFALLlogfile.txt'

		# run SFALL job
		job = ccp4Job('SFALL',self.commandInput1,self.commandInput2,self.outputDir,self.outputLogfile,self.outputMapFile)
		self.jobSuccess = job.checkJobSuccess()

	def provideFeedback(self):
		# provide some feedback
		print '--------------------------'
		print 'SFALL Summary:'
		print 'Input pdb file: {}'.format(self.inputPDBfile)
		print 'Output map file: {}'.format(self.outputMapFile)
		Map = mapTools(self.outputMapFile)
		Map.printMapInfo()
		print '--------------------------'



