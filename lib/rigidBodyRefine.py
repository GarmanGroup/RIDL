import os

class reRefine():

	# the aim of this class is to allow later datasets within a damage 
	# series to be rigid body refined (using REFMAC), using the coordinate
	# model from the first dataset coupled with the observed Fobs columns
	# from the higher dataset mtz files

	def __init__(self,
				 inputFile       = 'untitled.txt',
				 makeNewRIDLfile = True,
				 numCycles       = 10):

		self.inputFile = inputFile
		self.printPurpose()
		self.parseInputFile()

		newPDBs = []
		for mtz,cols in zip(self.mtzList,self.mtzCols):
			self.runREFMAC(pdbIn     = self.initialPDB,
						   mtzIn     = mtz,
						   mtzCols   = cols,
						   rFree     = self.rFree,
						   numCycles = numCycles)

			newPDBs.append(self.refmacPDBout)

		if makeNewRIDLfile is True:
			self.createNewRIDLfile(newPDBs)

	def printPurpose(self):

		# print purpose of script

		ln = 'Running refmac rigid body refinement to generate '+\
			 'higher dose pdb coordinate files specified in RIDL '+\
			 'input file "{}".\nInitial dataset coordiate model and '.format(self.inputFile)+\
			 'higher dose structure factor information have been used' 
		print ln

	def runREFMAC(self,
				  refineType = 'RIGID',
				  pdbIn      = 'untitled.pdb',
				  mtzIn      = 'untitled.mtz',
				  numCycles  = 1,
				  mtzCols    = 'P',
				  rFree      = 'FreeR_Flag',
				  inputDir   = './'):

		# run n = 'numCycles' cycles of 'refineType'
		# refinement in refmac to get phases

		print 'Running refmac ({} {} cycles)...'.format(refineType,numCycles)
		print '---> mtz: "{}"'.format(mtzIn.split('/')[-1])

		self.jobName = 'refmac'
		if refineType == 'RIGID':
			bref = 'over'
			numCycString = 'rigid ncycle {}'.format(numCycles) 
		elif refineType == 'REST':
			bref = 'ISOT'
			numCycString = 'ncyc {}'.format(numCycles) 
		else: 
			print 'Unreadable refinement type.. selecting 0 cycles of rigid body refinement'
			bref = 'over'
			numCycString = 'rigid ncycle 0'
			refineType = 'RIGID'
			numCycles = 0

		# make a refinement type id to append to file names from current job
		fileInd = '_{}_{}cycles'.format(refineType,numCycles)

		# check if input files exist
		if (self.checkFileExists('{}'.format(pdbIn)) is False or
			self.checkFileExists('{}'.format(mtzIn)) is False):			
			return	

		self.refmacPDBout = '{}_refmac{}.pdb'.format(mtzIn.replace('.mtz',''),fileInd)
		self.refmacMTZout = '{}_refmac{}.mtz'.format(mtzIn.replace('.mtz',''),fileInd)
		self.refmacLIBout = '{}_refmac{}.cif'.format(mtzIn.replace('.mtz',''),fileInd)

		self.commandInput1 = 'refmac5 '+\
							 'XYZIN {} '.format(pdbIn)+\
							 'XYZOUT {} '.format(self.refmacPDBout)+\
							 'HKLIN {} '.format(mtzIn)+\
						     'HKLOUT {} '.format(self.refmacMTZout)+\
							 'LIBOUT {} '.format(self.refmacLIBout)
		self.commandInput2   = 	'make check NONE\n'+\
								'make -\n'+\
								'    hydrogen ALL -\n'+\
								'     hout NO -\n'+\
								'     peptide NO -\n'+\
								'    cispeptide YES -\n'+\
								'    ssbridge YES -\n'+\
								'    symmetry YES -\n'+\
								'    sugar YES -\n'+\
								'    connectivity NO -\n'+\
								'    link NO\n'+\
								'refi -\n'+\
								'    type {} -\n'.format(refineType)+\
								'    resi MLKF -\n'+\
								'    meth CGMAT -\n'+\
								'    bref {}\n'.format(bref)+\
								'{}\n'.format(numCycString)+\
								'scal -\n'+\
								'    type SIMP -\n'+\
								'    LSSC -\n'+\
								'    ANISO -\n'+\
								'    EXPE\n'+\
								'solvent YES\n'+\
								'weight -\n'+\
								'    AUTO\n'+\
								'monitor MEDIUM -\n'+\
								'    torsion 10.0 -\n'+\
								'    distance 10.0 -\n'+\
								'    angle 10.0 -\n'+\
								'    plane 10.0 -\n'+\
								'    chiral 10.0 -\n'+\
								'    bfactor 10.0 -\n'+\
								'    bsphere 10.0 -\n'+\
								'    rbond 10.0 -\n'+\
								'    ncsr 10.0\n'+\
								'labin  FP=F{} SIGFP=SIGF{} -\n'.format(mtzCols,mtzCols)+\
								'   FREE={}\n'.format(rFree)+\
								'labout  FC=FC FWT=FWT PHIC=PHIC PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM\n'+\
								'PNAME {}\n'.format(pdbIn.strip('.pdb'))+\
								'DNAME 1\n'+\
								'RSIZE 80\n'+\
								'EXTERNAL WEIGHT SCALE 10.0\n'+\
								'EXTERNAL USE MAIN\n'+\
								'EXTERNAL DMAX 4.2\n'+\
								'END'
		self.outputLogfile = 'REFMAClogfile{}.txt'.format(fileInd)

		# run REFMAC job
		self.runCCP4program()

	def runCCP4program(self):

		# generic method to run a ccp4 program on command line

		textinput = open('{}{}inputfile.txt'.format(self.inputDir,self.jobName),'w')
		textinput.write(self.commandInput2)
		textinput.close()
		os.system('{} < {}{}inputfile.txt > {}{}'.format(self.commandInput1,
												   		   self.inputDir,
												   		   self.jobName,
												   		   self.inputDir,
												   		   self.outputLogfile))

	def checkFileExists(self,
						filename = 'untitled.pdb'):

		# method to check if file exists

		if os.path.isfile(filename) is False:
			ErrorString = 'File {} not found'.format(filename)
			print ErrorString
			return False
		else: 
			return True

	def parseInputFile(self):

		# information of mtz list and columns 
		# can be extracted from an input file

		print 'Reading input file...'

		fileIn = open(self.inputFile,'r')
		for l in fileIn.readlines():
			try:
				l.split()[0]
			except IndexError:
				continue
			if l.split()[0] == 'pdb1':
				self.initialPDB = l.split()[1]
			if l.split()[0] == 'mtz2':
				self.mtzList = (''.join(l.split()[1:])).split(',')
			if l.split()[0] == 'mtzlabels2':
				self.mtzCols = (''.join(l.split()[1:])).split(',')
			if l.split()[0] == 'RfreeFlag1':
				self.rFree = l.split()[1]
		fileIn.close()

		self.inputDir = '/'.join(self.initialPDB.replace('.pdb','').split('/')[0:-1])+'/'

	def createNewRIDLfile(self,
						  newPDBs = []):

		# create a new updated RIDL input file which will 
		# use the new coordinate models for the higher dose
		# datasets in the damage series

		print 'Writing updated RIDL input file...'

		newFile = self.inputFile.replace('.txt','-RigidBodyRefine.txt')

		fileIn = open(self.inputFile,'r')
		fileOut = open(newFile,'w')
		for l in fileIn.readlines():
			try:
				l.split()[0]
			except IndexError:
				continue
			if l.split()[0] == 'pdb2':
				ln = 'pdb2 {}\n'.format(','.join(newPDBs))
			elif l.split()[0] == 'dir':
				ln = 'dir {}_rigidBodyRefine/\n'.format(l.split()[1][0:-1])
			else:
				ln = l
			fileOut.write(ln)

		print '\nNew input file has been generated: "{}"'.format(newFile)

		self.newInputFile = newFile




