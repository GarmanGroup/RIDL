import os
import sys
FULL_PATH = os.path.abspath(os.path.dirname(sys.argv[0]))
sys.path.insert(0, FULL_PATH+'/lib')
from errors import error


class reRefine():

    # the aim of this class is to allow later datasets within a damage
    # series to be rigid body refined (using REFMAC), using the coordinate
    # model from the first dataset coupled with the observed Fobs columns
    # from the higher dataset mtz files

    def __init__(
        self, inputFile='untitled.txt', makeNewRIDLfile=True,
            numCycles=10):

        self.inputFile = inputFile
        self.printPurpose()
        self.parseInputFile()

        newPDBs = []
        for mtz, cols in zip(self.mtzList, self.mtzCols):
            self.runREFMAC(
                pdbIn=self.initialPDB, mtzIn=mtz,
                mtzCols=cols, rFree=self.rFree,
                numCycles=numCycles)

            newPDBs.append(self.refmacPDBout)

            # check the output log indicates success
            with open(self.outputLogfile, 'r') as log:
                for ln in log.readlines():
                    if 'Refmac:  Error' in ln:
                        error(text='Refmac did not run to completion. ' +
                                   'Please refer to the log file for ' +
                                   'Refmac "{}"'.format(self.outputLogfile))

            # check if required files from refinement are present
            self.checkFileExists(self.refmacPDBout)

        if makeNewRIDLfile:
            self.createNewRIDLfile(newPDBs)

    def printPurpose(self):

        # print purpose of script

        ln = 'Running refmac rigid body refinement to generate ' +\
             'higher dose pdb coordinate files specified in RIDL ' +\
             'input file "{}".\nInitial dataset coordiate model and '.format(
                self.inputFile) +\
             'higher dose structure factor information have been used'
        print(ln)

    def runREFMAC(
        self, refineType='RIGID', pdbIn='untitled.pdb',
        mtzIn='untitled.mtz', numCycles=1, mtzCols='FP',
            rFree='FreeR_Flag', inputDir='./'):

        # run n = 'numCycles' cycles of 'refineType'
        # refinement in refmac to get phases

        print('Running refmac ({} {} cycles)...'.format(refineType, numCycles))
        print('---> mtz: "{}"'.format(mtzIn.split('/')[-1]))

        # check if input files exist and stop if not
        print('--> checking that all files for rigid-body refinement are present')
        self.checkFileExists(str(pdbIn))
        self.checkFileExists(str(mtzIn))

        self.jobName = 'refmac'
        if refineType == 'RIGID':
            bref = 'over'
            numCycString = 'rigid ncycle {}'.format(numCycles)
        elif refineType == 'REST':
            bref = 'ISOT'
            numCycString = 'ncyc {}'.format(numCycles)
        else:
            print('Warning: unreadable refinement type.. ' +
                  'selecting 0 cycles of rigid body refinement')
            bref = 'over'
            numCycString = 'rigid ncycle 0'
            refineType = 'RIGID'
            numCycles = 0

        # make a refinement type id to append to file names from current job
        fileInd = '_{}_{}cycles'.format(refineType, numCycles)

        self.refmacPDBout = '{}_refmac{}.pdb'.format(
            mtzIn.replace('.mtz', ''), fileInd)
        self.refmacMTZout = '{}_refmac{}.mtz'.format(
            mtzIn.replace('.mtz', ''), fileInd)
        self.refmacLIBout = '{}_refmac{}.cif'.format(
            mtzIn.replace('.mtz', ''), fileInd)

        self.commandInput1 = 'refmac5 ' +\
                             'XYZIN {} '.format(pdbIn) +\
                             'XYZOUT {} '.format(self.refmacPDBout) +\
                             'HKLIN {} '.format(mtzIn) +\
                             'HKLOUT {} '.format(self.refmacMTZout) +\
                             'LIBOUT {} '.format(self.refmacLIBout)
        self.commandInput2 = 'make check NONE\n' +\
                             'make -\n' +\
                             '    hydrogen ALL -\n' +\
                             '     hout NO -\n' +\
                             '     peptide NO -\n' +\
                             '    cispeptide YES -\n' +\
                             '    ssbridge YES -\n' +\
                             '    symmetry YES -\n' +\
                             '    sugar YES -\n' +\
                             '    connectivity NO -\n' +\
                             '    link NO\n' +\
                             'refi -\n' +\
                             '    type {} -\n'.format(refineType) +\
                             '    resi MLKF -\n' +\
                             '    meth CGMAT -\n' +\
                             '    bref {}\n'.format(bref) +\
                             '{}\n'.format(numCycString) +\
                             'scal -\n' +\
                             '    type SIMP -\n' +\
                             '    LSSC -\n' +\
                             '    ANISO -\n' +\
                             '    EXPE\n' +\
                             'solvent YES\n' +\
                             'weight -\n' +\
                             '    AUTO\n' +\
                             'monitor MEDIUM -\n' +\
                             '    torsion 10.0 -\n' +\
                             '    distance 10.0 -\n' +\
                             '    angle 10.0 -\n' +\
                             '    plane 10.0 -\n' +\
                             '    chiral 10.0 -\n' +\
                             '    bfactor 10.0 -\n' +\
                             '    bsphere 10.0 -\n' +\
                             '    rbond 10.0 -\n' +\
                             '    ncsr 10.0\n' +\
                             'labin  FP={} SIGFP=SIG{} -\n'.format(mtzCols, mtzCols) +\
                             '   FREE={}\n'.format(rFree) +\
                             'labout  FC=FC FWT=FWT PHIC=PHIC PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM\n' +\
                             'PNAME {}\n'.format(pdbIn.strip('.pdb')) +\
                             'DNAME 1\n' +\
                             'RSIZE 80\n' +\
                             'EXTERNAL WEIGHT SCALE 10.0\n' +\
                             'EXTERNAL USE MAIN\n' +\
                             'EXTERNAL DMAX 4.2\n' +\
                             'END'
        self.outputLogfile = '{}REFMAClogfile{}.txt'.format(
            self.inputDir, fileInd)

        # run REFMAC job
        self.runCCP4program()

    def runCCP4program(self):

        # generic method to run a ccp4 program on command line

        textinput = open('{}{}inputfile.txt'.format(
            self.inputDir, self.jobName), 'w')
        textinput.write(self.commandInput2)
        textinput.close()
        os.system('{} < {}{}inputfile.txt > {}'.format(
            self.commandInput1, self.inputDir, self.jobName,
            self.outputLogfile))

    def checkFileExists(
            self, filename='untitled.pdb'):

        # method to check if file exists

        if not os.path.isfile(filename):
            error(text='file "{}" could not be located'.format(filename),
                  type='error')
            return False
        else:
            return True

    def parseInputFile(self):

        # information of mtz list and columns
        # can be extracted from an input file

        print('Reading input file...')

        fileIn = open(self.inputFile, 'r')
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

        self.inputDir = '/'.join(
            self.initialPDB.replace('.pdb', '').split('/')[0:-1])+'/'

        # perform basic checks here
        if len(self.mtzList) != len(self.mtzCols):
            error(text='Inconsistent number of comma-separated entries in ' +
                       'input file for line "mtz2" and "mtzlabels2"')

        try:
            self.rFree
        except AttributeError:
            error('"RfreeFlag1" has not been specified in input file')

        try:
            self.initialPDB
        except AttributeError:
            error('"pdb1" has not been specified in the input file')

        self.checkFileExists(self.initialPDB)
        for mtz in self.mtzList:
            self.checkFileExists(mtz)

    def createNewRIDLfile(
            self, newPDBs=[]):

        # create a new updated RIDL input file which will
        # use the new coordinate models for the higher dose
        # datasets in the damage series

        print('Writing updated RIDL input file...')

        newFile = self.inputFile.replace('.txt', '-RigidBodyRefine.txt')

        fileIn = open(self.inputFile, 'r')
        fileOut = open(newFile, 'w')
        for l in fileIn.readlines():
            try:
                l.split()[0]
            except IndexError:
                continue
            if l.split()[0] == 'pdb2':
                continue
            if l.split()[0] == 'mtzlabels2':
                ln = l+'pdb2 {}\n'.format(','.join(newPDBs))
            elif l.split()[0] == 'dir':
                ln = 'dir {}_rigidBodyRefine/\n'.format(l.split()[1][0:-1])
            else:
                ln = l
            fileOut.write(ln)
        fileIn.close()
        fileOut.close()

        print('\nNew input file has been generated: "{}"'.format(newFile))

        self.newInputFile = newFile
