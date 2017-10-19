from ccp4Job import ccp4Job, checkInputsExist


class PDBCURjob():

    def __init__(self,
                 inputPDBfile='', outputDir='./', runLog=''):

        self.inputPDBfile = inputPDBfile
        self.outputPDBfile = '{}{}_pdbcur.pdb'.format(
            outputDir, inputPDBfile.split('/')[-1].split('.pdb')[0])
        self.outputDir = outputDir
        self.runLog = runLog

    def run(self):

        inputFiles = [self.inputPDBfile]

        if not checkInputsExist(inputFiles, self.runLog):
            return False
        self.runPDBCUR()
        if self.jobSuccess:
            self.provideFeedback()
            return True
        else:
            self.runLog.writeToLog(
                'Job did not run successfully, see ' +
                'job log file "{}"'.format(self.outputLogfile))
            return False

    def runPDBCUR(self):

        # PDBcur job to specify to remove hydrogen atoms and pick
        # on the most probably conformation (for two conformations
        # with 0.5 occupancy, the first - A is chosen and occupancy
        # set to 1.00). Also remove all anisou info from file
        # - since it is not needed for current analysis

        self.printPurpose()
        self.jobName = 'PDBCUR'

        # run PDBCUR from command line

        self.commandInput1 = 'pdbcur ' +\
                             'XYZIN {} '.format(self.inputPDBfile) +\
                             'XYZOUT {} '.format(self.outputPDBfile)

        self.commandInput2 = 'delhydrogen\nmostprob\nnoanisou\nEND'

        self.outputLogfile = 'PDBCURlogfile.txt'

        # run PDBCUR job
        job = ccp4Job(jobName='PDBCUR', commandInput1=self.commandInput1,
                      commandInput2=self.commandInput2,
                      outputDir=self.outputDir, outputLog=self.outputLogfile,
                      outputFile=self.outputPDBfile)

        self.jobSuccess = job.checkJobSuccess()

    def provideFeedback(self,
                        includeDir=False):

        # provide some feedback

        txt = 'PDBCUR Summary:\n'

        files = {'Input': self.inputPDBfile, 'Output': self.outputPDBfile}

        for k in files.keys():
            if not includeDir:
                f = files[k].split('/')[-1]
            else:
                f = files[k]
            txt += '{} pdb file: {}\n'.format(k, f)

            # determine initial number of atoms in pdb file
            pdbin = open(files[k], 'r')
            counter = 0
            for line in pdbin.readlines():
                if 'ATOM' in line[0:5]:
                    counter += 1
            pdbin.close()
            txt += '\t# atoms in file: {}\n'.format(counter)

        self.runLog.writeToLog(txt)

    def printPurpose(self,
                     include=True):

        # provide a summary of what this does
        # (within RIDL) to the command line

        self.runLog.writeToLog(
            'Stripping hydrogens, removing secondary conformations, ' +
            'and removing anisotropic B-factors')
