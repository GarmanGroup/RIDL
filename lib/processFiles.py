from makeMapsFromMTZs import makeMapsFromMTZs
from cleanUpFiles import cleanUpFinalFiles
from errors import error
import difflib
import shutil
import os


class processFiles():

    def __init__(self,
                 inputFile='', proceedToMetricCalc=False,
                 skipToMetricCalc=False, outputGraphs='yes',
                 metCalcInput='metricCalc_inputfile.txt',
                 cleanFinalFiles=False, logFileObj='',
                 skipToSummaryFiles=False, writeSummaryFiles=False,
                 keepMapDir=True, includeSIGF=True):

        # class to read an input file and generate a set of density
        # and atom-tagged maps for a damage series. Can handle a single
        # low & high dose dataset pair, or a low dose dataset and a series
        # of high dose datasets. Can optionally proceed directly on to
        # the metric calculation stage of the pipeline with
        # 'proceedToMetricCalc' or skip directly to it (if the maps have
        # already been created) with 'skipToMetricCalc'.

        self.inputFile = inputFile
        self.outputGraphs = outputGraphs
        self.metCalcInput = metCalcInput
        self.logFile = logFileObj
        self.skipToSumFiles = skipToSummaryFiles
        self.writeSumFiles = writeSummaryFiles
        self.includeSIGF = includeSIGF

        self.props1 = ['name1', 'mtz1', 'mtzlabels1', 'pdb1', 'RfreeFlag1']

        self.props2 = ['name2', 'mtz2', 'mtzlabels2', 'pdb2']

        self.props3 = ['name3', 'mtz3', 'phaseLabel', 'FcalcLabel']

        if not skipToMetricCalc:
            success = self.runFileProcessing()
            if success:
                success = self.runMapGeneration()
                if success:
                    if proceedToMetricCalc:
                        self.runMetricCalcStep()

        else:
            success = self.skipToMetricCalc()

        if success and cleanFinalFiles:
            if skipToMetricCalc or proceedToMetricCalc:
                cleanUpFinalFiles(outputDir=self.dir,
                                  cleanMapDir=not skipToSummaryFiles,
                                  keepMapDir=keepMapDir)

        self.jobSuccess = success

    def runFileProcessing(self):

        self.logFile.writeToLog(str='\n**** INPUT FILE PROCESSING ****\n')

        success = self.readMainInputFile()
        if not success:
            return

        self.checkOutputDirsExists()
        self.findFilesInDir()
        self.checkForMultipleDatasets()

        # don't proceed if error in input file
        try:
            self.multiDatasets
        except AttributeError:
            return
        if self.multiDatasets:
            try:
                self.repeatedFile1InputsUsed
                self.repeatedPhaseInputsUsed
            except AttributeError:
                return

        success = self.checkCorrectInputFormats()
        if not success:
            return success

        self.checkForSeparateSIGFlabel()

        success = self.checkMtzLabelsExist()
        if not success:
            return success

        return success

    def runMapGeneration(self):

        self.logFile.writeToLog(str='\n\n**** DENSITY MAP GENERATION ****\n')

        if not self.multiDatasets:
            self.logFile.writeToLog(
                str='Generating suitable electron density maps for run.')
            self.getCurrentInputParams()
            success = self.runMapGenerationPipeline()

        else:
            self.logFile.writeToLog(
                str='Generating suitable electron density maps for each ' +
                    'higher dose dataset in turn.\n{} '.format(self.numDsets) +
                    'higher dose datasets located for this job.')

            for i in range(self.numDsets):
                self.logFile.writeToLog(
                    str='\n{}\nHigher dose '.format('-'*33) +
                        'dataset {} starts here'.format(i+1))
                self.getCurrentInputParams(jobNumber=i)
                success = self.runMapGenerationPipeline()
                if not success:
                    return success

        return success

    def runMapGenerationPipeline(self):

        # run the map generation pipeline for a single dataset

        self.setJobName()

        if self.sepSIGFPlabel1:
            sigFP1 = self.mtzSIGFPlabel1_current
        else:
            sigFP1 = self.mtzlabels1_current

        if self.sepSIGFPlabel2:
            sigFP2 = self.mtzSIGFPlabel2_current
        else:
            sigFP2 = self.mtzlabels2_current

        p = makeMapsFromMTZs(
            outputDir=self.mapProcessDir, jobName=self.jobName,
            log=self.logFile, mtzIn1=self.mtz1_current,
            Mtz1LabelName=self.mtzlabels1_current,
            RfreeFlag1=self.RfreeFlag1_current, Mtz1SIGFPlabel=sigFP1,
            Mtz1LabelRename=self.name1_current, mtzIn2=self.mtz2_current,
            Mtz2LabelName=self.mtzlabels2_current, Mtz2SIGFPlabel=sigFP2,
            Mtz2LabelRename=self.name2_current, mtzIn3=self.mtz3_current,
            Mtz3phaseLabel=self.phaseLabel_current,
            Mtz3FcalcLabel=self.FcalcLabel_current,
            Mtz3LabelRename=self.name3_current, inputPDBfile=self.pdb2_current,
            densMapType=self.densMapType, scaleType=self.scaleType,
            FOMweight=self.FFTmapWeight,
            deleteMtzs=self.deleteIntermediateFiles,
            sfall_VDWR=self.sfall_VDWR, mapResLimits=self.mapResLimits,
            includeFCmaps=self.calculateFCmaps,
            useLaterCellDims=self.useLaterCellDims)

        success = p.runPipeline()
        self.mtzToMapsPipelineLog = p.runLog.logFile
        if success:
            self.logFile.writeToLog(str='---> Subroutine ran to completion.')
            if self.deleteIntermediateFiles.lower() == 'true':
                success = self.cleanUpIntermediateFiles()
            else:
                success = self.cleanUpIntermediateFiles(
                        removeMtzs=False, removeMaps=False, removePdbs=False)
        else:
            self.writeError(text='Subroutine failed to run to completion')

        return success

    def runMetricCalcStep(self,
                          write=True, run=True, useImports=True,
                          skipToMetricCalc=False):

        # writes an input file for the run of the metric calculation
        # part of RIDL pipeline after the map generation pipeline has
        # completed. Allows the user to immediately run the rest of the
        # program, IF the user did specify all required datasets
        # within a damage series. If the user only processed a subset of
        # required datasets within a damage series, another metric
        # calculation input file will need to be manually created to
        # run the rest of the program with ALL datasets within the
        # damage series included.

        if useImports:
            from runRIDL_metricCalc import run as run_metricCalc

        if not skipToMetricCalc:
            if not write:
                return True

            r = run_metricCalc(calculate=False)
            seriesName = self.dir.split('/')[-2]

            if self.multiDatasets and not self.repeatedFile1InputsUsed:
                # currently if more than one initial dataset specified then
                # parts of further analysis are untested. Derivation of Dloss
                # values should not be affected

                self.writeError(
                    text='More than one INITIALDATASET input specified! ' +
                         'Taking only first coordinate file specified in ' +
                         'input file. This file will used for the x,y,z ' +
                         'coordinates of each input atom. If you are  using ' +
                         'outputs from RIDL other than Dloss, please contact' +
                         ' charles.bury@dtc.ox.ac.uk for concerns ' +
                         'over validity.',
                    type='warning')

            if self.dose2 == 'NOTCALCULATED':
                doses = ','.join(map(str, range(len(self.name2.split(',')))))
            else:
                doses = self.dose2

            r.writeInputFile(inDir=self.mapProcessDir, outDir=self.dir,
                             damSetName=seriesName, laterDatasets=self.name2,
                             initialPDB=self.name1, doses=doses,
                             outputGraphs=self.outputGraphs)

            shutil.move(r.inputFileName,
                        '{}/{}'.format(self.dir, r.inputFileName))

        if run:
            if self.calculateFCmaps.upper() == 'FALSE':
                inclFCmets = False
            else:
                inclFCmets = True

            r = run_metricCalc(inputFileLoc=self.dir,
                               calculate=not self.skipToSumFiles,
                               logFile=self.logFile,
                               skipToSumFiles=self.skipToSumFiles,
                               writeSumFiles=self.writeSumFiles,
                               inclFCmets=inclFCmets)

        return True

    def skipToMetricCalc(self):

        # do not generate maps for job, but proceed
        # directly to metric calculations, if correct
        # input file exists in working directory

            success = self.readMainInputFile()
            if not success:
                return False

            self.checkOutputDirsExists(makeProcessDir=False)
            self.findFilesInDir(mapProcessDir=False)

            if self.metCalcInput not in self.filesInDir:
                self.writeError(
                    text='Unable to find input file ' +
                         '"{}" in "{}"\n'.format(self.metCalcInput, self.dir) +
                         'Ensure "python runRIDL.py -i <input.txt> ' +
                         '-p" has been run prior to -c')
                return False

            success = self.runMetricCalcStep(run=True, useImports=True,
                                             skipToMetricCalc=True)
            return True

    def readMainInputFile(self):

        # split the input file into separate input
        # files for each section of the pipeline

        # if Input.txt not found, flag error
        if not os.path.isfile(self.inputFile):
            self.writeError(
                text='Input file "{}"" not found..'.format(self.inputFile))
            return False
        else:
            self.logFile.writeToLog(
                str='Reading input file "{}"'.format(self.inputFile))

        fileIn = open(self.inputFile, 'r')

        # parse input file and ignore blank/comment lines
        for line in fileIn.readlines():
            try:
                line.split()[1]
            except IndexError:
                continue
            if line.strip()[0] == '#':
                continue
            inputPart = ''.join(line.split()[1:])
            setattr(self, line.split()[0], inputPart)
        fileIn.close()

        self.checkNonNecessaryInputs()
        success = self.checkAllRequiredInputsFound()
        return success

    def checkAllRequiredInputsFound(self):

        # check that all required properties have been found

        requiredProps = self.props1+self.props2+self.props3
        requiredProps += ['dir', 'dose1', 'dose2']

        for prop in requiredProps:
            try:
                getattr(self, prop)
            except AttributeError:
                self.writeError(
                    text='Necessary input not found: {}'.format(prop))
                return False

        self.logFile.writeToLog(
            str='All necessary inputs found in input file.')

        return True

    def checkNonNecessaryInputs(self):

        # check whether non necessary inputs are present
        # and set to default values if not specified

        props = ['sfall_VDWR', 'mapResLimits', 'scaleType',
                 'densMapType', 'FFTmapWeight', 'calculateFCmaps',
                 'deleteIntermediateFiles', 'useLaterCellDims']

        defaults = [1, ',', 'ANISOTROPIC', 'DIFF',
                    'False', 'FALSE', 'TRUE', 'TRUE']

        for i, prop in enumerate(props):
            try:
                getattr(self, prop)
            except AttributeError:
                setattr(self, prop, defaults[i])
                self.logFile.writeToLog(
                    str='Input "{}" not found in input file. '.format(prop) +
                        'Setting to default value: "{}"'.format(defaults[i]))

    def checkCorrectInputFormats(self):

        # check that the required properties have correct formats

        props = ['mtz1', 'pdb1', 'mtz2', 'pdb2', 'mtz3']

        fType = ['.mtz', '.pdb']*2+['.mtz']

        for p, t in zip(props, fType):
            if not self.multiDatasets:

                name = getattr(self, p)
                success = self.checkCorrectFileExtension(
                    fileName=name, fileType=t, property=p)
                if not success:
                    return False

                success = self.checkFileExists(fileName=name)
                if not success:
                    return False
            else:
                for f in getattr(self, p).split(','):

                    success = self.checkCorrectFileExtension(
                        fileName=f, fileType=t, property=p)
                    if not success:
                        return False

                    success = self.checkFileExists(fileName=f)
                    if not success:
                        return False

        if not self.multiDatasets:
            if self.name1 == self.name2:
                self.writeError(
                    text='"name1" and "name2" inputs must be different, ' +
                         'otherwise CAD will fail. Both currently set ' +
                         'as "{}".'.format(self.name1))
                return False

            for n in ('name1', 'name2'):
                success = self.checkNameLength(
                    name=getattr(self, n), property=n)
                if not success:
                    return False
        else:
            isError = False
            if self.repeatedFile1InputsUsed:
                success = self.checkNameLength(
                    name=self.name1, property='name1')
                if not success:
                    return False

                for n in self.name2.split(','):
                    if self.name1 == n:
                        isError = True
                        sameName = n
                        break

                    success = self.checkNameLength(name=n, property='name2')
                    if not success:
                        return False

            else:
                nameZip = zip(self.name1.split(','), self.name2.split(','))
                for n1, n2 in nameZip:
                    if n1 == n2:
                        isError = True
                        sameName = n1
                        break

                    for n, m in zip([n1, n2], ['name1', 'name2']):
                        success = self.checkNameLength(name=n, property=m)
                        if not success:
                            return False

            if isError:
                self.writeError(
                    text='"name1" and "name2" inputs must be different ' +
                         'for each job in batch, otherwise CAD will fail. ' +
                         'Currently for one batch, the "initial" and "later"' +
                         ' datasets are both called "{}".'.format(sameName))
                return False

        if self.multiDatasets:
            if self.repeatedPhaseInputsUsed:
                names = [self.name3]
            else:
                names = self.name3.split(',')
        else:
            names = [self.name3]

        for name in names:
            success = self.checkNameLength(
                name=name, property='name3', maxLength=22)
            if not success:
                return False

        # add '-ph' to name3 to identify it from other datasets
        if not self.multiDatasets:
            self.name3 += '-ph'
        else:
            if self.repeatedPhaseInputsUsed:
                self.name3 += '-ph'
            else:
                self.name3 = ','.join([n+'-ph' for n in self.name3.split(',')])

        if self.densMapType not in ('DIFF', 'SIMPLE', '2FOFC', 'END'):
            self.writeError(
                text='"densMapType" input of incompatible format, ' +
                     '(default is "DIFF"). Currently set as ' +
                     '"{}" in input file.'.format(self.densMapType))
            return False

        if self.deleteIntermediateFiles.lower() not in ('true', 'false'):
            self.writeError(
                text='"deleteIntermediateFiles" input of incompatible format' +
                     ' ("true","false"), case insensitive. Currently set as ' +
                     '"{}" in input file'.format(self.deleteIntermediateFiles))
            return False
        try:
            float(self.sfall_VDWR)
        except ValueError:
            self.writeError(
                text='"sfall_VDWR" input must be a float. Currently set as ' +
                     '"{}" in input file.'.format(self.sfall_VDWR))
            return False

        if float(self.sfall_VDWR) <= 0:
            self.writeError(
                text='"sfall_VDWR" input must be a positive float. Currently' +
                     ' set as "{}" in input file.'.format(self.sfall_VDWR))
            return False

        if ('preset' not in self.FFTmapWeight and
                self.FFTmapWeight not in ('recalculate', 'False')):
            self.writeError(
                text='"FFTmapWeight" input must take either "recalculate",' +
                     '"False" or "preset,x" where "x" is the FOM weight ' +
                     'label (of the form FOMx) of a preset FOM column ' +
                     'within the input .mtz file. Currently set as ' +
                     '"{}" in input file'.format(self.FFTmapWeight))
            return False

        if self.FFTmapWeight.startswith('preset'):
            if self.FFTmapWeight.replace('preset', '')[0] != ',':
                self.writeError(
                    text='If "FFTmapWeight" input is specified by "preset", ' +
                         'then this must be followed a comma and then the ' +
                         'FOM label name within the input .mtz file (i.e. ' +
                         '"preset,x" for column label "FOMx"). If column ' +
                         'label is simply  "FOM" then use "preset,". ' +
                         'Currently set as "{}" in input file'.format(
                            self.FFTmapWeight))
                return False

        doseStr = 'Doses can be calculated using RADDOSE-3D (visit ' +\
                  'www.raddo.se for more details). If no doses have been' +\
                  'calculated but you wish to run this program, please' +\
                  ' set dose inputs to NOTCALCULATED within input file.\n'

        if self.dose1 != 'NOTCALCULATED':
            doses = self.dose1.split(',')
            for dose in doses:
                err = 'Each "dose1" input must be a positive float. A dose ' +\
                      'is currently set as "{}" in input file.\n{}'.format(
                        dose, doseStr)
                try:
                    float(dose)
                except ValueError:
                    self.writeError(text=err)
                    return False
                if float(dose) < 0:
                    self.writeError(text=err)
                    return False

        if self.dose2 != 'NOTCALCULATED':
            if not self.multiDatasets:
                err = '"dose2" input must be a positive float. ' +\
                      'Currently set as "{}" in input file.\n{}'.format(
                        self.dose2, doseStr)
                try:
                    float(self.dose2)
                except ValueError:
                    self.writeError(text=err)
                    return False
                if float(self.dose2) < 0:
                    self.writeError(text=err)
                    return False

            else:
                for dose in self.dose2.split(','):
                    err = 'All doses in "dose2" input must be positive ' +\
                          'floats. Currently set as "{}"'.format(self.dose2) +\
                          ' in input file.\n{}'.format(doseStr)
                    try:
                        float(dose)
                    except ValueError:
                        self.writeError(text=err)
                        return False
                    if float(dose) < 0:
                        self.writeError(text=err)
                        return False

        self.logFile.writeToLog(
            str='All input file parameters appear to be of suitable format.')

        return True

    def checkNameLength(self,
                        maxLength=24, name='', property='name1'):

        # check the length of an input name within input file

        if len(name) > maxLength:
            self.writeError(
                text='Each "{}" input must be less than '.format(property) +
                     '{} characters long. '.format(maxLength) +
                     'Currently {} characters long..'.format(len(name)))
            return False

        return True

    def checkForSeparateSIGFlabel(self):

        # check whether there is a separate SIGFP column label
        # specified within the RIDL input file and factor if so

        try:
            self.mtzSIGFPlabel1
            sep = True
        except AttributeError:
            self.sepSIGFPlabel1 = False
            sep = False

        if sep:
            self.sepSIGFPlabel1 = True
            self.props1.append('mtzSIGFPlabel1')

        try:
            self.mtzSIGFPlabel2
            sep = True
        except AttributeError:
            self.sepSIGFPlabel2 = False
            sep = False

        if sep:
            self.sepSIGFPlabel2 = True
            self.props2.append('mtzSIGFPlabel2')

    def checkMtzLabelsExist(self):

        # for each mtz file input, check that the correct labels as
        # specified within the txt input file are successfully found

        # initial dataset mtz files checked here
        if self.multiDatasets:
            if self.repeatedFile1InputsUsed:
                mtzFiles = [self.mtz1]
                if self.sepSIGFPlabel1:
                    mtzLabels = [[self.mtzlabels1, self.mtzSIGFPlabel1]]
                else:
                    mtzLabels = [self.mtzlabels1]
            else:
                mtzFiles = self.mtz1.split(',')
                mtzLabels = self.mtzlabels1.split(',')
                self.sepSIGFPlabel1 = False
        else:
            mtzFiles = [self.mtz1]
            if self.sepSIGFPlabel1:
                mtzLabels = [[self.mtzlabels1, self.mtzSIGFPlabel1]]
            else:
                mtzLabels = [self.mtzlabels1]

        for (f, lab) in zip(mtzFiles, mtzLabels):
            foundLabels = self.getLabelsFromMtz(fileName=f)

            if self.includeSIGF:
                if self.sepSIGFPlabel1:
                    labels = [lab[0], lab[1]]
                else:
                    labels = [lab, 'SIG'+lab]
            else:
                labels = [lab]

            for lab in labels:
                if lab not in foundLabels:
                    self.mtzLabelNotFound(
                        mtzFile=f, label=lab, labelList=foundLabels)
                    return False

        # later dataset mtz files checked here
        if not self.multiDatasets:
            mtzFiles = [self.mtz2]
            if self.sepSIGFPlabel2:
                mtzLabels = [[self.mtzlabels2, self.mtzSIGFPlabel2]]
            else:
                mtzLabels = [self.mtzlabels2]
        else:
            mtzFiles = self.mtz2.split(',')
            if self.sepSIGFPlabel2:
                zipLabels = zip(self.mtzlabels2.split(','),
                                self.mtzSIGFPlabel2.split(','))
                mtzLabels = [[a, b] for a, b in zipLabels]
            else:
                mtzLabels = self.mtzlabels2.split(',')

        for (f, lab) in zip(mtzFiles, mtzLabels):
            foundLabels = self.getLabelsFromMtz(fileName=f)

            if self.includeSIGF:
                if self.sepSIGFPlabel2:
                    labels = [lab[0], lab[1]]
                else:
                    labels = [lab, 'SIG'+lab]
            else:
                labels = [lab]

            for lab in labels:
                if lab not in foundLabels:
                    self.mtzLabelNotFound(mtzFile=f, label=lab,
                                          labelList=foundLabels)
                    return False

        # phase information dataset mtz file checked here
        if self.multiDatasets:
            if not self.repeatedPhaseInputsUsed:
                mtzFiles = self.mtz3.split(',')
                phaseLabels = self.phaseLabel.split(',')
                FcalcLabels = self.FcalcLabel.split(',')
            else:
                mtzFiles = [self.mtz3]
                phaseLabels = [self.phaseLabel]
                FcalcLabels = [self.FcalcLabel]
        else:
            mtzFiles = [self.mtz3]
            phaseLabels = [self.phaseLabel]
            FcalcLabels = [self.FcalcLabel]

        for (f, phaseLab, FcalcLab) in zip(mtzFiles, phaseLabels, FcalcLabels):
            foundLabels = self.getLabelsFromMtz(fileName=f)
            for lab in (phaseLab, FcalcLab):
                if lab not in foundLabels:
                    self.mtzLabelNotFound(mtzFile=f, label=lab,
                                          labelList=foundLabels)
                    return False

        return True

    def mtzLabelNotFound(self,
                         mtzFile='untitled.mtz', label='unspecified',
                         labelList=[]):

        # if an mtz label has not been found, print an error message
        err = 'Column "{}" not found in file '.format(label) +\
              '"{}".\n'.format(mtzFile.split('/')[-1]) +\
              'Please check the .txt input file used for the current' +\
              ' job.\n Also check the mtz file to ensure the ' +\
              'column name "{}" exists. '.format(label)

        if labelList != []:
            err += '\nLabels found in mtz file:\n' +\
                   '{}'.format(', '.join(labelList))
            closeMatch = difflib.get_close_matches(label, labelList)

            if len(closeMatch) == 1:
                err += '\nDid you mean "{}"?'.format(closeMatch[0])
            elif len(closeMatch) > 1:
                err += '\nDid you mean any of ({})?'.format(
                    ', '.join(['"{}"'.format(i) for i in closeMatch]))

        self.writeError(text=err)

    def getLabelsFromMtz(self,
                         fileName='untitled.mtz'):

        # get list of column names from an mtz file

        # write commands for mtzdump to run
        inFile = open('mtzdumpInput.txt', 'w')
        inFile.write('header\nend')
        inFile.close()

        # run mtzdump and parse output to find column labels
        os.system('mtzdump hklin {} < mtzdumpInput.txt > mtzdump.tmp'.format(
            fileName))
        output = open('mtzdump.tmp', 'r')
        labelsFound = False
        labels = ''
        for l in output.readlines():
            if l.replace('\n', '') == '':
                continue
            if '* Column Types :' in l:
                break
            if labelsFound:
                labels = l.split()
            if '* Column Labels :' in l:
                labelsFound = True
                continue
        output.close()
        os.remove('mtzdump.tmp')

        return labels

    def checkFileExists(self,
                        fileName=''):

        # check that file exists and return bool

        if os.path.isfile(fileName) and os.access(fileName, os.R_OK):
            return True
        else:
            self.writeError(
                text='File "{}" could not be located..'.format(fileName))
            return False

    def checkCorrectFileExtension(self,
                                  fileName='', fileType='.pdb',
                                  property='pdb1'):

        # check that a file has the required file extension

        if not fileName.endswith(fileType):
            self.writeError(
                text='"{}" input file input must end '.format(property) +
                     'with extension  "{}". '.format(fileType) +
                     'Currently file is set as "{}".'.format(fileName))
            return False
        return True

    def checkForMultipleDatasets(self):

        # multiple processing jobs can be run sequentially if input
        # file contains comma separated lists for particular files.
        # Determine whether this is the case and check that correctly
        # formatted.

        props = self.props2
        if self.dose2 != 'NOTCALCULATED':
            props.append('dose2')

        lengths = []
        for p in props:
            val = getattr(self, p)
            lengths.append(len(val.split(',')))

        self.numDsets = lengths[0]
        if lengths == [1]*len(props):
            self.logFile.writeToLog(
                str='Only single dataset located and will be processed')
            self.multiDatasets = False

        elif lengths[1:] != lengths[:-1]:
            self.writeError(
                text='Error! Input file properties ' +
                     '({}) must'.format(','.join(props)) +
                     'have same number of comma-separated inputs')

        else:
            self.logFile.writeToLog(
                str='Multiple higher dose datasets located in ' +
                    ' input file to be processed.')
            self.multiDatasets = True

            # check whether 1 single initial dataset present, or separate
            # 'initial' dataset for each corresponding higher dose dataset
            self.checkForOptionalMultiInputs(props=self.props1, type='initial')

            # check whether 1 single phase dataset present, or separate
            # phases dataset for each corresponding higher dose dataset
            self.checkForOptionalMultiInputs(props=self.props3, type='phase')

    def checkForOptionalMultiInputs(self,
                                    props=[], type='initial'):

        # check for optional multiple inputs in input file
        # (for initital dataset and phase datasets)

        lengths = []
        for p in props:
            val = getattr(self, p)
            lengths.append(len(val.split(',')))

        if lengths not in ([1]*len(props), [self.numDsets]*len(props)):
            self.writeError(
                text='Error! Input file properties ' +
                     '({})must all '.format(','.join(props)) +
                     'each have either 1 or {} '.format(self.numDsets) +
                     'comma-separated inputs')
        else:
            s = False
            if lengths == [1]*len(props):
                s = True
            if type == 'initial':
                self.repeatedFile1InputsUsed = s
            elif type == 'phase':
                self.repeatedPhaseInputsUsed = s

    def checkOutputDirsExists(self,
                              makeProcessDir=True):

        # check whether output directory exists and make if not

        if not self.dir.endswith('/'):
            self.logFile.writeToLog(
                str='Working directory specified in input ' +
                    'file must end in "/" - appending.')
            self.dir += '/'

        if not os.path.exists(self.dir):
            self.logFile.writeToLog(
                str='Output directory "{}"'.format(self.dir) +
                    ' not found, making directory')

            os.makedirs(self.dir)
        self.logFile.writeToLog(
            str='Working directory set to "{}"'.format(self.dir))

        self.mapProcessDir = self.dir + 'RIDL-maps/'
        if makeProcessDir:
            if 'RIDL-maps' not in os.listdir(self.dir):
                os.makedirs(self.mapProcessDir)

    def getCurrentInputParams(self,
                              jobNumber=0):

        # select correct input parameters if multiple jobs within
        # main input file. used when writing separate input files
        # for each part of pipeline.

        if not self.multiDatasets:
            for prop in self.props1+self.props2+self.props3:
                setattr(self, prop+'_current', getattr(self, prop))
            self.createDatasetName()
            return

        for prop in self.props2:
            setattr(self, prop+'_current',
                    getattr(self, prop).split(',')[jobNumber])

        for prop in self.props1:
            if self.repeatedFile1InputsUsed:
                setattr(self, prop+'_current', getattr(self, prop))
            else:
                setattr(self, prop+'_current',
                        getattr(self, prop).split(',')[jobNumber])

        for prop in self.props3:
            if self.repeatedPhaseInputsUsed:
                setattr(self, prop+'_current', getattr(self, prop))
            else:
                setattr(self, prop+'_current',
                        getattr(self, prop).split(',')[jobNumber])

        self.createDatasetName()

    def createDatasetName(self):

        # create dataset name here, used to distinguish
        # input file naming scheme (if name2 not the same
        # as pdb2 then this is required)

        self.dsetName = (self.pdb2_current).split('/')[-1].replace('.pdb', '')

    def setJobName(self):

        # set a name for the current job

        self.jobName = '{}-{}'.format(self.name2_current, self.name1_current)

    def findFilesInDir(self,
                       mapProcessDir=True):

        # find files initially in working directory

        if mapProcessDir:
            self.filesInDir = os.listdir(self.mapProcessDir)
        else:
            self.filesInDir = os.listdir(self.dir)

    def cleanUpIntermediateFiles(self,
                                 removeMtzs=True, removeMaps=True,
                                 removePdbs=True, includeFCmap=True):

        # after successful completion of the map processing
        # part of the pipeline for each SINGLE dataset
        # clean up working directory

        self.logFile.writeToLog(str='\nCleaning up working directory.')

        # distinguish between FFT and END map output formats
        # depending on program used (FFT/END shellscript)
        if self.densMapType == 'END':
            densMapProg = 'END_switchedAxes'
        else:
            densMapProg = 'fft'

        m = [self.mapProcessDir]
        params = m + [self.dsetName]
        renameParams = m + [self.name2_current]
        # renameParams2 = m + [self.name1_current]

        keyLogFiles = [self.mtzToMapsPipelineLog]

        if self.densMapType != 'END':
            mapFiles = ['{}{}-{}-{}_cropped_cropped.map'.format(
                self.mapProcessDir, self.dsetName,
                self.densMapType, densMapProg)]
        else:
            mapFiles = ['{}{}-{}_cropped_cropped.map'.format(
                self.mapProcessDir, self.dsetName, densMapProg)]
        mapFiles += ['{}{}_sfall_cropped.map'.format(*params)]

        renameMaps = ['{}{}_density.map'.format(*renameParams),
                      '{}{}_atoms.map'.format(*renameParams)]

        # cannot keep FC maps if they were never made
        if self.calculateFCmaps.upper() == 'FALSE':
            includeFCmap = False

        if includeFCmap:
            mapFiles += ['{}{}-FC-{}_cropped_cropped.map'.format(
                self.mapProcessDir, self.dsetName, densMapProg)]
            renameMaps += ['{}{}_FC.map'.format(*renameParams)]

        pdbFiles = ['{}{}_reordered.pdb'.format(*params)]
        renamePDB = ['{}{}.pdb'.format(*renameParams)]

        outputFiles = mapFiles + pdbFiles
        renameFiles = renameMaps + renamePDB

        subdir = '{}{}_additionalFiles/'.format(
            self.mapProcessDir, self.jobName)
        self.makeOutputDir(dirName=subdir)

        keyFiles = keyLogFiles + mapFiles + pdbFiles
        for file in os.listdir(self.mapProcessDir):
            if file.endswith('_additionalFiles'):
                continue
            if file not in self.filesInDir:
                fileName = '{}{}'.format(self.mapProcessDir, file)
                if fileName not in keyFiles:
                    shutil.move(fileName, '{}{}'.format(subdir, file))

        for file in os.listdir(subdir):
            remove = False
            if removeMtzs and file.endswith('.mtz'):
                remove = True
            if removeMaps and file.endswith('.map'):
                remove = True
            if removePdbs and file.endswith('.pdb'):
                remove = True
            if remove:
                os.remove(subdir + file)

        shutil.make_archive(subdir, 'zip', subdir)
        shutil.rmtree(subdir)

        # rename final map & pdb files
        for i in range(len(outputFiles)):
            shutil.move(outputFiles[i], renameFiles[i])

        # move initial dataset pdb files to working directory
        self.moveInitialPDBfile()

        # check that resulting files are found
        self.findFilesInDir()
        for f in renameFiles:
            if f.split('/')[-1] not in self.filesInDir:
                self.writeError(text='Not all key output files found')
                return False
        return True

    def moveInitialPDBfile(self):

        # copy the CURRENT initial dataset pdb file to the working
        # directory useful if further RIDL jobs to be run).
        # ONLY do this if file not found in directory already. This
        # avoids PDBCUR-curated files being overwritten from other
        # section of pipeline. This probably needs simplifying in
        # the long run, but avoids difficulties when multiple initial
        # dataset inputs have been specified in the input file

        f = '{}.pdb'.format(self.name1_current)
        d = self.mapProcessDir
        if f not in os.listdir(d):
            shutil.copy2(self.pdb1_current, d+f)

    def makeOutputDir(self,
                      dirName='./'):

        # if the above sub directory does not exist, make it

        if not os.path.exists(dirName):
            os.makedirs(dirName)
            self.logFile.writeToLog(
                str='New sub directory "{}"'.format(dirName) +
                    'created to contain output files.')

    def writeError(self,
                   text='', type='error'):

        # write error message to log file

        error(text=text, log=self.logFile, type=type)
