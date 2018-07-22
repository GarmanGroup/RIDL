from makeMapsFromMTZs import makeMapsFromMTZs
from calculateMetrics import calculateMetrics
from cleanUpFiles import cleanUpFinalFiles
from ridlFeedback import provideFeedback
from furtherOutput import furtherAnalysis
from savevariables import retrieveGenericObject
from errors import error
import difflib
import shutil
import os
import sys


class processFiles():

    def __init__(self,
                 inputFile='', makeMaps=True, makeMetrics=True,
                 metCalcInput='metricCalc_inputfile.txt',
                 cleanFinalFiles=False, logFileObj='',
                 makeSummaryFile=False,
                 keepMapDir=True, includeSIGF=True):

        # class to read an input file and generate a set of density
        # and atom-tagged maps for a damage series. Can handle a single
        # low & high dose dataset pair, or a low dose dataset and a series
        # of high dose datasets. Can optionally proceed directly on to
        # the metric calculation stage of the pipeline with
        # 'proceedToMetricCalc' or skip directly to it (if the maps have
        # already been created) with 'skipToMetricCalc'.

        self.inputFile = inputFile
        self.metCalcInput = metCalcInput
        self.logFile = logFileObj
        self.includeSIGF = includeSIGF

        self.runFileProcessing()

        if makeMaps:
            self.runMapGeneration()

        if makeMetrics:
            self.runMetricCalcStep()
            if not makeSummaryFile:
                self.writeSummaryFiles(csvOnly=True)

        if makeSummaryFile:
            self.writeSummaryFiles(csvOnly=False)
            if cleanFinalFiles:
                cleanUpFinalFiles(outputDir=self.dir,
                                  keepMapDir=keepMapDir)

    def runFileProcessing(self):

        self.logFile.writeToLog(str='\n**** INPUT FILE PROCESSING ****\n')

        self.readMainInputFile()
        self.checkForRelFilePaths()
        self.checkNonNecessaryInputs()
        self.checkAllRequiredInputsFound()
        self.checkWhetherParamsConsistent()
        self.checkOutputDirsExists()
        self.findFilesInDir()
        self.checkForMultipleDatasets()
        self.defineMetricNormSet()

        # for safety, at this point check for internal breakdown
        try:
            self.multiDatasets
        except AttributeError:
            print('Unexpected error')
            return
        if self.multiDatasets:
            try:
                if not self.highDsetOnly:
                    self.repeatedFile1InputsUsed
                self.repeatedPhaseInputsUsed
            except AttributeError:
                print('Unexpected error')
                return

        self.checkCorrectInputFormats()
        self.extractInfoFromMtzs()
        self.checkMtzLabelsExist()
        self.compareMtzFiles()

    def runMapGeneration(self):

        self.logFile.writeToLog(str='\n\n**** MAP GENERATION ****\n')

        if not self.multiDatasets:
            self.logFile.writeToLog(
                str='Generating suitable electron density maps for run.')
            self.getCurrentInputParams()
            self.runMapGenerationPipeline()

        else:
            self.logFile.writeToLog(
                str='Generating suitable electron density maps for each ' +
                    'higher dose dataset in turn.\n{} '.format(self.numDsets) +
                    'higher dose datasets located for this job.')

            for i in range(self.numDsets):
                if i == 0:
                    firstTimeRun = True
                else:
                    firstTimeRun = False
                self.logFile.writeToLog(
                    str='\n{}\nHigher dose '.format('-'*33) +
                        'dataset {} starts here'.format(i+1))
                self.getCurrentInputParams(jobNumber=i)
                self.runMapGenerationPipeline(firstTimeRun)


    def runMapGenerationPipeline(self,
                                 firstTimeRun=True):

        # run the map generation pipeline for a single dataset

        if not self.highDsetOnly:
            if self.sepSIGFPlabel1:
                sigFP1 = self.mtzSIGFPlabel1_current
            else:
                sigFP1 = 'SIG' + self.mtzlabels1_current

        if self.sepSIGFPlabel2:
            sigFP2 = self.mtzSIGFPlabel2_current
        else:
            sigFP2 = 'SIG' + self.mtzlabels2_current

        pam, sg, am, ma, gs = '', '', '', [], []
        fcMaps = self.includeFCmaps()
        ignoreSIGFs = self.whetherIgnoreSIGFs()
        if self.useSeparatePDBperDataset():
            pdb = self.pdb2_current
        else:
            # note that if self.highDsetOnly = True and this is
            # trigged, the program will crash, so ensure that
            # useLaterCellDims is set to "true" in input file
            pdb = self.pdb1_current
            am = self.name1_current
            if not firstTimeRun:
                pam = self.copySameRunInfo['atomMapName']
                sg = self.copySameRunInfo['spaceGroup']
                ma = self.copySameRunInfo['axes']
                gs = self.copySameRunInfo['gridsamp']
                fcMaps = False

        p = makeMapsFromMTZs(
            outputDir=self.mapProcessDir, densMapNaming=self.name2_current,
            atomMapNaming=am, log=self.logFile, mtzIn2=self.mtz2_current,
            Mtz2LabelName=self.mtzlabels2_current, Mtz2SIGFPlabel=sigFP2,
            Mtz2LabelRename=self.name2_current, mtzIn3=self.mtz3_current,
            Mtz3phaseLabel=self.phaseLabel_current,
            Mtz3FcalcLabel=self.FcalcLabel_current,
            Mtz3LabelRename=self.name3_current, inputPDBfile=pdb,
            densMapType=self.densMapType, scaleType=self.scaleType,
            FOMweight=self.FFTmapWeight, firstTimeRun=firstTimeRun,
            sfall_VDWR=self.sfall_VDWR, mapResLimits=self.mapResLimits,
            includeFCmaps=fcMaps, ignoreSIGFs=ignoreSIGFs,
            useLaterCellDims=self.useSeparatePDBperDataset(),
            deleteIntermediateFiles=self.deleteUnwantedFiles(),
            premadeAtomMap=pam, gridSampBeforeCropping=gs,
            mapAxisOrder=ma, spaceGroup=sg)

        if not self.highDsetOnly:
            p.mtzIn1 = self.mtz1_current
            p.Mtz1LabelName = self.mtzlabels1_current
            p.Mtz1SIGFPlabel = sigFP1
            p.Mtz1LabelRename = self.name1_current

            # only require an Rfree flag if SIGMAA will be used
            if self.FFTmapWeight.lower() == 'recalculate':
                p.RfreeFlag1 = self.RfreeFlag1_current

        success = p.runPipeline()
        self.mtzToMapsPipelineLog = p.runLog.logFile
        if success:
            self.logFile.writeToLog(str='---> Subroutine ran to completion.')

            # move initial dataset pdb files to working directory
            self.moveInitialPDBfile()

        else:
            self.writeError(text='Subroutine failed to run to completion')

        if firstTimeRun and not self.useSeparatePDBperDataset():
            # in the case where the same PDB file is used for each dataset
            # the atom-map is only generated once. Use this information
            # to properly define and density maps generated for each dataset
            self.copySameRunInfo = {'atomMapName': p.atomTaggedMap,
                                    'axes': p.axes,
                                    'gridsamp': p.gridSamps,
                                    'spaceGroup': p.spaceGroup}

    def runMetricCalcStep(self):

        # writes an input file for the run of the metric calculation part of
        # RIDL pipeline after the map generation pipeline has completed. Allows
        # the user to immediately run the rest of the program, IF the user did
        # specify all required datasets within a damage series. If the user
        # only processed a subset of required datasets within a damage series,
        # another metric calculation input file will need to be manually
        # created to run the rest of the program with ALL datasets within the
        # damage series included.

        seriesName = self.dir.split('/')[-2]

        if self.multiDatasets:
            if not self.repeatedFile1InputsUsed:
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

        self.logFile.writeToLog(str='\n\n**** METRIC CALCULATIONS ****\n')

        names2 = self.name2.split(',')

        if self.multiDatasets:
            if self.repeatedFile1InputsUsed:
                names1 = [self.name1]*len(names2)
            else:
                names1 = self.name1.split(',')
        else:
            names1 = [self.name1]

        densMapList = ['{}_density.map'.format(d) for d in names2]

        if self.useSeparatePDBperDataset():
            atomMapList = ['{}_atoms.map'.format(d) for d in names2]
            pdbFileList = ['{}.pdb'.format(d) for d in names2]
        else:
            atomMapList = ['{}_atoms.map'.format(d) for d in names1]
            pdbFileList = ['{}.pdb'.format(d) for d in names1]

        c = calculateMetrics(logFile=self.logFile, densMapList=densMapList,
                             atomMapList=atomMapList, outDir=self.dir,
                             mapDir=self.mapProcessDir, initialPDB=self.name1,
                             seriesName=seriesName, doses=self.getDoses(),
                             pklDataFile=self.pklDataFile,
                             pdbFileList=pdbFileList, normSet=self.normSet,
                             RIDLinputFile=self.inputFile,
                             sepPDBperDataset=self.useSeparatePDBperDataset())

        # decide whether Fcalc data is present and should be used
        c.inclFCmets = self.includeFCmaps()
        if self.includeFCmaps():
            if self.useSeparatePDBperDataset():
                c.FcMapList = ['{}_FC.map'.format(d) for d in names2]
            else:
                c.FcMapList = ['{}_FC.map'.format(d) for d in names1]

        c.runPipeline()

        self.pklDataFile = c.pklDataFile

    def writeSummaryFiles(self,
                          csvOnly=False, includeTests=False):

        # write feedback files for current RIDL job. if csvOnly is True then
        # ONLY csv files will be output from the run (i.e. no html summary
        # file and no plots)

        self.logFile.writeToLog(str='\n**** RIDL FEEDBACK ****\n')

        # retrieve list of atom objects from .pkl file
        self.logFile.writeToLog(
            str='Metric data retrieved from pkl file:\n' +
                '\t{}'.format(self.pklDataFile))

        # retrieve the combinedAtoms object from the pkl file
        combinedAtoms = retrieveGenericObject(fileName=self.pklDataFile)

        if not includeTests:
            outputDir = self.dir + 'RIDL-metrics/'
            provideFeedback(csvOnly=csvOnly, atmsObjs=combinedAtoms,
                            logFile=self.logFile, outputDir=outputDir,
                            doses=self.getDoses(), pklSeries=self.pklDataFile,
                            inputDir=self.mapProcessDir,
                            densMaps=self.name2.split(','),
                            initialPDB=self.name1, normSet=self.normSet,
                            inclFCmetrics=self.includeFCmaps())
        else:
            furtherAnalysis(csvOnly=csvOnly, atmsObjs=combinedAtoms,
                            logFile=self.logFile, outputDir=outputDir,
                            doses=self.getDoses(), normSet=self.normSet,
                            pklSeries=self.pklDataFile,
                            inputDir=self.mapProcessDir,
                            pdbNames=self.name2, initialPDB=self.name1,
                            inclFCmetrics=self.includeFCmaps())

    def readMainInputFile(self):

        # split the input file into separate input
        # files for each section of the pipeline

        # if Input.txt not found, flag error
        if not os.path.isfile(self.inputFile):
            self.writeError(
                text='Input file "{}"" not found..'.format(self.inputFile))
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

    def checkForRelFilePaths(self):

        #  convert any relative file paths (with ./) to
        # absolute file paths

        file_params = ['mtz1', 'pdb1', 'pdb1', 'mtz2', 'mtz3']
        for file_param in file_params:
            current_val = getattr(self, file_param)

            for c in current_val.split(','):
                if not os.path.isabs(c):
                    self.writeError(
                        text='Relative file paths have been found in input ' +
                        'file. Please convert to absolute paths to continue')

    def checkAllRequiredInputsFound(self):

        # check that all required properties have been found

        self.props1 = ['name1', 'pdb1', 'mtz1', 'mtzlabels1']
        self.props2 = ['name2', 'mtz2', 'mtzlabels2']
        self.props3 = ['name3', 'mtz3', 'phaseLabel', 'FcalcLabel']

        if self.useSeparatePDBperDataset():
            self.props2.append('pdb2')

        # only need Rfree set if SIGMAA will be used
        if self.FFTmapWeight.lower() == 'recalculate':
            self.props1 += ['RfreeFlag1']

        requiredProps = self.props2 + self.props3 + ['dir', 'dose2']

        # for the case where a map is generated ONLY using higher dataset
        # information, do not strictly require any INITIALDATASET parameters
        if self.densMapType != 'HIGHONLY':
            requiredProps += self.props1 + ['dose1']

        for prop in requiredProps:
            try:
                getattr(self, prop)
            except AttributeError:
                self.writeError(
                    text='Necessary input not found: {}'.format(prop))

        self.logFile.writeToLog(
            str='--> All necessary inputs found in input file.')

    def checkNonNecessaryInputs(self):

        # check whether non necessary inputs are present
        # and set to default values if not specified

        props = ['sfall_VDWR', 'mapResLimits', 'scaleType',
                 'densMapType', 'FFTmapWeight', 'calculateFCmaps',
                 'deleteIntermediateFiles', 'useLaterCellDims',
                 'pklDataFile', 'normSet', 'ignoreSIGFs']

        defaults = [1, ',', 'anisotropic', 'DIFF',
                    'false', 'true', 'true', 'false', '', 'CALPHA', 'false']

        for i, prop in enumerate(props):
            try:
                getattr(self, prop)
            except AttributeError:
                setattr(self, prop, defaults[i])
                self.logFile.writeToLog(
                    str='--> Input "{}" not found in input file'.format(prop) +
                        '. Setting to default value: "{}"'.format(defaults[i]))

        if self.densMapType == 'HIGHONLY':
            self.highDsetOnly = True
        else:
            self.highDsetOnly = False

    def checkWhetherParamsConsistent(self):

        # check that the specified parameters are compatible with each other

        if self.densMapType == 'HIGHONLY':
            if self.scaleType != 'NONE':
                self.writeError(
                    text='"scaleType" input parameter must take value "NONE"' +
                         ' when "densMapType" is set to "HIGHONLY"')
            if self.useLaterCellDims.lower() != 'true':
                self.writeError(
                    text='"useLaterCellDims" input parameter must take value' +
                         '"false" when "densMapType" is set to "HIGHONLY"')

            for prop in self.props1:
                foundProp = False
                try:
                    getattr(self, prop)
                    foundProp = True
                except AttributeError:
                    pass
                if foundProp:
                    self.writeError(
                        text='"{}" input parameter has been '.format(prop) +
                             'specified, yet since "densMapType" set to ' +
                             '"HIGHONLY", it will not be used.. check this ' +
                             'is intentional!',
                        type='warning')

    def checkCorrectInputFormats(self):

        # check that the required properties have correct formats

        props = ['mtz2', 'mtz3']
        fType = ['.mtz', '.mtz']

        if not self.highDsetOnly:
            props += ['mtz1', 'pdb1']
            fType += ['.mtz', '.pdb']

        if self.useSeparatePDBperDataset():
            props.append('pdb2')
            fType.append('.pdb')

        for p, t in zip(props, fType):
            if not self.multiDatasets:
                name = getattr(self, p)
                self.checkCorrectFileExtension(
                    fileName=name, fileType=t, property=p)
                self.checkFileExists(fileName=name)

            else:
                for f in getattr(self, p).split(','):
                    self.checkCorrectFileExtension(
                        fileName=f, fileType=t, property=p)
                    self.checkFileExists(fileName=f)

        if not self.highDsetOnly:
            # ensure that initial and later dataset names don't clash
            if not self.multiDatasets:
                if self.name1 == self.name2:
                    self.writeError(
                        text='"name1" and "name2" inputs must be different, ' +
                             'otherwise CAD will fail. Both currently set ' +
                             'as "{}".'.format(self.name1))
            else:
                isError = False
                if self.repeatedFile1InputsUsed:
                    for n in self.name2.split(','):
                        if self.name1 == n:
                            isError = True
                            sameName = n
                            break

                else:
                    nameZip = zip(self.name1.split(','), self.name2.split(','))
                    for n1, n2 in nameZip:
                        if n1 == n2:
                            isError = True
                            sameName = n1
                            break

                if isError:
                    self.writeError(
                        text='"name1" and "name2" inputs must be different ' +
                             'for each job in batch, otherwise CAD will fail. ' +
                             'Currently for one batch, the "initial" and "later"' +
                             ' datasets are both called "{}".'.format(sameName))

        # ensure that initial dataset names of suitable lengths
        if not self.highDsetOnly:
            for n in self.name1.split(','):
                self.checkNameLength(name=n, property='name1')

        # ensure that later dataset names of suitable lengths
        for n in self.name2.split(','):
            self.checkNameLength(name=n, property='name2')

        # ensure that multiple later dataset names don't clash
        if self.multiDatasets:
            names = self.name2.split(',')
            if len(list(set(names))) != len(names):
                self.writeError(
                    text='Comma-separated entries in "name2" must be unique')

        # ensure that phase dataset names of suitable lengths
        for n in self.name3.split(','):
            self.checkNameLength(
                name=n, property='name3', maxLength=22)

        # add '-ph' to name3 to identify it from other datasets
        if not self.multiDatasets:
            self.name3 += '-ph'
        else:
            if self.repeatedPhaseInputsUsed:
                self.name3 += '-ph'
            else:
                self.name3 = ','.join([n+'-ph' for n in self.name3.split(',')])

        doseStr = 'Doses can be calculated using RADDOSE-3D (visit ' +\
                  'www.raddo.se for more details). If no doses have been' +\
                  'calculated but you wish to run this program, please' +\
                  ' set dose inputs to NOTCALCULATED within input file.\n'

        if not self.highDsetOnly:
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
                    if float(dose) < 0:
                        self.writeError(text=err)

        if self.dose2 != 'NOTCALCULATED':
            if not self.multiDatasets:
                err = '"dose2" input must be a positive float. ' +\
                      'Currently set as "{}" in input file.\n{}'.format(
                        self.dose2, doseStr)
                try:
                    float(self.dose2)
                except ValueError:
                    self.writeError(text=err)
                if float(self.dose2) < 0:
                    self.writeError(text=err)

            else:
                for dose in self.dose2.split(','):
                    err = 'All doses in "dose2" input must be positive ' +\
                          'floats. Currently set as "{}"'.format(self.dose2) +\
                          ' in input file.\n{}'.format(doseStr)
                    try:
                        float(dose)
                    except ValueError:
                        self.writeError(text=err)
                    if float(dose) < 0:
                        self.writeError(text=err)

        # if a separate input exists for the SIGFP columns, check consistency
        self.checkForSeparateSIGFlabel()

        if not self.highDsetOnly:
            if self.sepSIGFPlabel1:
                l1 = len(self.mtzlabels1.split(','))
                l2 = len(self.mtzSIGFPlabel1.split(','))
                if l1 != l2:
                    self.writeError(
                        text='"mtzSIGFPlabel1" and "mtzlabels1" do not ' +
                             'contain same number of entries')

        if self.sepSIGFPlabel2:
            l1 = len(self.mtzlabels2.split(','))
            l2 = len(self.mtzSIGFPlabel2.split(','))
            if l1 != l2:
                self.writeError(
                    text='"mtzSIGFPlabel2" and "mtzlabels2" do not contain ' +
                         'same number of entries')

        # also check the optimal inputs that may be at the bottom of the file

        if self.densMapType not in ('DIFF', 'SIMPLE', '2FOFC', 'END', 'HIGHONLY'):
            self.writeError(
                text='"densMapType" input of incompatible format, ' +
                     '(default is "DIFF"). Currently set as ' +
                     '"{}" in input file.'.format(self.densMapType))

        if self.deleteIntermediateFiles.lower() not in ('true', 'false'):
            self.writeError(
                text='"deleteIntermediateFiles" input of incompatible format' +
                     ' ("true","false"), case insensitive. Currently set as ' +
                     '"{}" in input file'.format(self.deleteIntermediateFiles))

        try:
            float(self.sfall_VDWR)
        except ValueError:
            self.writeError(
                text='"sfall_VDWR" input must be a float. Currently set as ' +
                     '"{}" in input file.'.format(self.sfall_VDWR))

        if float(self.sfall_VDWR) <= 0:
            self.writeError(
                text='"sfall_VDWR" input must be a positive float. Currently' +
                     ' set as "{}" in input file.'.format(self.sfall_VDWR))

        if float(self.sfall_VDWR) > 3:
            self.writeError(
                text='"sfall_VDWR" input can only take values less than 3. ' +
                     'Larger values will be set to the upper limit of 3. ' +
                     'Currently "{}" in input file.'.format(self.sfall_VDWR),
                type='warning')

        if float(self.sfall_VDWR) < 1:
            self.writeError(
                text='A "sfall_VDWR" input of 1 or greater . ' +
                     'is recommended - consider increasing this..',
                type='warning')

        if ('preset' not in self.FFTmapWeight.lower() and
                self.FFTmapWeight.lower() not in ('recalculate', 'false')):
            self.writeError(
                text='"FFTmapWeight" input must take either "recalculate",' +
                     '"false" or "preset,x" where "x" is the FOM weight ' +
                     'label (of the form FOMx) of a preset FOM column ' +
                     'within the input .mtz file. Currently set as ' +
                     '"{}" in input file'.format(self.FFTmapWeight))

        if self.FFTmapWeight.lower().startswith('preset'):
            if self.FFTmapWeight.replace('preset', '')[0] != ',':
                self.writeError(
                    text='If "FFTmapWeight" input is specified by "preset", ' +
                         'then this must be followed a comma and then the ' +
                         'FOM label name within the input .mtz file (i.e. ' +
                         '"preset,x" for column label "FOMx"). If column ' +
                         'label is simply  "FOM" then use "preset,". ' +
                         'Currently set as "{}" in input file'.format(
                            self.FFTmapWeight))

        if self.useLaterCellDims.lower() not in ('true', 'false'):
            self.writeError(
                text='"useLaterCellDims" input of incompatible format' +
                     ' ("true","false"), case insensitive. Currently set as ' +
                     '"{}" in input file'.format(self.useLaterCellDims))

        if self.calculateFCmaps.lower() not in ('true', 'false'):
            self.writeError(
                text='"calculateFCmaps" input of incompatible format' +
                     ' ("true","false"), case insensitive. Currently set as ' +
                     '"{}" in input file'.format(self.calculateFCmaps))

        if self.ignoreSIGFs.lower() not in ('true', 'false'):
            self.writeError(
                text='"ignoreSIGFs" input of incompatible format' +
                     ' ("true","false"), case insensitive. Currently set as ' +
                     '"{}" in input file'.format(self.ignoreSIGFs))

        if self.scaleType.lower() not in ('anisotropic', 'isotropic', 'scale', 'none', 'phenix'):
            self.writeError(
                text='"scaleType" input of incompatible format, ' +
                     '(default is "anisotropic"). Currently set as ' +
                     '"{}" in input file.'.format(self.scaleType))

        if self.mapResLimits != ',':
            resLims = self.mapResLimits.split(',')
            if len(resLims) != 2:
                self.writeError(
                    text='"mapResLimits" input of incompatible format, ' +
                         ' must be of form "x,y" for upper resolution ' +
                         'cutoff "x" and lower cutoff "y". Currently set as ' +
                         '"{}" in input file.'.format(self.mapResLimits))
            else:
                try:
                    float(resLims[0])
                    float(resLims[1])
                except ValueError:
                    self.writeError(
                        text='"mapResLimits" input of incompatible format, ' +
                             ' must be of form "x,y" for floats "x" and "y" ' +
                             'for upper resolution cutoff "x" and lower ' +
                             'cutoff "y". Currently set as ' +
                             '"{}" in input file.'.format(self.mapResLimits))
                if float(resLims[0]) == float(resLims[1]):
                    self.writeError(
                        text='"mapResLimits" input of incompatible format, ' +
                             ' must be of form "x,y" for "x" and "y" ' +
                             'for upper resolution cutoff "x" and lower ' +
                             'cutoff "y", BUT DIFFERENT. Currently set as ' +
                             '"{}" in input file.'.format(self.mapResLimits))

        self.logFile.writeToLog(
            str='--> All input file parameters of suitable format.')

    def checkNameLength(self,
                        maxLength=24, name='', property='name1'):

        # check the length of an input name within input file

        if len(name) > maxLength:
            self.writeError(
                text='Each "{}" input must be less than '.format(property) +
                     '{} characters long. '.format(maxLength) +
                     'Currently {} characters long..'.format(len(name)))

    def checkForSeparateSIGFlabel(self):

        # check whether there is a separate SIGFP column label
        # specified within the RIDL input file and factor if so

        if not self.highDsetOnly:
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
        if not self.highDsetOnly:
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

                if self.includeSIGF:
                    if self.sepSIGFPlabel1:
                        labels = [lab[0], lab[1]]
                    else:
                        labels = [lab, 'SIG'+lab]
                else:
                    labels = [lab]

                for lab in labels:
                    if lab not in self.mtzContents[f]['labels']:
                        self.mtzLabelNotFound(mtzFile=f, label=lab)

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

            if self.includeSIGF:
                if self.sepSIGFPlabel2:
                    labels = [lab[0], lab[1]]
                else:
                    labels = [lab, 'SIG'+lab]
            else:
                labels = [lab]

            for lab in labels:
                if lab not in self.mtzContents[f]['labels']:
                    self.mtzLabelNotFound(mtzFile=f, label=lab)

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
            for lab in (phaseLab, FcalcLab):
                if lab not in self.mtzContents[f]['labels']:
                    self.mtzLabelNotFound(mtzFile=f, label=lab)

    def mtzLabelNotFound(self,
                         mtzFile='untitled.mtz', label='unspecified'):

        labelList = self.mtzContents[mtzFile]['labels']

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

    def extractInfoFromMtzs(self):

        # get list of column names from an mtz file

        files = self.mtz1.split(',') + self.mtz2.split(',') + self.mtz3.split(',')

        self.mtzContents = {}

        for f in files:
            # write commands for mtzdump to run
            inFile = open('mtzdumpInput.txt', 'w')
            inFile.write('stats nbin 1\nend')
            inFile.close()

            # run mtzdump and parse output to find column labels
            os.system('mtzdump hklin {} < mtzdumpInput.txt > mtzdump.tmp'.format(f))
            output = open('mtzdump.tmp', 'r')
            labelsFound, overallStatsFound = False, False
            labels, overallStats = '', ''
            for l in output.readlines():
                if l.replace('\n', '') == '':
                    continue

                if '* Column Types :' in l:
                    labelsFound = False
                if labelsFound:
                    labels = l.split()
                if '* Column Labels :' in l:
                    labelsFound = True
                    continue

                if 'OVERALL FILE STATISTICS' in l:
                    overallStatsFound = True
                if 'LIST OF REFLECTIONS' in l:
                    break
                if overallStatsFound:
                    overallStats += l

            output.close()
            os.remove('mtzdump.tmp')

            self.mtzContents[f] = {'labels': labels, 'stats': overallStats}

    def compareMtzFiles(self):

        # compare the initial dataset with later datasets to check that
        # all mtz files are unique. If a low and high dose mtz are identical
        # this will cause FFT to crash later when calculating a Fo-Fo map

        # currently the check is only made for difference maps
        if self.densMapType != 'DIFF':
            return

        for f in self.mtz1.split(','):
            for f2 in self.mtz2.split(','):
                if self.mtzContents[f]['stats'] == self.mtzContents[f2]['stats']:
                    self.writeError(
                        text='Initial mtz "{}" and later mtz "{}" '.format(f, f2) +
                        'appear to contain identical sets of reflections. ' +
                        'This not allowed since it will cause FFT to fail ' +
                        'when calculating a difference map. '
                        'This decision was made by inspecting the overall ' +
                        'statistics output of MTZDUMP.')

    def checkFileExists(self,
                        fileName=''):

        # check that file exists and write a fatal error if not

        if not os.path.isfile(fileName) or not os.access(fileName, os.R_OK):
            self.writeError(
                text='File "{}" could not be located..'.format(fileName))

    def checkCorrectFileExtension(self,
                                  fileName='', fileType='.pdb',
                                  property='pdb1'):

        # check that a file has the required file extension

        if not fileName.endswith(fileType):
            self.writeError(
                text='"{}" input file input must end '.format(property) +
                     'with extension  "{}". '.format(fileType) +
                     'Currently file is set as "{}".'.format(fileName))

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
                str='--> Only single dataset located and will be processed')
            self.multiDatasets = False

        elif lengths[1:] != lengths[:-1]:
            self.writeError(
                text='Error! Input file properties ' +
                     '({}) must'.format(', '.join(props)) +
                     'have same number of comma-separated inputs')

        else:
            self.logFile.writeToLog(
                str='--> Multiple higher dose datasets located in ' +
                    ' input file to be processed.')
            self.multiDatasets = True

            # check whether 1 single initial dataset present, or separate
            # 'initial' dataset for each corresponding higher dose dataset
            if not self.highDsetOnly:
                self.checkForOptionalMultiInputs(props=self.props1,
                                                 type='initial')

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
                     '({}) must all '.format(', '.join(props)) +
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
                str='--> Working directory specified in input ' +
                    'file must end in "/" - appending.')
            self.dir += '/'

        if not os.path.exists(self.dir):
            self.logFile.writeToLog(
                str='--> Output directory "{}"'.format(self.dir) +
                    ' not found, making directory')

            os.makedirs(self.dir)
        self.logFile.writeToLog(
            str='--> Working directory set to "{}"'.format(self.dir))

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
            for prop in self.props2+self.props3:
                setattr(self, prop+'_current', getattr(self, prop))
            if not self.highDsetOnly:
                for prop in self.props1:
                    setattr(self, prop+'_current', getattr(self, prop))
            self.createDatasetName()
            return

        for prop in self.props2:
            setattr(self, prop+'_current',
                    getattr(self, prop).split(',')[jobNumber])

        if not self.highDsetOnly:
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

    def includeFCmaps(self):

        # interpret from input file whether to calculate FCALC maps

        if self.calculateFCmaps.lower() == 'false':
            return False
        else:
            return True

    def whetherIgnoreSIGFs(self):

        # interpret from input file whether to use SIGF coefficients

        if self.ignoreSIGFs.lower() == 'false':
            return False
        else:
            return True

    def defineMetricNormSet(self):

        # define how metrics should be normalised from input file

        n = self.normSet.upper()
        if n == 'CALPHA':
            self.normSet = [['', 'CA']]
        elif n == 'NONE':
            self.normSet = [[]]
        elif n.split(',')[0] == 'CUSTOM':
            s = []
            for n2 in n.split(',')[1:]:
                s.append(n2.split('='))
            self.normSet = s
        else:
            self.normSet = [[]]

    def getDoses(self):

        # interpret set of doses from the input file

        if self.dose2 == 'NOTCALCULATED':
            doses = ','.join(map(str, range(len(self.name2.split(',')))))
        else:
            doses = self.dose2

        return doses

    def deleteUnwantedFiles(self):

        # interpret from input file whether to delete intermediate files

        if self.deleteIntermediateFiles.lower() == 'true':
            return True
        else:
            return False

    def useSeparatePDBperDataset(self):

        # decide whether a separate pdb file should be used per dataset
        # (using the PDB2 line). This is currently controlled by the
        # 'useLaterCellDims' input, since maps must be generated over
        # later unit cell dimensions if the pdb file changes with dataset
        # number

        if self.useLaterCellDims.lower() == 'true':
            return True
        else:
            return False

    def createDatasetName(self,
                          fmt=2):

        # create dataset name here, used to distinguish
        # input file naming scheme (if name2 not the same
        # as pdb2 then this is required)

        if fmt == 1:
            # this will not work if self.highDsetOnly = True
            d = '{}-{}'.format(self.name2_current, self.name1_current)
        else:
            # this is the recommended default
            d = self.name2_current

        self.dsetName = d

    def setJobName(self,
                   fmt=2):

        # set a name for the current job

        if fmt == 1:
            # this will not work if self.highDsetOnly = True
            j = '{}-{}'.format(self.name2_current, self.name1_current)
        else:
            # this is the recommended default
            j = self.name2_current

        self.jobName = j

    def findFilesInDir(self,
                       mapProcessDir=True):

        # find files initially in working directory

        if mapProcessDir:
            self.filesInDir = os.listdir(self.mapProcessDir)
        else:
            self.filesInDir = os.listdir(self.dir)

    def moveInitialPDBfile(self):

        # copy the CURRENT initial dataset pdb file to the working
        # directory useful if further RIDL jobs to be run).
        # ONLY do this if file not found in directory already. This
        # avoids PDBCUR-curated files being overwritten from other
        # section of pipeline. This probably needs simplifying in
        # the long run, but avoids difficulties when multiple initial
        # dataset inputs have been specified in the input file

        # skip if no initial dataset information provided
        if self.highDsetOnly:
            return

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
