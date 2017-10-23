from MAPMASKjob import MAPMASKjob
from PDBCURjob import PDBCURjob
from SFALLjob import SFALLjob
from mapTools import mapTools
from logFile import logFile
from FFTjob import FFTjob
from ENDjob import ENDjob
from errors import error
import shutil
import os


class pipeline():

    # class to run SFALL job to generate atom-tagged map, and FFT to
    # generate a corresponding density map (typically a Fobs(n)-Fobs(1)
    # Fourier difference map). The maps are then both restricted to the
    # crystal asymmetric unit, with identical grid sampling properties
    # (i.e. they are 'compatible'). This subroutine runs after the
    # CAD-SCALEIT subroutine (see processFiles.py for how these are called).

    def __init__(self,
                 outputDir='', jobName='untitled-job', log='',
                 inputPDBfile='', Mtz1LabelName='', Mtz2LabelName='',
                 phaseDataset='', sfall_VDWR=1, scaleType='ANISOTROPIC',
                 mapResLimits=',', inputMtzFile='', densMapType='DIFF',
                 FOMweight='NONE', includeFCmaps=True,
                 useLaterCellDims=True, sfallGRIDdims=[]):

        self.outputDir = outputDir
        self.jobName = jobName
        self.inputPDBfile = inputPDBfile
        self.Mtz1LabelName = Mtz1LabelName
        self.Mtz2LabelName = Mtz2LabelName
        self.phaseDataset = phaseDataset
        self.sfall_VDWR = sfall_VDWR
        self.scaleType = scaleType
        self.mapResLimits = mapResLimits
        self.inputMtzFile = inputMtzFile
        self.densMapType = densMapType
        self.FOMweight = FOMweight
        self.includeFCmaps = includeFCmaps
        self.useLaterCellDims = useLaterCellDims
        self.sfallGRIDdims = sfallGRIDdims
        self.findFilesInDir()

        # create log file
        if log == '':
            f = '{}{}_runLog1.log'.format(self.outputDir+'RIDL-log', jobName)
            self.runLog = logFile(fileName=f, fileDir=self.outputDir)
        else:
            self.runLog = log

    def runPipeline(self):

        # run the current subroutine within this class

        if self.scaleType in ('NONE', 'PHENIX'):
            self.inputMtzFile = self.inputMtzFile.replace('_SCALEIT', '_CAD')

        success = self.curatePdbFile()
        if not success:
            return False

        self.renumberPDBFile()

        success = self.getSpaceGroup()
        if not success:
            return False

        success = self.getAtomTaggedMap()
        if not success:
            return False

        makeFoFoMapWithPhenix = False
        if self.scaleType == 'PHENIX':
            makeFoFoMapWithPhenix = True
        if makeFoFoMapWithPhenix:
            self.generatePhenixDensMap()
        else:
            self.mapTag = ''
            self.fcalcMtz = self.inputMtzFile
            self.densMapMtz = self.inputMtzFile
            success = self.generateCCP4DensMap()

        if not success:
            return False

        if self.includeFCmaps:
            success = self.generateFcalcMap()
            if not success:
                return False

        success = self.cropAtmTaggedMapToAsymUnit()
        if not success:
            return False

        if self.densMapType == 'END':
            success = self.ensureSameMapAxesOrder()
            if not success:
                return False

        success = self.cropMapToAtomTaggedMap(densMap=self.densityMap)
        if not success:
            return False

        if self.includeFCmaps:
            success = self.cropMapToAtomTaggedMap(densMap=self.FcalcMap)
            if not success:
                return False

        self.reportCroppedMapInfo()
        success = self.mapConsistencyCheck()
        if not success:
            return False

        self.cleanUpDir()
        return True

    def curatePdbFile(self):

        # run pdbcur job

        self.printStepNumber()
        pdbcur = PDBCURjob(inputPDBfile=self.inputPDBfile,
                           outputDir=self.outputDir, runLog=self.runLog)
        success = pdbcur.run()
        self.PDBCURoutputFile = pdbcur.outputPDBfile
        return success

    def getAtomTaggedMap(self):

        # run SFALL job to generate atom tagged map

        self.printStepNumber()
        sfall = SFALLjob(inputPDBfile=self.reorderedPDBFile,
                         outputDir=self.outputDir, VDWR=self.sfall_VDWR,
                         symmetrygroup=self.spaceGroup,
                         gridDimensions=self.sfallGRIDdims, runLog=self.runLog)
        success = sfall.run()

        sfallMap = mapTools(mapName=sfall.outputMapFile)
        self.axes = [sfallMap.fastaxis, sfallMap.medaxis, sfallMap.slowaxis]
        self.gridSamps = [sfallMap.gridsamp1, sfallMap.gridsamp2,
                          sfallMap.gridsamp3]

        self.atomTaggedMap = sfall.outputMapFile

        return success

    def generatePhenixDensMap(self):

        # run phenix.fobs_minus_fobs to generate difference map

        cmd = 'phenix.fobs_minus_fobs_map ' +\
              'f_obs_1_file_name={} '.format(self.inputMtzFile) +\
              'f_obs_2_file_name={} '.format(self.inputMtzFile) +\
              'phase_source={} '.format(self.reorderedPDBFile) +\
              'f_obs_1_label=FP_{} '.format(self.Mtz2LabelName) +\
              'f_obs_2_label=FP_{}'.format(self.Mtz1LabelName)

        os.system(cmd+'>phenix.log')
        self.fcalcMtz = self.inputMtzFile
        self.densMapMtz = self.inputMtzFile.replace(
            '_CADcombined', 'phenixFoFo')
        shutil.move('FoFoPHFc.mtz', self.densMapMtz)
        shutil.move('phenix.log', '{}phenix.log'.format(self.outputDir))
        self.densMapType = 'FC'
        labelsLater = ['FoFo', '', '', 'PHFc']
        labelsInit = ['', '', '', '']
        self.mapTag = 'DIFF'
        success = self.generateCCP4DensMap(labelsInit, labelsLater)

        return success

    def generateCCP4DensMap(self,
                            labelsInit=[], labelsLater=[]):

        # run FFT job to generate density map

        self.printStepNumber()
        if self.densMapType in ('DIFF', 'SIMPLE'):
            tags = ['FP_', 'SIGFP_', 'FOM_']
            labelsInit = [i+self.Mtz1LabelName for i in tags] +\
                         ['PHIC_'+self.phaseDataset]
            labelsLater = [i+self.Mtz2LabelName for i in tags] +\
                          ['PHIC_'+self.phaseDataset]

        if self.densMapType == '2FOFC':
            if self.useLaterCellDims.upper() == 'TRUE':
                labelsInit = ['']*4
                labelsLater = ['FWT{}'.format(self.Mtz2LabelName),
                               '', '', 'PHIC']
            else:
                labelsLater = ['']*4
                labelsInit = ['FWT{}'.format(self.Mtz2LabelName),
                              '', '', 'PHIC']

        if self.densMapType != 'END':
            if self.useLaterCellDims.upper() == 'TRUE':
                fft = FFTjob(mapType=self.densMapType, mapTag=self.mapTag,
                             FOMweight=self.FOMweight, runLog=self.runLog,
                             pdbFile=self.reorderedPDBFile,
                             mtzFile=self.densMapMtz, outputDir=self.outputDir,
                             axes=self.axes, gridSamps=self.gridSamps,
                             labels1=labelsLater, labels2=labelsInit,
                             lowResCutoff=self.mapResLimits.split(',')[1],
                             highResCutoff=self.mapResLimits.split(',')[0])
            else:
                fft = FFTjob(mapType=self.densMapType, mapTag=self.mapTag,
                             FOMweight=self.FOMweight, runLog=self.runLog,
                             pdbFile=self.reorderedPDBFile,
                             mtzFile=self.densMapMtz, outputDir=self.outputDir,
                             axes=self.axes, gridSamps=self.gridSamps,
                             labels1=labelsInit, labels2=labelsLater,
                             lowResCutoff=self.mapResLimits.split(',')[1],
                             highResCutoff=self.mapResLimits.split(',')[0])
            success = fft.run()
            self.densityMap = fft.outputMapFile

            if success and self.useLaterCellDims.upper() != 'TRUE':
                m = self.multipleMapByFactor(map=self.densityMap)
                self.densityMap = m

        else:
            # run END job if required (may take time to run!!)
            endInputPDB = self.inputPDBfile
            endInputMTZ = ''.join(endInputPDB.split('.')[:-1]+['.mtz'])
            endInputEFF = ''.join(endInputPDB.split('.')[:-1]+['.eff'])

            end = ENDjob(pdbFile=endInputPDB, mtzFile=endInputMTZ,
                         effFile=endInputEFF, outputDir=self.outputDir,
                         gridSamps=self.gridSamps, runLog=self.runLog)
            success = end.run()
            self.densityMap = end.outputMapFile
        return success

    def generateFcalcMap(self,
                         method='FFT2'):

        # generate FC map using FFT

        self.printStepNumber()
        if method == 'FFT':
            fcLabels = ['FC_{}'.format(self.phaseDataset), '', '',
                        'PHIC_'+self.phaseDataset]
            fft_FC = FFTjob(mapType='FC', FOMweight=self.FOMweight,
                            pdbFile=self.reorderedPDBFile,
                            mtzFile=self.fcalcMtz, outputDir=self.outputDir,
                            axes=self.axes, gridSamps=self.gridSamps,
                            labels1=fcLabels,
                            runLog=self.runLog)
            success = fft_FC.run()
            self.FcalcMap = fft_FC.outputMapFile

        elif method == 'SFALL':
            sfall = SFALLjob(inputPDBfile=self.reorderedPDBFile,
                             outputDir=self.outputDir, VDWR=self.sfall_VDWR,
                             symmetrygroup=self.spaceGroup, runLog=self.runLog,
                             gridDimensions=self.sfallGRIDdims,
                             task='fcalc mtz')
            success = sfall.run()
            FcalcMtz = sfall.outputMtzFile

            fft_FC = FFTjob(mapType='FC', FOMweight=self.FOMweight,
                            pdbFile=self.reorderedPDBFile, mtzFile=FcalcMtz,
                            outputDir=self.outputDir, axes=self.axes,
                            gridSamps=self.gridSamps, runLog=self.runLog,
                            labels1=['FCalc', '', '', 'PHICalc'])

            success = fft_FC.run()
            self.FcalcMap = fft_FC.outputMapFile

        elif method == 'FFT2':
            tags = ['FP_', 'SIGFP_', 'FOM_']
            labels1 = [i+self.Mtz2LabelName for i in tags] +\
                      ['PHIC_'+self.phaseDataset]
            labels2 = ['FC_{}'.format(self.phaseDataset), '', '',
                       'PHIC_'+self.phaseDataset]
            fft_FC = FFTjob(mapType='DIFF', mapTag='FC',
                            FOMweight=self.FOMweight, runLog=self.runLog,
                            pdbFile=self.reorderedPDBFile,
                            mtzFile=self.fcalcMtz, outputDir=self.outputDir,
                            axes=self.axes, gridSamps=self.gridSamps,
                            F1Scale=0.0, labels1=labels1, labels2=labels2)
            success = fft_FC.run()

            tmpMap = fft_FC.outputMapFile

            if success:
                m = self.multipleMapByFactor(map=tmpMap)
                self.FcalcMap = m

        return success

    def cropAtmTaggedMapToAsymUnit(self):

        # crop atom-tagged map to asymmetric unit:

        self.printStepNumber()
        mapmask1 = MAPMASKjob(mapFile1=self.atomTaggedMap,
                              outputDir=self.outputDir, runLog=self.runLog)
        success = mapmask1.crop2AsymUnit()
        self.atomTaggedMap = mapmask1.outputMapFile

        return success

    def ensureSameMapAxesOrder(self):

        # switch map axes to match SFALL atom-tagged map if
        # required(only typically required for END maps)

        mapmask = MAPMASKjob(mapFile1=self.densityMap,
                             outputDir=self.outputDir, runLog=self.runLog)
        success = mapmask.switchAxisOrder(order=self.axes,
                                          symGroup=self.spaceGroup)

        self.densityMap = mapmask.outputMapFile

        return success

    def multipleMapByFactor(self,
                            factor=-1.0, map='./untitled.map'):

        # multiple all points in a density map by a value.
        # useful to switch positive and negative in a map

        mapmask = MAPMASKjob(mapFile1=map, outputDir=self.outputDir,
                             runLog=self.runLog)
        mapmask.multipleByFactor(factor=factor, symGroup=self.spaceGroup)

        return mapmask.outputMapFile

    def renumberPDBFile(self):

        # reorder atoms in pdb file since some may be
        # missing now after pdbcur has been run

        self.runLog.writeToLog(
            str='Renumbering input pdb file: {}'.format(self.PDBCURoutputFile))

        self.reorderedPDBFile = self.PDBCURoutputFile.split(
            '_pdbcur.pdb')[0]+'_reordered.pdb'

        pdbin = open(self.PDBCURoutputFile, 'r')
        pdbout = open(self.reorderedPDBFile, 'w')

        counter = 0
        for line in pdbin.readlines():
            if ('ATOM' not in line[0:4] and 'HETATM' not in line[0:6]):
                pdbout.write(line)
            else:
                counter += 1
                pdbout.write(line[0:6])
                new_atomnum = " "*(5-len(str(counter))) + str(counter)
                pdbout.write(new_atomnum)
                pdbout.write(line[11:80]+'\n')
        pdbin.close()
        pdbout.close()
        self.runLog.writeToLog(
            str='Output pdb file: {}'.format(self.reorderedPDBFile))

    def getSpaceGroup(self):

        # parse the space group from the input pdb file

        pdbin = open(self.reorderedPDBFile, 'r')
        for line in pdbin.readlines():
            if line.split()[0] == 'CRYST1':
                self.spaceGroup = line[55:66].replace(' ', '')
                self.runLog.writeToLog(
                    str='Retrieving space group from file:' +
                        '{}.\nSpace group determined to be {}'.format(
                            self.PDBCURoutputFile, self.spaceGroup))

        try:
            self.spaceGroup
        except AttributeError:
            error(
                text='Unable to find space group from file: {}'.format(
                    self.PDBCURoutputFile),
                log=self.runLog, type='warning')
            return False
        return True

    def cropMapToAtomTaggedMap(self,
                               densMap='untitled.map'):

        # crop the density map to exact same dimensions
        # as SFALL atom-tagged map

        # run MAPMASK job to crop fft density map to asym unit
        mapmask2 = MAPMASKjob(mapFile1=densMap, outputDir=self.outputDir,
                              runLog=self.runLog)
        success = mapmask2.crop2AsymUnit()
        if not success:
            return False

        # run MAPMASK job to crop fft density map to same
        # grid sampling dimensions as SFALL atom map
        mapmask3 = MAPMASKjob(mapFile1=mapmask2.outputMapFile,
                              mapFile2=self.atomTaggedMap,
                              outputDir=self.outputDir, runLog=self.runLog)
        success = mapmask3.cropMap2Map()
        self.croppedDensityMap = mapmask3.outputMapFile

        return success

    def reportCroppedMapInfo(self):

        # print the map info (grid size, num voxels)

        self.runLog.writeToLog(
            str='All generated maps should now be restricted to asym unit ' +
                'with properties: ')
        Map = mapTools(mapName=self.atomTaggedMap, logFile=self.runLog)
        Map.printMapInfo()

    def mapConsistencyCheck(self):

        # this function determines whether the atom map and density
        # map calculated using SFALL and FFT are compatible - meaning
        # the grid dimensions/filtering are the same and the ordering
        # of the fast, medium, and slow axes are identical.

        fftMap = mapTools(self.croppedDensityMap)
        fftMap.readHeader()
        sfallMap = mapTools(self.atomTaggedMap)
        sfallMap.readHeader()

        self.runLog.writeToLog(
            str='Checking that atom map (SFALL) and density ' +
                'map (FFT) are compatible...')

        if (sfallMap.gridsamp1 != fftMap.gridsamp1 or
            sfallMap.gridsamp2 != fftMap.gridsamp2 or
                sfallMap.gridsamp3 != fftMap.gridsamp3):
            error(text='Incompatible grid sampling found...',
                  log=self.runLog, type='error')
            return False

        if (sfallMap.fastaxis != fftMap.fastaxis or
            sfallMap.medaxis != fftMap.medaxis or
                sfallMap.slowaxis != fftMap.slowaxis):
            error(text='Incompatible fast,med,slow axes ordering found...',
                  log=self.runLog, type='error')
            return False

        if (sfallMap.numCols != fftMap.numCols or
            sfallMap.numRows != fftMap.numRows or
                sfallMap.numSecs != fftMap.numSecs):
            error(text='Incompatible number of rows, columns and sections...',
                  log=self.runLog, type='error')
            return False

        if sfallMap.getMapSize() != fftMap.getMapSize():
            error(text='Incompatible map file sizes',
                  log=self.runLog, type='error')
            return False

        self.runLog.writeToLog(str='---> success!')
        return True

    def cleanUpDir(self):

        # give option to clean up working directory

        self.runLog.writeToLog(str='\nCleaning up working directory...')

        # move txt files to subdir
        self.makeOutputDir(dirName='{}txtFiles/'.format(self.outputDir))

        for file in os.listdir(self.outputDir):
            if file.endswith('.txt') and file not in self.filesInDir:
                args = [self.outputDir, file]
                shutil.move('{}{}'.format(*args),
                            '{}txtFiles/{}'.format(*args))

    def findFilesInDir(self):

        # find files initially in working directory

        self.filesInDir = os.listdir(self.outputDir)

    def makeOutputDir(self,
                      dirName='./'):

        # if the above sub directory does not exist, make it

        if not os.path.exists(dirName):
            os.makedirs(dirName)
            self.runLog.writeToLog(
                str='New sub directory "{}" '.format(dirName) +
                    'created to contain output files.')

    def printStepNumber(self):

        # print a string indicating the current pipeline
        # step number directory to the command line

        try:
            self.stepNumber
        except AttributeError:
            self.stepNumber = 3

        self.runLog.writeToLog(
            str='\n_______\nSTEP {})'.format(self.stepNumber))

        self.stepNumber += 1
