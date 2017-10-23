from CADjob import CADjob
from SCALEITjob import SCALEITjob
from SIGMAAjob import SIGMAAjob
from MAPMASKjob import MAPMASKjob
from PDBCURjob import PDBCURjob
from SFALLjob import SFALLjob
from mapTools import mapTools
from logFile import logFile
from FFTjob import FFTjob
from ENDjob import ENDjob
from errors import error
import os.path
import shutil
import os


class makeMapsFromMTZs():

    # this class handles the generation of suitable MAP-format files
    # from input MTZ-format files. It uses a range of CCP4 programs to
    # achieve this task. In short:
    # SIGMAA performs any figure of merit recalculation (non default)
    # CAD combines columns from mtz files corresponding to different datasets
    # SCALEIT scales FP columns from later dataset mtz to first dataset
    # PDBCUR strips hydrogens, anisotropic B-factors and 2nd conformations
    # from pdb file
    # SFALL generates a tagged atom-map from new pdb file
    # FFT generates a density map (either Fobs-Fobs, 2Fobs-Fcalc etc depending)
    # on the input flags with exact same grid dimensions as atom-map
    # MAPMASK used to crop both maps to an identical volume

    def __init__(self,
                 outputDir='', jobName='untitled-job', log='',
                 mtzIn1='./untitled.mtz', Mtz1LabelName='FP',
                 Mtz1SIGFPlabel='SIGFP', RfreeFlag1='FREE',
                 Mtz1LabelRename='D1', mtzIn2='./untitled2.mtz',
                 Mtz2LabelName='FP',  Mtz2SIGFPlabel='SIGFP',
                 Mtz2LabelRename='D2', mtzIn3='./untitled3.mtz',
                 Mtz3phaseLabel='PHIC', Mtz3FcalcLabel='FC',
                 Mtz3LabelRename='DP', inputPDBfile='./untitled.pdb',
                 densMapType='DIFF', scaleType='ANISOTROPIC',
                 deleteMtzs=True, FOMweight='NONE', sfall_VDWR=1,
                 mapResLimits=',', includeFCmaps=True,
                 useLaterCellDims=True, sfallGRIDdims=[]):

        # specify where output files should be written
        self.outputDir = outputDir
        self.makeOutputDir(dirName=self.outputDir)
        self.findFilesInDir()

        self.jobName = jobName
        self.mtzIn1 = mtzIn1
        self.Mtz1LabelName = Mtz1LabelName
        self.Mtz1SIGFPlabel = Mtz1SIGFPlabel
        self.RfreeFlag1 = RfreeFlag1
        self.Mtz1LabelRename = Mtz1LabelRename
        self.mtzIn2 = mtzIn2
        self.Mtz2LabelName = Mtz2LabelName
        self.Mtz2SIGFPlabel = Mtz2SIGFPlabel
        self.Mtz2LabelRename = Mtz2LabelRename
        self.mtzIn3 = mtzIn3
        self.Mtz3phaseLabel = Mtz3phaseLabel
        self.Mtz3FcalcLabel = Mtz3FcalcLabel
        self.Mtz3LabelRename = Mtz3LabelRename
        self.inputPDBfile = inputPDBfile
        self.densMapType = densMapType
        self.scaleType = scaleType
        self.deleteMtzs = deleteMtzs
        self.FOMweight = FOMweight
        self.sfall_VDWR = sfall_VDWR
        self.mapResLimits = mapResLimits
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

        # copy input mtz files to working directory and rename
        self.moveInputMtzs()

        skipStep = False
        if self.FOMweight == 'recalculate':
            success = self.generateNewFOMcolumn()
            if not success:
                return False

            # if 2FO-FC map required, use FWT column from
            # sigmaa-output mtz (we are done here)
            if self.densMapType == '2FOFC':
                skipStep = True
                self.mtzForMaps = '{}{}_sigmaa.mtz'.format(
                    self.mapProcessDir, self.name2_current)
            else:
                pass

        else:
            self.CADinputMtz1 = self.SIGMAAinputMtz

        if not skipStep:

            success = self.combineMTZcolumns()
            if not success:
                return False

            if self.scaleType != 'NONE':
                success = self.scaleFPcolumnsTogether()
                if not success:
                    return False

            if self.scaleType in ('NONE', 'PHENIX'):
                self.mtzForMaps = self.CADoutputMtz
            else:
                self.mtzForMaps = self.SCALEIToutputMtz

        self.cleanUpDir()

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

        # phenix maps are not currently a tested option
        makeFoFoMapWithPhenix = False
        if self.scaleType == 'PHENIX':
            makeFoFoMapWithPhenix = True
        if makeFoFoMapWithPhenix:
            self.generatePhenixDensMap()
        else:
            self.mapTag = ''
            self.fcalcMtz = self.mtzForMaps
            self.densMapMtz = self.mtzForMaps
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

    def moveInputMtzs(self):

        # move input mtz files to working directory and rename as suitable

        if self.densMapType == '2FOFC':
            self.SIGMAAinputMtz = '{}{}.mtz'.format(
                self.outputDir, self.Mtz2LabelRename.strip())
            shutil.copy2(self.mtzIn2, self.SIGMAAinputMtz)

        else:
            self.SIGMAAinputMtz = '{}{}.mtz'.format(
                self.outputDir, self.Mtz1LabelRename.strip())
            self.CADinputMtz2 = '{}{}.mtz'.format(
                self.outputDir, self.Mtz2LabelRename.strip())
            self.CADinputMtz3 = '{}{}.mtz'.format(
                self.outputDir, self.Mtz3LabelRename.strip())
            shutil.copy2(self.mtzIn1, self.SIGMAAinputMtz)
            shutil.copy2(self.mtzIn2, self.CADinputMtz2)
            shutil.copy2(self.mtzIn3, self.CADinputMtz3)

    def generateNewFOMcolumn(self):

        # run SIGMAA job if required to generate a new FOM weight column

        self.printStepNumber()

        if self.densMapType == '2FOFC':
            mtzLbls_in = self.Mtz2LabelName
            mtzLbls_out = self.Mtz2LabelRename

        else:
            mtzLbls_in = self.Mtz1LabelName
            mtzLbls_out = self.Mtz1LabelName

        sigmaa = SIGMAAjob(inputMtz=self.SIGMAAinputMtz,
                           MtzLabelNameIn=mtzLbls_in,
                           MtzLabelNameOut=mtzLbls_out,
                           RfreeFlag=self.RfreeFlag1,
                           inputPDB=self.inputPDBfile,
                           outputDir=self.outputDir,
                           runLog=self.runLog)
        success = sigmaa.run()
        self.CADinputMtz1 = sigmaa.outputMtz

        return success

    def combineMTZcolumns(self):

        # run CAD to combine necessary columns from input mtz files

        self.printStepNumber()

        self.CADoutputMtz = '{}{}_CADcombined.mtz'.format(
            self.outputDir, self.jobName)

        cad = CADjob(inputMtz1=self.CADinputMtz1, inputMtz2=self.CADinputMtz2,
                     inputMtz3=self.CADinputMtz3,
                     Mtz1LabelName=self.Mtz1LabelName,
                     Mtz2LabelName=self.Mtz2LabelName,
                     Mtz3phaseLabel=self.Mtz3phaseLabel,
                     Mtz3FcalcLabel=self.Mtz3FcalcLabel,
                     Mtz1LabelRename=self.Mtz1LabelRename,
                     Mtz2LabelRename=self.Mtz2LabelRename,
                     Mtz3LabelRename=self.Mtz3LabelRename,
                     outputMtz=self.CADoutputMtz,
                     outputDir=self.outputDir,
                     runLog=self.runLog,
                     FOMWeight=self.FOMweight)
        success = cad.run()

        return success

    def scaleFPcolumnsTogether(self):

        # run SCALEIT to scale later dataset FP columns to first dataset

        self.printStepNumber()

        # specify output files for parts of pipeline

        self.SCALEIToutputMtz = '{}{}_SCALEITcombined.mtz'.format(
            self.outputDir, self.jobName)

        scaleit = SCALEITjob(inputMtz=self.CADoutputMtz,
                             outputMtz=self.SCALEIToutputMtz,
                             Mtz1Label=self.Mtz1LabelRename,
                             Mtz2Label=self.Mtz2LabelRename,
                             outputDir=self.outputDir,
                             scaling=self.scaleType,
                             runLog=self.runLog)
        success = scaleit.run()

        return success

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
              'f_obs_1_file_name={} '.format(self.mtzForMaps) +\
              'f_obs_2_file_name={} '.format(self.mtzForMaps) +\
              'phase_source={} '.format(self.reorderedPDBFile) +\
              'f_obs_1_label=FP_{} '.format(self.Mtz2LabelRename) +\
              'f_obs_2_label=FP_{}'.format(self.Mtz1LabelRename)

        os.system(cmd+'>phenix.log')
        self.fcalcMtz = self.mtzForMaps
        self.densMapMtz = self.mtzForMaps.replace(
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
            labelsInit = [i+self.Mtz1LabelRename for i in tags] +\
                         ['PHIC_'+self.Mtz3LabelRename]
            labelsLater = [i+self.Mtz2LabelRename for i in tags] +\
                          ['PHIC_'+self.Mtz3LabelRename]

        if self.densMapType == '2FOFC':
            if self.useLaterCellDims.upper() == 'TRUE':
                labelsInit = ['']*4
                labelsLater = ['FWT{}'.format(self.Mtz2LabelRename),
                               '', '', 'PHIC']
            else:
                labelsLater = ['']*4
                labelsInit = ['FWT{}'.format(self.Mtz2LabelRename),
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
            fcLabels = ['FC_{}'.format(self.Mtz3LabelRename), '', '',
                        'PHIC_'+self.Mtz3LabelRename]
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
            labels1 = [i+self.Mtz2LabelRename for i in tags] +\
                      ['PHIC_'+self.Mtz3LabelRename]
            labels2 = ['FC_{}'.format(self.Mtz3LabelRename), '', '',
                       'PHIC_'+self.Mtz3LabelRename]
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

        self.deleteNonFinalMtzs()

        # move txt files to subdir
        self.makeOutputDir(dirName='{}txtFiles/'.format(self.outputDir))

        for file in os.listdir(self.outputDir):
            if file.endswith('.txt') and file not in self.filesInDir:
                args = [self.outputDir, file]
                shutil.move('{}{}'.format(*args),
                            '{}txtFiles/{}'.format(*args))

    def deleteNonFinalMtzs(self):

        # delete all non-final mtz files within run
        # (such as those output by CAD before SCALEIT etc).
        # Not current used at runtime!
        return

        # give option to delete all mtz files within
        # output directory except the final resulting
        # mtz for job - used to save room if necessary
        if self.deleteMtzs.lower() != 'true':
            return
        if self.densMapType == '2FOFC':
            fileEnd = 'sigmaa.mtz'
        else:
            fileEnd = 'SCALEITcombined.mtz'
        for f in os.listdir(self.outputDir):
            if ((f.endswith('.mtz') and not f.endswith(fileEnd)) or
               f.endswith('.tmp')):
                os.remove(self.outputDir+f)

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
            self.stepNumber = 1
        self.runLog.writeToLog(
            str='\n_______\nSTEP {})'.format(self.stepNumber))
        self.stepNumber += 1
