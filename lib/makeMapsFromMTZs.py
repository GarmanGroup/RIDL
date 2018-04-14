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
                 outputDir='', log='', densMapNaming='', atomMapNaming='',
                 mtzIn1='./untitled.mtz', Mtz1LabelName='FP',
                 Mtz1SIGFPlabel='SIGFP', RfreeFlag1='FREE',
                 Mtz1LabelRename='D1', mtzIn2='./untitled2.mtz',
                 Mtz2LabelName='FP',  Mtz2SIGFPlabel='SIGFP',
                 Mtz2LabelRename='D2', mtzIn3='./untitled3.mtz',
                 Mtz3phaseLabel='PHIC', Mtz3FcalcLabel='FC', ignoreSIGFs=False,
                 Mtz3LabelRename='DP', inputPDBfile='./untitled.pdb',
                 densMapType='DIFF', scaleType='anisotropic',
                 deleteIntermediateFiles=True, FOMweight='none',
                 sfall_VDWR=1, mapResLimits=',', includeFCmaps=True,
                 useLaterCellDims=True, sfallGRIDdims=[], mapAxisOrder=[],
                 firstTimeRun=True, premadeAtomMap='', spaceGroup='',
                 gridSampBeforeCropping=[], cropToModel=False):

        # specify where output files should be written
        self.outputDir = outputDir
        self.makeOutputDir(dirName=self.outputDir)
        self.findFilesInDir()

        self.densMapNaming = densMapNaming
        if atomMapNaming == '':
            atomMapNaming = densMapNaming
        self.atomMapNaming = atomMapNaming
        self.scaleType = scaleType

        self.mtzIn1 = mtzIn1
        self.Mtz1LabelName = Mtz1LabelName
        self.RfreeFlag1 = RfreeFlag1
        self.Mtz1LabelRename = Mtz1LabelRename
        self.Mtz1SIGFPlabel = Mtz1SIGFPlabel

        self.mtzIn2 = mtzIn2
        self.Mtz2LabelName = Mtz2LabelName
        self.Mtz2LabelRename = Mtz2LabelRename
        self.Mtz2SIGFPlabel = Mtz2SIGFPlabel

        # Provide option to ignore SIGFs during map calculations.
        # Must ensure scaling set to NONE if this chosen
        self.ignoreSIGFs = ignoreSIGFs
        if ignoreSIGFs:
            if scaleType.lower() != 'none':
                error(text='Cross-dataset scaling specified, but user has ' +
                      'specified to ignore all SIGF terms. Incompatible!',
                      log=self.runLog, type='error')

        if densMapType == 'HIGHONLY':
            if scaleType.lower() != 'none':
                error(text='Cross-dataset scaling specified, but user has ' +
                      'specified densMapType as "HIGHONLY". Incompatible!',
                      log=self.runLog, type='error')

        self.mtzIn3 = mtzIn3
        self.Mtz3phaseLabel = Mtz3phaseLabel
        self.Mtz3FcalcLabel = Mtz3FcalcLabel
        self.Mtz3LabelRename = Mtz3LabelRename
        self.inputPDBfile = inputPDBfile

        self.densMapType = densMapType
        self.deleteIntermediateFiles = deleteIntermediateFiles
        self.FOMweight = FOMweight
        self.sfall_VDWR = sfall_VDWR
        self.mapResLimits = mapResLimits
        self.includeFCmaps = includeFCmaps
        self.useLaterCellDims = useLaterCellDims
        self.sfallGRIDdims = sfallGRIDdims
        self.firstTimeRun = firstTimeRun
        self.spaceGroup = spaceGroup
        self.axes = mapAxisOrder
        self.gridSamps = gridSampBeforeCropping

        # usually not crop to model, but crop to asym unit
        # to get a smaller output map file
        self.cropToModel = cropToModel

        self.findFilesInDir()

        # create log file
        if log == '':
            f = '{}{}_runLog1.log'.format(
              self.outputDir+'RIDL-log', densMapNaming)
            self.runLog = logFile(fileName=f, fileDir=self.outputDir)
        else:
            self.runLog = log

        self.CADoutputMtz = '{}{}_CADcombined.mtz'.format(
            self.outputDir, densMapNaming)

        self.SCALEIToutputMtz = '{}{}_SCALEITcombined.mtz'.format(
            self.outputDir, densMapNaming)

        self.PDBCURoutputFile = '{}{}_PDBCUR.pdb'.format(
            self.outputDir, atomMapNaming)

        self.reorderedPDBFile = '{}{}_curated.pdb'.format(
            self.outputDir, atomMapNaming)

        if premadeAtomMap == '':
            self.atomTaggedMap = '{}{}_SFALL.map'.format(
                self.outputDir, atomMapNaming)
        else:
            self.atomTaggedMap = premadeAtomMap

    def runPipeline(self):

        # run the current subroutine within this class

        self.moveInputMtzs()

        skipStep = False
        if self.FOMweight.lower() == 'recalculate':
            self.generateNewFOMcolumn()

            # if 2FO-FC map required, use FWT column from
            # sigmaa-output mtz (we are done here)
            if self.densMapType == '2FOFC':
                skipStep = True
                self.mtzForMaps = '{}{}_sigmaa.mtz'.format(
                    self.mapProcessDir, self.name2_current)
            else:
                pass

        else:
            if not self.densMapType == 'HIGHONLY':
                self.CADinputMtz1 = self.SIGMAAinputMtz

        if not skipStep:
            self.combineMTZcolumns()

            if self.scaleType.lower() != 'none':
                self.scaleFPcolumnsTogether()

            if self.scaleType.lower() in ('none', 'phenix'):
                self.mtzForMaps = self.CADoutputMtz
            else:
                self.mtzForMaps = self.SCALEIToutputMtz

        self.cleanUpDir()

        if self.firstTimeRun or self.useLaterCellDims:

            self.curatePdbFile()
            self.renumberPDBFile()

            if self.spaceGroup == '':
                self.getSpaceGroup()

            self.getAtomTaggedMap()

            if self.cropToModel:
                self.atomTaggedMap = self.cropMapToModel(self.atomTaggedMap)
            else:
                self.cropAtmTaggedMapToAsymUnit()

        # phenix maps are not currently a tested option
        if self.scaleType == 'PHENIX':
            self.generatePhenixDensMap()
        else:
            self.mapTag = ''
            self.fcalcMtz = self.mtzForMaps
            self.densMapMtz = self.mtzForMaps
            self.generateCCP4DensMap()

        if self.includeFCmaps:
            if self.firstTimeRun or self.useLaterCellDims:
                self.generateFcalcMap()

        if self.densMapType == 'END':
            self.ensureSameMapAxesOrder()

        if self.cropToModel:
            self.densityMap = self.cropMapToModel(self.densityMap)
            if self.includeFCmaps:
                if self.firstTimeRun or self.useLaterCellDims:
                    self.FcalcMap = self.cropMapToModel(self.FcalcMap)
        else:
            mapOut = self.cropMapToAtomTaggedMap(densMap=self.densityMap)
            self.densityMap = mapOut
            if self.includeFCmaps:
                if self.firstTimeRun or self.useLaterCellDims:
                    mapOut = self.cropMapToAtomTaggedMap(densMap=self.FcalcMap)
                    self.FcalcMap = mapOut

        self.reportCroppedMapInfo()
        self.mapConsistencyCheck()

        self.renameFinalMapFiles()
        if self.firstTimeRun or self.useLaterCellDims:
            self.renameFinalPDBFile()

        if self.deleteIntermediateFiles:
            self.cleanUpDir(deleteUnwantedMapFiles=True,
                            deleteUnwantedMtzFiles=True,
                            deleteUnwantedPdbFiles=True)
        else:
            self.cleanUpDir()

        return True

    def moveInputMtzs(self):

        # move input mtz files to working directory and rename as suitable

        if self.densMapType == '2FOFC':
            self.SIGMAAinputMtz = '{}{}.mtz'.format(
                self.outputDir, self.Mtz2LabelRename.strip())
            shutil.copy2(self.mtzIn2, self.SIGMAAinputMtz)

        else:
            if self.densMapType != 'HIGHONLY':
                self.SIGMAAinputMtz = '{}{}.mtz'.format(
                    self.outputDir, self.Mtz1LabelRename.strip())
                shutil.copy2(self.mtzIn1, self.SIGMAAinputMtz)
            self.CADinputMtz2 = '{}{}.mtz'.format(
                self.outputDir, self.Mtz2LabelRename.strip())
            self.CADinputMtz3 = '{}{}.mtz'.format(
                self.outputDir, self.Mtz3LabelRename.strip())
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

        if not success:
            error(text='Failure to successfully generate new FOM column',
                  log=self.runLog)

    def combineMTZcolumns(self):

        # run CAD to combine necessary columns from input mtz files

        self.printStepNumber()

        cad = CADjob(inputMtz2=self.CADinputMtz2,
                     inputMtz3=self.CADinputMtz3,
                     Mtz2FPlabel=self.Mtz2LabelName,
                     Mtz2SIGFPlabel=self.Mtz2SIGFPlabel,
                     Mtz2LabelName=self.Mtz2LabelName,
                     Mtz3phaseLabel=self.Mtz3phaseLabel,
                     Mtz3FcalcLabel=self.Mtz3FcalcLabel,
                     Mtz2LabelRename=self.Mtz2LabelRename,
                     Mtz3LabelRename=self.Mtz3LabelRename,
                     outputMtz=self.CADoutputMtz,
                     outputDir=self.outputDir,
                     runLog=self.runLog,
                     FOMWeight=self.FOMweight,
                     ignoreSIGFs=self.ignoreSIGFs)

        if self.densMapType != 'HIGHONLY':
            cad.inputMtz1 = self.CADinputMtz1
            cad.Mtz1FPlabel = self.Mtz1LabelName
            cad.Mtz1SIGFPlabel = self.Mtz1SIGFPlabel
            cad.Mtz1LabelRename = self.Mtz1LabelRename
        else:
            cad.ignoreDset1 = True

        success = cad.run()

        if not success:
            error(text='Failure to successfully combine mtz files',
                  log=self.runLog)

    def scaleFPcolumnsTogether(self):

        # run SCALEIT to scale later dataset FP columns to first dataset

        self.printStepNumber()

        # specify output files for parts of pipeline

        scaleit = SCALEITjob(inputMtz=self.CADoutputMtz,
                             outputMtz=self.SCALEIToutputMtz,
                             Mtz1Label=self.Mtz1LabelRename,
                             Mtz2Label=self.Mtz2LabelRename,
                             outputDir=self.outputDir,
                             scaling=self.scaleType,
                             runLog=self.runLog)
        success = scaleit.run()

        if not success:
            error(text='Failure to successfully scale Fobs columns together',
                  log=self.runLog)

    def curatePdbFile(self):

        # run pdbcur job

        self.printStepNumber()

        pdbcur = PDBCURjob(inputPDBfile=self.inputPDBfile,
                           outputPDBfile=self.PDBCURoutputFile,
                           outputDir=self.outputDir, runLog=self.runLog)
        success = pdbcur.run()

        if not success:
            error(text='Failure to successfully run PDBCUR',
                  log=self.runLog)

    def renumberPDBFile(self):

        # reorder atoms in pdb file since some may be
        # missing now after pdbcur has been run

        self.runLog.writeToLog(
            str='Renumbering input pdb file: {}'.format(self.PDBCURoutputFile))

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

    def getAtomTaggedMap(self):

        # run SFALL job to generate atom tagged map

        self.printStepNumber()

        sfall = SFALLjob(inputPDBfile=self.reorderedPDBFile,
                         outputDir=self.outputDir, VDWR=self.sfall_VDWR,
                         symmetrygroup=self.spaceGroup, runLog=self.runLog,
                         outputMapFile=self.atomTaggedMap,
                         gridDimensions=self.sfallGRIDdims)
        success = sfall.run()

        sfallMap = mapTools(mapName=self.atomTaggedMap)
        self.axes = [sfallMap.fastaxis, sfallMap.medaxis, sfallMap.slowaxis]
        self.gridSamps = [sfallMap.gridsamp1, sfallMap.gridsamp2,
                          sfallMap.gridsamp3]

        if not success:
            error(text='Failure to successfully generate tagged atom map',
                  log=self.runLog)

    def generatePhenixDensMap(self):

        # run phenix.fobs_minus_fobs to generate difference map
        # NOTE: there may be an unfixed issue with how PHENIX
        # generates this map, in particular it may take the same
        # mtz label set for both f_obs_1_label and f_obs_2_label,
        # despite being separate below. I think this is to do with
        # f_obs_1_file_name = f_obs_2_file_name here.
        # Not currently recommended to use unless carefully checked

        error(text='The current map calculation using PHENIX may not ' +
              'work. It is suggested that this option is not used',
              log=self.runLog, type='warning')

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
        if self.useLaterCellDims:
            success = self.generateCCP4DensMap(labelsInit, labelsLater)
        else:
            success = self.generateCCP4DensMap(labelsLater, labelsInit)

        if not success:
            error(text='Failure to successfully generate map using PHENIX',
                  log=self.runLog)

    def generateCCP4DensMap(self,
                            labelsInit=[], labelsLater=[]):

        # run FFT job to generate density map

        self.printStepNumber()

        densMap = '{}{}_FFT.map'.format(self.outputDir, self.densMapNaming)

        if self.densMapType in ('DIFF', 'HIGHONLY'):
            # note that the distinction in these options only comes
            # apparent in FFTjob.py. Even though HIGHONLY'
            # calculate maps using only the later dataset info, dummy
            # information for labelsInit is provided anyway for simplicity
            tags = ['FP_', 'SIGFP_', 'FOM_']
            labelsInit = [i+self.Mtz1LabelRename for i in tags] +\
                         ['PHIC_'+self.Mtz3LabelRename]
            labelsLater = [i+self.Mtz2LabelRename for i in tags] +\
                          ['PHIC_'+self.Mtz3LabelRename]

        if self.densMapType == '2FOFC':
            if self.useLaterCellDims:
                labelsInit = ['']*4
                labelsLater = ['FWT{}'.format(self.Mtz2LabelRename),
                               '', '', 'PHIC']
            else:
                labelsLater = ['']*4
                labelsInit = ['FWT{}'.format(self.Mtz2LabelRename),
                              '', '', 'PHIC']

        if self.densMapType != 'END':

            fft = FFTjob(
                mapType=self.densMapType, mapTag=self.mapTag,
                FOMweight=self.FOMweight, runLog=self.runLog,
                pdbFile=self.reorderedPDBFile, outputMapFile=densMap,
                mtzFile=self.densMapMtz, useSigLabs=not self.ignoreSIGFs,
                outputDir=self.outputDir, axes=self.axes,
                gridSamps=self.gridSamps, labels1=labelsLater,
                labels2=labelsInit,
                lowResCutoff=self.mapResLimits.split(',')[1],
                highResCutoff=self.mapResLimits.split(',')[0])

            if self.useLaterCellDims:
                fft.labels1 = labelsLater
                fft.labels2 = labelsInit
            else:
                fft.labels1 = labelsInit
                fft.labels2 = labelsLater

            success = fft.run()
            self.densityMap = fft.outputMapFile

            if success and not self.useLaterCellDims:
                m = self.multiplyMapByFactor(map=self.densityMap)
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

        if not success:
            error(text='Failure to successfully generate density map',
                  log=self.runLog)

    def generateFcalcMap(self,
                         method=3):

        # generate FC map. Three methods have currently been proposed:
        # In 1, the Fcalc map is calculated using the FC and PHIC cols
        # in the supplied phase dataset.. due to the order of columns
        # in the input mtz file, the map will take the unit cell dims
        # of the later dataset.
        # In 2, the Fcalc map is also calculated from the phase dataset
        # mtz, however is forced to take the cell dims of dataset 1.
        # In 3, the Fcalc and PHIC cols are instead generated from a
        # pdb file, and so take the cell dims of that pdb file also.
        # The important thing here is that the Fcalc map dimensions
        # must match those of the atom-map at the same dose. If the PDB
        # file is varied between datasets, then only methods 1 and 3
        # will allow this. It may be wise to use method 3 as default
        # since this does not depend on a separate mtz file for phases

        self.printStepNumber()

        fft_FC = FFTjob(
            FOMweight=self.FOMweight, pdbFile=self.reorderedPDBFile,
            outputDir=self.outputDir, axes=self.axes,
            gridSamps=self.gridSamps, runLog=self.runLog)

        if method == 1:
            fcLabels = ['FC_{}'.format(self.Mtz3LabelRename), '', '',
                        'PHIC_'+self.Mtz3LabelRename]

            fft_FC.inputMtzFile = self.fcalcMtz
            fft_FC.mapType = 'FC'
            fft_FC.labels1 = fcLabels

            success = fft_FC.run()
            self.FcalcMap = fft_FC.outputMapFile

        elif method == 2:
            tags = ['FP_', 'SIGFP_', 'FOM_']
            labels1 = [i+self.Mtz2LabelRename for i in tags] +\
                      ['PHIC_'+self.Mtz3LabelRename]
            labels2 = ['FC_{}'.format(self.Mtz3LabelRename), '', '',
                       'PHIC_'+self.Mtz3LabelRename]

            fft_FC.mapType = 'DIFF'
            fft_FC.mapTag = 'FC'
            fft_FC.inputMtzFile = self.fcalcMtz
            fft_FC.F1Scale = 0.0
            fft_FC.labels1 = labels1
            fft_FC.labels2 = labels2

            success = fft_FC.run()

            tmpMap = fft_FC.outputMapFile

            if success:
                m = self.multiplyMapByFactor(map=tmpMap)
                self.FcalcMap = m

        elif method == 3:
            sfall = SFALLjob(inputPDBfile=self.reorderedPDBFile,
                             outputDir=self.outputDir, VDWR=self.sfall_VDWR,
                             symmetrygroup=self.spaceGroup, runLog=self.runLog,
                             gridDimensions=self.sfallGRIDdims,
                             task='mtz from pdb')
            success = sfall.run()
            FcalcMtz = sfall.outputMtzFile

            fft_FC.inputMtzFile = FcalcMtz
            fft_FC.mapType = 'FC'
            fft_FC.labels1 = ['FCalc', '', '', 'PHICalc']

            success = fft_FC.run()
            self.FcalcMap = fft_FC.outputMapFile

        if not success:
            error(text='Failure to successfully generate Fcalc map',
                  log=self.runLog)

    def cropMapToModel(self,
                       map=''):

        # crop a map to the input coordinate model

        self.printStepNumber()
        mapmask1 = MAPMASKjob(
          mapFile1=map, pdbFile=self.inputPDBfile,
          outputDir=self.outputDir, runLog=self.runLog)
        success = mapmask1.crop2model(spaceGroup=self.spaceGroup)

        if not success:
            error(text='Failure to crop map in MAPMASK',
                  log=self.runLog)

        croppedMap = mapmask1.outputMapFile

        return croppedMap

    def cropAtmTaggedMapToAsymUnit(self):

        # crop atom-tagged map to asymmetric unit:

        self.printStepNumber()
        mapmask1 = MAPMASKjob(mapFile1=self.atomTaggedMap,
                              outputDir=self.outputDir, runLog=self.runLog)
        success = mapmask1.crop2AsymUnit()

        if not success:
            error(text='Failure to crop map in MAPMASK',
                  log=self.runLog)

        self.atomTaggedMap = mapmask1.outputMapFile

    def ensureSameMapAxesOrder(self):

        # switch map axes to match SFALL atom-tagged map if
        # required(only typically required for END maps)

        mapmask = MAPMASKjob(mapFile1=self.densityMap,
                             outputDir=self.outputDir, runLog=self.runLog)
        success = mapmask.switchAxisOrder(order=self.axes,
                                          symGroup=self.spaceGroup)

        self.densityMap = mapmask.outputMapFile

        if not success:
            error(text='Failure to ensure correct map axes order in MAPMASK',
                  log=self.runLog)

    def multiplyMapByFactor(self,
                            factor=-1.0, map='./untitled.map'):

        # multiply all points in a density map by a value.
        # useful to switch positive and negative in a map

        mapmask = MAPMASKjob(mapFile1=map, outputDir=self.outputDir,
                             runLog=self.runLog)
        mapmask.multiplyByFactor(factor=factor, symGroup=self.spaceGroup)

        return mapmask.outputMapFile

    def getSpaceGroup(self):

        # parse the space group from the input pdb file

        pdbin = open(self.reorderedPDBFile, 'r')
        for line in pdbin.readlines():
            if line.startswith('CRYST1'):
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
        croppedDensityMap = mapmask3.outputMapFile

        if not success:
            error(text='Failure to successfully crop atom map',
                  log=self.runLog)
        else:
            return croppedDensityMap

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

        fftMap = mapTools(self.densityMap)
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

        if (sfallMap.fastaxis != fftMap.fastaxis or
            sfallMap.medaxis != fftMap.medaxis or
                sfallMap.slowaxis != fftMap.slowaxis):
            error(text='Incompatible fast,med,slow axes ordering found...',
                  log=self.runLog, type='error')

        if (sfallMap.numCols != fftMap.numCols or
            sfallMap.numRows != fftMap.numRows or
                sfallMap.numSecs != fftMap.numSecs):
            error(text='Incompatible number of rows, columns and sections...',
                  log=self.runLog, type='error')

        if sfallMap.getMapSize() != fftMap.getMapSize():
            error(text='Incompatible map file sizes',
                  log=self.runLog, type='error')

        self.runLog.writeToLog(str='---> success!')

    def renameFinalMapFiles(self):

        # rename the final map files

        dm = '{}{}_density.map'.format(self.outputDir, self.densMapNaming)
        am = '{}{}_atoms.map'.format(self.outputDir, self.atomMapNaming)

        if self.useLaterCellDims:
            fm = '{}{}_FC.map'.format(self.outputDir, self.densMapNaming)
        else:
            fm = '{}{}_FC.map'.format(self.outputDir, self.atomMapNaming)

        shutil.move(self.densityMap, dm)
        shutil.move(self.atomTaggedMap, am)

        self.densityMap = dm
        self.atomTaggedMap = am

        if self.includeFCmaps:
            shutil.move(self.FcalcMap, fm)
            self.FcalcMap = fm

    def renameFinalPDBFile(self):

        # rename the final pdb file

        shutil.move(self.reorderedPDBFile,
                    '{}{}.pdb'.format(self.outputDir, self.atomMapNaming))

    def cleanUpDir(self,
                   deleteUnwantedMapFiles=False,
                   deleteUnwantedMtzFiles=False,
                   deleteUnwantedPdbFiles=False):

        # give option to clean up working directory

        self.runLog.writeToLog(str='\nCleaning up working directory...')

        # move txt files to subdir
        self.makeOutputDir(dirName='{}txtFiles/'.format(self.outputDir))

        if deleteUnwantedMapFiles:
            self.runLog.writeToLog(str='\tDeleting unwanted .map files')
        if deleteUnwantedMtzFiles:
            self.runLog.writeToLog(str='\tDeleting unwanted .mtz files')
        if deleteUnwantedPdbFiles:
            self.runLog.writeToLog(str='\tDeleting unwanted .pdb files')

        for f in os.listdir(self.outputDir):
            if f.endswith('.txt') and f not in self.filesInDir:
                args = [self.outputDir, f]
                shutil.move('{}{}'.format(*args),
                            '{}txtFiles/{}'.format(*args))

            if deleteUnwantedMapFiles:
                if f.endswith('.map'):
                    if not (f.endswith('_density.map') or
                            f.endswith('_atoms.map') or
                            f.endswith('_FC.map')):
                        os.remove(self.outputDir+f)

            if deleteUnwantedMtzFiles:
                if f.endswith('.mtz'):
                    os.remove(self.outputDir+f)

            if deleteUnwantedPdbFiles:
                if f.endswith('_PDBCUR.pdb'):
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
