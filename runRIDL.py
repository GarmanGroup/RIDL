import sys
sys.path.insert(0, './lib')
from checkDependencies import checkDependencies
import argparse
import sys


# an outer layer for the pipeline scripts. This allows
# the pipeline to be run from the command line by simply
# calling:
#
#   python runRIDL.py -i [inputfilename.txt] -pc
#
# or:
#
#   python runRIDL.py -i [inputfilename.txt] -pco
#
# This is the recommended run mode for the scripts.

ln = 'Run the RIDL pipeline from the command line.'
parser = argparse.ArgumentParser(description=ln)

h = 'Check RIDL dependencies are accessible and flag if not'
parser.add_argument('--dependencies',
                    dest='checkDependencies', action='store_const',
                    default=False, const=True,
                    help='Check RIDL dependencies are ' +
                         'accessible and flag if not')

parser.add_argument('--rigid',
                    dest='performRigidBodyRefine', action='store_const',
                    default=False, const=True,
                    help='Perform REFMAC rigid body refinement using initial' +
                         ' coordinate model to generate higher dose models')

parser.add_argument('-i',
                    type=str, dest='inputFile', action='store',
                    help='An input file containing pdb and mtz file ' +
                         'names, as well as relevant mtz column labels.')

parser.add_argument('-p',
                    dest='makeMaps', action='store_const',
                    default=False, const=True,
                    help='Generate electron density difference maps.')

parser.add_argument('-c',
                    dest='calcMetrics', action='store_const',
                    default=False, const=True,
                    help='Calculate damage metrics per atom. Will output ' +
                         'csv-format files of metric per atom upon completion')

parser.add_argument('-t',
                    type=int, dest='template',
                    action='store', default=0,
                    help='Create a template input file to be ' +
                         'completed by the user.')

parser.add_argument('-j',
                    dest='inputFileHelp', action='store_const',
                    default=False, const=True,
                    help='Provide information to help user ' +
                         'complete input file.')

parser.add_argument('-r',
                    dest='cleanUpFinalFiles', action='store_const',
                    default=True, const=False,
                    help='Do not clean up output directory at end of full ' +
                         'run. If not included, intermediate map files (e.g.' +
                         ' atom-tagged maps) will be included for each ' +
                         'dataset within damage series.')

parser.add_argument('--remove_maps',
                    dest='removeMaps', action='store_const',
                    default=False, const=True,
                    help='Remove the RIDL-maps/ subdirectory that is ' +
                         'created following a standard run.')

parser.add_argument('-o',
                    dest='output', action='store_const',
                    default=False, const=True,
                    help='Create output hmtl file and graphical analysis. ' +
                         'If included, it is assumed that -c step of RIDL ' +
                         'has already been performed and that a _data.pkl ' +
                         'file is available to retrieve the data from ' +
                         'this step.')

parser.add_argument('-s',
                    dest='suppressOutput', action='store_const',
                    default=True, const=False,
                    help='Suppress all output to command line. All output ' +
                         'text will be printed directly to log file and ' +
                         'not command line.')

parser.add_argument('-v',
                    dest='verboseOutput', action='store_const',
                    default=False, const=True,
                    help='Print verbose output to command line.')

args = parser.parse_args()

# check RIDL dependencies are present
if args.checkDependencies:
    checkDependencies(checkAll=True)
    sys.exit()

# these are here so that they don't interfere with dependency checks
from rigidBodyRefine import reRefine
from runRIDL_class import process

# create a template input file to be filled in manually by the user
if args.template != 0:
    p = process(inputFile='templateInputFile.txt')
    p.writeTemplateInputFile(numHigherDoseDatasets=args.template)
    print('Can use -j command for help on how to ' +
          'complete the generated input file')


# call the help information
if args.inputFileHelp:
    p = process()
    p.howToWriteInputFile()

inputFileToUse = args.inputFile

# decide whether rigid body refinement must be run
# first to generate higher dose coordinate models
if args.performRigidBodyRefine:
    if inputFileToUse is None:
        print('\nRIDL run error:\nMust specify input file with ' +
              '-i tag to perform rigid body refinement job')
        sys.exit()
    else:
        r = reRefine(inputFile=inputFileToUse)
        inputFileToUse = r.newInputFile
        if args.makeMaps or args.calcMetrics or args.output:
            print('\nUsing newly generated RIDL input file ' +
                  'for remainder of pipeline')
        else:
            print('\nUse this newly generated RIDL input ' +
                  'file on next RIDL run')

if args.makeMaps or args.calcMetrics or args.output:
    # run the pipeline, including generating density and atom-tagged
    # maps from input pdb and mtz files (in a damage series), as
    # specified within an input file

    if inputFileToUse is None:
        print('\nRIDL run error:\nMust specify input file with ' +
              '-i tag to when using either -p or -c flags')
        sys.exit()

    p = process(inputFile=inputFileToUse,
                makeMaps=args.makeMaps,
                makeMetrics=args.calcMetrics,
                cleanUpFinalFiles=args.cleanUpFinalFiles,
                printverboseOutput=args.verboseOutput,
                printOutput=args.suppressOutput,
                makeSummaryFile=args.output,
                keepMapDir=not args.removeMaps)
    p.run()
