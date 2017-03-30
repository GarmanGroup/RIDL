from runRIDL_class import process
import sys
sys.path.insert(0,'./lib')
from checkDependencies import checkDependencies
from rigidBodyRefine import reRefine
import argparse

# an outer layer for the pipeline scripts. This allows 
# the pipeline to be run from the command line by simply 
# calling:
#
# 	python runRIDL.py -i [inputfilename.txt] -pc
#
# This is the recommended run mode for the scripts.

parser = argparse.ArgumentParser(description = 'Run the RIDL pipeline from the command line.')

parser.add_argument('--dependencies',
					dest    = 'checkDependencies',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Check RIDL dependencies are accessible and flag if not')

parser.add_argument('--rigid',
					dest    = 'performRigidBodyRefine',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Perform REFMAC rigid body refinement using initial '+\
                    		 'coordinate model to generate higher dose models')

parser.add_argument('-i',
					type   = str,
					dest   = 'inputFile',
					action = 'store',
                    help   = 'An input file containing pdb and mtz file names, '+\
                    		 'as well as relevant mtz column labels.')

parser.add_argument('-p',
					dest    = 'process',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Generate electron density difference maps.')

parser.add_argument('-c',
					dest    = 'calculate',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Calculate damage metrics per atom. Will output '+\
                    		  'csv-format files of metric per atom upon completion')

parser.add_argument('-t',
					type    = int,
					dest    = 'template',
					action  = 'store',
					default = 0,
                    help    = 'Create a template input file to be completed by the user.')

parser.add_argument('-j',
					dest    = 'inputFileHelp',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Provide information to help user complete input file.')

parser.add_argument('-g',
					dest    = 'noGraphs',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Include if no graphs should be output (for case where '+\
                    		  'seaborn plotting library not accessible, for example).')

parser.add_argument('--slim',
					dest    = 'slim',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Include only key graphs in the generated output. Per-atom per-dataset '+\
                    		  'heatmaps are not produced. Recommended if many datasets on a large '+\
                    		  'structure are to be input, since the resulting heatmap files will be large')

parser.add_argument('-r',
					dest    = 'cleanUpFinalFiles',
					action  = 'store_const',
					default = True,
					const   = False,
                    help    = 'Do not clean up output directory at end of full run. '+\
                    		  'If not included, intermediate map files (e.g. '+\
                    		  'atom-tagged maps) will be included for each '+\
							  'dataset within damage series.')

parser.add_argument('--remove_maps',
					dest    = 'removeMaps',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Remove the RIDL-maps/ subdirectory that is created '+\
                    		  'following a standard run.')

parser.add_argument('-o',
					dest    = 'output',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Create output hmtl file and graphical analysis. If included, it '+\
                    		  'is assumed that -c step of RIDL has already been performed and that '+\
                    		  'a _data.pkl file is available to retrieve the data from this step.')

parser.add_argument('-s',
					dest    = 'suppressOutput',
					action  = 'store_const',
					default = True,
					const   = False,
                    help    = 'Suppress all output to command line. All output text '+\
                    		  'will be printed directly to log file and not command line.')

args = parser.parse_args()

# check RIDL dependencies are present
if args.checkDependencies:
	checkDependencies(checkAll = True)

# create a template input file to be filled in manually by the user
if args.template != 0:
	p = process(run       = False,
				inputFile = 'templateInputFile.txt')
	p.writeTemplateInputFile(numHigherDoseDatasets = args.template)
	print 'Can use -j command for help on how to complete the generated input file'

# call the help information
if args.inputFileHelp:
	p = process(run = False)
	p.howToWriteInputFile()

# decide whether to include no graphs (no),
# key graphs only (slim) or all graphs (yes)
if args.noGraphs:
	plot = 'no'
else:
	if args.slim:
		plot = 'slim'
	else:
		plot = 'yes'

inputFileToUse = args.inputFile

# decide whether rigid body refinement must be run
# first to generate higher dose coordinate models
if args.performRigidBodyRefine:
	if inputFileToUse == None:
		print '\nRun error:\nMust specify input file with -i tag to perform rigid body refinement job'
	else:
		r=reRefine(inputFile = inputFileToUse)
		inputFileToUse = r.newInputFile
		if args.process or args.calculate or args.output:
			print '\nUsing newly generated RIDL input file for remainder of pipeline'
		else:
			print '\nUse this newly generated RIDL input file on next RIDL run'

if args.process or args.calculate or args.output:
	# run the pipeline, including generating density and atom-tagged
	# maps from input pdb and mtz files (in a damage series), as 
	# specified within an input file
	if args.process:
		p = process(inputFile           = inputFileToUse,
					proceedToMetricCalc = args.calculate,
					outputGraphs        = plot,
					cleanUpFinalFiles   = args.cleanUpFinalFiles,
					printOutput         = args.suppressOutput,
					writeSummaryFiles   = args.output,
					keepMapDir          = not args.removeMaps)

	# do not run code to create atom-tagged and density maps, but 
	# proceed directly to the code to calculate per-atom damage 
	# metrics. This works provided that the map generation step
	# has been performed beforehand

	else:
		p = process(inputFile          = inputFileToUse,
					skipToMetricCalc   = True,
					outputGraphs       = plot,
					cleanUpFinalFiles  = args.cleanUpFinalFiles,
					printOutput        = args.suppressOutput,
					skipToSummaryFiles = not args.calculate,
					writeSummaryFiles  = args.output,
					keepMapDir         = not args.removeMaps)



