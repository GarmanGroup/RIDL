from runRIDL_class import process
import argparse

# an outer layer for the pipeline scripts. This allows 
# the pipeline to be run from the command line by simply 
# calling:
#
# 	python runRIDL.py -i [inputfilename.txt] -pc
#
# This is the recommended run mode for the scripts.

parser = argparse.ArgumentParser(description = 'Run the RIDL pipeline from the command line.')

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
                    help    = 'Calculate damage metrics per atom.')

parser.add_argument('-t',
					type    = int,
					dest    = 'integer',
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

parser.add_argument('--replot',
					dest    = 'replot',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Remake output files. If included, is assumed that a previous run of '+\
                    		  'RIDL has been performed and that a _data.pkl file is available to '+\
                    		  'retrieve the data of the previous run.')

parser.add_argument('-s',
					dest    = 'suppressOutput',
					action  = 'store_const',
					default = True,
					const   = False,
                    help    = 'Suppress all output to command line. All output text '+\
                    		  'will be printed directly to log file and not command line.')

args = parser.parse_args()

# create a template input file to be filled in manually by the user
if args.integer != 0:
	p = process(run       = False,
				inputFile = 'templateInputFile.txt')
	p.writeTemplateInputFile(numHigherDoseDatasets = args.template)

# call the help information
if args.inputFileHelp is True:
	p = process(run = False)
	p.howToWriteInputFile()

# decide whether to include no graphs (no),
# key graphs only (slim) or all graphs (yes)
if args.noGraphs is True:
	plot = 'no'
else:
	if args.slim is True:
		plot = 'slim'
	else:
		plot = 'yes'

# run the pipeline, including generating density and atom-tagged
# maps from input pdb and mtz files (in a damage series), as 
# specified within an input file
if args.process is True:
	p = process(inputFile           = args.inputFile,
				proceedToMetricCalc = args.calculate,
				outputGraphs        = plot,
				cleanUpFinalFiles   = args.cleanUpFinalFiles,
				printOutput         = args.suppressOutput)

# do not run code to create atom-tagged and density maps, but 
# proceed directly to the code to calculate per-atom damage 
# metrics. This works provided that the map generation step
# has been performed beforehand
else:
	if args.calculate is True:
		p = process(inputFile         = args.inputFile,
					skipToMetricCalc  = True,
					outputGraphs      = plot,
					cleanUpFinalFiles = args.cleanUpFinalFiles,
					printOutput       = args.suppressOutput,
					skipToSummaryFile = args.replot)
