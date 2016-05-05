import argparse
from runProcessFiles import process

# an outer layer for the pipeline scripts. This allows the pipeline
# to be run from the command line by simply calling:
# python ETRACK.py -i [inputfilename.txt] -pc
# This is the recommended run mode for the scripts

parser = argparse.ArgumentParser(description = 'Run the ETRACK pipeline from the command line')

parser.add_argument('-i',
					type   = str,
					dest   = 'inputFile',
					action = 'store',
                    help   = 'An input file containing pdb and mtz file names, '+\
                    		 'as well as relevant mtz column labels')

parser.add_argument('-p',
					dest    = 'process',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Generate electron density difference maps')

parser.add_argument('-c',
					dest    = 'calculate',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Calculate damage metrics per atom')

parser.add_argument('-t',
					type    = int,
					dest    = 'integer',
					action  = 'store',
					default = 0,
                    help    = 'Create a template input file to be completed by the user')

parser.add_argument('-j',
					dest    = 'inputFileHelp',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Provide information to help user complete input file')

parser.add_argument('-g',
					dest    = 'noGraphs',
					action  = 'store_const',
					default = True,
					const   = False,
                    help    = 'Include if no graphs should be output (for case where '+\
                    		  'seaborn plotting library not accessible, for example)')

parser.add_argument('-r',
					dest    = 'cleanUpFinalFiles',
					action  = 'store_const',
					default = False,
					const   = True,
                    help    = 'Clean up output directory at end of full run. '+\
                    		  'If not included, intermediate map files (e.g. '+\
                    		  'atom-tagged maps) will be included for each '+\
							  'dataset within damage series')

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

# run the pipeline, including generating density and atom-tagged
# maps from input pdb and mtz files (in a damage series), as 
# specified within an input file
if args.process is True:
	p = process(inputFile         = args.inputFile,
				proceedToETRACK   = args.calculate,
				outputGraphs      = args.noGraphs,
				cleanUpFinalFiles = args.cleanUpFinalFiles)

# do not run code to create atom-tagged and density maps, but 
# proceed directly to the code to calculate per-atom damage 
# metrics. This works provided that the map generation step
# has been performed beforehand
else:
	if args.calculate is True:
		p = process(inputFile         = args.inputFile,
					skipToETRACK      = True,
					outputGraphs      = args.noGraphs,
					cleanUpFinalFiles = args.cleanUpFinalFiles)
