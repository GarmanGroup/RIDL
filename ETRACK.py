import argparse
from runProcessFiles import process

parser = argparse.ArgumentParser(description='Run the ETRACK pipeline from the command line')

parser.add_argument('-i',
					type=str,
					dest='inputFile',
					action='store',
                    help='An input file containing pdb and mtz file names, as well as relevant mtz column labels')

parser.add_argument('-p',
					dest='process',
					action='store_const',
					default=False,
					const=True,
                    help='Generate electron density difference maps')

parser.add_argument('-c',
					dest='calculate',
					action='store_const',
					default=False,
					const=True,
                    help='Calculate damage metrics per atom')

parser.add_argument('-t',
					type=int,
					dest='template',
					action='store',
					default=0,
                    help='Create a template input file to be completed by the user')

parser.add_argument('-j',
					dest='inputFileHelp',
					action='store_const',
					default=False,
					const=True,
                    help='Provide information to help user complete input file')

args = parser.parse_args()

if args.template != 0:
	p = process(run=False,inputFile='templateInputFile.txt')
	p.writeTemplateInputFile(numHigherDoseDatasets=args.template)

if args.inputFileHelp is True:
	p = process(run=False)
	p.howToWriteInputFile()

if args.process is True:
	p = process(inputFile=args.inputFile,proceedToETRACK=args.calculate)
else:
	if args.calculate is True:
		from runETRACK import run as runETRACK
		r = runETRACK()
