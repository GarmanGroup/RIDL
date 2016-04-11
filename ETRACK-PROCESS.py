import argparse
from runProcessFiles import process

parser = argparse.ArgumentParser(description='Run the ETRACK pipeline from the command line')

parser.add_argument('-i',type=str,dest='inputFile',action='store',
                   help='An input file containing pdb and mtz file names, as well as relevant mtz column labels')

parser.add_argument('-p',dest='process',action='store_const',default=False,
                   help='Generate electron density difference maps',const=True)

parser.add_argument('-c',dest='calculate',action='store_const',default=False,
                   help='Calculate damage metrics per atom',const=True)

args = parser.parse_args()

if args.process is True:
	p = process(inputFile=args.inputFile,proceedToETRACK=args.calculate)
