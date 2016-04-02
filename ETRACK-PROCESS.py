import argparse
from runProcessFiles import process

parser = argparse.ArgumentParser(description='Process input files for ETRACK job')
parser.add_argument('-i', type=str, dest='inputFile', action='store',
                   help='An input file containing pdb and mtz file names, as well as relevant mtz column labels')

args = parser.parse_args()
p = process(inputFile=args.inputFile)
