from PDBFileManipulation import PDBtoCLASSARRAY_v2,multiARRAY_diffatomnumbers
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import string
from scipy import stats

def batchRun():
	for i in range(1,11):
		readTwoFiles('TRAP_phenix1.pdb','TRAP_phenix{}.pdb'.format(i),i)

def readTwoFiles(lowDoseFile,highDoseFile,index):
	# read in low dose pdb file
	# lowDoseFile = 'TRAP_phenix1.pdb'
	lowDoseAtoms = PDBtoCLASSARRAY_v2(lowDoseFile,[])

	# read in high dose pdb file
	# highDoseFile = 'TRAP_phenix3.pdb'
	highDoseAtoms = PDBtoCLASSARRAY_v2(highDoseFile,[])

	# combine to one list of lowDoseAtoms
	PDBmulti = multiARRAY_diffatomnumbers([lowDoseAtoms,highDoseAtoms])

	# need to calculate Bfactor change between 2 datasets
	for atom in PDBmulti:
		atom.Bfactorchange = float(atom.Bfactor[1]) - float(atom.Bfactor[0])

	printBfactorPerChain(PDBmulti)
	DlossHistogram(PDBmulti,index)

	return PDBmulti

def printBfactorPerChain(PDBmulti):

	chainDict = {}
	for atm in PDBmulti:
		if atm.chaintype not in chainDict.keys():
			chainDict[atm.chaintype] = [atm.Bfactorchange]
		else:
			chainDict[atm.chaintype].append(atm.Bfactorchange)

	for key in chainDict.keys():
		avBfactorchange = np.mean(chainDict[key])
		print 'Chain {} --> {}'.format(key,avBfactorchange)


def DlossHistogram(atoms,index):
	# plot a histogram of the distribution of Dloss values within a structure
	# for a given dataset number. Dmetric specifies density metric

	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (10, 6)})
	fig = plt.figure()
	ax = plt.subplot(1,1,1)

	counter = -1 
	colors = ['red','blue','green']
	for boundType in ['unbound protein','bound protein','rna']:
		counter += 1
		colorType = colors[counter]

		if boundType == 'unbound protein':
			chains = list(string.ascii_lowercase.upper())[:11]
		elif boundType == 'bound protein':
			chains = list(string.ascii_lowercase.upper())[11:22]
		elif boundType == 'rna':
			chains = list(string.ascii_lowercase.upper())[22:23]

		datax = [atm.Bfactorchange for atm in atoms if atm.chaintype in chains]
		plt.hist(datax, 300, histtype="stepfilled", alpha=.7,color=colorType,label=boundType)

		print "Component '{}' length: {}".format(boundType,len(datax))
		# test for normality in data
		print 'Testing for normality in dataset: {}'.format(boundType)
		W,critVal,SigLevel = stats.anderson(datax,dist='norm')
		print 'Test statistic: {}'.format(W)
		print 'critical value: {}'.format(critVal)
		print 'SigLevel: {}'.format(SigLevel)

		plt.xlabel('Bfactor change value')
		plt.ylabel('Frequency of atoms')
		ax.legend(loc='best')
		plt.title('Histrogram of Bfactor change per atom')
		fig.savefig('BfactorchangePerAtom_{}.png'.format(index))



