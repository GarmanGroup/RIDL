# -*- coding: utf-8 -*-
import numpy as np
from progbar import progress
import sys

def find_Bchange(initialPDB,PDBmulti,Bmetric):
	# function to determine the Bfactor/Bdamage (specified by Bmetric)
	# change between the initial and later datasets --> becomes an 
	# object attribute for the later datasets

	# check that valid metric specified
	if Bmetric not in ('Bfactor','Bdamage'):
		print 'Unrecognised metric (choose between Bfactor and Bdamage)'
		print '---> terminating script...'
		sys.exit()

	print '------------------------------------------------------------'
	print 'Determining {} change between initial and later datasets'.format(str(Bmetric))
	num_atoms = len(PDBmulti)
	counter = 0

	# ensure atom list ordered by number of atom in structure (atomnum)
	PDBmulti.sort(key=lambda x: x.atomnum)
	initialPDB.sort(key=lambda x: x.atomnum)

	initpdbindices = range(0,len(initialPDB))
	for atom in PDBmulti:

		# unessential loading bar add-in
		counter += 1
		progress(counter, num_atoms, suffix='')

		k = -1
		for atomindex in initpdbindices:
			k += 1
			otheratom = initialPDB[atomindex]
			if (atom.atomtype == otheratom.atomtype and
			   atom.basetype == otheratom.basetype and
			   atom.chaintype == otheratom.chaintype and
			   atom.residuenum == otheratom.residuenum):

				# determine the Bmetric change between all later datasets
				# and initial dataset
				if Bmetric == 'Bdamage':
					atom.densMetric['Bdamagechange'] = list(np.array(map(float, atom.bdam)) - 
									 np.array([float(otheratom.bdam)]*len(PDBmulti[0].meandensity)))
				elif Bmetric == 'Bfactor':
					atom.densMetric['Bfactorchange'] = list(np.array(map(float, atom.Bfactor)) - 
									 np.array([float(otheratom.Bfactor)]*len(PDBmulti[0].meandensity)))
				break

		initpdbindices.pop(k)        

	print '\n---> success...'


