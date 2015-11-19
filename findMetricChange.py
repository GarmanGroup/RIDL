# -*- coding: utf-8 -*-
import numpy as np
from progbar import progress
import sys

def findBchange(initialPDB,multiDoseList,Bmetric):
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
	num_atoms = len(multiDoseList)
	counter = 0

	# ensure atom list ordered by number of atom in structure (atomnum)
	multiDoseList.sort(key=lambda x: x.atomnum)
	initialPDB.sort(key=lambda x: x.atomnum)

	initpdbindices = range(0,len(initialPDB))
	numDatasets = len(multiDoseList[0].densMetric[Bmetric]['Standard']['values'])
	for atom in multiDoseList:

		# unessential loading bar add-in
		counter += 1
		progress(counter, num_atoms, suffix='')

		Inds = ('residuenum','atomtype','basetype','chaintype')
		atomIndentifier = [getattr(atom, attr) for attr in Inds]
		k = -1
		for atomindex in initpdbindices:
			k += 1
			otheratom = initialPDB[atomindex]

	        if atomIndentifier == [getattr(otheratom, att) for att in Inds]:        
				# determine the Bmetric change between all later datasets and initial dataset
				BmetricChange = Bmetric+'Change'
				laterVals = np.array(map(float, atom.densMetric[Bmetric]['Standard']['values']))
				initialVal = np.array([float(otheratom.densMetric[Bmetric]['Standard']['values'])]*numDatasets)
				atom.densMetric[BmetricChange] = list(laterVals - initialVal)
				break

		initpdbindices.pop(k)        
	print '\n---> success...'


