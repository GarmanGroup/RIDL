# -*- coding: utf-8 -*-
import numpy as np
from progbar import progress
import sys


def findBchange(initialPDB, multiDoseList, Bmetric, relative=True):
    # function to determine the Bfactor/Bdamage (specified by Bmetric)
    # change between the initial and later datasets --> becomes an
    # object attribute for the later datasets

    # check that valid metric specified
    if Bmetric not in ('Bfactor', 'Bdamage'):
        print('Unrecognised metric (choose between Bfactor and Bdamage)')
        print('---> terminating script...')
        sys.exit()

    print('------------------------------------------------------------')
    print('Finding {} change between first and later datasets'.format(Bmetric))
    num_atoms = len(multiDoseList)

    # ensure atom list ordered by number of atom in structure (atomnum)
    multiDoseList.sort(key=lambda x: x.atomnum)
    initialPDB.sort(key=lambda x: x.atomnum)

    BmetDic = {}
    initBfacDic = {a.getAtomID(): getattr(a, Bmetric) for a in initialPDB}

    for c, atom in enumerate(multiDoseList):

        # unessential loading bar add-in
        progress(c+1, num_atoms, suffix='')

        atmID = atom.getAtomID()
        try:
            initB = initBfacDic[atmID]
        except KeyError:
            print('Error!! Atom "{}" not present in dataset 1'.format(atmID))
            initB = np.nan
        laterBs = np.array(
            map(float, atom.densMetric[Bmetric]['Standard']['values']))

        if not relative:
            metric = list(laterBs - initB)
        else:
            metric = list((laterBs - initB)/initB)
        BmetDic[atom.getAtomID()] = metric

    print('\n---> success...')
    return BmetDic
