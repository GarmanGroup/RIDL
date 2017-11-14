import numpy as np


class CalphaWeight(object):
    # class to determine the Calpha weighting per dataset,
    # to monitor the overall degregation of the electron
    # density map with increasing dose (due to global effects)

    def __init__(self,
                 atomList=[]):

        # list of atoms over multiple doses
        self.atomList = atomList

        # dictionary of weightings per metric type (mean, gain etc.)
        self.meanweight = {}

        # # dictionary of weightings per metric type (mean, gain etc.)
        self.stdweight = {}

    def calculateWeights(self,
                         metric='loss'):
        # collect the Calpha atoms and for each dataset
        # number determine the Calpha weight

        # don't run function if no list of atoms included
        if len(self.atomList) == 0:
            print 'Need to add list of atoms first'
            return

        # collect the Calpha atoms here
        Cas = []
        for atom in self.atomList:
            if atom.atomtype == 'CA':
                Cas.append(atom)

        # calculate the weighting for each dataset and
        # for each density metric individually here
        vals = [atom.densMetric[metric]['Standard']['values'] for atom in Cas]
        self.meanweight[metric] = np.nanmean(vals, 0)
        self.stdweight[metric] = np.nanstd(vals, 0)

    def printWeights(self,
                     metric='loss'):
        # print the Calpha weights per dataset to command line.
        # But don't run function if weights have not been calculated yet

        try:
            self.meanweight[metric]
        except AttributeError:
            print 'Need to calculate the weights for metric ' +\
                  '"{}" first (use calculateWeights method)'.format(metric)
            return

        print '-------------------------------------------------------'
        print 'Calpha weights as follows:'
        print '\tDataset\tD{}'.format(metric)
        vals = self.atomList[0].densMetric[metric]['Standard']['values']
        for i in range(0, len(vals)):
            print '\t{}:\t{}'.format(i+1, str(self.weight[metric][i])[:5])
