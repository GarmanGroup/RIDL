from errors import error
import numpy as np
import sys


class metricNormalisation(object):
    # class to determine a normalisation scheme per dataset.
    # If normalising to Calpha protein backbone, can then
    # monitor the overall degregation of the electron
    # density map with increasing dose (due to global effects).
    # But can also normalise to anything atom type in PDB file now.

    def __init__(self,
                 atomList=[],
                 normaliseTo=[['', 'CA']],
                 logFile=''):

        # list of atoms over multiple doses
        self.atomList = atomList

        # logfile for RIDL if it exists
        self.logFile = logFile

        # specify what to normalise to. Of the form
        # of [[atom type, res type], [atom type 2, res type 2],...]
        self.normaliseTo = normaliseTo

        # if in doubt, normalise to Calpha backbone
        if normaliseTo == []:
            self.normaliseTo = [['', 'CA']]

        # dictionary of weightings per metric type (mean, gain etc.)
        self.meanweight = {}

        # # dictionary of weightings per metric type (mean, gain etc.)
        self.stdweight = {}

    def calculateWeights(self,
                         metric='loss'):

        # calculate the normalisation weights here

        # don't run function if no list of atoms included
        if len(self.atomList) == 0:
            print 'Need to add list of atoms first'
            return

        # collect the atoms to normalise to for each dataset
        normSet = []
        for a in self.atomList:
            for n in self.normaliseTo:
                if a.atomtype == n[1] and (a.basetype == n[0] or n[0] == ''):
                    normSet.append(a)

        # check if a non empty set of atoms found
        if len(normSet) == 0:
            sys.exit('No suitable atoms found for metric normalisation! ' +
                     '\nCheck they exist')
            error(text='Failed to find the specified set of ' +
                       'atoms to normalise metrics', log=self.logFile)

        # calculate the weighting for each dataset and
        # for each density metric individually here
        vals = [a.densMetric[metric]['Standard']['values'] for a in normSet]
        self.meanweight[metric] = np.nanmean(vals, 0)
        self.stdweight[metric] = np.nanstd(vals, 0)

    def printWeights(self,
                     metric='loss'):
        # print the normalisation weights per dataset to command line.
        # But don't run function if weights have not been calculated yet

        try:
            self.meanweight[metric]
        except AttributeError:
            error(text='Need to calculate the weights for metric ' +
                  '"{}" first (use calculateWeights method)'.format(metric),
                  log=self.logFile)

        print '-------------------------------------------------------'
        print 'Normalisation weights as follows:'
        print '\tDataset\tD{}'.format(metric)
        vals = self.atomList[0].densMetric[metric]['Standard']['values']
        for i in range(0, len(vals)):
            print '\t{}:\t{}'.format(i+1, str(self.weight[metric][i])[:5])
