# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import math

try:
    import seaborn as sns
    seabornFound = True
except ImportError:
    seabornFound = False

def plotVxlsPerAtm(pdbName     = 'untitled',
                   where       = '',
                   vxlsPerAtom = {},
                   plotType    = 'histogram',
                   saveFig     = True,
                   printText   = True):

    # histogram/kde plot of number of voxels per atom
    # plotType in ('histogram','kde','both')

    if printText is True:
        print 'Plotting {} plot of number of voxels per atom...'.format(plotType)
        
    if seabornFound is True:    
        sns.set_palette("deep", desat=.6)
        sns.set_context(rc={"figure.figsize": (10, 6)})

    fig = plt.figure()

    datax = [len(vxlsPerAtom[key]) for key in vxlsPerAtom.keys()]

    if seabornFound is False:
        plotType = 'histogram'

    if plotType == 'histogram':
        plt.hist(datax, 300, histtype="stepfilled", alpha=.7)
        yTitle = 'Frequency'

    elif plotType == 'kde':
        sns.kdeplot(np.array(datax), shade=True)
        yTitle = 'Normed-Frequency'

    elif plotType == 'both':
        sns.distplot(np.array(datax))
        yTitle = 'Normed-frequency'

    else:
        print 'Unknown plotting type selected.. cannot plot..'
        return
    plt.xlabel('# voxels per atom')
    plt.ylabel(yTitle)
    plt.title('Plot of # voxels assigned per atom'.format(plotType))
    if saveFig is True:
        fig.savefig('{}plots/{}-vxlsperatm_{}.png'.format(where,pdbName,plotType))
    
