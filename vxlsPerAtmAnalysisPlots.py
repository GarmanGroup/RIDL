# -*- coding: utf-8 -*-

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import math

def plotVxlsPerAtm(pdbname,where,vxlsperatom,plotType):
    # histogram/kde plot of number of voxels per atom
    # plotType is 'histogram' or 'kde'
    print '\n------------------------------------------------'
    print 'Plotting {} plot of number of voxels per atom...'.format(plotType)
    sns.set_palette("deep", desat=.6)
    sns.set_context(rc={"figure.figsize": (10, 6)})
    fig = plt.figure()

    datax = [len(vxlsperatom[key]) for key in vxlsperatom.keys()]

    if plotType == 'histogram':
        plt.hist(datax, 300, histtype="stepfilled", alpha=.7)
    elif plotType == 'kde':
        sns.kdeplot(np.array(datax), shade=True);
    else:
        print 'Unknown plotting type selected.. cannot plot..'
        return
    plt.xlabel('Voxels per atom')
    plt.ylabel('Frequency')
    plt.title('{} plot of voxels per atom'.format(plotType))
    fig.savefig('{}output/plots/{}vxlsperatm_{}.png'.format(where,pdbname,plotType))
    