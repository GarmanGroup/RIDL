import matplotlib.pyplot as plt
import numpy as np
import os

try:
    import seaborn as sns
    seabornFound = True
except ImportError:
    seabornFound = False


def plotHist(datax=[],  plotType='histogram', xlabel='', title='',
             outDir='./', saveFig=True, fName=''):

    # plot histogram or kde plot

    if seabornFound:
        sns.set_style("ticks")
        sns.set_context("talk", rc={"figure.figsize": (10, 6)})

    fig = plt.figure()

    if not seabornFound:
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

    if seabornFound:
        sns.despine(offset=0, trim=True)

    plt.xlabel(xlabel)
    plt.ylabel(yTitle)
    plt.title(title)

    if saveFig:
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        fig.savefig(outDir + fName)


def plotVxlsPerAtm(pdbName='untitled', where='', vxlsPerAtom={},
                   plotType='histogram', saveFig=True, saveType='.svg',
                   printText=False, returnStats=False):

    # histogram/kde plot of number of voxels per atom
    # plotType in ('histogram','kde','both')

    if printText:
        print 'Plotting {} number of voxels per atom...'.format(plotType)

    datax = [len(vxlsPerAtom[key]) for key in vxlsPerAtom.keys()]

    plotHist(datax=datax, plotType=plotType, xlabel='# voxels per atom',
             title='Plot of # voxels assigned per atom',
             outDir=where + 'plots/voxelsPerAtom_allAtoms/', saveFig=saveFig,
             fName='{}-vxlsperatm_{}{}'.format(pdbName, plotType, saveType))

    if returnStats:
        meanVal = np.mean(datax)
        stdVal = np.std(datax)
        maxVal = max(datax)
        minVal = min(datax)
        return (meanVal, stdVal, maxVal, minVal)
    else:
        return ''


def plotDensForAtm(pdbName='untitled', where='', vxlsPerAtom={},
                   plotType='histogram', saveFig=True, saveType='.svg',
                   printText=False, PDBarray=[]):

    # histogram/kde plot of per-voxel density distribution per atom
    # plotType in ('histogram','kde','both')

    if printText:
        print 'Plotting {} density values per-atom...'.format(plotType)

    for atm in PDBarray:
        datax = vxlsPerAtom[atm.atomnum]

        plotHist(datax=datax, plotType=plotType,
                 xlabel='Map density (e per cubic angstrom)',
                 title='Plot of raw densities assigned to atom: {}'.format(
                    atm.getAtomID()),
                 outDir=where + 'plots/rawDensitiesForAtoms/', saveFig=saveFig,
                 fName='{}-densForAtm-{}_{}{}'.format(
                    pdbName, atm.getAtomID(), plotType, saveType))
