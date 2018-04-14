from matplotlib import pyplot as plt
from scipy import stats
import numpy as np
import seaborn as sns
import sys
import os


def edens_scatter(outputDir='./', metrics=['meandensity', 'mindensity'],
                  PDBarray=[], pdbName='untitled', fileType='.svg',
                  printText=False, savefig=True, titleFont=20, axesFont=18,
                  edgeColor='#FFFFFF', spotColor='#a93b1c', ignoreNan=True):

    # plot scatter plot of two selected metrics. 'metrics' of
    # form ['meandensity','mediandensity'] to choose two metrics
    # of electron density per atom to plot against each other
    # in a scatter plot.

    valsPerParam = []
    for metric in metrics:
        valPerAtom = [getattr(atom, metric) for atom in PDBarray]
        valsPerParam.append(valPerAtom)

    # ignore any points for which either metric is not defined
    if ignoreNan:
        plotData = [[], []]
        for v1, v2 in zip(valsPerParam[0], valsPerParam[1]):
            if np.isnan(v1) or np.isnan(v2):
                continue
            else:
                plotData[0].append(v1)
                plotData[1].append(v2)
    else:
        plotData = valsPerParam

    # check two generated lists of same length:
    if len(plotData[0]) != len(plotData[1]):
        print('Error: lists not same length for scatter plot')
        sys.exit()

    scatter1 = plt.figure()
    plt.scatter(x=plotData[0], y=plotData[1], color=spotColor,
                edgecolors=edgeColor, marker='o', s=100)

    # calculate linear regression and R-squared value
    slope, intercept, r_value, p_value, std_err = stats.linregress(*plotData)

    infoStr = '------------------------------------------------------\n' +\
              'Scatter plot: {} density vs {} density\n'.format(*metrics) +\
              'r-squared: {}\n'.format(r_value**2) +\
              'p-value: {}\n'.format(p_value)

    if printText:
        print(infoStr)

    scatter1.suptitle('{} vs {} density'.format(*metrics), fontsize=titleFont)
    plt.xlabel('{} density'.format(metrics[0]), fontsize=axesFont)
    plt.ylabel('{} density'.format(metrics[1]), fontsize=axesFont)

    sns.despine()
    subdir = 'plots'
    subsubdir = 'metricVmetric-scatterplots'

    if subdir not in os.listdir(outputDir):
        os.mkdir(outputDir+subdir+'/')
    if subsubdir not in os.listdir(outputDir+'/'+subdir+'/'):
        os.mkdir('{}/{}/{}/'.format(outputDir, subdir, subsubdir))

    figName = '{}/{}/{}/{}_{}_vs_{}{}'.format(outputDir, subdir, subsubdir,
                                              pdbName, metrics[0],
                                              metrics[1], fileType)
    if savefig:
        scatter1.savefig(figName)

    return infoStr
