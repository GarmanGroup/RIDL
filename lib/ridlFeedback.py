from checkDependencies import checkDependencies
from os import path, makedirs, system
from PDBFileManipulation import writePDBline
from time import gmtime, strftime
from errors import error
import numpy as np
import sys
import os


class provideFeedback(object):
    # create series of output feedback files and graphs
    def __init__(self,
                 csvOnly=False, atmsObjs=[], outputDir='./',
                 csvExtent='simple', logFile=[], pklSeries='',
                 doses=[], densMaps=[], inputDir='./', autoRun=True,
                 initialPDB='untitled.pdb', inclFCmetrics=False,
                 normSet=[['', 'CA']]):

        # object containing per atom metric information across all doses
        self.atmsObjs = atmsObjs

        # the output directory
        self.outputDir = os.path.abspath(outputDir)+'/'

        # the log file used for the whole RIDL run
        self.logFile = logFile

        # the pkl file where the per-atom metric info is stored
        # (for reference only)
        self.pklSeries = pklSeries

        # set of atoms to metrics have been normalised against
        self.normSet = normSet

        # decide whether to plot a per-atom metric
        # heatmap to visualise damage sites (very slow)
        self.plotHeatMaps = False

        # list of doses for higher dose datasets
        self.doses = list(map(float, doses.split(',')))

        # list of density map names (for reference only)
        self.densMaps = [d + '_density.map' for d in densMaps]

        # input file directory
        self.inDir = inputDir

        # the input coordinate model file
        self.initialPDB = initialPDB

        # whether to include density-weighted metrics
        self.inclFCmetrics = inclFCmetrics

        # whether to write CSV-format files
        self.writeCsvs = True

        # whether to write HTML-format summary file
        self.writeSumFile = False

        # whether to write PDB-format file of top dam sites
        self.writeTopSites = False

        # whether to include only key metrics (='simple') or all metrics
        self.csvExtent = csvExtent

        if not csvOnly:
            self.writeSumFile = True
            self.writeTopSites = True
            self.outputPlotDir = '{}data/plots/'.format(self.outputDir)
            self.makeOutputDir(dirName=self.outputPlotDir)

        if autoRun:
            self.run()

    def run(self):

        # main procedure for writing output files

        self.checkToIncludeAtomNormSet()

        # no plotting (and no html summary) if seaborn not found
        self.checkSeaborn()

        if self.writeCsvs:
            self.writeCsvFiles()

        if self.writeTopSites:
            self.writeDamSitesToFile()

        if self.writeSumFile:

            if not self.inclNormSet:
                n = 'Standard'
            else:
                n = self.normMetName

            self.summaryHTML(
                primaryMetric='density_weighted_mean_negOnly', primaryNorm=n)

        # create heatmap plots (large files) if requested
        if self.plotHeatMaps:
            subdir = 'metric_heatmap/'
            outDir = self.makeNewPlotSubdir(subdir=subdir)
            for norm in ('Standard', self.normMetName):
                self.atmsObjs.densityMetricHeatMap(
                    saveFig=True, metric='loss', normType=norm,
                    outputDir=outDir)

    def writeCsvFiles(self,
                      inclGroupby=False, inclGainMet=False,
                      inclMeanMet=False, inclBfactor=False, numDP=4):

        # write atom numbers and density metrics to simple
        # csv files,one for each density metric separately.
        # inclGroupby = True also writes csv files for atoms
        # grouped by residue type and atom type

        self.fillerLine()
        self.logFile.writeToLog(
            str='Writing .csv files for per-atom density metrics')

        csvDir = self.outputDir+'csvFiles/'
        self.makeOutputDir(dirName=csvDir)
        self.makeOutputDir(dirName=csvDir+'Standard/')
        if self.inclNormSet:
            self.makeOutputDir(dirName=csvDir+'Normalised/')

        if self.csvExtent == 'simple':
            n = 'Standard'
            metrics = [['loss', n]]

            if inclBfactor:
                metrics += [['Bfactor', n]]

            if inclMeanMet:
                metrics += [['mean', n], ['mean_negOnly', n], ['median', n]]

            if inclGainMet:
                metrics += [['gain', n]]

            if self.inclFCmetrics:

                metrics += [['density_weighted_loss', n],
                            ['density_weighted_mean_negOnly', n]]

                if inclGainMet:
                    metrics += [['density_weighted_gain', n]]

                if inclMeanMet:
                    metrics += [['density_weighted_mean', n]]

            if self.inclNormSet:
                n = self.normMetName
                metrics += [['loss', n]]

                if inclMeanMet:
                    metrics += [['mean', n]]

                if self.inclFCmetrics:
                    metrics += [['density_weighted_loss', n],
                                ['density_weighted_mean_negOnly', n]]

                    if inclMeanMet:
                        metrics += [['density_weighted_mean', n]]
        else:
            metrics = self.atmsObjs.getDensMetrics()

        for densMet in metrics:
            self.logFile.writeToLog(
                str='\tmetric: {}, normalisation: {}'.format(*densMet),
                priority='minor')

            if densMet[1] == 'Standard':
                dr = csvDir+'Standard/'
            else:
                dr = csvDir+'Normalised/'

            self.atmsObjs.writeMetric2File(
                where=dr, metric=densMet[0],
                normType=densMet[1], numDP=numDP)

        if inclGroupby:
            for m in ['loss']:
                for groupBy in ('atomtype', 'residue'):
                    self.atmsObjs.writeMetric2File(
                        where=csvDir+'Standard/', metric=m, groupBy=groupBy)
                    if self.inclNormSet:
                        self.atmsObjs.writeMetric2File(
                            where=csvDir+'Normalised/', groupBy=groupBy,
                            metric=m, normType=self.normMetName)

    def summaryHTML(self,
                    primaryMetric='loss', primaryNorm='X-normalised'):

        # produce a selection of per-dataset summary statistics
        # and write in a HTML summary format

        # perform some initial checks to ensure metric info is present
        if primaryNorm == self.normMetName and not self.inclNormSet:
            error(
                text='Error! Tried to write summary file including ' +
                     'normalised metric but no such metric exists',
                log=self.logFile, type='error')

        self.logFile.writeToLog(
            str='Writing HTML summary output file for metric analysis')

        # create a formatted version of the primary metric
        primMetricName = self.atmsObjs.getFormattedmetricName(
            primaryMetric, primaryNorm)
        primMetricNameTEX = self.atmsObjs.getFormattedmetricName(
            primaryMetric, primaryNorm, form='TEX')

        # decide which metrics to feature in the HTML summary file
        if primaryNorm == self.normMetName:
            norms = ['Standard', self.normMetName]
        else:
            norms = [primaryNorm]
        plotNorm = primaryNorm

        numDsets = self.getNumDatasets()
        summaryFile = open('{}summaryFile.html'.format(self.outputDir), 'w')
        summaryFile.write('<!DOCTYPE html>\n<html>\n')

        # create html head
        summaryFile.write(
            '<head>\n' +
            '<meta name="viewport" content="width=device-width, initial-' +
            'scale=1">\n<link rel="stylesheet" href="http://maxcdn.' +
            'bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css">\n' +
            '<script src="https://ajax.googleapis.com/ajax/libs/jquery/' +
            '1.12.2/jquery.min.js"></script>\n<script src="http://maxcdn' +
            '.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js">' +
            '</script>\n<title>RIDL summary file</title>\n<style>\ntable,' +
            ' th, td {\nborder: 1px solid black;\nborder-collapse: ' +
            'collapse;\n}\nth, td {\npadding: 5px;\ntext-align: center;' +
            '\n}\n</style>\n<script type="text/javascript" async\n' +
            'src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/' +
            '2.7.2/MathJax.js?config=TeX-MML-AM_CHTML"></script>' +
            ' </head>\n')

        summaryFile.write(
            '<body>\n<div class="container">\n<h1>RIDL summary file</h1>\n' +
            'Created: {}<br>\n'.format(
                strftime("%Y-%m-%d %H:%M:%S", gmtime())) +
            'Summary information derived from {}<br>'.format(self.pklSeries) +
            '\nEmail charles.bury@dtc.ox.ac.uk for queries<br>\n' +
            'Number of datasets reported in file: {}<br>\n'.format(numDsets))

        # provide some links to useful output files
        bodyString = '<h3>Links to useful output files</h3><ul>\n'

        # definition of metrics present in the file
        bodyString += r"""
<button type="button" class="panel panel-default" data-toggle="collapse" data-target="#metricDefs">Click for definitions of metrics</button>
<div id="metricDefs" class="collapse">
<h3>Indicators of radiation damage:</h3>
For a given atom, let \(V_\text{atom}\) be the set of density map voxels assigned to the atom.
The \(D_\text{loss}\text{(atom)}\) metric is defined as:
$$D_\text{loss}\text{(atom)} = \max_{v\in V_\text{(atom)}} [ -\rho_\Delta(v)]$$
where \(\rho_\Delta(v)\) is the density map value at voxel \(v\in V_\text{atom}\). This metric provides a worst case scenario for the electron density lost by an atom, and provides a useful indicator to screen for potential damage sites.
<p><p>
The \(D_\text{loss}^{\rho}\text{(a)}\) is defined as:
$$D_\text{loss}^{\rho}\text{(atom)} = \frac{\max_{v \in V_\text{atom}} - \rho_\Delta (v)\cdot \rho_\text{calc} (v) }{\max_{v \in V_\text{atom}} \rho_\text{calc} (v)}$$
where \(\rho_\text{calc}(v)\) is the value of the calculated electron density at voxel \(v\), as derived from the input coordinate model. It provides a weighted-version of the \(D_\text{loss}\text{(atom)}\) metric, where voxels contributing more to the electron density of the atom are up-weighted, and noisy background in the vicinity of each atom should be suppressed.
<p><p>
<h3>Quantifying the magnitude of electron density loss:</h3>
Whereas the \(D_\text{loss}\text{(atom)}\) and \(D_\text{loss}^{\rho}\text{(atom)}\) metrics provide useful indicators for potential damage sites, they does not provide a suitable description of the overall magnitude of electron density lost at an atomic site. It is non-trivial to define such metrics, since metrics that require any averaging or integration of voxel values inside a specified search radius are particularly sensitive to the search radius used, due to the inclusion/exclusion of voxels as the search radius is varied. A possible solution is to apply a voxel-weighting scheme, in which the each voxel is weighted by how much of an atom's electron density contributed to that voxel. The \(D_\text{neg}\text{(atom)}\) is defined as:
$$D_\text{neg}\text{(atom)} = \frac{\sum_{v \in V_\text{atom}^{-}} - \rho_\Delta (v)\cdot \rho_\text{calc} (v) }{\sum_{v \in V_\text{atom}^{-}} \rho_\text{calc}(v) }$$
as a weighted average over all voxels \(V_\text{atom}^{-}\) in the vicinity of the atom that attain \(\rho_\Delta(v) < 0\). 
<p><p>
<h3>The \(C_\alpha\)-normalisation scheme:</h3>
It is anticipated that even for non-radiation sensitive atoms, the above metrics will increase with dose, due mainly to the inevitable increase in density map noise as a result of diminishing diffraction data quality with increasing dose. As such the metric values are insufficient alone to distinguish whether an atom has experienced significant electron density loss, without reference to the rest of the structure. This is particularly relevant if comparing damage between independent crystals. 
A normalisation scheme can be constructed in which metric values are compared relative to the reference set of \(C_\alpha\) atoms of the protein backbone. The normalised version of the \(D_\text{mean}^{-,\,\rho}\text{(atom)}\) metric is:
$$C_\alpha\text{-normalised}\ D_\text{neg}\text{(atom)} = \frac{D_\text{neg}\text{(atom)}-\langle D_\text{neg}\text{(a)}\rangle_{a\in C_\alpha} }{\sigma_{a\in C_\alpha} [D_\text{neg}\text{(a)}]}$$
where \(\langle D_\text{neg}\text{(a)}\rangle_{a\in C_\alpha}\) and \(\sigma_{a\in C_\alpha} [D_\text{neg}\text{(a)}]\) represent the average and standard deviation of \(D_\text{neg}\text{(atom)}\) values attained by the set of \(C_\alpha\) atoms. The \(C_\alpha\)-normalised versions of the other metrics are defined similarly.
</div>
"""

        # provide links to csv files
        bodyString += '<h4>CSV-format metric files:</h4><ul>\n'

        metToDisplay = [['loss', 'Standard']]
        if self.inclNormSet:
            metToDisplay += [['loss', self.normMetName]]
        if self.inclFCmetrics:
            metToDisplay += [['density_weighted_loss', 'Standard'],
                             ['density_weighted_mean_negOnly', 'Standard']]
            if self.inclNormSet:
                metToDisplay += [['density_weighted_loss', self.normMetName],
                                 ['density_weighted_mean_negOnly', self.normMetName]]
        for m in metToDisplay:
            if m[1] == self.normMetName:
                subDir = 'Normalised'
            else:
                subDir = 'Standard'
            t = self.atmsObjs.getFormattedmetricName(m[0], m[1])
            bodyString += '<li><a href = "{}csvFiles/{}/{}-'.format(
                self.outputDir, subDir, m[0]) +\
                '{}.csv">{}</a></li>\n'.format(m[1].replace(' ', ''), t)

        # provide links to top 25 damage site information
        topDamSitesMetric = 'loss'
        metName = self.atmsObjs.getFormattedmetricName(topDamSitesMetric)
        subdir = 'topDamSites/'
        plotDir = self.makeNewPlotSubdir(subdir=subdir)
        figName = self.atmsObjs.getTopAtomsStackedBarplot(
            outputDir=plotDir, metric=topDamSitesMetric)
        bodyString += '</ul><h4>Top 25 damage sites:</h4><ul>\n'
        bodyString += '<li><a href = "{}{}">'.format(
            plotDir, figName.split('/')[-1]) +\
            'Per residue/nucleotide barplot ' +\
            '(detected by {} metric) </a></li>'.format(metName)

        if self.writeTopSites:
            bodyString += '<li><a href = "{}">'.format(self.damSitesPDB) +\
                'PDB-format file (detected by ' +\
                '{} </a></li>'.format(
                    self.atmsObjs.getFormattedmetricName('loss'))

        # create heatmap plots (large files) if requested
        # NOT DEFAULT - TAKES A LONG TIME TO RUN
        if self.plotHeatMaps:
            for m, norm in zip([primaryMetric]*len(norms), norms):
                bodyString += '<li><a href="RIDL-metrics/plots/metric_heatm' +\
                              'ap/heatmap_metric-{}_normalisation'.format(m) +\
                              '-{}.svg">{} D<sub>{}'.format(
                                norm.replace(' ', ''), norm, m) +\
                              '</sub> per-atom heat map</a></li></ul>'
        else:
            bodyString += '</ul>'

        summaryFile.write(bodyString)

        # find top overall damaged atoms and plot line plots
        # versus dataset. Only plot if more than 1 high dose
        # dataset is present.
        if numDsets == 1:
            pass

        else:

            # Distn plots over entire structure to indicate
            # the role of any metric normalisation applied
            subdir = 'metricDistn_allAtoms/'
            plotDir = self.makeNewPlotSubdir(subdir=subdir)
            info = '<h3>Metric distribution over all present atoms</h3>\n<div class = "row">\n'

            for n in norms:
                m = self.atmsObjs.getFormattedmetricName(
                    primaryMetric, n)
                info += '<div class = "col-sm-6">{} metric:</div>\n'.format(m)

            for n in norms:

                m = self.atmsObjs.getFormattedmetricName(
                    primaryMetric, n, form='TEX')

                data, saveName = self.makeDistnPlots(
                    densMet=primaryMetric, normType=n, plotSet=4,
                    outputDir=plotDir, dataset='all', xLabel=m)

                f = '<img class="img-responsive" src="{}{}">'.format(
                    plotDir, saveName.split('/')[-1])

                info += '<div class = "col-sm-6">{}</div>\n'.format(f)

            info += '</div>\n'
            summaryFile.write(info)

            # Line plots of metric values versus dose provided for top and
            # bottom 10 atoms located in structure (top and bottom determined
            # over whole dataset range)

            subdir = 'topDamSites/'
            plotDir = self.makeNewPlotSubdir(subdir=subdir)
            numDamSites = 10

            info = '<h3>Top 10 density loss sites with dose:</h3>\n' +\
                   '<div class = "row">\n'

            for n in norms:
                m = self.atmsObjs.getFormattedmetricName(
                    primaryMetric, n)
                info += '<div class = "col-sm-6">{} metric:</div>\n'.format(m)

            for i, n in enumerate(norms):
                topAtoms = self.atmsObjs.getTopNAtoms(
                    dataset='all', metric=primaryMetric,
                    normType=n, n=numDamSites)

                m = self.atmsObjs.getFormattedmetricName(
                    primaryMetric, n, form='TEX')

                keys = ['chain', 'num', 'res', 'atm']
                atmInfo = {k: [] for k in keys}

                for atm in topAtoms:
                    aInfo = atm.split('-')
                    for i, j in zip(keys, aInfo):
                        atmInfo[i].append(j)

                figName = self.atmsObjs.graphMetric(
                    atomType=atmInfo['atm'], resType=atmInfo['res'],
                    chainType=atmInfo['chain'], resiNum=atmInfo['num'],
                    densMet=primaryMetric, normType=n, saveFig=True,
                    outputDir=self.outputPlotDir+subdir, figTitle=' ',
                    saveName='Lineplot_Metric-D{}'.format(primaryMetric) +
                             '_Normalisation-{}_topDamSites'.format(n),
                    yLabel=m)

                f = '<img class="img-responsive" src="{}{}">'.format(
                        plotDir, figName.split('/')[-1])

                info += '<div class = "col-sm-6">{}</div>\n'.format(f)
            info += '</div>\n'
            summaryFile.write(info)

        # make a set of tabs for easy navigation to each dataset
        summaryFile.write('<h2>Per-dataset analysis</h2>\n')
        tabs = '<ul class="nav nav-tabs"><li class="active">' +\
               '<a data-toggle="tab" href="#dset-tab0">Dataset 1</a></li>'

        for i in range(1, numDsets):
            tabs += '<li><a data-toggle="tab" href="#dset-tab' +\
                    '{}">Dataset {}</a></li>'.format(i, i+1)
        tabs += '</ul>'

        summaryFile.write(tabs)
        summaryFile.write('<div class="tab-content">')

        #######################################################################
        # Per-dataset breakdown analysis begins at this point

        for i in range(numDsets):

            if i != 0:
                str = '<div id="dset-tab{}" class="tab-pane fade">'.format(i)
            else:
                str = '<div id="dset-tab{}"'.format(i) +\
                      ' class="tab-pane fade in active">'
            summaryFile.write(str + '<hr>')

            # start block of collapsible panels here
            summaryFile.write(
                '<div class = "panel-group" id = "datasetinfo{}">'.format(i))

            ###################################################################
            # General information regarding current dataset provided here
            t = 'Dataset info'
            c = 'Number in series : {}<br>\n'.format(i+1) +\
                'Dose (MGy)       : {}<br>\n'.format(self.doses[i]) +\
                'Number of atoms  : {}<br>\n'.format(
                    self.atmsObjs.getNumAtoms()) +\
                'Fourier diff map : <a href ="../RIDL-maps/' +\
                '{}">Download</a><br>\n'.format(self.densMaps[i])

            if self.inclNormSet:
                c += 'Normalisation weight for this dataset: {}<br>\n'.format(
                    round(self.atmsObjs.metricNormWeights.meanweight[
                        primaryMetric][i], 3))

            self.writeHtmlDropDownPanel(
                title=t, content=c, dataset=i, sumFile=summaryFile)

            # # TO MAKE METRIC VERSUS METRIC SCATTER PLOTS, UNCOMMENT THIS BIT
            # subdir = 'metricVmetric-scatterplots/'
            # outDir = self.makeNewPlotSubdir(subdir=subdir)
            # m1s = ['loss', 'loss']
            # m2s = ['density_weighted_loss', 'density_weighted_mean_negOnly']
            # for m1, m2 in zip(m1s, m2s):
            #     self.atmsObjs.compareMetrics(
            #         metric1=m1, metric2=m2, atomtype='',
            #         dSet=i, outputDir=outDir, fileType='.svg')

            ###################################################################
            # Create distribution plots for metric values over whole structure
            # and key residue types that are typically damaged

            t = '{} Distribution Plots'.format(primMetricName)

            subdir = 'metricDistn_allAtoms/'
            plotDir = self.makeNewPlotSubdir(subdir=subdir)

            failedPlots = self.makeDistnPlots(
                densMet=primaryMetric, normType=plotNorm, plotSet=4,
                outputDir=plotDir, dataset=i, xLabel=primMetricNameTEX)

            figCall = '<img class="img-responsive" src="' +\
                      '{}DistnPlot_Residues-all_Metric-D'.format(plotDir) +\
                      '{}_Normalisation-{}_Dataset-{}.svg">'.format(
                        primaryMetric, plotNorm.replace(' ', ''), i)
            figInfo = '<div class = "row">\n' +\
                      '<div class = "col-sm-6">{}</div>\n'.format(figCall)

            params = [['Residues', 1, 'GLUASPCYSMETTYR', 'residue'],
                      ['nucleotides', 3, 'DADCDGDT', 'nucleotide']]

            for paramSet in params:
                subdir = 'metricDistn_key{}/'.format(paramSet[0])
                plotDir = self.makeNewPlotSubdir(subdir=subdir)

                failedPlots, saveName = self.makeDistnPlots(
                    densMet=primaryMetric, normType=plotNorm,
                    plotSet=paramSet[1], outputDir=plotDir, dataset=i,
                    xLabel=primMetricNameTEX)

                if not failedPlots[paramSet[2]]:
                    figCall = '<img class="img-responsive" src="' +\
                              '{}{}"><br>'.format(
                                plotDir, saveName.split('/')[-1])
                    figInfo += '<div class = "col-sm-6">{}</div>\n'.format(
                        figCall)

            figInfo += '</div>\n'

            self.writeHtmlDropDownPanel(
                title=t, content=figInfo, dataset=i, sumFile=summaryFile)

            ###################################################################
            # Determination of top N=25 damage sites for
            # current dataset given here

            subdir = 'topNatoms/'
            plotDir = self.makeNewPlotSubdir(subdir=subdir)

            # choose the order to plot the metrics for the dot plot
            # (including normalisation types). atoms will be ordered
            # by the first metric in the list

            if self.inclFCmetrics:
                metsToPlot = [['density_weighted_mean_negOnly', primaryNorm],
                              ['loss', primaryNorm],
                              ['density_weighted_loss', primaryNorm],
                              ['Bfactor', 'Standard']]
            else:
                metsToPlot = [['loss', primaryNorm],
                              ['Bfactor', 'Standard']]

            firstToPlot = [primaryMetric, plotNorm]
            if firstToPlot in metsToPlot:
                metsToPlot.remove(firstToPlot)
            metsToPlot = [firstToPlot] + metsToPlot
            metsToPlot = [[mt[j] for mt in metsToPlot] for j in range(2)]

            saveName = self.atmsObjs.getTopNAtomsDotPlot(
                dataset=i, metrics=metsToPlot[0], normTypes=metsToPlot[1],
                numHits=25)

            t = 'Top 25 damage sites'
            c = '<img class="img-responsive" src="{}{}">'.format(
                plotDir, saveName.split('/')[-1])

            self.writeHtmlDropDownPanel(
                title=t, content=c, dataset=i, sumFile=summaryFile)

            ###################################################################
            # statistics breakdown (per-atom, per-chain,
            # per-residue, full structure) starts here

            c = 'Key:\n' +\
                '<ul><li>mean: average of metric calculated over all atoms ' +\
                'of a specified type</li>\n<li>std: standard deviation of ' +\
                'metric calculated over all atoms of a specified type</li>' +\
                '\n<li>#atoms: total number of atoms of a specified type' +\
                '</li>\n<li>outliers: assuming a symmetric distn around ' +\
                'the mode, number of atoms that fall outside this domain' +\
                '</li>\n<li>skew: skewness of metric distribution for atoms' +\
                'of specified type</li>\n<li>kurtosis: kurtosis of metric ' +\
                'distribution for atoms of specified type</li>\n' +\
                '<li>asymmetry score: A measure of the asymmetry of ' +\
                'distribution of metric values about 0, taken as the ' +\
                'ratio of the sum of all metric values >0 divided by that ' +\
                'for metric values <0 (needs both >0 and <0 values ' +\
                'present, otherwise N/A</li></ul>'

            c += '<ul class="nav nav-tabs">' +\
                 '<li class="active"><a data-toggle="tab" href="#byatom' +\
                 '{}">By Atom</a></li><li><a data-toggle="tab" '.format(i) +\
                 'href="#byresidue{}">By Residue</a></li>'.format(i) +\
                 '<li><a data-toggle="tab" href="#bychain{}">'.format(i) +\
                 'By Chain</a></li><li><a data-toggle="tab" href=' +\
                 '"#bystructure{}">Full Structure</a></li></ul>'.format(i)

            c += '<div class="tab-content">'

            # per-atom statistics
            c += '<div id="byatom{}" class="tab-pane '.format(i) +\
                 'fade in active">'
            statsOut = self.atmsObjs.getPerAtmtypeStats(
                metric=primaryMetric, normType=primaryNorm, dataset=i)
            c += '{} '.format(primMetricName) +\
                 'metric ranked by mean value:<br><br>\n'
            c += self.convertPlainTxtTable2html(statsOut[0])
            c += '</div>'

            # per-residue statistics
            c += '<div id="byresidue{}" class="tab-pane fade">'.format(i)
            statsOut = self.atmsObjs.getPerResidueStats(
                metric=primaryMetric, normType=primaryNorm, dataset=i)
            c += '{} '.format(primMetricName) +\
                 'ranked by mean value:<br><br>\n'
            c += self.convertPlainTxtTable2html(statsOut[0])
            c += '</div>'

            # per-chain statistics
            c += '<div id="bychain{}" class="tab-pane fade">'.format(i)
            statsOut = self.atmsObjs.getPerChainStats(
                metric=primaryMetric, normType=primaryNorm, dataset=i, n='all')
            c += '{} '.format(primMetricName) +\
                 'ranked by mean value:<br><br>\n'
            c += self.convertPlainTxtTable2html(statsOut[0])
            c += '</div>'

            # full-structure statistics
            c += '<div id="bystructure{}" class="tab-pane fade">'.format(i)
            statsOut = self.atmsObjs.getStructureStats(
                metric=primaryMetric, normType=primaryNorm, dataset=i)
            c += '{} '.format(primMetricName) +\
                 'ranked by mean value:<br><br>\n'
            c += self.convertPlainTxtTable2html(statsOut[0])
            c += '</div>'
            c += '</div>'

            self.writeHtmlDropDownPanel(
                title='Statistics Breakdown', content=c,
                dataset=i, sumFile=summaryFile)

            ###################################################################
            # Detection of any atoms that behaviour unusually
            # (unlike rest of that atomtype) is presented here

            infoString = 'Atoms with unusually <font color="red">high' +\
                         '</font> or <font color="blue">low</font> D<sub>' +\
                         '{}</sub> values relative '.format(primaryMetric) +\
                         'to the mean D<sub>{}'.format(primaryMetric) +\
                         '</sub> value for that specific ' +\
                         'atom type will be reported.<br><br>'

            lastInfoString = ''
            suspAtomLens = []
            notDoneBefore = True

            for t in [6, 5, 4]:
                suspAtoms, highOrLow = self.atmsObjs.detectSuspiciousAtoms(
                    dataset=i, metric=primaryMetric, threshold=t)

                tmpInfoString = '{} atoms found with'.format(len(suspAtoms)) +\
                                ' unusually high/low metric values compared' +\
                                ' to average for that atom type (over ' +\
                                '{} standard devs from average)'.format(t)

                if len(suspAtoms) > 0:

                    if notDoneBefore:
                        infoString += lastInfoString + tmpInfoString
                        notDoneBefore = False
                    else:
                        infoString += tmpInfoString

                    infoString += ':<br>'
                    suspAtomLens.append(len(suspAtoms))
                    for s, h in zip(suspAtoms, highOrLow):
                        if h == 'high':
                            c = 'red'
                        else:
                            c = 'blue'
                        infoString += '<font color="{}">{}'.format(c, s) +\
                                      '</font><br>'
                    infoString += '<br>'

                if len(suspAtomLens) > 1:
                    if suspAtomLens[-1] > 0 and suspAtomLens[-2] > 0:
                        break

                lastInfoString = tmpInfoString + '<br>'

            t = 'Suspicious Atoms'
            c = infoString
            self.writeHtmlDropDownPanel(
                title=t, content=c, dataset=i, sumFile=summaryFile)

            # end block of collapsible panels here
            summaryFile.write('</div>')
            summaryFile.write('</div>')

        summaryFile.write('<div class="tab-content">')

        endString = '</div>\n</body>\n</html>'

        summaryFile.write(endString)
        summaryFile.close()

    def writeHtmlDropDownPanel(self,
                               title='untitled', header='h3',
                               content='no content',  panelColor='alternate',
                               dataset=0, sumFile='untitled.txt'):

        # write the info for a html collapsible panel

        try:
            self.panelIndex
        except AttributeError:
            self.panelIndex = 1

        if panelColor == 'alternate':
            if self.panelIndex % 2 == 0:
                panelColor = 'link'
            else:
                panelColor = 'info'

        if self.panelIndex == 1:
            expanded = 'in'
        else:
            expanded = ''
        i = dataset
        txt = '<div class="panel panel-{}">\n'.format(panelColor) +\
              '<div class="panel-heading">\n' +\
              '<{} class="panel-title">\n'.format(header) +\
              '<a data-toggle="collapse" data-parent="#datasetinfo' +\
              '{}" href="#collapse{}">{}</a>\n'.format(
                i, self.panelIndex, title) +\
              '</{}></div>\n'.format(header) +\
              '<div id="collapse{}" '.format(self.panelIndex) +\
              'class="panel-collapse collapse {}">'.format(expanded) +\
              '<div class="panel-body">{}'.format(content) +\
              '</div>\n</div>\n</div>\n'

        self.panelIndex += 1

        if sumFile != '':
            sumFile.write(txt)
        else:
            return txt

    def convertPlainTxtTable2html(self,
                                  plainText='', numHeaderLines=1,
                                  splitBy='\t\t'):

        # convert a plain text table of values into a html format table

        htmlTable = '<table class="table table-striped">\n'
        for i, l in enumerate(plainText.split('\n')):
            if i < numHeaderLines:
                newStr = '<tr><th>'+'</th><th>'.join(l.split(splitBy)) +\
                         '</th></tr>'
            else:
                newStr = '<tr><td>'+'</td><td>'.join(l.split(splitBy)) +\
                         '</td></tr>'
            htmlTable += newStr+'\n'
        htmlTable += '</table>\n'
        return htmlTable

    def getNumDatasets(self):

        # get number of datasets in damage series

        numDsets = self.atmsObjs.atomList[0].getNumDatasets()
        return numDsets

    def writeDamSitesToFile(self,
                            metric='loss', normType='Standard',
                            numDamSites=25):

        # write top damage sites to .pdb file for each dataset
        # 'numDamSites' takes 'all' or integer

        pdbTemplate = self.get1stDsetPDB()
        if not os.path.exists(pdbTemplate):
            error(
                text='Could not coordinate file "{}" '.format(pdbTemplate) +
                     '--> could not write top damage sites to .pdb file..',
                log=self.logFile, type='warning')
            return

        damPDB = self.atmsObjs.getTopNAtomsPDBfile(
            metric=metric, normType=normType, dataset='all',
            n=numDamSites, pdbFile=pdbTemplate)

        self.damSitesPDB = os.path.abspath(damPDB)

    def colorByMetric(self,
                      metric='loss', normType='Standard',
                      dataset=0, singleRes=''):

        # for the initial pdb file in damage series, convert the Bfactor
        # column to values for the specified metric. If 'singleRes' is
        # specified (use 3-letter residue code) then the average density
        # metric for each atom within this residue is calculated within
        # the structure and only this is output to resulting PDB file,
        # otherwise use ''

        if normType == self.normMetName:
            self.atmsObjs.calcAdditionalMetrics(metric=metric,
                                                normType=normType)

        pdbIn = open(self.get1stDsetPDB(), 'r')
        pdbTemplate = self.get1stDsetPDB()
        fileOut = pdbTemplate.replace('.pdb', '') +\
            '_{}D{}_{}.pdb'.format(normType.replace(" ", ""), metric, dataset)
        if singleRes != '':
            fileOut = fileOut.replace('.pdb', '-{}.pdb'.format(singleRes))

        pdbOut = open(fileOut, 'w')
        pdbOut.write(
            'REMARK\tBfactor column replaced by {} D{} metric values\n'.format(
                normType, metric))
        for l in pdbIn.readlines():
            if l.split()[0] in ('CRYST1', 'SCALE1', 'SCALE2', 'SCALE3'):
                pdbOut.write(l)
            elif 'ATOM' in l.split()[0]:
                break
        pdbIn.close()

        if singleRes == '':
            for atm in self.atmsObjs.atomList:
                dens = atm.densMetric[metric][normType]['values'][dataset]
                if not np.isnan(dens):
                    l = writePDBline(atm, dens)
                    pdbOut.write(l+'\n')
        else:
            atmDic = self.atmsObjs.getAvMetricPerAtmInRes(
                singleRes, metric, normType, dataset)
            for atm in self.atmsObjs.atomList:
                if atm.basetype == singleRes:
                    resNum = atm.residuenum
                    break
            for atm in self.atmsObjs.atomList:
                if atm.residuenum == resNum:
                    dens = atmDic[atm.atomtype]
                    l = writePDBline(atm, dens)
                    pdbOut.write(l+'\n')
        pdbOut.write('END')
        pdbOut.close()

    def visualiseDamSites(self,
                          dataset=0, metric='loss', software='pymol',
                          size=1, savePic=True):

        # open coot/pymol to view top damage sites. size is the
        # density metric scale factor used when visualising damage
        # sites as spheres of vdw = size*{density metric} within pymol

        if software not in ('coot', 'pymol'):
            print('Damage sites can be visualised in either "coot" or "pymol"')
            print('Please make sure the paths to these programs are ' +
                  'correctly set up before proceeding')
            return
        try:
            self.damSitesPDB
        except AttributeError:
            return 'Must run .writeDamSitesToFile() before damage ' +\
                   'sites can be read in {}'.format(software)
        if software == 'coot':
            system('coot -pdb {} -pdb2 {}'.format(
                self.get1stDsetPDB(), self.damSitesPDB[dataset]))
        else:
            # need to write script for pymol to run with
            damSitesTag = (self.damSitesPDB[dataset].split('/')[-1]).replace(
                '.pdb', '')
            pdbTemplate = self.get1stDsetPDB()
            structureTag = pdbTemplate.replace('.pdb', '')
            scriptName = self.outputDir + 'runPymolScript.pml'
            pymolScript = open(scriptName, 'w')

            pymolStr = 'load {}\n'.format(self.get1stDsetPDB()) +\
                       'load {}\n'.format(self.damSitesPDB[dataset]) +\
                       'hide lines\nshow cartoon\n' +\
                       'set antialias, 1\n' +\
                       'set ray_trace_mode, 0\n' +\
                       'set cartoon_fancy_helices, 1\n' +\
                       'set cartoon_side_chain_helper, on\n' +\
                       'set ray_opaque_background, 0\n' +\
                       'show sphere, {}\n'.format(damSitesTag) +\
                       'alter {}, vdw=b*{}\n'.format(damSitesTag, size) +\
                       'rebuild\n' +\
                       'select nearDamage, {} w. 4 of {}\n'.format(
                        structureTag, damSitesTag) +\
                       'show sticks, nearDamage\n' +\
                       'set stick_transparency, 0.6\n' +\
                       'color red, {}\n'.format(damSitesTag) +\
                       'color white, {}\n'.format(structureTag)

            if savePic:
                pass
                # DOESN'T WORK CURRENTLY
                pymolStr += 'ray\n png {}\nquit'.format(
                    self.damSitesPDB[dataset].replace('.pdb', '.png'))

            pymolScript.write(pymolStr)
            pymolScript.close()
            system('pymol -q {}'.format(scriptName))

    def makeDistnPlots(self,
                       densMet='loss', normType='Standard', plotType='both',
                       inclTitle=False, plotSet=1, calcKSstat=False, xLabel='',
                       calcADstat=False, sideOnly=False, resSplit=False,
                       outputDir='', dataset=0, requireAll=False):

        # retrieve the distribution for damage metric values for all
        # atoms of specified residues (see 'resList' below), and then
        # plot a kde plot for each.
        # If 'calcKSstat' is True, perform the 2-sample Kolmogorov-Smirnov
        # test (will only calculate when two histograms plotted on 1 axis)

        if outputDir == '':
            outputDir = self.outputPlotDir

        title = ''
        if plotSet == 1:
            resList = [['GLU', 'ASP', 'CYS', 'MET', 'TYR']]
            title = 'Predicted susceptible residue types'
        elif plotSet == 2:
            resList = [['GLU', 'GLN'], ['ASP', 'ASN'], ['ILE', 'LEU'],
                       ['TYR', 'PHE'], ['GLU', 'ASP'], ['GLU', 'GLY'],
                       ['ASP', 'GLY'], ['ALA', 'GLY'], ['TYR', 'GLY'],
                       ['CYS', 'GLY'], ['MET', 'GLY'], ['ARG', 'LYS'],
                       ['A', 'GLY'], ['G', 'GLY'], ['U', 'GLY']]
        elif plotSet == 3:
            resList = [['DA', 'DC', 'DG', 'DT']]

        elif plotSet == 4:
            resList = ['all']
            resSplit = True
            plotType = 'both'

        elif plotSet == 5:
            resList = [['CYS'], ['GLU'], ['GLN'], ['TYR'], ['PHE'],
                       ['ASP'], ['ASN'], ['LYS'], ['ARG'], ['SER'],
                       ['GLY'], ['DA'], ['DC'], ['DG'], ['DT'],
                       ['A'], ['C'], ['G'], ['U']]
            resSplit = True

        for resGroup in resList:
            k = ''.join(resGroup)
            failedPlots = {k: []}
            if len(resGroup) > 6:
                plotType = 'kde'
            else:
                plotType = 'both'
            data, saveName = self.atmsObjs.graphMetricDistn(
                metric=densMet, normType=normType, valType=dataset,
                plotType=plotType, resiType=resGroup, resSplit=resSplit,
                sideOnly=sideOnly, outputDir=outputDir, printText=False,
                plotTitle=title, inclTitle=inclTitle, calcKSstat=calcKSstat,
                calcADstat=calcADstat, requireAll=requireAll, xLabel=xLabel)
            if data == {}:
                failedPlots[k] = True
            else:
                failedPlots[k] = False

        return failedPlots, saveName

    def checkToIncludeAtomNormSet(self):

        # check whether structure contains the atom
        # set against which metric normalisation will be done

        self.normMetName = 'Calpha normalised'

        if self.normSet == [[]]:
            # don't include norm set if not calculated
            self.inclNormSet = False

        else:
            # check atoms were found, if not ignore the set
            self.inclNormSet = self.atmsObjs.checkSpecificAtomsExist(
                self.normSet)
            if self.normSet != [['', 'CA']]:
                self.normMetName = 'X-normalised'

    def get1stDsetPDB(self):

        # retrieve name of first dataset pdb coordinate file.

        if isinstance(self.initialPDB, list):
            a = self.initialPDB[0]
        elif isinstance(self.initialPDB, str):
            a = self.initialPDB

        if not a.endswith('.pdb'):
            a += '.pdb'

        pdbFile = self.inDir + a

        return pdbFile

    def getSpaceGroup(self):

        # parse the first dataset pdb file and retrieve the space group

        pdbFile = self.get1stDsetPDB()
        pdbin = open(pdbFile, 'r')

        for line in pdbin.readlines():
            if line.split()[0] == 'CRYST1':
                self.spaceGroup = line[55:66].replace(' ', '')
        pdbin.close()

        try:
            self.spaceGroup
        except AttributeError:
            self.logFile.writeToLog(
                str='Unable to find space group from file: {}'.format(pdbFile))
            return False
        return self.spaceGroup

    def fillerLine(self,
                   blank=False):

        # print a filler line to command line

        if not blank:
            ln = '\n----------------------------------------------------------'
        else:
            ln = '\n'
        self.logFile.writeToLog(str=ln)

    def checkSeaborn(self):

        # force no plotting if seaborn library not found.
        # Since the html summary file seriously depends on seaborn
        # also do not allow this to be made

        c = checkDependencies()
        if not c.checkPyPackage(logFile=self.logFile, packageName='seaborn'):
            self.writeSumFile = False

    def makeOutputDir(self,
                      dirName='./'):

        # if the above sub directory does not exist, make it

        if not path.exists(dirName):
            makedirs(dirName)
            self.logFile.writeToLog(
                str='New sub directory "{}" '.format(
                    dirName.replace(self.outputDir, '')) +
                'created to contain output files')

    def makeNewPlotSubdir(self,
                          subdir='untitled/'):

        # make a new sub directory for plots within plots/ subdir

        outDir = self.outputPlotDir + subdir
        self.makeOutputDir(dirName=outDir)

        return outDir

    def endAfter(self):

        # force quit of the class

        sys.exit()
