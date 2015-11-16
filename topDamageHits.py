# -*- coding: utf-8 -*-

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from operator import truediv
from PDBFileManipulation import convertPDBobject_toPDBline_sudoWater
from confidenceIntervalCalculator import mean_confidence_interval

def topNdamsites_printer(PDBmulti,N,where,densmet,initialPDBname):
	# determine the top N atoms with respect to density 
	# metric and print the top N atom results to an 
	# output file

	# for each dataset, determine top N atoms
	for i in range(0,len(PDBmulti[0].meandensity)):

		# create a new pdb file to show the N atoms as sudo-waters
		initialPDBin = open(initialPDBname,'r')
		PDBout = open(where+'topNby'+str(densmet)+'_'+str(i+1)+'.pdb','w')
		for line in initialPDBin.readlines():
			if line.split()[0] in ('CRYST1','SCALE1','SCALE2','SCALE3'):
				PDBout.write(line)
		initialPDBin.close()

		topNfile = open(where+'topNby'+str(densmet)+'_'+str(i+1)+'.txt','w')
		topNfile.write('Counter\tatom type\tres type\tres num\tchain\n')

		# order atoms by density metric
		if densmet == 'mean':
			PDBmulti.sort(key=lambda x: x.meandensity[i])
		elif densmet == 'median':
			PDBmulti.sort(key=lambda x: x.mediandensity[i])
		elif densmet == 'min':
			PDBmulti.sort(key=lambda x: x.mindensity[i])
		elif densmet == 'max':
			PDBmulti.sort(key=lambda x: x.maxdensity[i])
		else:
			print 'Unknown density metric... '
			print 'Choose out of mean,median,max,min'
			print '---> Terminating script...'
			sys.exit()

		for j in range(0,int(N)):
			topNfile.write('%s\t%s\t%s\t%s\t%s\n'\
			%(str(j+1),str(PDBmulti[j].atomtype),str(PDBmulti[j].basetype),
			  str(PDBmulti[j].residuenum),str(PDBmulti[j].chaintype)))

			# write coordinates of sudo-water to new pdb file
			if densmet == 'mean':
				sudoWaterLine = convertPDBobject_toPDBline_sudoWater(PDBmulti[j],PDBmulti[j].meandensity[i])
			elif densmet == 'median':
				sudoWaterLine = convertPDBobject_toPDBline_sudoWater(PDBmulti[j],PDBmulti[j].mediandensity[i])
			elif densmet == 'min':
				sudoWaterLine = convertPDBobject_toPDBline_sudoWater(PDBmulti[j],PDBmulti[j].mindensity[i])
			elif densmet == 'max':
				sudoWaterLine = convertPDBobject_toPDBline_sudoWater(PDBmulti[j],PDBmulti[j].maxdensity[i])
			
			PDBout.write(sudoWaterLine + '\n')
		PDBout.write('END')
		topNfile.close()

def topNdamsites_resibarplotter(PDBmulti,N,where,densmet,plotstyle):
	# determine the top N atoms with respect to density 
	# metric and plot bar plot to show frequency of each
	# residue type present.
	sns.set(style="white", context="talk")
	f, axes = plt.subplots(len(PDBmulti[0].meandensity), 1, figsize=(16, 50), sharex=True)

	# determine full list of residue types present
	uniq_resis = []
	for atom in PDBmulti:
		if atom.basetype not in uniq_resis:
			uniq_resis.append(atom.basetype)
	uniq_resis.sort()
	x = np.array(uniq_resis)

	# determine how many atoms of each res type present
	totfreq_resis = []
	allresis = [atom.basetype for atom in PDBmulti]
	for res in uniq_resis:
		totfreq_resis.append(allresis.count(res))
	print 'Total number of residues of each type:'
	print totfreq_resis

	# for each dataset, determine top N atoms
	csvSummary = open('normalisedFreqSummary.csv','w')

	for i in range(0,len(PDBmulti[0].meandensity)):

		csvSummary.write('Dataset {} currently analysed...\n'.format(i+2))

		# order atoms by density metric
		if densmet == 'mean':
			PDBmulti.sort(key=lambda x: x.meandensity[i])
		elif densmet == 'median':
			PDBmulti.sort(key=lambda x: x.mediandensity[i])
		elif densmet == 'min':
			PDBmulti.sort(key=lambda x: x.mindensity[i])
		elif densmet == 'max':
			PDBmulti.sort(key=lambda x: x.maxdensity[i])
		else:
			print 'Unknown density metric... '
			print 'Choose out of mean,median,max,min'
			print '---> Terminating script...'
			sys.exit()

		# determine for each dataset i, the top N atoms 
		topNatms = []
		counter = 0
		for atom in PDBmulti:
			counter += 1
			if counter <= int(N):
				topNatms.append(atom)
			else:
				print 'Top %s atoms selected' %(str(N))
				break

		resi_list = [atom.basetype for atom in topNatms]

		# calculate frequency of each residue/base type present
		uniq_reslist = []
		resfreq_list = []
		print '\nFor dataset %s:' %(str(i+2))
		for res in resi_list:
			if res not in uniq_reslist:
				uniq_reslist.append(res)
				freq_res = resi_list.count(res)
				print '\t%s ---> %s' %(str(res),str(freq_res))
				resfreq_list.append(freq_res)

		# for non occurring res/base types in top N, want to
		# 'bulk out' the lists to include zero of these (for
		# bar plotting purposes later)
		for res in uniq_resis:
			if res not in uniq_reslist:
				uniq_reslist.append(res)
				resfreq_list.append(0)
		# now order both the res type and res freq lists by res type
		uniq_reslist, resfreq_list = (list(t) for t in zip(*sorted(zip(uniq_reslist,resfreq_list))))

		# determine the fraction of each residue in the top N atoms
		resfrac_list = map(truediv, resfreq_list, totfreq_resis)

		for element in uniq_reslist:
			csvSummary.write('{},'.format(element))
		csvSummary.write('\n')
		for element in resfrac_list:
			csvSummary.write('{},'.format(element))
		csvSummary.write('\n')
	
		# add barplot for this dataset 
		if plotstyle == 'normalised':
			y = np.array(resfrac_list)
			hval = 0
		else:
			y = np.array(resfreq_list)
			hval = 0.000001

		sns.barplot(x, y, ci=None, palette="Paired", hline=hval, ax=axes[i])
		axes[i].set_ylabel('data '+str(i+1))

	plt.xlabel('Residue type', fontsize=18)
	f.suptitle(str(densmet)+' density: top '+str(N)+ ' hits',fontsize=20)
	sns.despine(bottom=True)
	plt.setp(f.axes)
	plt.tight_layout(h_pad=3)
	plt.subplots_adjust(top=0.95)
	if plotstyle == 'normalised':
		f.savefig(str(where) + 'normalised_top'+str(N)+'damsites_residues_'+str(densmet)+'.png')
	else:
		f.savefig(str(where) + 'top'+str(N)+'damsites_residues_'+str(densmet)+'.png')


def topNdamsites_resibarplotter_Single(PDBmulti,N,where,densmet,plotstyle,datasetnum):
	# determine the top N atoms with respect to density 
	# metric and plot bar plot to show frequency of each
	# residue type present. This function plots for a 
	# single specified dataset number
	sns.set(style="white", context="talk")
	f, axes = plt.subplots(1, 1)

	# determine full list of residue types present
	uniq_resis = []
	for atom in PDBmulti:
		if atom.basetype not in uniq_resis:
			uniq_resis.append(atom.basetype)
	uniq_resis.sort()
	x = np.array(uniq_resis)

	# determine how many atoms of each res type present
	totfreq_resis = []
	allresis = [atom.basetype for atom in PDBmulti]
	for res in uniq_resis:
		totfreq_resis.append(allresis.count(res))
	print 'Total number of residues of each type:'
	print totfreq_resis

	# for specified dataset number 'datasetnum', determine top N atoms
	i  = datasetnum

	# order atoms by density metric
	if densmet == 'mean':
		PDBmulti.sort(key=lambda x: x.meandensity[i])
	elif densmet == 'median':
		PDBmulti.sort(key=lambda x: x.mediandensity[i])
	elif densmet == 'min':
		PDBmulti.sort(key=lambda x: x.mindensity[i])
	elif densmet == 'max':
		PDBmulti.sort(key=lambda x: x.maxdensity[i])
	else:
		print 'Unknown density metric... '
		print 'Choose out of mean,median,max,min'
		print '---> Terminating script...'
		sys.exit()

	# determine for dataset i, the top N atoms 
	topNatms = []
	counter = 0
	for atom in PDBmulti:
		counter += 1
		if counter <= int(N):
			topNatms.append(atom)
		else:
			print 'Top %s atoms selected' %(str(N))
			break

	resi_list = [atom.basetype for atom in topNatms]

	# calculate frequency of each residue/base type present
	uniq_reslist = []
	resfreq_list = []
	print '\nFor dataset %s:' %(str(i+1))
	for res in resi_list:
		if res not in uniq_reslist:
			uniq_reslist.append(res)
			freq_res = resi_list.count(res)
			print '\t%s ---> %s' %(str(res),str(freq_res))
			resfreq_list.append(freq_res)

	# for non occurring res/base types in top N, want to
	# 'bulk out' the lists to include zero of these (for
	# bar plotting purposes later)
	for res in uniq_resis:
		if res not in uniq_reslist:
			uniq_reslist.append(res)
			resfreq_list.append(0)
	# now order both the res type and res freq lists by res type
	uniq_reslist, resfreq_list = (list(t) for t in zip(*sorted(zip(uniq_reslist,resfreq_list))))

	print uniq_reslist
	print resfreq_list

	# determine the fraction of each residue in the top N atoms
	resfrac_list = map(truediv, resfreq_list, totfreq_resis)

	# create barplot for this dataset 
	if plotstyle == 'normalised':
		y = np.array(resfrac_list)
	else:
		y = np.array(resfreq_list)

	sns.barplot(x, y, ci=None, palette="Paired", hline=0.000001, ax=axes)
	axes.set_ylabel('data '+str(i+1))

	plt.xlabel('Residue type', fontsize=18)
	f.suptitle(str(densmet)+' density: top '+str(N)+ ' hits',fontsize=20)
	sns.despine(bottom=True)
	plt.setp(f.axes)
	plt.tight_layout(h_pad=3)
	plt.subplots_adjust(top=0.95)
	if plotstyle == 'normalised':
		f.savefig(str(where) + 'normalised_top'+str(N)+'damsites_residues_'+str(densmet)+'-'+str(datasetnum)+'.png')
	else:
		f.savefig(str(where) + 'top'+str(N)+'damsites_residues_'+str(densmet)+'-'+str(datasetnum)+'.png')

def damSitesAboveThreshold_resibarplotter(PDBmulti,where,plotstyle,threshold,perDatasetFraction,multipleOfMean):
	# determine all damage sites above a threshold (calculated a % increase from mean density 
	# metric map value) with respect to density metric and plot bar plot to show frequency of
	# each residue type present.
	# If 'perDatasetFraction' is true, then 'threshold' is a per-dataset multiple of the average Dloss throughout
	# the structure which atoms are tested, otherwise 'threshold' is an overall multiple of average Dloss within
	# first difference map dataset tested (i=0).
	# If 'multipleOfMean' is True then a per dataset threshold is chosen as a multiple of the mean Dloss, else if 
	# 'multipleOfMean' is False then a per dataset threshold is chosen as a multiple of standard deviations away 
	# from the Dloss mean
	numDatasets = min([len(PDBmulti[0].meandensity),9])
	sns.set(style="white", context="talk")
	f, axes = plt.subplots(numDatasets, 1, figsize=(16, 50), sharex=True)

	# determine full list of residue types present
	uniq_resis = []
	for atom in PDBmulti:
		if atom.basetype not in uniq_resis:
			uniq_resis.append(atom.basetype)
	uniq_resis.sort()
	x = np.array(uniq_resis)

	# determine how many atoms of each res type present
	totfreq_resis = []
	allresis = [atom.basetype for atom in PDBmulti]
	for res in uniq_resis:
		totfreq_resis.append(allresis.count(res))
	print 'Total number of residues of each type:'
	print totfreq_resis

	# for each dataset, determine top N atoms
	csvSummary = open('normalisedFreqSummary.csv','w')

	# calculate overall threshold for each dataset based on average Dloss 
	# metric calculated within first dataset
	if perDatasetFraction == False:
		meanOfMetric = np.mean([atom.mindensity[0] for atom in PDBmulti])
		densityThreshold = meanOfMetric*threshold
		print 'Overall threshold set at {} e/A^3'.format(densityThreshold)

	else:
		pass

	for i in range(0,numDatasets):

		csvSummary.write('Dataset {} currently analysed...\n'.format(i+2))

		# determine for each dataset i,the mean of the Dloss density 
		# metric over all atoms in structure
		meanOfMetric = np.mean([atom.mindensity[i] for atom in PDBmulti])
		stdOfMetric = np.std([atom.mindensity[i] for atom in PDBmulti])
		print 'Mean of Dloss metric: {}'.format(meanOfMetric)
		print 'Std of Dloss metric: {}'.format(stdOfMetric)

		# calculate a dataset-specific threshold here
		if perDatasetFraction == True:
			if multipleOfMean == True:
				densityThreshold = meanOfMetric*threshold
			else:
				densityThreshold = meanOfMetric - threshold*stdOfMetric
			print 'Dataset threshold set at {} e/A^3'.format(densityThreshold) 

		else:
			pass

		# determine for each dataset i, the atoms above threshold 
		keptAtoms = []
		for atom in PDBmulti:
			if atom.mindensity[i] <= densityThreshold:
				keptAtoms.append(atom)

		resi_list = [atom.basetype for atom in keptAtoms]

		# calculate frequency of each residue/base type present
		uniq_reslist = []
		resfreq_list = []
		print '\nFor dataset %s:' %(str(i+2))
		for res in resi_list:
			if res not in uniq_reslist:
				uniq_reslist.append(res)
				freq_res = resi_list.count(res)
				print '\t%s ---> %s' %(str(res),str(freq_res))
				resfreq_list.append(freq_res)

		# for non occurring res/base types in kept atoms, want to
		# 'bulk out' the lists to include zero of these (for
		# bar plotting purposes later)
		for res in uniq_resis:
			if res not in uniq_reslist:
				uniq_reslist.append(res)
				resfreq_list.append(0)
		# now order both the res type and res freq lists by res type
		uniq_reslist, resfreq_list = (list(t) for t in zip(*sorted(zip(uniq_reslist,resfreq_list))))

		# determine the fraction of each residue in the top N atoms
		resfrac_list = map(truediv, resfreq_list, totfreq_resis)

		for element in uniq_reslist:
			csvSummary.write('{},'.format(element))
		csvSummary.write('\n')
		for element in resfrac_list:
			csvSummary.write('{},'.format(element))
		csvSummary.write('\n')
	
		# add barplot for this dataset 
		if plotstyle == 'normalised':
			y = np.array(resfrac_list)
			hval = 0
		else:
			y = np.array(resfreq_list)
			hval = 0.000001

		sns.barplot(x, y, ci=None, palette="Paired", hline=hval, ax=axes[i])
		axes[i].set_ylabel('data '+str(i+1))

	plt.xlabel('Residue type', fontsize=18)
	f.suptitle('Top damage sites per residue type. threshold: {}'.format(threshold),fontsize=20)
	sns.despine(bottom=True)
	plt.setp(f.axes)
	plt.tight_layout(h_pad=3)
	plt.subplots_adjust(top=0.95)
	if plotstyle == 'normalised':
		f.savefig('{}DlossNormalisedDamSites_residues_{}.png'.format(where,threshold))
	else:
		f.savefig('{}DlossDamSites_residues_.png'.format(where,threshold))


def topNdamsites_chainbarplotter(PDBmulti,N,where,densmet):
	# determine the top N atoms with respect to density 
	# metric and plot bar plot to show frequency of each
	# chain type present
	sns.set(style="white", context="talk")
	f, axes = plt.subplots(len(PDBmulti[0].meandensity), 1, figsize=(16, 50), sharex=True)

	# determine full list of chain types present
	uniq_chains = []
	for atom in PDBmulti:
		if atom.chaintype not in uniq_chains:
			uniq_chains.append(atom.chaintype)
	uniq_chains.sort()
	x = np.array(uniq_chains)


	# for each dataset, determine top N atoms
	for i in range(0,len(PDBmulti[0].meandensity)):

		# order atoms by density metric
		if densmet == 'mean':
			PDBmulti.sort(key=lambda x: x.meandensity[i])
		elif densmet == 'median':
			PDBmulti.sort(key=lambda x: x.mediandensity[i])
		elif densmet == 'min':
			PDBmulti.sort(key=lambda x: x.mindensity[i])
		elif densmet == 'max':
			PDBmulti.sort(key=lambda x: x.maxdensity[i])
		else:
			print 'Unknown density metric... '
			print 'Choose out of mean,median,max,min'
			print '---> Terminating script...'
			sys.exit()

		# determine for each dataset i, the top N atoms 
		topNatms = []
		counter = 0
		for atom in PDBmulti:
			counter += 1
			if counter <= int(N):
				topNatms.append(atom)
			else:
				print 'Top %s atoms selected' %(str(N))
				break

		chain_list = [atom.chaintype for atom in topNatms]

		# calculate frequency of each chain type present
		uniq_chainlist = []
		chainfreq_list = []
		print '\nFor dataset %s:' %(str(i+1))
		for chain in chain_list:
			if chain not in uniq_chainlist:
				uniq_chainlist.append(chain)
				freq_chain = chain_list.count(chain)
				print '\t%s ---> %s' %(str(chain),str(freq_chain))
				chainfreq_list.append(freq_chain)

		# for non occurring chain types in top N, want to
		# 'bulk out' the lists to include zero of these (for
		# bar plotting purposes later)
		for chain in uniq_chains:
			if chain not in uniq_chainlist:
				uniq_chainlist.append(chain)
				chainfreq_list.append(0)
		# now order both the chain type and chain freq lists by chain type
		uniq_chainlist, chainfreq_list = (list(t) for t in zip(*sorted(zip(uniq_chainlist,chainfreq_list))))

		print uniq_chainlist
		print chainfreq_list

		# add barplot for this dataset 
		y = np.array(chainfreq_list)
		sns.barplot(x, y, ci=None, palette="Paired", hline=.1, ax=axes[i])
		axes[i].set_ylabel('data '+str(i+1))

	plt.xlabel('Chain type', fontsize=18)
	f.suptitle(str(densmet)+' density: top '+str(N)+ ' hits',fontsize=20)
	sns.despine(bottom=True)
	plt.setp(f.axes)
	plt.tight_layout(h_pad=3)
	plt.subplots_adjust(top=0.95)
	f.savefig(str(where) + 'top'+str(N)+'damsites_chains_'+str(densmet)+'.png')


def topNdamsites_chainbarplotter_Single(PDBmulti,N,where,densmet,plotstyle,datasetnum):
	# determine the top N atoms with respect to density 
	# metric and plot bar plot to show frequency of each
	# chain type present. This function plots for a 
	# single specified dataset number
	sns.set(style="white", context="talk")
	f, axes = plt.subplots(1, 1)

	# determine full list of chain types present
	uniq_chains = []
	for atom in PDBmulti:
		if atom.chaintype not in uniq_chains:
			uniq_chains.append(atom.chaintype)
	uniq_chains.sort()
	x = np.array(uniq_chains)

	# determine how many atoms of each chain type present
	totfreq_chains = []
	allchains = [atom.chaintype for atom in PDBmulti]
	for chain in uniq_chains:
		totfreq_chains.append(allchains.count(chain))
	print 'Total number of chains of each type:'
	print totfreq_chains

	# for specified dataset number 'datasetnum', determine top N atoms
	i  = datasetnum

	# order atoms by density metric
	if densmet == 'mean':
		PDBmulti.sort(key=lambda x: x.meandensity[i])
	elif densmet == 'median':
		PDBmulti.sort(key=lambda x: x.mediandensity[i])
	elif densmet == 'min':
		PDBmulti.sort(key=lambda x: x.mindensity[i])
	elif densmet == 'max':
		PDBmulti.sort(key=lambda x: x.maxdensity[i])
	else:
		print 'Unknown density metric... '
		print 'Choose out of mean,median,max,min'
		print '---> Terminating script...'
		sys.exit()

	# determine for dataset i, the top N atoms 
	topNatms = []
	counter = 0
	for atom in PDBmulti:
		counter += 1
		if counter <= int(N):
			topNatms.append(atom)
		else:
			print 'Top %s atoms selected' %(str(N))
			break

	chain_list = [atom.chaintype for atom in topNatms]

	# calculate frequency of each chain type present
	uniq_chainlist = []
	chainfreq_list = []
	print '\nFor dataset %s:' %(str(i+1))
	for chain in chain_list:
		if chain not in uniq_chainlist:
			uniq_chainlist.append(chain)
			freq_chain = chain_list.count(chain)
			print '\t%s ---> %s' %(str(chain),str(freq_chain))
			chainfreq_list.append(freq_chain)

	# for non occurring res/base types in top N, want to
	# 'bulk out' the lists to include zero of these (for
	# bar plotting purposes later)
	for chain in uniq_chains:
		if chain not in uniq_chainlist:
			uniq_chainlist.append(chain)
			chainfreq_list.append(0)
	# now order both the chain type and chain freq lists by chain type
	uniq_chainlist, chainfreq_list = (list(t) for t in zip(*sorted(zip(uniq_chainlist,chainfreq_list))))

	print uniq_chainlist
	print chainfreq_list

	# determine the fraction of each chain in the top N atoms
	chainfrac_list = map(truediv, chainfreq_list, totfreq_chains)

	# create barplot for this dataset 
	if plotstyle == 'normalised':
		y = np.array(chainfrac_list)
	else:
		y = np.array(chainfreq_list)

	sns.barplot(x, y, ci=None, palette="Paired", hline=0.000001, ax=axes)
	axes.set_ylabel('data '+str(i+1))

	plt.xlabel('Chain type', fontsize=18)
	f.suptitle(str(densmet)+' density: top '+str(N)+ ' hits',fontsize=20)
	sns.despine(bottom=True)
	plt.setp(f.axes)
	plt.tight_layout(h_pad=3)
	plt.subplots_adjust(top=0.95)
	if plotstyle == 'normalised':
		f.savefig(str(where) + 'normalised_top'+str(N)+'damsites_chains_'+str(densmet)+'-'+str(datasetnum)+'.png')
	else:
		f.savefig(str(where) + 'top'+str(N)+'damsites_chains_'+str(densmet)+'-'+str(datasetnum)+'.png')



def calculateAvDlossPerDose(PDBmulti,n):
	# calculated the average Dloss metric over all atoms of each residue/nucleotide type within
	# TRAP complex

	# determine full list of residue types present
	uniq_resis = []
	for atom in PDBmulti:
		if atom.basetype not in uniq_resis:
			uniq_resis.append(atom.basetype)
	uniq_resis.sort()
	x = np.array(uniq_resis)

	# for each residue type, find average Dloss at each dose
	# first group atoms by residue type:
	atomsByRes = {} 
	for res in uniq_resis:
		atomsByRes[res] = []
	for atom in PDBmulti:
		atomsByRes[atom.basetype].append(atom.mindensity[0:9])

	# calculate average Dloss over structure for each dose
	avDlossOverall = np.mean([atom.mindensity[0:9] for atom in PDBmulti],0)
	print 'Average Dloss for structure per dose as following:'
	print ','.join([str(val) for val in avDlossOverall])

	# for each residue type, calculate mean Dloss at each dose
	avDlossByRes = {}
	stdDlossByRes = {}
	nthTileDlossByRes = {}
	confIntDlossByRes = {}
	avDlossByResNorm = {}
	for key in atomsByRes.keys():
		avDloss = np.mean(atomsByRes[key],0)
		stdDloss = np.std(atomsByRes[key],0)
		nthTileDloss = np.percentile(atomsByRes[key],n,0)

		avDlossNorm = avDloss - avDlossOverall

		# calculate 95% confidence interval
		confIntDloss = []
		for i in range(0,9):
			confIntDloss.append(mean_confidence_interval([val[i] for val in atomsByRes[key]]))

		avDlossByRes[key] = avDloss
		stdDlossByRes[key] = stdDloss
		nthTileDlossByRes[key] = nthTileDloss
		confIntDlossByRes[key] = confIntDloss
		avDlossByResNorm[key] = avDlossNorm

	# get results
	csvOutput_stdError = open('avDlossPerDataset_stdError.csv','w')
	csvOutput_ConfIntError = open('avDlossPerDataset_confIntError.csv','w')
	csvOutput_nthTile = open('{}thTile_DlossPerDataset.csv'.format(n),'w')
	csvOutput_DlossNormalised = open('avDlossPerDataset_DlossNormalised.csv','w')

	for key in avDlossByRes.keys():
		avDlossList = list([str(element) for element in avDlossByRes[key]])
		stdDlossList = list([str(element) for element in stdDlossByRes[key]])
		nthTileDlossList = list([str(element) for element in nthTileDlossByRes[key]])
		confIntDlossList = list([str(element) for element in confIntDlossByRes[key]])
		avDlossByResNormList = list([str(element) for element in avDlossByResNorm[key]])

		print '\n***\n{}: {}'.format(key,'-->'.join(avDlossList))
		print '{}: {}'.format(key,'-->'.join(stdDlossList))
		print '{}: {}'.format(key,'-->'.join(confIntDlossList))

		csvOutput_stdError.write('{},{}\n'.format(key,','.join(avDlossList)))
		csvOutput_stdError.write('{},{}\n'.format(key,','.join(stdDlossList)))

		csvOutput_nthTile.write('{},{}\n'.format(key,','.join(nthTileDlossList)))

		csvOutput_DlossNormalised.write('{},{}\n'.format(key,','.join(avDlossByResNormList)))

		csvOutput_ConfIntError.write('{},{}\n'.format(key,','.join(avDlossList)))
		csvOutput_ConfIntError.write('{},{}\n'.format(key,','.join(confIntDlossList)))

	csvOutput_stdError.close()
	csvOutput_ConfIntError.close()
	csvOutput_nthTile.close()
	csvOutput_DlossNormalised.close()


def topDlossPerResiType(PDBmulti):
	# for each dose and each residue/nucleotide type, find the top Dloss value and print it

	# print for each dose separately
	for i in range(0,9):
		print 'For dataset {}:'.format(i+2)
		avDloss = np.mean([atom.mindensity[i] for atom in PDBmulti])
		print 'average Dloss: {}'.format(avDloss)
		# get unique residue/nucleotide name keys
		topDloss = {}
		for atom in PDBmulti:
			if atom.basetype not in topDloss.keys():
				topDloss[atom.basetype] = [atom.mindensity[i]]
			else:
				topDloss[atom.basetype].append(atom.mindensity[i])


		for key in topDloss.keys():
			print '{} --> {}'.format(key,(min(topDloss[key])-avDloss))




