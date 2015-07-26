# -*- coding: utf-8 -*-

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from operator import truediv
from PDBFileManipulation import convertPDBobject_toPDBline_sudoWater

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


