# -*- coding: utf-8 -*-

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy import stats
from matplotlib.mlab import PCA
from mpl_toolkits.mplot3d import Axes3D
from retrieveAtomList import retrieveAtomList
from confidenceIntervalCalculator import mean_confidence_interval

def TRAPRNA_heterostudy(PDBmulti,datasetnum,atomtypes,resitypes):
	# function to plot violin plots for the TRAP-RNA case study
	# to determine whether there is any heterogeneity between 
	# the spread of protein damage between the RNA-bound and non-RNA
	# protein TRAP macromolecules

	# specify which chains are TRAP RNA bound, not RNA bound 
	# and also RNA
	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']
	RNA = ['W']

	# determine the atom types, residue types to be included
	PDBmulti_included = []
	for atom in PDBmulti:
		if atomtypes != 'all' and resitypes != 'all':
			if (atom.atomtype in atomtypes and atom.basetype in resitypes):
				PDBmulti_included.append(atom)
		elif atomtypes == 'all' and resitypes != 'all':
			if atom.basetype in resitypes:
				PDBmulti_included.append(atom)
		elif atomtypes != 'all' and resitypes == 'all':
			if atom.atomtype in atomtypes:
				PDBmulti_included.append(atom)
		else:
			PDBmulti_included.append(atom)

	# loop to only include atoms in RNA bound TRAP closer than
	# 4.5 Angstroms to RNA. 
	distancelist = []
	PDBmulti_included2 = []
	for atom in PDBmulti_included:
		if atom.chaintype in RNAbound:
			for otheratom in PDBmulti:
				if otheratom.chaintype in RNA:
					distance = np.sqrt(np.square(atom.X_coord - otheratom.X_coord) +
	                           		   np.square(atom.Y_coord - otheratom.Y_coord) + 
	                                   np.square(atom.Z_coord - otheratom.Z_coord))
					distancelist.append(distance)
					if distance < 4.5:
						PDBmulti_included2.append(atom)
						break
	# additional loop to include respective atoms in RNA unbound TRAP
	# that would be closer than 4.5 Angstroms to RNA (if bound)
	for atom in PDBmulti_included:
		if atom.chaintype in nonRNAbound:
			for otheratom in PDBmulti_included2:
				if (atom.atomtype == otheratom.atomtype and
				   atom.residuenum == otheratom.residuenum and
				   atom.basetype == otheratom.basetype):
					PDBmulti_included2.append(atom)
					break


	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (8, 4)})
	fig = plt.figure()

	# determine the dataset number i to be chosen:
	i = datasetnum


	# determine mean density for RNA bound and RNA unbound TRAP
	meandens_nonRNAbound = []
	meandens_RNAbound = []
	for atom in PDBmulti_included2:
		if atom.chaintype in nonRNAbound:
			meandens_nonRNAbound.append(atom.meandensity[i])
		elif atom.chaintype in RNAbound:
			meandens_RNAbound.append(atom.meandensity[i])
		else:
			pass

	data = [meandens_nonRNAbound,meandens_RNAbound]

	# Create an axes instance
	plt.subplot(2, 2, 1)
	sns.violinplot(data, names=["nonRNA", "RNAbound"], color="muted", lw=2)

	# ## Remove top axes and right axes ticks
	# ax_sub.get_xaxis().tick_bottom()
	# ax_sub.get_yaxis().tick_left()

	## Custom x-label,y-label             
	plt.xlabel('TRAP components', fontsize=18)
	plt.ylabel('Mean density change', fontsize=12)



	# determine min density change for RNA bound and RNA unbound TRAP
	meandens_nonRNAbound = []
	meandens_RNAbound = []
	for atom in PDBmulti_included2:
		if atom.chaintype in nonRNAbound:
			meandens_nonRNAbound.append(atom.mindensity[i])
		elif atom.chaintype in RNAbound:
			meandens_RNAbound.append(atom.mindensity[i])
		else:
			pass

	data = [meandens_nonRNAbound,meandens_RNAbound]

	# Create an axes instance
	plt.subplot(2, 2, 2)
	sns.violinplot(data, names=["nonRNA", "RNAbound"], color="deep", lw=2)

	# ## Remove top axes and right axes ticks
	# ax_sub.get_xaxis().tick_bottom()
	# ax_sub.get_yaxis().tick_left()

	## Custom x-label,y-label            
	plt.xlabel('TRAP components', fontsize=18)
	plt.ylabel('min density change', fontsize=12)



	# determine max density change for RNA bound and RNA unbound TRAP
	meandens_nonRNAbound = []
	meandens_RNAbound = []
	for atom in PDBmulti_included2:
		if atom.chaintype in nonRNAbound:
			meandens_nonRNAbound.append(atom.maxdensity[i])
		elif atom.chaintype in RNAbound:
			meandens_RNAbound.append(atom.maxdensity[i])
		else:
			pass

	data = [meandens_nonRNAbound,meandens_RNAbound]

	# Create an axes instance
	plt.subplot(2, 2, 3)
	sns.violinplot(data, names=["nonRNA", "RNAbound"], color="pastel", lw=2)

	# ## Remove top axes and right axes ticks
	# ax_sub.get_xaxis().tick_bottom()
	# ax_sub.get_yaxis().tick_left()

	## Custom x-label,y-label            
	plt.xlabel('TRAP components', fontsize=18)
	plt.ylabel('max density change', fontsize=12)



	# determine Bfactor for RNA bound and RNA unbound TRAP
	meandens_nonRNAbound = []
	meandens_RNAbound = []
	for atom in PDBmulti_included2:
		if atom.chaintype in nonRNAbound:
			meandens_nonRNAbound.append(atom.Bfactor[i])
		elif atom.chaintype in RNAbound:
			meandens_RNAbound.append(atom.Bfactor[i])
		else:
			pass

	data = [meandens_nonRNAbound,meandens_RNAbound]

	# Create an axes instance
	plt.subplot(2, 2, 4)
	sns.violinplot(data, names=["nonRNA", "RNAbound"], color="bright", lw=2)

	# ## Remove top axes and right axes ticks
	# ax_sub.get_xaxis().tick_bottom()
	# ax_sub.get_yaxis().tick_left()

	## Custom x-label,y-label            
	plt.xlabel('TRAP components', fontsize=18)
	plt.ylabel('Bfactor', fontsize=12)


	fig.suptitle('TRAPRNA binding/non-binding statistics',
	             fontsize=20)   
	# Save the figure
	fig.savefig('TRAPRNA_meandensVSproteintype.png',bbox_inches='tight')

	return distancelist


def densitychange_v_atomnum(where,PDBmulti):
	# function to plot mean density change and min density change 
	# as function of atom number on same plot

	sns.set(style="white", context="talk")
	f, axes = plt.subplots(2, 1, figsize=(100, 20), sharex=True)

	PDBmulti.sort(key=lambda x: x.atomnum)

	y = [atom.meandensity for atom in PDBmulti]
	x = [atom.atomnum for atom in PDBmulti]

	# add line plot for this dataset 
	axes[0].plot(x, y)
	axes[0].set_ylabel('mean change')

	y = [atom.mindensity for atom in PDBmulti]

	# add line plot for this dataset 
	axes[1].plot(x, y)
	axes[1].set_ylabel('min change')


	plt.xlabel('atom number', fontsize=18)
	f.suptitle('density change vs atom number',fontsize=20)
	sns.despine(bottom=True)
	plt.setp(f.axes)
	plt.tight_layout(h_pad=3)
	plt.subplots_adjust(top=0.95)
	f.savefig(where+'densitychange_vs_atom_number.png')



def TRAPringcompare_format(PDBmulti):
	# function to output a formatted multi-dim list with each element
	# as a sublist containing the equivalent atom for each of the 22
	# TRAP protein chains present in the structure. Also formats a separate
	# list of atoms included in every 4 nucleotide repeat in the RNA sequence

	PDBmulti.sort(key=lambda x: x.atomnum)
	# specify which chains are TRAP RNA bound, not RNA bound 
	# and also RNA
	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']
	RNA = ['W']

	# determine the corresponding atom in each chain, starting with chain A
	# only considering protein TRAP chains here
	atomsbychain = []
	for atom in PDBmulti:
		if atom.chaintype in nonRNAbound[0]:
			atomperchain = []
			atomperchain.append(atom)
			for otheratom in PDBmulti:
				if otheratom.chaintype not in nonRNAbound[0]:
					if (atom.atomtype == otheratom.atomtype and
					    atom.residuenum == otheratom.residuenum and
					    atom.basetype == otheratom.basetype):
						atomperchain.append(otheratom)
			atomsbychain.append(atomperchain)

	# only want to include atoms which are present in all refined chains
	# since it was observed that some atoms only appeared in particular 
	# chains in the refined model
	atomsbychain_present_PRO = []
	for j in range(0,len(atomsbychain)):
		if len(atomsbychain[j]) != (len(nonRNAbound) + len(RNAbound)):
			pass
		else:
			atomsbychain_present_PRO.append(atomsbychain[j])

	# next repeat the same process as above for the RNA chain. Here every 4 refined
	# RNA nucleotides are considered as 1 'chain' (to correspond with the protein notation
	# above) and later sets of 4 nucleotides are considered equivalent 'chains'
	atomsby4repeats = []
	for atom in PDBmulti:
		if atom.chaintype in RNA and atom.residuenum in [101,102,103,104]:
			atomper4repeats = [atom]
			for otheratom in PDBmulti:
				if otheratom.chaintype in RNA and otheratom.residuenum not in [101,102,103,104]:
					if (atom.atomtype == otheratom.atomtype and
					    atom.basetype == otheratom.basetype and
					    atom.residuenum%5 == otheratom.residuenum%5):
						atomper4repeats.append(otheratom)
			atomsby4repeats.append(atomper4repeats)

	# only want to include atoms which are present in all refined 4 nucleotide
	# repeats (so 11 copies of each!)
	atomsby4repeats_present_RNA = []
	for j in range(0,len(atomsby4repeats)):
		if len(atomsby4repeats[j]) != 11:
			pass
		else:
			atomsby4repeats_present_RNA.append(atomsby4repeats[j])

	return atomsbychain_present_PRO,atomsby4repeats_present_RNA



def TRAPringcompare_plot(atomsbychain_present,atomtype,resitype,resinum,densmet):
	# following on from formatting function above, plot density metric for corresponding 
	# atom found in all chains against chain number, with a separate subplot for each
	# dataset number

	# determine the desired atom by atom type, residue type and residue number
	counter = -1
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype and
		   atom[0].basetype == resitype and
		   atom[0].residuenum == resinum):
			break

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	x = range(1,len(atomsbychain_present[0])+1)

	for d_num in range(1,10):
		ax = plt.subplot(3, 3, d_num)
		if densmet == 'mindensitychange':
			y = [atom.mindensity[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'meandensitychange':
			y = [atom.meandensity[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'maxdensitychange':
			y = [atom.maxdensity[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'Bfactor':
			y = [atom.Bfactor[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'Bfactorchange':
			y = [atom.Bfactorchange[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'Bdamage':
			y = [atom.bdam[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'bdamagechange':
			y = [atom.bdamchange[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'stdchange':
			y = [atom.stddensity[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'Datm':
			y = [(float(atom.meandensity[d_num-1])/atom.stddensity[d_num-1]) for atom in atomsbychain_present[counter]]

		else:
			print 'unrecognised density metric\n---> terminating script'
			sys.exit()

		ax.plot(x, y,linewidth=3)
		plt.axvline(11.5, color='#d64d4d', linestyle='dashed', linewidth=3)
		ax.text(0.05, 0.90, '# '+str(d_num), color='black',transform=ax.transAxes,fontsize=18)
	
	f.text(0.5, 0.04, 'chain number', ha='center',fontsize=18)
	f.text(0.04, 0.5, str(densmet)+' density loss', va='center', rotation='vertical',fontsize=18)

	f.suptitle(str(densmet)+' vs chain number: '+\
		       str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		       str(atomsbychain_present[counter][0].residuenum) + ' ' +\
		       str(atomsbychain_present[counter][0].basetype),fontsize=24)
	plt.setp(f.axes)

	f.savefig('%s_vs_chain_number_%s_%s_%s.png' %(str(densmet),str(atomsbychain_present[counter][0].atomtype),
		       										        str(atomsbychain_present[counter][0].residuenum),
		       											    str(atomsbychain_present[counter][0].basetype)))




def TRAPringcompare_vioplot(atomsbychain_present,atomtype,resitype,resinum,densmet):
	# following on from formatting function above, plots a violin plot for density or 
	# Bfactor/Bdamage change for corresponding atom found in all chains against chain 
	# number, with a separate subplot for each dataset number. Separate violin for 
	# RNA-bound and not bound.

	# determine the desired atom by atom type, residue type and residue number
	counter = -1
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype and
		   atom[0].basetype == resitype and
		   atom[0].residuenum == resinum):
			break

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	for d_num in range(1,10):
		ax = plt.subplot(3, 3, d_num)
		if densmet == 'mindensitychange':
			y = [atom.mindensity[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'meandensitychange':
			y = [atom.meandensity[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'maxdensitychange':
			y = [atom.maxdensity[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'Bfactor':
			y = [atom.Bfactor[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'Bfactorchange':
			y = [atom.Bfactorchange[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'Bdamage':
			y = [atom.bdam[d_num-1] for atom in atomsbychain_present[counter]]
		elif densmet == 'bdamagechange':
			y = [atom.bdamchange[d_num-1] for atom in atomsbychain_present[counter]]
		else:
			print 'unrecognised density metric\n---> terminating script'
			sys.exit()

		sns.violinplot([y[0:11],y[11:]], names=["nonRNA", "RNAbound"], color="deep", lw=2)

		ax.text(0.05, 0.90, '# '+str(d_num), color='black',transform=ax.transAxes,fontsize=18)
	
	f.text(0.5, 0.04, 'Chain type', ha='center',fontsize=18)
	f.text(0.04, 0.5, str(densmet)+' density loss', va='center', rotation='vertical',fontsize=18)

	f.suptitle(str(densmet)+' vs chain number: '+\
		       str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		       str(atomsbychain_present[counter][0].residuenum) + ' ' +\
		       str(atomsbychain_present[counter][0].basetype),fontsize=24)
	plt.setp(f.axes)

	f.savefig('%s_vs_chain_number_%s_%s_%s_violinplot.png' %(str(densmet),str(atomsbychain_present[counter][0].atomtype),
		       										        str(atomsbychain_present[counter][0].residuenum),
		       											    str(atomsbychain_present[counter][0].basetype)))




def numsurroundatoms_perchainplot(atomsbychain_present,atomtype,resitype,resinum):
	# following on from formatting function above, plots a plot for number of surrounding
	# atoms calculated for corresponding atom found in all chains against chain 
	# number

	# determine the desired atom by atom type, residue type and residue number
	counter = -1
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype and
		   atom[0].basetype == resitype and
		   atom[0].residuenum == resinum):
			break

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	x = range(1,len(atomsbychain_present[0])+1)
	y = [atom.numsurroundatoms for atom in atomsbychain_present[counter]]


	plt.plot(x, y,linewidth=3)
	plt.axvline(11.5, color='#d64d4d', linestyle='dashed', linewidth=3)

	
	f.text(0.5, 0.04, 'Chain number', ha='center',fontsize=18)
	f.text(0.04, 0.5, '# surrounding atoms', va='center', rotation='vertical',fontsize=18)

	f.suptitle('Num neighbours vs chain number: '+\
		       str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		       str(atomsbychain_present[counter][0].residuenum) + ' ' +\
		       str(atomsbychain_present[counter][0].basetype),fontsize=24)
	plt.setp(f.axes)

	f.savefig('numsurroundatoms_vs_chain_number_%s_%s_%s_plot.png' %(str(atomsbychain_present[counter][0].atomtype),
		       										        str(atomsbychain_present[counter][0].residuenum),
		       											    str(atomsbychain_present[counter][0].basetype)))




def TRAPringcompare_metricVnumsurroundatms(atomsbychain_present,atomtype,resitype,resinum,densmet):
	# following on from formatting function above, scatter plot atomic or density 
	# metric for corresponding atom found in all chains against number of surrounding
	# atoms, with a separate subplot for each dataset number

	# determine the desired atom by atom type, residue type and residue number
	counter = -1
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype and
		   atom[0].basetype == resitype and
		   atom[0].residuenum == resinum):
			break

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	x = [atom.numsurroundatoms for atom in atomsbychain_present[counter][0:11]]

	# NOTE! EDITTED FOR TEST TO ONLY USE UNBOUND RING CHAINS HERE
	for d_num in range(1,10):
		ax = plt.subplot(3, 3, d_num)
		if densmet == 'mindensitychange':
			y = [atom.mindensity[d_num-1] for atom in atomsbychain_present[counter][0:11]]
		elif densmet == 'meandensitychange':
			y = [atom.meandensity[d_num-1] for atom in atomsbychain_present[counter][0:11]]
		elif densmet == 'maxdensitychange':
			y = [atom.maxdensity[d_num-1] for atom in atomsbychain_present[counter][0:11]]
		elif densmet == 'Bfactor':
			y = [atom.Bfactor[d_num-1] for atom in atomsbychain_present[counter][0:11]]
		elif densmet == 'Bfactorchange':
			y = [atom.Bfactorchange[d_num-1] for atom in atomsbychain_present[counter][0:11]]
		elif densmet == 'Bdamage':
			y = [atom.bdam[d_num-1] for atom in atomsbychain_present[counter][0:11]]
		elif densmet == 'bdamagechange':
			y = [atom.bdamchange[d_num-1] for atom in atomsbychain_present[counter][0:11]]
		else:
			print 'unrecognised density metric\n---> terminating script'
			sys.exit()

		plt.scatter(x, y, s=30, c='#d64d4d', alpha=0.5)
		ax.text(0.05, 0.90, '# '+str(d_num), color='black',transform=ax.transAxes,fontsize=18)
	
	f.text(0.5, 0.04, '# surrounding atoms', ha='center',fontsize=18)
	f.text(0.04, 0.5, str(densmet), va='center', rotation='vertical',fontsize=18)

	f.suptitle(str(densmet)+' vs # surrounding atoms: '+\
		       str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		       str(atomsbychain_present[counter][0].residuenum) + ' ' +\
		       str(atomsbychain_present[counter][0].basetype),fontsize=24)
	plt.setp(f.axes)

	f.savefig('%s_vs_numsurroundatms_%s_%s_%s.png' %(str(densmet),str(atomsbychain_present[counter][0].atomtype),
		       										        str(atomsbychain_present[counter][0].residuenum),
		       											    str(atomsbychain_present[counter][0].basetype)))




def densitychange_v_dataset(atomsbychain_present,atomtype,resitype,
							resinum,densmet):
	# function to plot density change as function of dataset number
	# for a specific atom in the structure, with all the inclusion
	# of equivalent atoms in other chains also (in the case of the 
	# RNA this will include all symmetrically related atoms, by 
	# 11-fold symmetry)

	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']
	# determine the desired atom by atom type, residue type and residue number.
	# for RNA only the first 4 nucleotides are chosen, so the residue number
	# must be selected out of 101,102,103,104
	counter = -1
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype and
		   atom[0].basetype == resitype and
		   atom[0].residuenum == resinum):
			break

	sns.set(style="white", context="talk")
	f, axes = plt.subplots(1, 1, figsize=(12, 12), sharex=True)

	# define x range here (damage set numbers)
	x = range(2,len(atom[0].meandensity)+2)

	# colours for plot defined here
	jet = cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(atom))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

	for i in range(0,len(atom)):
		if densmet == 'mindensitychange':
			y = atom[i].mindensity
		elif densmet == 'meandensitychange':
			y = atom[i].meandensity
		elif densmet == 'maxdensitychange':
			y = atom[i].maxdensity
		elif densmet == 'Bfactor':
			y = atom[i].Bfactor
		elif densmet == 'Bfactorchange':
			y = atom[i].Bfactorchange
		elif densmet == 'Bdamage':
			y = atom[i].bdam
		elif densmet == 'Bdamagechange':
			y = atom[i].bdamchange
		else:
			print 'unrecognised density metric\n---> terminating script'
			sys.exit()

		# set colour of current line here
		colorVal = scalarMap.to_rgba(i)

		# if protein atom label by chain, if RNA atom label by nucleotide number
		if atom[0].chaintype != 'W':
			if atom[i].chaintype in nonRNAbound:
				plt.plot(x,y,color='#99ccff',label='chain: '+str(atom[i].chaintype))
			elif atom[i].chaintype in RNAbound:
				plt.plot(x,y,color='#f47835',label='chain: '+str(atom[i].chaintype))
		else:
			plt.plot(x,y,color=colorVal,label='nucleotide: '+str(atom[i].residuenum))

	plt.xlabel('damage set', fontsize=18)
	plt.ylabel(str(densmet), fontsize=18)
	
	if atom[0].chaintype != 'W':
		f.suptitle(str(densmet)+' vs damage set: '+\
		  	       str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		    	   str(atomsbychain_present[counter][0].residuenum) + ' ' +\
		     	   str(atomsbychain_present[counter][0].basetype),fontsize=24)
	else:
		f.suptitle(str(densmet)+' vs damage set: '+\
		  	       str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		    	   str(atomsbychain_present[counter][0].chaintype) + ' ' +\
		     	   str(atomsbychain_present[counter][0].basetype),fontsize=24)
	plt.setp(f.axes)
	plt.subplots_adjust(top=0.90)

	# Shrink current axis by 20%
	box = axes.get_position()
	axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	# Put a legend to the right of the current axis
	axes.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=18)

	f.savefig('%s_vs_damageset_%s_%s_%s.png' %(str(densmet),str(atomsbychain_present[counter][0].atomtype),
		       										        str(atomsbychain_present[counter][0].residuenum),
		       											    str(atomsbychain_present[counter][0].basetype)))



def densitychange_v_dataset_errorbars(atomsbychain_present,atomtype,resitype,
							resinum,densmet):
	# function to plot density change as function of dataset number
	# for a specific atom in the structure, with the mean value over 
	# all protein chains plotted, along with error bars for the 22 
	# equivalent atoms present

	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']

	# determine the desired atom by atom type, residue type and residue number.
	# for RNA only the first 4 nucleotides are chosen, so the residue number
	# must be selected out of 101,102,103,104
	counter = -1
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype and
		   atom[0].basetype == resitype and
		   atom[0].residuenum == resinum):
			break

	# only can handle protein here so check it is a protein atom that 
	# has been selected
	if atom[0].chaintype not in (nonRNAbound + RNAbound):
		print 'Script cannot handle anything but protein atoms'
		print '---> terminating script'
		sys.exit()

	sns.set(style="white", context="talk")
	f, axes = plt.subplots(1, 1, figsize=(12, 12), sharex=True)

	# define x range here (damage set numbers)
	x = range(2,len(atom[0].meandensity)+2)

	y_mean_unboundTRAP = []
	y_mean_boundTRAP = []
	y_std_unboundTRAP = []
	y_std_boundTRAP = []
	for j in range(len(x)):
		y_equivatms = []
		for i in range(0,len(atom)):
			if densmet == 'mindensitychange':
				y = atom[i].mindensity[j]
			elif densmet == 'meandensitychange':
				y = atom[i].meandensity[j]
			elif densmet == 'maxdensitychange':
				y = atom[i].maxdensity[j]
			elif densmet == 'Bfactor':
				y = atom[i].Bfactor[j]
			elif densmet == 'Bfactorchange':
				y = atom[i].Bfactorchange[j]
			elif densmet == 'Bdamage':
				y = atom[i].bdam[j]
			elif densmet == 'bdamagechange':
				y = atom[i].bdamchange[j]
			else:
				print 'unrecognised density metric\n---> terminating script'
				sys.exit()
			y_equivatms.append(y)

		y_mean_unboundTRAP.append(np.mean(y_equivatms[0:11]))
		y_mean_boundTRAP.append(np.mean(y_equivatms[11:22]))
		y_std_unboundTRAP.append(np.std(y_equivatms[0:11]))
		y_std_boundTRAP.append(np.std(y_equivatms[11:22]))

	plt.errorbar(x,y_mean_unboundTRAP,yerr=y_std_unboundTRAP, fmt='-o',color='#99ccff',label='unbound')
	plt.errorbar(x,y_mean_boundTRAP,yerr=y_std_boundTRAP,fmt='-o',color='#f47835',label='bound')

	plt.xlabel('damage set', fontsize=18)
	plt.ylabel(str(densmet), fontsize=18)
	
	f.suptitle(str(densmet)+' vs damage set: '+\
		  	   str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		       str(atomsbychain_present[counter][0].residuenum) + ' ' +\
		       str(atomsbychain_present[counter][0].basetype),fontsize=24)

	plt.setp(f.axes)
	plt.subplots_adjust(top=0.90)

	# # Shrink current axis by 20%
	# box = axes.get_position()
	# axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	# Put a legend to the right of the current axis
	# axes.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=18)
	axes.legend(loc='best',fontsize=18)
	f.savefig('%s_vs_damageset_%s_%s_%s_errorbars.png' %(str(densmet),str(atomsbychain_present[counter][0].atomtype),
		       										        str(atomsbychain_present[counter][0].residuenum),
		       											    str(atomsbychain_present[counter][0].basetype)))






def densitychange_v_dataset_errorbars_4subplots(atomsbychain_present,atomtype,resitype,
							resinum):
	# function to plot density change as function of dataset number
	# for a specific atom in the structure, with the mean value over 
	# all protein chains plotted, along with error bars for the 22 
	# equivalent atoms present

	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']

	# determine the desired atom by atom type, residue type and residue number.
	# for RNA only the first 4 nucleotides are chosen, so the residue number
	# must be selected out of 101,102,103,104
	counter = -1
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype and
		   atom[0].basetype == resitype and
		   atom[0].residuenum == resinum):
			break

	# only can handle protein here so check it is a protein atom that 
	# has been selected
	if atom[0].chaintype not in (nonRNAbound + RNAbound):
		print 'Script cannot handle anything but protein atoms'
		print '---> terminating script'
		sys.exit()

	sns.set(style="white", context="talk")
	f = plt.figure()
	# plt.subplots_adjust(hspace=0.4,)

	# f = plt.subplots_adjust(hspace=0.4)

	# define x range here (damage set numbers)
	x = range(2,len(atom[0].meandensity)+2)[0:9]

	# mean density axis here
	y_mean_unboundTRAP = []
	y_mean_boundTRAP = []
	y_std_unboundTRAP = []
	y_std_boundTRAP = []
	for j in range(len(x)):
		y_equivatms = []
		for i in range(0,len(atom)):
			y_equivatms.append(atom[i].meandensity[j])
		y_mean_unboundTRAP.append(np.mean(y_equivatms[0:11]))
		y_mean_boundTRAP.append(np.mean(y_equivatms[11:22]))
		y_std_unboundTRAP.append(np.std(y_equivatms[0:11]))
		y_std_boundTRAP.append(np.std(y_equivatms[11:22]))

	ax = plt.subplot(221)
	plt.errorbar(x,y_mean_unboundTRAP,yerr=y_std_unboundTRAP, fmt='-o',capthick=2,color='#99ccff',label='unbound')
	plt.errorbar(x,y_mean_boundTRAP,yerr=y_std_boundTRAP,fmt='-o',capthick=2,color='#f47835',label='bound')
	ax.set_xlim([1, 11])
	ax.legend(loc='best')

	plt.xlabel('Damage set')
	plt.ylabel('mean density change')

	# min density axis here
	y_mean_unboundTRAP = []
	y_mean_boundTRAP = []
	y_std_unboundTRAP = []
	y_std_boundTRAP = []
	for j in range(len(x)):
		y_equivatms = []
		for i in range(0,len(atom)):
			y_equivatms.append(atom[i].mindensity[j])
		y_mean_unboundTRAP.append(np.mean(y_equivatms[0:11]))
		y_mean_boundTRAP.append(np.mean(y_equivatms[11:22]))
		y_std_unboundTRAP.append(np.std(y_equivatms[0:11]))
		y_std_boundTRAP.append(np.std(y_equivatms[11:22]))

	ax = plt.subplot(222)
	plt.errorbar(x,y_mean_unboundTRAP,yerr=y_std_unboundTRAP, fmt='-o',capthick=2,color='#99ccff',label='unbound')
	plt.errorbar(x,y_mean_boundTRAP,yerr=y_std_boundTRAP,fmt='-o',capthick=2,color='#f47835',label='bound')
	ax.set_xlim([1, 11])
	ax.legend(loc='best')

	plt.xlabel('Damage set')
	plt.ylabel('min density change')

	# Bdamage change axis here
	y_mean_unboundTRAP = []
	y_mean_boundTRAP = []
	y_std_unboundTRAP = []
	y_std_boundTRAP = []
	for j in range(len(x)):
		y_equivatms = []
		for i in range(0,len(atom)):
			y_equivatms.append(atom[i].bdamchange[j])
		y_mean_unboundTRAP.append(np.mean(y_equivatms[0:11]))
		y_mean_boundTRAP.append(np.mean(y_equivatms[11:22]))
		y_std_unboundTRAP.append(np.std(y_equivatms[0:11]))
		y_std_boundTRAP.append(np.std(y_equivatms[11:22]))

	ax = plt.subplot(223)
	plt.errorbar(x,y_mean_unboundTRAP,yerr=y_std_unboundTRAP, fmt='-o',capthick=2,color='#99ccff',label='unbound')
	plt.errorbar(x,y_mean_boundTRAP,yerr=y_std_boundTRAP,fmt='-o',capthick=2,color='#f47835',label='bound')
	ax.set_xlim([1, 11])
	ax.legend(loc='best')

	plt.xlabel('Damage set')
	plt.ylabel('Bdamage change')

	# Bfactor axis here
	y_mean_unboundTRAP = []
	y_mean_boundTRAP = []
	y_std_unboundTRAP = []
	y_std_boundTRAP = []
	for j in range(len(x)):
		y_equivatms = []
		for i in range(0,len(atom)):
			y_equivatms.append(atom[i].Bfactorchange[j])
		y_mean_unboundTRAP.append(np.mean(y_equivatms[0:11]))
		y_mean_boundTRAP.append(np.mean(y_equivatms[11:22]))
		y_std_unboundTRAP.append(np.std(y_equivatms[0:11]))
		y_std_boundTRAP.append(np.std(y_equivatms[11:22]))

	ax = plt.subplot(224)
	plt.errorbar(x,y_mean_unboundTRAP,yerr=y_std_unboundTRAP, fmt='-o',capthick=2,color='#99ccff',label='unbound')
	plt.errorbar(x,y_mean_boundTRAP,yerr=y_std_boundTRAP,fmt='-o',capthick=2,color='#f47835',label='bound')
	ax.set_xlim([1, 11])
	ax.legend(loc='best')

	plt.xlabel('Damage set')
	plt.ylabel('Bfactor change')

	plt.subplots_adjust(top=0.90)
	f.subplots_adjust(hspace=.3)
	f.subplots_adjust(wspace=.4)

	f.suptitle('damage metrics vs damage set: '+\
		  	   str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		       str(atomsbychain_present[counter][0].residuenum) + ' ' +\
		       str(atomsbychain_present[counter][0].basetype),fontsize=20)

	f.savefig('4metrics_vs_damageset_%s_%s_%s_errorbars_4subplots.png' %(str(atomsbychain_present[counter][0].atomtype),
		       										        str(atomsbychain_present[counter][0].residuenum),
		       											    str(atomsbychain_present[counter][0].basetype)))



def student_t_test(atomsbychain_present,atomtype,resitype,resinum,densmet):
	# function to calculate student t test between two sets of observations (ie unbound 
	# and bound equivalent atoms) for each damage set individually
	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']
	# determine the desired atom by atom type, residue type and residue number.
	# for RNA only the first 4 nucleotides are chosen, so the residue number
	# must be selected out of 101,102,103,104
	counter = -1
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype and
		   atom[0].basetype == resitype and
		   atom[0].residuenum == resinum):
			break

	# for each dataset number, calculate t test statistic
	t_perdataset,p_perdataset = [],[]
	for d_num in range(1,10):
		vals = [atom.meandensity[d_num-1] for atom in atomsbychain_present[counter]]
		vals_unbound = vals[0:11]
		vals_bound = vals[11:22]
		# calculates t test statistic and two-tailed p-value
		t,p = stats.ttest_rel(vals_unbound,vals_bound)
		t_perdataset.append(t)
		p_perdataset.append(p)

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f, ax1 = plt.subplots()

	# x axis of plot is number of damage sets 
	x = range(2,len(atomsbychain_present[counter][0].meandensity)+2)

	ax1.plot(x, t_perdataset,linewidth=3)
	ax1.set_xlabel('dataset')
	ax1.set_ylabel('t test statistic', color='b')
	for tl in ax1.get_yticklabels():
		tl.set_color('b')

	ax2 = ax1.twinx()
	ax2.plot(x, p_perdataset, 'r')
	ax2.set_ylabel('2-tailed p-value', color='r')
	for tl in ax2.get_yticklabels():
		tl.set_color('r')

	f.suptitle('t-test/2-tailed p-value vs dataset number: '+\
		       str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		       str(atomsbychain_present[counter][0].residuenum) + ' ' +\
		       str(atomsbychain_present[counter][0].basetype),fontsize=24)
	plt.setp(f.axes)

	f.savefig('t_test_vs_datasetnum_%s_%s_%s_plot.png' %(str(atomsbychain_present[counter][0].atomtype),
		       										        str(atomsbychain_present[counter][0].residuenum),
		       											    str(atomsbychain_present[counter][0].basetype)))



def GLUandASPs(PDBmulti,resinum):
	# function finds all the GLU and ASP residues in the structure, and groups
	# by equivalent atoms (using 11-fold ring symmetry). Here a specific residue
	# number is specified and PCA restricted to this residue type. The PCA input 
	# dimensions are mean density and min density change (one dataset = one 
	# separate dimension)

	keptatoms = []
	# find the GLU and ASP carboxyl oxygens here
	for atom in PDBmulti:
		if (atom.basetype in ('GLU','ASP') and
			atom.atomtype in ('OD1','OD2','OE1','OE2') and
			atom.residuenum in [int(resinum)]):
			keptatoms.append(atom)


	# linearly fit to each density metric (over the number of damage sets present)
	x = np.array(range(len(keptatoms[0].meandensity)))
	dataperatom = []
	for atom in keptatoms:
		print str(atom.atomtype) + str(atom.basetype) + str(atom.chaintype) + str(atom.residuenum)

		atom.meandensity + atom.mindensity
		dataperatom.append(atom.meandensity + atom.mindensity)

	# perform PCA on dataperatom data specified above
	data = np.array(dataperatom)
	result = PCA(data) 

	x = []
	y = []
	z = []
	for item in result.Y:
	 x.append(item[0])
	 y.append(item[1])
	 z.append(item[2])

	plt.close('all') # close all latent plotting windows
	fig = plt.figure() # Make a plotting figure
	fig, axes = plt.subplots(1, 1, figsize=(16, 16))

	pltData = [x,y,z] 
	colors = ['r']*22+['g']*22
	plt.scatter(pltData[0], pltData[1],c=colors,s=50)
	 
	# label the axes 
	axes.set_xlabel("component 1",fontsize = 24) 
	axes.set_ylabel("component 2",fontsize = 24)
	axes.set_title("PCA analysis: %s %s carboxyl oxygens\n red: unbound; green: bound" %(str(keptatoms[0].basetype),str(keptatoms[0].residuenum)),
				fontsize = 28)


	fig.savefig("PCAanalysis:%s%s_carboxyl_oxys_red-unbound_green-bound.png" %(str(keptatoms[0].basetype),str(keptatoms[0].residuenum)))


def GLUandASPs_all(PDBmulti,colorscheme,boundmarker,unboundmarker):
	# function finds all the GLU and ASP residues in the structure, and groups
	# by equivalent atoms (using 11-fold ring symmetry). Function then calculates PCA 
	# over all GLU and ASP residues found. The PCA input dimensions are mean density 
	# and min density change (one dataset = one separate dimension)

	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']

	keptatoms = []
	# find the GLU and ASP carboxyl oxygens here
	for atom in PDBmulti:
		if (atom.basetype in ('GLU','ASP') and
			atom.atomtype in ('OD1','OD2','OE1','OE2')):
			keptatoms.append(atom)


	# linearly fit to each density metric (over the number of damage sets present)
	x = np.array(range(len(keptatoms[0].meandensity)))
	dataperatom = []
	for atom in keptatoms:
		dataperatom.append(atom.meandensity + atom.mindensity)# + atom.bdamchange + atom.Bfactorchange)

	# perform PCA on dataperatom data specified above
	data = np.array(dataperatom)
	result = PCA(data) 

	x = []
	y = []
	for item in result.Y:
	 x.append(item[0])
	 y.append(item[1])

	plt.close('all') # close all latent plotting windows
	fig = plt.figure() # Make a plotting figure
	fig, axes = plt.subplots(1, 1, figsize=(16, 16))

	pltData = [x,y] 

	# colours for plot defined here
	jet = cm = plt.get_cmap(str(colorscheme)) 
	cNorm  = colors.Normalize(vmin=0, vmax=9)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

	# determine list of unique residue numbers present here
	uniq_resnums = []
	for atom in keptatoms:
		if atom.residuenum not in uniq_resnums:
			uniq_resnums.append(atom.residuenum)

	# make a list of seen residue numbers so only included once in legend
	seen_resnums_bound = []
	seen_resnums_unbound = []

	counter = -1
	for atom in keptatoms:
		counter += 1
		ind = uniq_resnums.index(atom.residuenum)
		color = scalarMap.to_rgba(ind)

		if atom.chaintype in RNAbound:
			if atom.residuenum in seen_resnums_bound:
				plt.scatter(pltData[0][counter], pltData[1][counter],
					   		c=color,marker=str(boundmarker),s=50)
			else:
				plt.scatter(pltData[0][counter], pltData[1][counter],
					        c=color,marker=str(boundmarker),s=50,label=str(atom.basetype)+str(atom.residuenum)+' bound')
				seen_resnums_bound.append(atom.residuenum)

		else:
			if atom.residuenum in seen_resnums_unbound:
				plt.scatter(pltData[0][counter], pltData[1][counter],
					   		c=color,marker=str(unboundmarker),s=50)
			else:
				plt.scatter(pltData[0][counter], pltData[1][counter],
					        c=color,marker=str(unboundmarker),s=50,label=str(atom.basetype)+str(atom.residuenum)+' unbound')
				seen_resnums_unbound.append(atom.residuenum)


	# label the axes 
	axes.set_xlabel("component 1",fontsize = 24) 
	axes.set_ylabel("component 2",fontsize = 24)
	axes.set_title("PCA analysis: ALL carboxyl oxygens",
				   fontsize = 28)
 
	# Shrink current axis by 20%
	box = axes.get_position()
	axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	# Put a legend to the right of the current axis
	axes.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=18)

	fig.savefig("PCAanalysis:ALL_carboxyl_oxys.png")




def disttoRNA(PDBmulti):
	# this function calculates them minimum distance between the RNA and
	# GLU and ASP carboxyl oxygens. Outputs as a list of distances with 
	# same order as atom order in structure. This function then takes 
	# the list of distances of GLU and ASP carboxyl oxygens from the RNA 
	# and creates a scatter plot of distance from RNA against density loss 

	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound = ['L','M','N','O','P','Q','R','S','T','U','V']
	RNA = ['W']

	keptatoms = []
	# find the GLU and ASP carboxyl oxygens here
	for atom in PDBmulti:
		if (atom.basetype in ('GLU','ASP') and
			atom.atomtype in ('OD1','OD2','OE1','OE2')):
			keptatoms.append(atom)

	# calculating min distance from RNA for each protein atom in keptatoms
	distancelist = []
	for atom in keptatoms:
		distancelist4atm = []
		for otheratom in PDBmulti:
			if otheratom.chaintype in RNA:
				distance = np.sqrt(np.square(atom.X_coord - otheratom.X_coord) +
	                           	   np.square(atom.Y_coord - otheratom.Y_coord) + 
	                               np.square(atom.Z_coord - otheratom.Z_coord))
				distancelist4atm.append(distance)
		distancelist.append(min(distancelist4atm))

	# # # The code here plots a graph of distance from RNA for atoms in keptatoms.
	# # Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()
	# # x axis of plot is number of atoms in keptatoms
	# x = range(len(keptatoms))
	# ax.plot(x, distancelist,linewidth=3)
	# ax.set_xlabel('counter',fontsize = 20)
	# ax.set_ylabel('Distance from RNA',fontsize = 20)
	# ax.set_title("Distance from RNA: Carboxyl Oxygens",fontsize = 24)
	# f.savefig("Distance_from_RNA_Carboxyl_Oxygens.png")

	# The code here plots a scatter plot of distance from RNA against min density
	# change, for atoms in keptatoms
	for d_num in range(len(keptatoms[0].mindensity)):
		ax = plt.subplot(3, 3, d_num+1)
		y = [atom.mindensity[d_num] for atom in keptatoms]
		ax.scatter(distancelist,y,c='g')
		ax.text(0.05, 0.90, '# '+str(d_num+2), color='black',transform=ax.transAxes,fontsize=18)
	
	f.text(0.5, 0.04, 'Distance from RNA (Angstrom)', ha='center',fontsize=22)
	f.text(0.04, 0.5, 'Min density change (e/cubicA)', va='center', rotation='vertical',fontsize=22)

	f.suptitle("Distance from RNA Vs Min density change:\n Carboxyl Oxygens",fontsize=24)
	plt.setp(f.axes)
	f.savefig("RNAdistVSmindens_Carboxyl_Oxygens.png")


def disttoRNA_1atmtype(PDBmulti,atomsbychain_present,atomtype,resitype,resinum):

	# determine the desired atom by atom type, residue type and residue number
	counter = -1
	for atomset in atomsbychain_present:
		counter += 1
		if (atomset[0].atomtype == atomtype and
		   atomset[0].basetype == resitype and
		   atomset[0].residuenum == resinum):
			break

	# calculating min distance from RNA for each protein atom in atomset (the 
	# located set of equivalent atoms in each chain on a specific type)
	distancelist = []
	for atom in atomset:
		distancelist4atm = []
		for otheratom in PDBmulti:
			if otheratom.chaintype == 'W':
				distance = np.sqrt(np.square(atom.X_coord - otheratom.X_coord) +
	                           	   np.square(atom.Y_coord - otheratom.Y_coord) + 
	                               np.square(atom.Z_coord - otheratom.Z_coord))
				distancelist4atm.append(distance)
		distancelist.append(min(distancelist4atm))

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	# The code here plots a scatter plot of distance from RNA against min density
	# change, for atoms in atomset
	for d_num in range(len(atomset[0].mindensity)):
		ax = plt.subplot(3, 3, d_num+1)
		y = [atom.mindensity[d_num] for atom in atomset]
		ax.scatter(distancelist,y,c='r')
		ax.text(0.05, 0.90, '# '+str(d_num+2), color='black',transform=ax.transAxes,fontsize=18)
	
	f.text(0.5, 0.04, 'Distance from RNA (Angstrom)', ha='center',fontsize=22)
	f.text(0.04, 0.5, 'Min density change (e/cubicA)', va='center', rotation='vertical',fontsize=22)

	f.suptitle("Distance from RNA Vs min density change:\n%s %s %s" %(str(atomset[0].basetype),str(atomset[0].residuenum),str(atomset[0].atomtype)),fontsize=24)
	plt.setp(f.axes)
	f.savefig("RNAdistVSmindens_%s%s%s.png" %(str(atomset[0].basetype),str(atomset[0].residuenum),str(atomset[0].atomtype)))




def numWaters_1atmtype(PDBmulti,atomsbychain_present,atomtype,resitype,resinum,threshold):
	# This function determines the number of waters within a specified threshold from a 
	# a given atom (and equivalent atoms) and plots a scatter plot to determine any 
	# possible correlation between hydration level and damage to specific residue types

	# determine the desired atom by atom type, residue type and residue number
	counter = -1
	for atomset in atomsbychain_present:
		counter += 1
		if (atomset[0].atomtype == atomtype and
		   atomset[0].basetype == resitype and
		   atomset[0].residuenum == resinum):
			break

	# calculating number of waters within a threshold distance atom the atoms chosen
	numWaters = []
	for atom in atomset:
		watercounter = 0 
		for otheratom in PDBmulti:
			if otheratom.chaintype == 'Y':
				distance = np.sqrt(np.square(atom.X_coord - otheratom.X_coord) +
	                           	   np.square(atom.Y_coord - otheratom.Y_coord) + 
	                               np.square(atom.Z_coord - otheratom.Z_coord))
				if distance < threshold:
					watercounter += 1

		numWaters.append(watercounter)

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	# The code here plots a scatter plot of distance from RNA against min density
	# change, for atoms in atomset
	for d_num in range(len(atomset[0].mindensity)):
		ax = plt.subplot(3, 3, d_num+1)
		y = [atom.mindensity[d_num] for atom in atomset]
		ax.scatter(numWaters,y,c='r')	

		slope, intercept, r_value, p_value, std_err = stats.linregress(numWaters, y)
		ax.text(0.05, 0.85, '# '+str(d_num+2)+'\nR2: '+str("{0:.2f}".format(r_value**2)),
			    color='black',transform=ax.transAxes,fontsize=18)
		x2 = range(min(numWaters),max(numWaters)+1)
		y2 = slope*np.array(x2) + np.array(intercept)
		plt.plot(x2, y2, '-',color='#d11141',linewidth=3)
	
	f.text(0.5, 0.04, '# waters', ha='center',fontsize=22)
	f.text(0.04, 0.5, 'Min density change (e/A^3)', va='center', rotation='vertical',fontsize=22)

	f.suptitle("Proximal waters (< %s Angstrom) vs min density change:\n%s %s %s" %(str(threshold),str(atomset[0].basetype),str(atomset[0].residuenum),str(atomset[0].atomtype)),fontsize=24)
	plt.setp(f.axes)
	f.savefig("waters%sAngVSmindens_%s%s%s.png" %(str(threshold),str(atomset[0].basetype),str(atomset[0].residuenum),str(atomset[0].atomtype)))




def refbarplot_variabilityperrestype(PDBmulti,atomsbychain_present,threshold):

	# determine the locations of GLU and ASP residues in the atomsbychain_present
	# list (sorted into sublists of equivalent chain atoms)
	counter = -1
	counter_list = []
	for atomset in atomsbychain_present:
		counter += 1
		if (atomset[0].atomtype in ('OD1','OD2','OE1','OE2') and
		   atomset[0].basetype in ('GLU','ASP')):
			counter_list.append(counter)

	# loop over each equivalent atom type and determine mean and std of metrics 
	# specified in loop over each atom in an equivalent atom type
	# FIRST FOR UNBOUND RING CASE:
	mean_Bfactors_unbound,mean_numWaters_unbound,mean_packdens_unbound = [],[],[]
	mean_packdens_protons_unbound,xTickMarks_unbound = [],[]
	std_Bfactors_unbound,std_numWaters_unbound = [],[]
	std_packdens_unbound,std_packdens_protons_unbound = [],[]
	for counter in counter_list:
		atomgroup = atomsbychain_present[counter][0:11]
		# find list of Bfactors for these equivalent atoms
		Bfactors = [float(atom.Bfactor[0]) for atom in atomgroup]

		# calculating number of waters within a threshold distance atom 
		# of the atoms chosen
		numWaters = []
		for atom in atomgroup:
			watercounter = 0 
			for otheratom in PDBmulti:
				if otheratom.chaintype == 'Y':
					distance = np.sqrt(np.square(atom.X_coord - otheratom.X_coord) +
		                           	   np.square(atom.Y_coord - otheratom.Y_coord) + 
		                               np.square(atom.Z_coord - otheratom.Z_coord))
					if distance < threshold:
						watercounter += 1

			numWaters.append(watercounter)

		# find list of packing density for these equivalent atoms
		packdens = [int(atom.numsurroundatoms) for atom in atomgroup]

		# find list of # proton neighbours for equivalent atoms
		packdens_protons = [int(atom.numsurroundprotons) for atom in atomgroup]

		# calculate means of above lists
		mean_Bfactors_unbound.append(np.mean(Bfactors))
		mean_numWaters_unbound.append(np.mean(numWaters))
		mean_packdens_unbound.append(np.mean(packdens))
		mean_packdens_protons_unbound.append(np.mean(packdens_protons))

		# calculate standard deviations of above lists
		std_Bfactors_unbound.append(np.std(Bfactors))
		std_numWaters_unbound.append(np.std(numWaters))
		std_packdens_unbound.append(np.std(packdens))
		std_packdens_protons_unbound.append(np.std(packdens_protons))

		# determine the reference name of the current residue type,
		# and append this to list of x-axis labels
		resitype = str(atomgroup[0].basetype)+str(atomgroup[0].residuenum)
		xTickMarks_unbound.append(resitype)

	# NEXT FOR BOUND RING CASE:
	mean_Bfactors_bound,mean_numWaters_bound,mean_packdens_bound = [],[],[]
	mean_packdens_protons_bound,xTickMarks_bound = [],[]
	std_Bfactors_bound,std_numWaters_bound = [],[]
	std_packdens_bound,std_packdens_protons_bound = [],[]
	for counter in counter_list:
		atomgroup = atomsbychain_present[counter][11:22]
		# find list of Bfactors for these equivalent atoms
		Bfactors = [float(atom.Bfactor[0]) for atom in atomgroup]

		# calculating number of waters within a threshold distance atom 
		# of the atoms chosen
		numWaters = []
		for atom in atomgroup:
			watercounter = 0 
			for otheratom in PDBmulti:
				if otheratom.chaintype == 'Y':
					distance = np.sqrt(np.square(atom.X_coord - otheratom.X_coord) +
		                           	   np.square(atom.Y_coord - otheratom.Y_coord) + 
		                               np.square(atom.Z_coord - otheratom.Z_coord))
					if distance < threshold:
						watercounter += 1

			numWaters.append(watercounter)

		# find list of packing density for these equivalent atoms
		packdens = [int(atom.numsurroundatoms) for atom in atomgroup]

		# find list of # proton neighbours for equivalent atoms
		packdens_protons = [int(atom.numsurroundprotons) for atom in atomgroup]

		# calculate means of above lists
		mean_Bfactors_bound.append(np.mean(Bfactors))
		mean_numWaters_bound.append(np.mean(numWaters))
		mean_packdens_bound.append(np.mean(packdens))
		mean_packdens_protons_bound.append(np.mean(packdens_protons))

		# calculate standard deviations of above lists
		std_Bfactors_bound.append(np.std(Bfactors))
		std_numWaters_bound.append(np.std(numWaters))
		std_packdens_bound.append(np.std(packdens))
		std_packdens_protons_bound.append(np.std(packdens_protons))

		# determine the reference name of the current residue type,
		# and append this to list of x-axis labels
		resitype = str(atomgroup[0].basetype)+str(atomgroup[0].residuenum)+str(atomgroup[0].atomtype)
		xTickMarks_bound.append(resitype)

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (40, 12)})
	fig = plt.figure()

	# width of bars in plot
	width = 0.35 

	# Create an axes instance
	ax = plt.subplot(2, 2, 1)
	# x axis is discrete number of different equivalent atom groups
	x = np.arange(len(xTickMarks_bound))
	unboundplt = plt.bar(x, mean_Bfactors_unbound,width,
						 yerr=std_Bfactors_unbound,color='#6897bb',
						 error_kw=dict(elinewidth=2,ecolor='#31698a'))
	boundplt = plt.bar(x+width, mean_Bfactors_bound,width,
		               yerr=std_Bfactors_bound,color='#c0d6e4',
		               error_kw=dict(elinewidth=2,ecolor='#31698a'))

	## Custom x-label,y-label    
	ax.set_xlim(-width,len(x)+width)                    
	plt.xlabel('Residue', fontsize=18)
	plt.ylabel('Bfactor', fontsize=18)
	ax.set_xticks(x+width)
	xtickNames = ax.set_xticklabels(xTickMarks_bound)
	plt.setp(xtickNames, rotation=90, fontsize=10)
	# add a legend
	ax.legend( (unboundplt[0], boundplt[0]), ('unbound', 'bound'),loc='best')
	
	# Create an axes instance
	ax = plt.subplot(2, 2, 2)
	unboundplt = plt.bar(x, mean_numWaters_unbound,width,
					     yerr=std_numWaters_unbound,color='#6897bb',
					     error_kw=dict(elinewidth=2,ecolor='#31698a'))
	boundplt = plt.bar(x+width, mean_numWaters_bound,width,
					   yerr=std_numWaters_bound,color='#c0d6e4',
					   error_kw=dict(elinewidth=2,ecolor='#31698a'))

	## Custom x-label,y-label
	ax.set_xlim(-width,len(x)+width)                       
	plt.xlabel('Residue', fontsize=18)
	plt.ylabel('# waters', fontsize=18)
	ax.set_xticks(x+width)
	xtickNames = ax.set_xticklabels(xTickMarks_bound)
	plt.setp(xtickNames, rotation=90, fontsize=10)
	# add a legend
	ax.legend( (unboundplt[0], boundplt[0]), ('unbound', 'bound'),loc='best')
	
	# Create an axes instance
	ax = plt.subplot(2, 2, 3)
	unboundplt = plt.bar(x, mean_packdens_unbound,width,
					     yerr=std_packdens_unbound,color='#6897bb',
					     error_kw=dict(elinewidth=2,ecolor='#31698a'))
	boundplt = plt.bar(x+width,mean_packdens_bound,width,
					   yerr=std_packdens_bound,color='#c0d6e4',
					   error_kw=dict(elinewidth=2,ecolor='#31698a'))

	## Custom x-label,y-label   
	ax.set_xlim(-width,len(x)+width)                 
	plt.xlabel('Residue', fontsize=18)
	plt.ylabel('# atom neighbours', fontsize=18)
	ax.set_xticks(x+width)
	xtickNames = ax.set_xticklabels(xTickMarks_bound)
	plt.setp(xtickNames, rotation=90, fontsize=10)
	# add a legend
	ax.legend( (unboundplt[0], boundplt[0]), ('unbound', 'bound'),loc='best')
	
	# Create an axes instance
	ax = plt.subplot(2, 2, 4)
	unboundplt = plt.bar(x, mean_packdens_protons_unbound,width,
						 yerr=std_packdens_protons_unbound,color='#6897bb',
						 error_kw=dict(elinewidth=2,ecolor='#31698a'))
	boundplt = plt.bar(x+width,mean_packdens_protons_bound,width,
					   yerr=std_packdens_protons_bound,color='#c0d6e4',
					   error_kw=dict(elinewidth=2,ecolor='#31698a'))

	## Custom x-label,y-label  
	ax.set_xlim(-width,len(x)+width)           
	plt.xlabel('Residue', fontsize=18)
	plt.ylabel('# proton neighbours', fontsize=18)
	ax.set_xticks(x+width)
	xtickNames = ax.set_xticklabels(xTickMarks_bound)
	plt.setp(xtickNames, rotation=90, fontsize=10)
	# add a legend
	ax.legend( (unboundplt[0], boundplt[0]), ('unbound', 'bound'),loc='best')
	
	plt.subplots_adjust(top=0.90)
	fig.subplots_adjust(hspace=.4)

	fig.suptitle('GLU and ASP variability per chain',
	             fontsize=24)   
	# Save the figure
	fig.savefig('GLU_and_ASP_variability_per_chain.png',bbox_inches='tight')



def RNAdamage(PDBmulti,densmet,atomtype):
	# function to directly analyse the damage to RNA and any heterogeneity damage to 
	# the RNA bases around the ring.

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (12, 12)})
	fig = plt.figure()

	# find RNA atoms and group into lists corresponding to each basetype
	for d_num in range(1,10):
		ax = plt.subplot(3, 3, d_num)
		RNA_atoms_A,RNA_atoms_G,RNA_atoms_U = [],[],[]
		for atom in PDBmulti:
			# filter to only keep base or backbone atoms (specified by sideormain)

			if atom.chaintype == 'W' and atom.atomtype == atomtype:
				# first determine the correct metric
				if densmet == 'mindensitychange':
					y = atom.mindensity[d_num-1] 
				elif densmet == 'meandensitychange':
					y = atom.meandensity[d_num-1] 
				elif densmet == 'maxdensitychange':
					y = atom.maxdensity[d_num-1] 
				elif densmet == 'Bfactor':
					y = atom.Bfactor[d_num-1] 
				elif densmet == 'Bfactorchange':
					y = atom.Bfactorchange[d_num-1] 
				elif densmet == 'Bdamage':
					y = atom.bdam[d_num-1] 
				elif densmet == 'bdamagechange':
					y = atom.bdamchange[d_num-1] 
				elif densmet == 'stdchange':
					y = atom.stddensity[d_num-1] 
				else:
					print 'unrecognised metric\n---> terminating script'
					sys.exit()
				# now assign to correct nucleotide list
				if atom.basetype == 'A':
					RNA_atoms_A.append(y)
				elif atom.basetype == 'G':
					RNA_atoms_G.append(y)
				elif atom.basetype == 'U':
					RNA_atoms_U.append(y)


		sns.violinplot([RNA_atoms_A,RNA_atoms_G,RNA_atoms_U], names=["A","G","U"], color="muted", lw=2)
		ax.text(0.05, 0.90, '# '+str(d_num), color='black',transform=ax.transAxes,fontsize=18)

	fig.text(0.5, 0.04, 'RNA base type', ha='center',fontsize=18)
	fig.text(0.04, 0.5, str(densmet), va='center', rotation='vertical',fontsize=18)
	fig.subplots_adjust(hspace=.4)
	fig.subplots_adjust(wspace=.4)

	fig.suptitle("RNA base type vs %s" %(str(densmet)),fontsize = 24)
	plt.setp(fig.axes)
	fig.savefig("RNAbasetype_vs_%s.png" %(str(densmet)))



def RNAdamage_vsatomnum(PDBmulti,densmet):
	# function to directly analyse the damage to RNA and any heterogeneity damage to 
	# the RNA bases around the ring.

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (100, 100)})
	fig = plt.figure()

	# find RNA atoms and group into lists corresponding to each basetype
	for d_num in range(1,10):
		ax = plt.subplot(3, 3, d_num)
		RNA_atomvals = []
		RNA_atomnums = []
		for atom in PDBmulti:
			# filter to only keep base or backbone atoms (specified by sideormain)

			if atom.chaintype == 'W':
				# first determine the correct metric
				if densmet == 'mindensitychange':
					y = atom.mindensity[d_num-1] 
				elif densmet == 'meandensitychange':
					y = atom.meandensity[d_num-1] 
				elif densmet == 'maxdensitychange':
					y = atom.maxdensity[d_num-1] 
				elif densmet == 'Bfactor':
					y = atom.Bfactor[d_num-1] 
				elif densmet == 'Bfactorchange':
					y = atom.Bfactorchange[d_num-1] 
				elif densmet == 'Bdamage':
					y = atom.bdam[d_num-1] 
				elif densmet == 'bdamagechange':
					y = atom.bdamchange[d_num-1] 
				elif densmet == 'stdchange':
					y = atom.stddensity[d_num-1] 
				else:
					print 'unrecognised metric\n---> terminating script'
					sys.exit()
				RNA_atomvals.append(y)
				RNA_atomnums.append(atom.atomnum)

		ax.plot(RNA_atomnums,RNA_atomvals)
		ax.text(0.05, 0.90, '# '+str(d_num), color='black',transform=ax.transAxes,fontsize=18)

	fig.text(0.5, 0.04, 'Atom Number', ha='center',fontsize=18)
	fig.text(0.04, 0.5, str(densmet), va='center', rotation='vertical',fontsize=18)
	fig.subplots_adjust(hspace=.4)
	fig.subplots_adjust(wspace=.4)

	fig.suptitle("RNA atom number vs %s" %(str(densmet)),fontsize = 24)
	plt.setp(fig.axes)
	fig.savefig("RNAatomnum_vs_%s.png" %(str(densmet)))




def barplot_densitychangeperrestype(PDBmulti,atomsbychain_present,i,atomtype,residuetype,residuenums):
	# this function plots a series of bar plots for a given damage level i, for all atoms
	# present in the structure specified by atomtype and resitype (which can be strings or 
	# lists of strings). A figure is produced with separate subplots for mean density, min 
	# density, bdamage change, and bfactor change. The barplots show separate bars for 
	# unbound and bound TRAP rings, to distinguish the different dynamics between residues
	# in the 2 different rings

	# determine the locations of residues in the atomsbychain_present list (sorted into 
	# sublists of equivalent chain atoms) corresponding to atomtype and resitype
	counter = -1
	counter_list = []
	for atomset in atomsbychain_present:
		counter += 1
		if (atomset[0].atomtype in atomtype and
		   atomset[0].basetype in residuetype and
		   atomset[0].residuenum in residuenums):
			counter_list.append(counter)

	# loop over each equivalent atom type and determine mean and std of metrics 
	# specified in loop over each atom in an equivalent atom type
	# FIRST FOR UNBOUND RING CASE:
	mean_meanchange_unbound,mean_minchange_unbound,mean_bfacchange_unbound = [],[],[]
	mean_bdamchange_unbound,xTickMarks_unbound = [],[]
	std_meanchange_unbound,std_minchange_unbound = [],[]
	std_bfacchange_unbound,std_bdamchange_unbound = [],[]
	for counter in counter_list:
		atomgroup = atomsbychain_present[counter][0:11]

		# find list of mean density changes for these equivalent atoms --> find gradient wrt dataset number
		meandensity = [atom.meandensity[i] for atom in atomgroup]
		# find list of min density changes for these equivalent atoms --> find gradient wrt dataset number
		mindensity = [atom.mindensity[i] for atom in atomgroup]
		# find list of B factor changes for these equivalent atoms --> find gradient wrt dataset number
		bfacchange = [atom.Bfactorchange[i] for atom in atomgroup]
		# find list of Bdamage changes for these equivalent atoms --> find gradient wrt dataset number
		bdamchange = [atom.bdamchange[i] for atom in atomgroup]

		# calculate means of above lists
		mean_meanchange_unbound.append(np.mean(meandensity))
		mean_minchange_unbound.append(np.mean(mindensity))
		mean_bfacchange_unbound.append(np.mean(bfacchange))
		mean_bdamchange_unbound.append(np.mean(bdamchange))

		# calculate standard deviations of above lists
		std_meanchange_unbound.append(np.std(meandensity))
		std_minchange_unbound.append(np.std(mindensity))
		std_bfacchange_unbound.append(np.std(bfacchange))
		std_bdamchange_unbound.append(np.std(bdamchange))

		# determine the reference name of the current residue type,
		# and append this to list of x-axis labels
		resitype = str(atomgroup[0].basetype)+str(atomgroup[0].residuenum)+str(atomgroup[0].atomtype)
		xTickMarks_unbound.append(resitype)

	# NEXT FOR BOUND RING CASE:
	# loop over each equivalent atom type and determine mean and std of metrics 
	# specified in loop over each atom in an equivalent atom type
	mean_meanchange_bound,mean_minchange_bound,mean_bfacchange_bound = [],[],[]
	mean_bdamchange_bound,xTickMarks_bound = [],[]
	std_meanchange_bound,std_minchange_bound = [],[]
	std_bfacchange_bound,std_bdamchange_bound = [],[]
	for counter in counter_list:
		atomgroup = atomsbychain_present[counter][11:22]

		# find list of mean density changes for these equivalent atoms --> find gradient wrt dataset number
		meandensity = [atom.meandensity[i] for atom in atomgroup]
		# find list of min density changes for these equivalent atoms --> find gradient wrt dataset number
		mindensity = [atom.mindensity[i] for atom in atomgroup]
		# find list of B factor changes for these equivalent atoms --> find gradient wrt dataset number
		bfacchange = [atom.Bfactorchange[i] for atom in atomgroup]
		# find list of Bdamage changes for these equivalent atoms --> find gradient wrt dataset number
		bdamchange = [atom.bdamchange[i] for atom in atomgroup]

		# calculate means of above lists
		mean_meanchange_bound.append(np.mean(meandensity))
		mean_minchange_bound.append(np.mean(mindensity))
		mean_bfacchange_bound.append(np.mean(bfacchange))
		mean_bdamchange_bound.append(np.mean(bdamchange))

		# calculate standard deviations of above lists
		std_meanchange_bound.append(np.std(meandensity))
		std_minchange_bound.append(np.std(mindensity))
		std_bfacchange_bound.append(np.std(bfacchange))
		std_bdamchange_bound.append(np.std(bdamchange))

		# determine the reference name of the current residue type,
		# and append this to list of x-axis labels
		resitype = str(atomgroup[0].basetype)+str(atomgroup[0].residuenum)+str(atomgroup[0].atomtype)
		xTickMarks_bound.append(resitype)

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (40, 12)})
	fig = plt.figure()

	# width of bars in plot
	width = 0.35 

	# Create an axes instance
	ax = plt.subplot(2, 2, 1)
	# x axis is discrete number of different equivalent atom groups
	x = np.arange(len(xTickMarks_bound))
	unboundplt = plt.bar(x, mean_meanchange_unbound,width,
						 yerr=std_meanchange_unbound,color='#6897bb',
						 error_kw=dict(elinewidth=2,ecolor='#31698a'))
	boundplt = plt.bar(x+width, mean_meanchange_bound,width,
		               yerr=std_meanchange_bound,color='#c0d6e4',
		               error_kw=dict(elinewidth=2,ecolor='#31698a'))

	## Custom x-label,y-label    
	ax.set_xlim(-width,len(x)+width)                    
	plt.xlabel('Residue', fontsize=18)
	plt.ylabel('Mean density change', fontsize=18)
	ax.set_xticks(x+width)
	xtickNames = ax.set_xticklabels(xTickMarks_bound)
	plt.setp(xtickNames, rotation=90, fontsize=10)
	# add a legend
	ax.legend( (unboundplt[0], boundplt[0]), ('unbound', 'bound'),loc='best')
	
	# Create an axes instance
	ax = plt.subplot(2, 2, 2)
	unboundplt = plt.bar(x, mean_minchange_unbound,width,
					     yerr=std_minchange_unbound,color='#6897bb',
					     error_kw=dict(elinewidth=2,ecolor='#31698a'))
	boundplt = plt.bar(x+width, mean_minchange_bound,width,
					   yerr=std_minchange_bound,color='#c0d6e4',
					   error_kw=dict(elinewidth=2,ecolor='#31698a'))

	## Custom x-label,y-label
	ax.set_xlim(-width,len(x)+width)                       
	plt.xlabel('Residue', fontsize=18)
	plt.ylabel('Min density change', fontsize=18)
	ax.set_xticks(x+width)
	xtickNames = ax.set_xticklabels(xTickMarks_bound)
	plt.setp(xtickNames, rotation=90, fontsize=10)
	# add a legend
	ax.legend( (unboundplt[0], boundplt[0]), ('unbound', 'bound'),loc='best')
	
	# Create an axes instance
	ax = plt.subplot(2, 2, 3)
	unboundplt = plt.bar(x, mean_bfacchange_unbound,width,
					     yerr=std_bfacchange_unbound,color='#6897bb',
					     error_kw=dict(elinewidth=2,ecolor='#31698a'))
	boundplt = plt.bar(x+width,mean_bfacchange_bound,width,
					   yerr=std_bfacchange_bound,color='#c0d6e4',
					   error_kw=dict(elinewidth=2,ecolor='#31698a'))

	## Custom x-label,y-label   
	ax.set_xlim(-width,len(x)+width)                 
	plt.xlabel('Residue', fontsize=18)
	plt.ylabel('Bfactor change', fontsize=18)
	ax.set_xticks(x+width)
	xtickNames = ax.set_xticklabels(xTickMarks_bound)
	plt.setp(xtickNames, rotation=90, fontsize=10)
	# add a legend
	ax.legend( (unboundplt[0], boundplt[0]), ('unbound', 'bound'),loc='best')
	
	# Create an axes instance
	ax = plt.subplot(2, 2, 4)
	unboundplt = plt.bar(x, mean_bdamchange_unbound,width,
						 yerr=std_bdamchange_unbound,color='#6897bb',
						 error_kw=dict(elinewidth=2,ecolor='#31698a'))
	boundplt = plt.bar(x+width,mean_bdamchange_bound,width,
					   yerr=std_bdamchange_bound,color='#c0d6e4',
					   error_kw=dict(elinewidth=2,ecolor='#31698a'))

	## Custom x-label,y-label  
	ax.set_xlim(-width,len(x)+width)           
	plt.xlabel('Residue', fontsize=18)
	plt.ylabel('Bdamage change', fontsize=18)
	ax.set_xticks(x+width)
	xtickNames = ax.set_xticklabels(xTickMarks_bound)
	plt.setp(xtickNames, rotation=90, fontsize=10)
	# add a legend
	ax.legend( (unboundplt[0], boundplt[0]), ('unbound', 'bound'),loc='best')
	
	plt.subplots_adjust(top=0.90)
	fig.subplots_adjust(hspace=.4)

	fig.suptitle('%s %s damage variability per chain\n:damage set %s'%(str(residuetype),str(atomtype),str(i)),
	             fontsize=24)   
	# Save the figure
	fig.savefig('%s%s_DAMAGEvariability_per_chain_damage%s.png' %(str(residuetype),str(atomtype),str(i)),bbox_inches='tight')



def chiSquaredTest_GLU50unbound(PDBmulti):
	# this little function is designed to determine whether the observed
	# density values for the GLU50 residues on the unbound ring are likely 
	# to arise from the same distribution. In order to test whether the 
	# differing interface between chains in the unbound ring is playing a 
	# statistically significant role on the dynamics of GLU50 damage
	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']

	# first find all unbound chain GLU50 OE1 and OE2 atoms
	GLU50_OE1s,GLU50_OE2s = [],[]
	for atom in PDBmulti:
		if atom.basetype == 'GLU' and atom.residuenum == 50:
			if atom.chaintype in nonRNAbound:
				if atom.atomtype == 'OE1':
					GLU50_OE1s.append(atom)
				elif atom.atomtype == 'OE2':
					GLU50_OE2s.append(atom)
				else:
					pass

	# also make list of all GLU50 OE1 and OE2 atoms together
	GLU50_OEs = GLU50_OE1s + GLU50_OE2s

	# determine lists of GLU50 OE1 mean density change. Here I will treat 
	# each damage set observation for a specific atom as a dimension of the
	# observation
	refval = GLU50_OE1s[0].mindensity[0]
	meandensities_OE1s = [atom.mindensity for atom in GLU50_OE1s]
	meandensities_OE2s = [atom.mindensity for atom in GLU50_OE2s]
	meandensities_OEs = [atom.mindensity for atom in GLU50_OEs]

	meandensities_OE1s_norm,meandensities_OE2s_norm = [],[]
	meandensities_OEs_norm = []
	for atomdamlevels in meandensities_OE1s:
		meandensities_OE1s_norm.append([dens/refval for dens in atomdamlevels])
	print meandensities_OE1s_norm

	for atomdamlevels in meandensities_OE2s:
		meandensities_OE2s_norm.append([dens/refval for dens in atomdamlevels])

	for atomdamlevels in meandensities_OEs:
		meandensities_OEs_norm.append([dens/refval for dens in atomdamlevels])


	# perform chi squared test here
	chi_OE1 = stats.chisquare(meandensities_OE1s_norm)
	chi_OE2 = stats.chisquare(meandensities_OE2s_norm)
	chi_OEs = stats.chisquare(meandensities_OEs_norm)

	print '\nchi squared on GLU50 OE1 atoms:\n'
	print chi_OE1
	print '\nchi squared on GLU50 OE2 atoms:\n'
	print chi_OE2
	print '\nchi squared on GLU50 OE1 and OE2 atoms:\n'
	print chi_OEs


def GLU50unbound_densVmetricCorrelations(PDBmulti,densmet,atommet):

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']

	# first find all unbound chain GLU50 OE1 and OE2 atoms
	GLU50_OE1s,GLU50_OE2s = [],[]
	for atom in PDBmulti:
		if atom.basetype == 'GLU' and atom.residuenum == 50:
			if atom.chaintype in nonRNAbound:
				if atom.atomtype == 'OE1':
					GLU50_OE1s.append(atom)
				elif atom.atomtype == 'OE2':
					GLU50_OE2s.append(atom)
				else:
					pass

	# also make list of all GLU50 OE1 and OE2 atoms together
	GLU50_OEs = GLU50_OE1s + GLU50_OE2s

	# find atomic metrics each atom
	numneighbours = [atom.numsurroundatoms for atom in GLU50_OEs]

	if atommet == 'numneighbours':
		x = [atom.numsurroundatoms for atom in GLU50_OEs]
	elif atommet == 'Bfactor':
		x = [atom.Bfactor[0] for atom in GLU50_OEs]
	else:
		print 'unrecognised density metric\n---> terminating script'
		sys.exit()

	for d_num in range(1,10):
		ax = plt.subplot(3, 3, d_num)
		if densmet == 'mindensitychange':
			y = [atom.mindensity[d_num-1] for atom in GLU50_OEs]
		elif densmet == 'meandensitychange':
			y = [atom.meandensity[d_num-1] for atom in GLU50_OEs]
		elif densmet == 'maxdensitychange':
			y = [atom.maxdensity[d_num-1] for atom in GLU50_OEs]
		elif densmet == 'Bfactor':
			y = [atom.Bfactor[d_num-1] for atom in GLU50_OEs]
		elif densmet == 'Bfactorchange':
			y = [atom.Bfactorchange[d_num-1] for atom in GLU50_OEs]
		elif densmet == 'Bdamage':
			y = [atom.bdam[d_num-1] for atom in GLU50_OEs]
		elif densmet == 'bdamagechange':
			y = [atom.bdamchange[d_num-1] for atom in GLU50_OEs]
		else:
			print 'unrecognised density metric\n---> terminating script'
			sys.exit()

		plt.scatter(x, y, s=30, c='#008080', alpha=0.5)

		# calculate linear regression for each subplot and state r-squared on plot
		slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
		ax.text(0.05, 0.85, '# '+str(d_num+1)+'\nR2: '+str("{0:.2f}".format(r_value**2)),
			    color='black',transform=ax.transAxes,fontsize=18)
		x2 = range(min(x),max(x)+1)
		y2 = slope*np.array(x2) + np.array(intercept)
		plt.plot(x2, y2, '-',color='#d11141',linewidth=3)
	
	f.text(0.5, 0.04, str(atommet), ha='center',fontsize=18)
	f.text(0.04, 0.5, str(densmet), va='center', rotation='vertical',fontsize=18)

	f.suptitle('unbound ring GLU50 OEs:'+str(densmet)+' vs '+str(atommet),fontsize=24)
	plt.setp(f.axes)

	f.savefig('unboundringGLU50OEs_'+str(densmet)+'V'+str(atommet)+'.png')



def GLUASPoxyCorrel(PDBmulti,i,xWhat,yWhat,densmet):
	# this function is designed to determine a correlation between the extent of damage for
	# both carboxyl group side chain oxygens and carbon. The function plots a scatter plot of 
	# xWhat damage vs yWhat damage for each dataset damage level. xWhat and yWhat specify what
	# x and y scatter axes should be - takes values O1,O2,C (corresponding to GLUOE1/ASPOD1,
	# GLUOE2/ASPOE2, GLUCD/ASPCG). densmet specifies the density metric to use and currently 
	# takes values of 'min' and 'mean'

	# first find all GLU/ASP carboxyl side chain oxy atoms
	keptAtoms = []
	for atom in PDBmulti:
		if atom.basetype in ('GLU','ASP'):
			if atom.atomtype in ('OE1','OE2','OD1','OD2'):
				keptAtoms.append(atom)

	# next find all GLU/ASP carboxyl side chain end carbon atoms
	for atom in PDBmulti:
		if atom.basetype in ('GLU') and atom.atomtype == 'CD':
			keptAtoms.append(atom)
		elif atom.basetype in ('ASP') and atom.atomtype == 'CG':
			keptAtoms.append(atom)

	# next create a gluOE1, gluOE2 and gluCD list such that GLU OE1 atoms in gluOE1 list
	# match to GLU OE2 atoms in gluOE2 list and GLU CD atoms in gluCD list
	gluOE1,gluOE2,gluCD=[],[],[]
	for atom in keptAtoms:
		if atom.basetype in 'GLU':
			if atom.atomtype in ('OE1'):
				gluOE1.append(atom)
				for otheratom in keptAtoms:
					if (otheratom.chaintype == atom.chaintype and
						otheratom.residuenum == atom.residuenum): 
						if otheratom.atomtype == 'OE2':
							gluOE2.append(otheratom)
						elif otheratom.atomtype == 'CD':
							gluCD.append(otheratom)
	# check that three lists are same length 
	if (len(gluOE2) != len(gluOE1)) or (len(gluOE2) != len(gluCD)):
		print 'something went wrong!'
		sys.exit()

	# same again for the ASP OD1, OD2 and CG atoms
	aspOD1,aspOD2,aspCG=[],[],[]
	for atom in keptAtoms:
		if atom.basetype in 'ASP':
			if atom.atomtype in ('OD1'):
				aspOD1.append(atom)
				for otheratom in keptAtoms:
					if (otheratom.chaintype == atom.chaintype and
						otheratom.residuenum == atom.residuenum): 
						if otheratom.atomtype == 'OD2':
							aspOD2.append(otheratom)
						elif otheratom.atomtype == 'CG':
							aspCG.append(otheratom)
	# check that three lists are same length 
	if (len(aspOD2) != len(aspOD1)) or (len(aspOD2) != len(aspCG)):
		print 'something went wrong!'
		sys.exit()

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	# want to create x and y lists to plot in scatter plot, GLU first
	if densmet == 'min':
		if xWhat == 'O1':
			x = [atom.mindensity[i] for atom in gluOE1]
		elif xWhat == 'O2':
			x = [atom.mindensity[i] for atom in gluOE2]
		elif xWhat == 'C':
			x = [atom.mindensity[i] for atom in gluCD]
		if yWhat == 'O1':
			y = [atom.mindensity[i] for atom in gluOE1]
		elif yWhat == 'O2':
			y = [atom.mindensity[i] for atom in gluOE2]
		elif yWhat == 'C':
			y = [atom.mindensity[i] for atom in gluCD]
	elif densmet == 'mean':
		if xWhat == 'O1':
			x = [atom.meandensity[i] for atom in gluOE1]
		elif xWhat == 'O2':
			x = [atom.meandensity[i] for atom in gluOE2]
		elif xWhat == 'C':
			x = [atom.meandensity[i] for atom in gluCD]
		if yWhat == 'O1':
			y = [atom.meandensity[i] for atom in gluOE1]
		elif yWhat == 'O2':
			y = [atom.meandensity[i] for atom in gluOE2]
		elif yWhat == 'C':
			y = [atom.meandensity[i] for atom in gluCD]

	gluScatter = plt.scatter(x, y, s=100, c='#d64525', alpha=0.5)

	# calculate linear regression for GLU points here
	slope, intercept, r_valueGLU, p_value, std_err = stats.linregress(x, y)
	# x2 = range(np.floor(min(x)),np.ceil(max(x))+1)
	x2 = [min(x), min(x)+ (max(x) - min(x))/2, max(x)]
	y2 = slope*np.array(x2) + np.array(intercept)
	plt.plot(x2, y2, '-',color='#d64525',linewidth=3)

	# want to create x and y lists to plot in scatter plot, ASP next
	if densmet == 'min':
		if xWhat == 'O1':
			x = [atom.mindensity[i] for atom in aspOD1]
		elif xWhat == 'O2':
			x = [atom.mindensity[i] for atom in aspOD2]
		elif xWhat == 'C':
			x = [atom.mindensity[i] for atom in aspCG]
		if yWhat == 'O1':
			y = [atom.mindensity[i] for atom in aspOD1]
		elif yWhat == 'O2':
			y = [atom.mindensity[i] for atom in aspOD2]
		elif yWhat == 'C':
			y = [atom.mindensity[i] for atom in aspCG]
	elif densmet == 'mean':
		if xWhat == 'O1':
			x = [atom.meandensity[i] for atom in aspOD1]
		elif xWhat == 'O2':
			x = [atom.meandensity[i] for atom in aspOD2]
		elif xWhat == 'C':
			x = [atom.meandensity[i] for atom in aspCG]
		if yWhat == 'O1':
			y = [atom.meandensity[i] for atom in aspOD1]
		elif yWhat == 'O2':
			y = [atom.meandensity[i] for atom in aspOD2]
		elif yWhat == 'C':
			y = [atom.meandensity[i] for atom in aspCG]

	aspScatter = plt.scatter(x, y, s=100, c='#35b9ac', alpha=0.5,marker='d')

	# calculate linear regression for ASP points here
	slope, intercept, r_valueASP, p_value, std_err = stats.linregress(x, y)
	x2 = [min(x), min(x)+ (max(x) - min(x))/2, max(x)]
	y2 = slope*np.array(x2) + np.array(intercept)
	plt.plot(x2, y2, '-',color='#35b9ac',linewidth=3)

	plt.legend([gluScatter,aspScatter],['GLU, R^2: '+str("{0:.2f}".format(r_valueGLU**2)),
										'ASP, R^2: '+str("{0:.2f}".format(r_valueASP**2))],
			   loc=4, borderaxespad=0.,fontsize=24)

	plt.xlabel('%s damage: %s density' %(xWhat,densmet),fontsize=24)
	plt.ylabel('%s damage: %s density' %(yWhat,densmet),fontsize=24)
	f.suptitle('GLU/ASP %s vs %s damage %s' %(xWhat,yWhat,str(i)),fontsize=24)
	plt.setp(f.axes)
	f.savefig('GLUASP_%s_vs_%s_%sdamage%s.png'%(xWhat,yWhat,densmet,str(i)))


def GLUASPdensHist(atomsbychain_present,j,densmet):
	# this function is designed to create a histogram plot of the density loss (specified by the
	# density metric: mean density loss) for different GLU or ASP residues (grouped by unbound and bound
	# TRAP rings)

	# for atoms grouped by equivalency over rings (in atomsbychain_present list), find GLU CD atoms and
	# also ASP CG atoms
	keptAtoms = []
	counter = -1
	for atomGroup in atomsbychain_present:
		counter += 1
		if (atomGroup[0].atomtype == 'CD' and
		    atomGroup[0].basetype == 'GLU' and
		    atomGroup[0].residuenum in j):
			keptAtoms.append(atomGroup)		
		elif (atomGroup[0].atomtype == 'CG' and
		    atomGroup[0].basetype == 'ASP' and
		    atomGroup[0].residuenum in j):
			keptAtoms.append(atomGroup)

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()	

	for d_num in range(1,10):
		ax = plt.subplot(3, 3, d_num)
		for atomGroup in keptAtoms:

			if densmet == 'min':
				data_unbound = [atom.mindensity[d_num-1] for atom in atomGroup[0:11]]
				data_bound = [atom.mindensity[d_num-1] for atom in atomGroup[11:22]]
			elif densmet == 'mean':
				data_unbound = [atom.meandensity[d_num-1] for atom in atomGroup[0:11]]
				data_bound = [atom.meandensity[d_num-1] for atom in atomGroup[11:22]]

			sns.kdeplot(np.array(data_unbound), shade=True, label = str(atomGroup[0].residuenum)+ ' unbound')
			sns.kdeplot(np.array(data_bound), shade=True, label = str(atomGroup[0].residuenum)+ ' bound')
		# make a legend	for kde plots
		plt.legend(loc='best')

	varAtomGroups,var2AtomGroups,varAtomNames = [],[],[]
	for atomGroup in keptAtoms:
		varAtomGroup_unbound,varAtomGroup_bound = [],[]
		var2AtomGroup_unbound,var2AtomGroup_bound = [],[]
		for d_num in range(1,10):

			if densmet == 'min':
				data_unbound = [atom.mindensity[d_num-1] for atom in atomGroup[0:11]]
				data_bound = [atom.mindensity[d_num-1] for atom in atomGroup[11:22]]
			elif densmet == 'mean':
				data_unbound = [atom.meandensity[d_num-1] for atom in atomGroup[0:11]]
				data_bound = [atom.meandensity[d_num-1] for atom in atomGroup[11:22]]

			# determine the variance in each set of grouped equivalent atoms
			varAtomGroup_unbound.append(np.var(data_unbound))
			varAtomGroup_bound.append(np.var(data_bound))
			# also determine the variance normalised by (mean)^2 in each set of grouped
			# equivalent atoms - this is in an attempt to create a metric involving
			# variance that is not dependent on the mean level of damage (since more
			# damage will imply more variance between equivalent atoms)
			var2AtomGroup_unbound.append(np.var(data_unbound)/(np.mean(data_unbound)**2))
			var2AtomGroup_bound.append(np.var(data_bound)/(np.mean(data_unbound)**2))
		
		varAtomGroups.append(varAtomGroup_unbound)
		varAtomGroups.append(varAtomGroup_bound)
		var2AtomGroups.append(var2AtomGroup_unbound)
		var2AtomGroups.append(var2AtomGroup_bound)
		varAtomNames.append(str(atomGroup[0].residuenum)+ ' unbound')
		varAtomNames.append(str(atomGroup[0].residuenum)+ ' bound')

	f.text(0.5, 0.04, 'GLU CD/ASP CG damage: min density change', ha='center',fontsize=18)
	f.text(0.04, 0.5, 'Freq', va='center', rotation='vertical',fontsize=18)

	f.suptitle('GLU CD/ASP CG %s density change distribution' %(densmet),fontsize=24)
	plt.setp(f.axes)
	f.savefig('GLUCD_ASPCG_%sDensDistn_%s.png' %(densmet,str(j)))

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	g = plt.figure()	

	counter = -1
	colorlist = ['#fbf96d','#f57d62','#e15b64','#abbd83','#83abbd','#9e4a9e',
	  			 '#573f94','#ff0000','#014a01','#06357a','#ec891d','#507e4e']
	for atomGroup in varAtomGroups:
		counter += 1
		plt.plot(range(len(atomGroup)),atomGroup,label=varAtomNames[counter],
			     color=colorlist[counter],marker='o',markersize=10,linestyle='None')
		m,c = np.polyfit(range(len(atomGroup)), atomGroup, 1)
		y = np.array(range(len(atomGroup)))*m + c
		plt.plot(range(len(atomGroup)),y,color=colorlist[counter],linewidth=3)

	# make a legend	for line plot
	plt.legend(loc='best',fontsize=24)
	plt.xlabel('dataset #',fontsize=24)
	plt.ylabel('variance: %s density change ' %(densmet),fontsize=24)
	g.suptitle('variance vs dataset #: GLU CD/ASP CG %s' %(str(j)),fontsize=24)
	plt.setp(g.axes)
	g.savefig('varianceVdatasetGLUCD_ASPCG_%s_%s.png' %(str(j),densmet))


	# Create a figure instance for the variance metric normalised by mean^2
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	h = plt.figure()	

	counter = -1
	colorlist = ['#fbf96d','#f57d62','#e15b64','#abbd83','#83abbd','#9e4a9e',
	  			 '#573f94','#ff0000','#014a01','#06357a','#ec891d','#507e4e']
	for atomGroup in var2AtomGroups:
		counter += 1
		plt.plot(range(len(atomGroup)),atomGroup,label=varAtomNames[counter],
			     color=colorlist[counter],marker='o',markersize=10)

	# make a legend	for line plot
	plt.legend(loc='best',fontsize=24)
	plt.xlabel('dataset #',fontsize=24)
	plt.ylabel('var/mean^2: %s density change ' %(densmet),fontsize=24)
	h.suptitle('var/mean^2 vs dataset #: GLU CD/ASP CG %s' %(str(j)),fontsize=24)
	plt.setp(h.axes)
	h.savefig('var2MetricVdatasetGLUCD_ASPCG_%s_%s.png' %(str(j),densmet))


def densitychange_v_dataset_linfit(atomsbychain_present,atomtype,resitype,
							resinum,densmet):
	# function to plot density change as function of dataset number
	# for a specific atom in the structure, with all the inclusion
	# of equivalent atoms in other chains also (in the case of the 
	# RNA this will include all symmetrically related atoms, by 
	# 11-fold symmetry). Linearly fits (using polyfit) and outputs gradients
	# of each linear fitting. Designed specifically to show relationship between
	# damage and crystal contact types for Glu50 on unbound TRAP ring.

	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	# RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']
	RNAbound = []
	# determine the desired atom by atom type, residue type and residue number.
	# for RNA only the first 4 nucleotides are chosen, so the residue number
	# must be selected out of 101,102,103,104
	counter = -1
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype and
		   atom[0].basetype == resitype and
		   atom[0].residuenum == resinum):
			break

	sns.set(style="white", context="talk")
	f, axes = plt.subplots(1, 1, figsize=(12, 12), sharex=True)

	# define x range here (damage set numbers)
	x = range(2,len(atom[0].meandensity)+2)

	# colours for plot defined here
	jet = cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(atom)-11)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

	for i in range(0,len(atom)):
		if densmet == 'mindensitychange':
			y = atom[i].mindensity
		elif densmet == 'meandensitychange':
			y = atom[i].meandensity
		elif densmet == 'maxdensitychange':
			y = atom[i].maxdensity
		elif densmet == 'Bfactor':
			y = atom[i].Bfactor
		elif densmet == 'Bfactorchange':
			y = atom[i].Bfactorchange
		elif densmet == 'Bdamage':
			y = atom[i].bdam
		elif densmet == 'Bdamagechange':
			y = atom[i].bdamchange
		else:
			print 'unrecognised density metric\n---> terminating script'
			sys.exit()

		# set colour of current line here
		colorVal = scalarMap.to_rgba(i)

		# if protein atom label by chain, if RNA atom label by nucleotide number
		if atom[0].chaintype != 'W':
			if atom[i].chaintype in nonRNAbound:
				plt.plot(x[0:9],y[0:9],color=colorVal,label='chain: '+str(atom[i].chaintype),linestyle='none',marker='o')
				z = np.polyfit(x[0:9], y[0:9], 1)
				plt.plot(x[0:9],z[0]*np.array(x[0:9]) + z[1],color=colorVal)
				print 'Chain: '+str(atom[i].chaintype)+', gradient: '+str(z[0])
			elif atom[i].chaintype in RNAbound:
				plt.plot(x[0:9],y[0:9],color=colorVal,label='chain: '+str(atom[i].chaintype),linestyle='none',marker='o')
				z = np.polyfit(x[0:9], y[0:9], 1)
				plt.plot(x[0:9],z[0]*np.array(x[0:9]) + z[1],color=colorVal)

		else:
			plt.plot(x,y,color=colorVal,label='nucleotide: '+str(atom[i].residuenum))

	plt.xlabel('damage set', fontsize=18)
	plt.ylabel(str(densmet), fontsize=18)
	
	if atom[0].chaintype != 'W':
		f.suptitle(str(densmet)+' vs damage set: '+\
		  	       str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		    	   str(atomsbychain_present[counter][0].residuenum) + ' ' +\
		     	   str(atomsbychain_present[counter][0].basetype),fontsize=24)
	else:
		f.suptitle(str(densmet)+' vs damage set: '+\
		  	       str(atomsbychain_present[counter][0].atomtype) + ' ' +\
		    	   str(atomsbychain_present[counter][0].chaintype) + ' ' +\
		     	   str(atomsbychain_present[counter][0].basetype),fontsize=24)
	plt.setp(f.axes)
	plt.subplots_adjust(top=0.90)

	# Shrink current axis by 20%
	box = axes.get_position()
	axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	# Put a legend to the right of the current axis
	axes.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=18)

	f.savefig('%s_vs_damageset_%s_%s_%s.png' %(str(densmet),str(atomsbychain_present[counter][0].atomtype),
		       										        str(atomsbychain_present[counter][0].residuenum),
		       											    str(atomsbychain_present[counter][0].basetype)))


def hotellingTsquareTest(pdbmulti,atomtypes,basetypes,residuenums):
	# this function is designed to calculate Hotelling's T-square test statistic
	# in order to compare bound vs bound sets of atoms over multiple datasets 
	# simultaneously (each dose level is a different dimension to the data)
	keptAtoms_unbound = []
	keptAtoms_bound = []

	print '---------------------------------------------------'
	print "Hotelling's T-square Test"
	print 'Residue type: %s' %(str(basetypes))
	print 'Residue number: %s' %(str(residuenums))
	print 'Atom type: %s' %(str(atomtypes))
	print 'Max density loss metric used...'

	# set the number of datasets needed to be included
	p = 9
	print 'First %s datasets included...' %(str(p))
	print 'Running test...'

	# find the correct atoms to be included from each chain
	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']
	for atom in pdbmulti:
		if (atom.atomtype in atomtypes and
			atom.basetype in basetypes and
			atom.residuenum in residuenums):
			col = []
			for dens in atom.mindensity[0:p]:
				col.append([dens])
			if atom.chaintype in nonRNAbound:
				keptAtoms_unbound.append(col)
			elif atom.chaintype in RNAbound:
				keptAtoms_bound.append(col)		

	# check that atom found in all chains
	if (len(keptAtoms_unbound) != 11 or len(keptAtoms_bound) != 11):
		print 'error --> terminating script'
		sys.exit()

	# calculate means for unbound and bound vectors
	nx,ny = 11,11
	unbound_mean = np.sum(np.array(keptAtoms_unbound),axis=0)/nx
	bound_mean = np.sum(np.array(keptAtoms_bound),axis=0)/ny

	# calculate the unbiased pooled covariance matrix estimate:
	W = 0
	for i in range(0,nx):
		unboundPart = (np.array(keptAtoms_unbound)[i] - unbound_mean)*np.transpose((np.array(keptAtoms_unbound)[i] - unbound_mean))
		boundPart = (np.array(keptAtoms_bound)[i] - bound_mean)*np.transpose((np.array(keptAtoms_bound)[i] - bound_mean))
		combinedParts = unboundPart + boundPart
		W = W + combinedParts
	W = W/(nx+ny-2)

	# now calculate the Hotelling's two-sample T-squared statistic:
	tSquared = (float(nx*ny)/(nx+ny))*np.dot(np.transpose(unbound_mean-bound_mean),np.dot(np.linalg.inv(W),(unbound_mean-bound_mean)))
	print 'tSquared is as follows: {}'.format(tSquared)

	# now calculate the related F test statistic from F distribution:
	F = (float(nx+ny-p-1)/((nx+ny-2)*p))*tSquared
	print 'F-value is as follows: {}'.format(F)

	# generate a p-value for the given F test statistic
	alpha = 0.05 
	df1 = p
	df2 = nx+ny-1-p
	p_value = 1-stats.f.cdf(F, df1, df2)
	print 'p-value is as follows: {}'.format(p_value)

	reject = 'notsureyet'
	if p_value < alpha:
		print 'p_value < %s' %(str(alpha))
		print '---> reject null hypothesis that means are equal'
		reject = 'YES'
	else:
		print 'p_value > %s' %(str(alpha))
		print '---> cannot reject null hypothesis that means are equal'
		reject = 'NO'
	print '---------------------------------------------------'
	return F,p_value,reject

def hotellingTsquareTest_batch(pdbmulti):
	# run multiple hotellingTsquareTests above in one go

	# for Glu atoms, use the following:
	# atomtypes = [['CA'],['CB'],['CG'],['CD'],['OE1'],['OE2']]
	# basetypes = [['GLU']]
	# residuenums = [[16],[36],[42],[50],[71],[73]]

	# for Asp atoms, use the following:
	atomtypes = [['CA'],['CB'],['CG'],['OD1'],['OD2']]
	basetypes = [['ASP']]
	residuenums = [[8],[17],[29],[39]]

	for residuenum in residuenums:
		for basetype in basetypes:
			for atomtype in atomtypes:
				F,p_value,reject = hotellingTsquareTest(pdbmulti,atomtype,basetype,residuenum)
				print '| %s-%s-%s | F-value: %s | p-value: %s | reject?: %s |' %(str(residuenum[0]),basetype[0],atomtype[0],str(F[0][0]),str(p_value[0][0]),reject)

def RNAdensitychange_v_dataset_errorbars_4subplots(atomsbychain_present,atomtype):
	# function to plot density change as function of dataset number
	# for a specific atom in the structure, with the mean value over 
	# all RNA 4-nucleotide repeats plotted, along with error bars calculated over the
	# 11 instances of the 4-nucleotide repeat

	# determine the desired atom by atom type, residue type and residue number.
	# for RNA only the first 4 nucleotides are chosen, so the residue number
	# must be selected out of 101,102,103,104
	counter = -1
	counter_list = []
	print 'Locating atoms of interest around RNA strand:'
	for atom in atomsbychain_present:
		counter += 1
		if (atom[0].atomtype == atomtype):
			counter_list.append(counter)
			print 'Chain: %s - Atom: %s - Nucleotide: %s' %(str(atom[0].atomtype),str(atom[0].chaintype),str(atom[0].basetype))

	# check that atom found in all 4 nucleotides:
	if len(counter_list) != 4:
		print len(counter_list)
		print 'error --> terminating script'
		sys.exit()

	# only can handle RNA here so check it is a RNA atom that 
	# has been selected
	if atom[0].chaintype not in ('W'):
		print 'Script cannot handle anything but RNA atoms'
		print '---> terminating script'
		sys.exit()

	sns.set(style="white", context="talk")
	f = plt.figure()

	# define x range here (damage set numbers)
	x = range(2,len(atom[0].meandensity)+2)

	# plot for G1, A2, G3, U4 nucleotides separately
	atom_G1 = atomsbychain_present[counter_list[0]]
	atom_A2 = atomsbychain_present[counter_list[1]]
	atom_G3 = atomsbychain_present[counter_list[2]]
	atom_U4 = atomsbychain_present[counter_list[3]]

	# mean density axis here
	y_mean_G1,y_mean_A2,y_mean_G3,y_mean_U4 = [],[],[],[]
	y_std_G1,y_std_A2,y_std_G3,y_std_U4 = [],[],[],[]
	for j in range(len(x)):
		y_equivatms_G1,y_equivatms_A2 = [],[]
		y_equivatms_G3,y_equivatms_U4 = [],[]
		for i in range(0,len(atom)):
			y_equivatms_G1.append(atom_G1[i].meandensity[j])
			y_equivatms_A2.append(atom_A2[i].meandensity[j])
			y_equivatms_G3.append(atom_G3[i].meandensity[j])
			y_equivatms_U4.append(atom_U4[i].meandensity[j])
		y_mean_G1.append(np.mean(y_equivatms_G1))
		y_mean_A2.append(np.mean(y_equivatms_A2))
		y_mean_G3.append(np.mean(y_equivatms_G3))
		y_mean_U4.append(np.mean(y_equivatms_U4))

		y_std_G1.append(np.std(y_equivatms_G1))
		y_std_A2.append(np.std(y_equivatms_A2))
		y_std_G3.append(np.std(y_equivatms_G3))
		y_std_U4.append(np.std(y_equivatms_U4))

	ax = plt.subplot(221)
	plt.errorbar(x,y_mean_G1,yerr=y_std_G1, fmt='-o',capthick=2,color='#00aedb',label='G1')
	plt.errorbar(x,y_mean_A2,yerr=y_std_A2,fmt='-o',capthick=2,color='#a200ff',label='A2')
	plt.errorbar(x,y_mean_G3,yerr=y_std_G3,fmt='-o',capthick=2,color='#f47835',label='G3')
	plt.errorbar(x,y_mean_U4,yerr=y_std_U4,fmt='-o',capthick=2,color='#d41243',label='U4')

	ax.set_xlim([1, 11])
	ax.legend(loc='best')

	plt.xlabel('Damage set')
	plt.ylabel('Mean density change')

	# min density axis here
	y_mean_G1,y_mean_A2,y_mean_G3,y_mean_U4 = [],[],[],[]
	y_std_G1,y_std_A2,y_std_G3,y_std_U4 = [],[],[],[]
	for j in range(len(x)):
		y_equivatms_G1,y_equivatms_A2 = [],[]
		y_equivatms_G3,y_equivatms_U4 = [],[]
		for i in range(0,len(atom)):
			y_equivatms_G1.append(atom_G1[i].mindensity[j])
			y_equivatms_A2.append(atom_A2[i].mindensity[j])
			y_equivatms_G3.append(atom_G3[i].mindensity[j])
			y_equivatms_U4.append(atom_U4[i].mindensity[j])
		y_mean_G1.append(np.mean(y_equivatms_G1))
		y_mean_A2.append(np.mean(y_equivatms_A2))
		y_mean_G3.append(np.mean(y_equivatms_G3))
		y_mean_U4.append(np.mean(y_equivatms_U4))

		y_std_G1.append(np.std(y_equivatms_G1))
		y_std_A2.append(np.std(y_equivatms_A2))
		y_std_G3.append(np.std(y_equivatms_G3))
		y_std_U4.append(np.std(y_equivatms_U4))

	ax = plt.subplot(222)
	plt.errorbar(x,y_mean_G1,yerr=y_std_G1, fmt='-o',capthick=2,color='#00aedb',label='G1')
	plt.errorbar(x,y_mean_A2,yerr=y_std_A2,fmt='-o',capthick=2,color='#a200ff',label='A2')
	plt.errorbar(x,y_mean_G3,yerr=y_std_G3,fmt='-o',capthick=2,color='#f47835',label='G3')
	plt.errorbar(x,y_mean_U4,yerr=y_std_U4,fmt='-o',capthick=2,color='#d41243',label='U4')

	ax.set_xlim([1, 11])
	ax.legend(loc='best')

	plt.xlabel('Damage set')
	plt.ylabel('Max density loss')

	# Max density gain axis here	
	y_mean_G1,y_mean_A2,y_mean_G3,y_mean_U4 = [],[],[],[]
	y_std_G1,y_std_A2,y_std_G3,y_std_U4 = [],[],[],[]
	for j in range(len(x)):
		y_equivatms_G1,y_equivatms_A2 = [],[]
		y_equivatms_G3,y_equivatms_U4 = [],[]
		for i in range(0,len(atom)):
			y_equivatms_G1.append(atom_G1[i].maxdensity[j])
			y_equivatms_A2.append(atom_A2[i].maxdensity[j])
			y_equivatms_G3.append(atom_G3[i].maxdensity[j])
			y_equivatms_U4.append(atom_U4[i].maxdensity[j])
		y_mean_G1.append(np.mean(y_equivatms_G1))
		y_mean_A2.append(np.mean(y_equivatms_A2))
		y_mean_G3.append(np.mean(y_equivatms_G3))
		y_mean_U4.append(np.mean(y_equivatms_U4))

		y_std_G1.append(np.std(y_equivatms_G1))
		y_std_A2.append(np.std(y_equivatms_A2))
		y_std_G3.append(np.std(y_equivatms_G3))
		y_std_U4.append(np.std(y_equivatms_U4))

	ax = plt.subplot(223)
	plt.errorbar(x,y_mean_G1,yerr=y_std_G1, fmt='-o',capthick=2,color='#00aedb',label='G1')
	plt.errorbar(x,y_mean_A2,yerr=y_std_A2,fmt='-o',capthick=2,color='#a200ff',label='A2')
	plt.errorbar(x,y_mean_G3,yerr=y_std_G3,fmt='-o',capthick=2,color='#f47835',label='G3')
	plt.errorbar(x,y_mean_U4,yerr=y_std_U4,fmt='-o',capthick=2,color='#d41243',label='U4')

	ax.set_xlim([1, 11])
	ax.legend(loc='best')

	plt.xlabel('Damage set')
	plt.ylabel('Max density gain')

	# Bfactor change axis here
	y_mean_G1,y_mean_A2,y_mean_G3,y_mean_U4 = [],[],[],[]
	y_std_G1,y_std_A2,y_std_G3,y_std_U4 = [],[],[],[]
	for j in range(len(x)):
		y_equivatms_G1,y_equivatms_A2 = [],[]
		y_equivatms_G3,y_equivatms_U4 = [],[]
		for i in range(0,len(atom)):
			y_equivatms_G1.append(atom_G1[i].Bfactorchange[j])
			y_equivatms_A2.append(atom_A2[i].Bfactorchange[j])
			y_equivatms_G3.append(atom_G3[i].Bfactorchange[j])
			y_equivatms_U4.append(atom_U4[i].Bfactorchange[j])
		y_mean_G1.append(np.mean(y_equivatms_G1))
		y_mean_A2.append(np.mean(y_equivatms_A2))
		y_mean_G3.append(np.mean(y_equivatms_G3))
		y_mean_U4.append(np.mean(y_equivatms_U4))

		y_std_G1.append(np.std(y_equivatms_G1))
		y_std_A2.append(np.std(y_equivatms_A2))
		y_std_G3.append(np.std(y_equivatms_G3))
		y_std_U4.append(np.std(y_equivatms_U4))

	ax = plt.subplot(224)
	plt.errorbar(x,y_mean_G1,yerr=y_std_G1, fmt='-o',capthick=2,color='#00aedb',label='G1')
	plt.errorbar(x,y_mean_A2,yerr=y_std_A2,fmt='-o',capthick=2,color='#a200ff',label='A2')
	plt.errorbar(x,y_mean_G3,yerr=y_std_G3,fmt='-o',capthick=2,color='#f47835',label='G3')
	plt.errorbar(x,y_mean_U4,yerr=y_std_U4,fmt='-o',capthick=2,color='#d41243',label='U4')

	ax.set_xlim([1, 11])
	ax.legend(loc='best')

	plt.xlabel('Damage set')
	plt.ylabel('Bfactor change')

	plt.subplots_adjust(top=0.90)
	f.subplots_adjust(hspace=.3)
	f.subplots_adjust(wspace=.4)

	f.suptitle('damage metrics vs damage set: '+\
		  	   'atom ' + str(atom_G1[0].atomtype) + ' ' +\
		       'chain ' + str(atomsbychain_present[counter][0].basetype),fontsize=20)

	f.savefig('4metrics_vs_damageset_RNA_%s_errorbars_4subplots.png' %(str(atom_G1[0].atomtype)))


def damageRel2CA(pdbmulti,i,threshold):
	# this little function calculates the damage for each atom relative to the CA of the same residue
	for atom in pdbmulti:
		if atom.chaintype not in ('W','Y'):
			for otheratom in pdbmulti:
				if (otheratom.residuenum == atom.residuenum and
					otheratom.chaintype == atom.chaintype and
					otheratom.atomtype == 'CA'):
					damRel2CA = float(atom.maxdensity[i])/otheratom.maxdensity[i]
					if damRel2CA > threshold:
						print '%s %s %s %s --> %s' %(str(atom.chaintype),str(atom.residuenum),str(atom.basetype),str(atom.atomtype),str(damRel2CA))
					break


def RNArankdamage(pdbmulti,N,basetype,basenumber):
	# this function ranks the atoms of the RNA in terms of density metric
	phosphate_counterList,deoxyribose_counterList,base_counterList = [],[],[]
	phosphate_densMeans,deoxyribose_densMeans,base_densMeans = [],[],[]
	phosphate_densStd,deoxyribose_densStd,base_densStd = [],[],[]
	GLU_densMeans = []
	GLU_densStd = []

	for i in range(len(pdbmulti[0].mindensity)):
		pdbmulti.sort(key=lambda x: (x.mindensity[i]))
		RNA_counter,phosphate_counter,deoxyribose_counter,base_counter = 0,0,0,0 
		phosphate_dens,deoxyribose_dens,base_dens = [],[],[]
		GLU_dens = []
		for atom in pdbmulti:
			if (atom.chaintype == 'W' and atom.basetype == basetype and atom.residuenum%5 == basenumber):
				RNA_counter += 1
				if atom.atomtype in ['P','OP1','OP2']:
					RNApart = 'Phosphate'
					phosphate_counter += 1
					phosphate_dens.append(atom.mindensity[i])
				elif atom.atomtype in ["C5'","C4'","C3'","C2'","O2'","C1'","O4'","O5'","O3'"]:
					RNApart = 'Deoxyribose'
					deoxyribose_counter += 1
					deoxyribose_dens.append(atom.mindensity[i])
				else:
					RNApart = 'Base'
					base_counter += 1
					base_dens.append(atom.mindensity[i])
			else:
				if atom.basetype == 'GLU' and atom.atomtype in ('CD','OE1','OE2'):
					GLU_dens.append(atom.mindensity[i])

				# print str(atom.basetype) + '---' + str(atom.residuenum) + '---' + str(atom.atomtype) + ' ---> ' + str(atom.mindensity[i]) + ' ---> ' + RNApart
			if RNA_counter > N:
				break

		print '----------------------------------'
		print 'Dataset %s' %(str(i+2))
		print '# phosphate: %s' %(str(phosphate_counter))
		print '# deoxyribose: %s' %(str(deoxyribose_counter))
		print '# base: %s' %(str(base_counter))
		phosphate_counterList.append(phosphate_counter)
		deoxyribose_counterList.append(deoxyribose_counter)
		base_counterList.append(base_counter)
		phosphate_densMeans.append(np.mean(phosphate_dens))
		deoxyribose_densMeans.append(np.mean(deoxyribose_dens))
		base_densMeans.append(np.mean(base_dens))
		phosphate_densStd.append(np.std(phosphate_dens))
		deoxyribose_densStd.append(np.std(deoxyribose_dens))
		base_densStd.append(np.std(base_dens))

		GLU_densMeans.append(np.mean(GLU_dens))
		GLU_densStd.append(np.std(GLU_dens))

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	#plt.plot(range(2,len(pdbmulti[0].mindensity)+2), phosphate_counterList,linewidth=3,label='phosphate')
	#plt.plot(range(2,len(pdbmulti[0].mindensity)+2), deoxyribose_counterList,linewidth=3,label='deoxyribose')
	#plt.plot(range(2,len(pdbmulti[0].mindensity)+2), base_counterList,linewidth=3,label='base')
	plt.errorbar(range(2,len(pdbmulti[0].mindensity)+2), phosphate_densMeans,yerr=phosphate_densStd, fmt='-o',capthick=2,label='phosphate')
	plt.errorbar(range(2,len(pdbmulti[0].mindensity)+2), deoxyribose_densMeans,yerr=deoxyribose_densStd, fmt='-o',capthick=2,label='deoxyribose')
	plt.errorbar(range(2,len(pdbmulti[0].mindensity)+2), base_densMeans,yerr=base_densStd, fmt='-o',capthick=2,label='base')
	plt.errorbar(range(2,len(pdbmulti[0].mindensity)+2), GLU_densMeans,yerr=GLU_densStd, fmt='-o',capthick=2,label='glu')

	# make a legend	for line plot
	plt.legend(loc='best',fontsize=24)
	plt.xlabel('Dataset #',fontsize=24)
	plt.ylabel('Max density loss',fontsize=24)
	f.suptitle('Max Density Loss vs Dataset',fontsize=24)
	plt.setp(f.axes)
	f.savefig('maxdenslossVdatasetRNA_%s%s.png'%(str(basetype),str(basenumber)))


def RNArankdamage_RNAsubsection(pdbmulti,RNApart):
	# this function ranks the atoms of the RNA in terms of density metric. Here either the base, sugar, or phosphate
	# part of the RNA for each nucleotide type are compared
	G1_densMeans,A2_densMeans,G3_densMeans,U4_densMeans = [],[],[],[]
	G1_densStd,A2_densStd,G3_densStd,U4_densStd = [],[],[],[]

	# find the correct atoms depending on the part of the RNA to compare. For the base, the complement of the 
	# atomtypes will be taken later
	if RNApart == 'phosphate':
		atomtypes = ['P','OP1','OP2']
	elif RNApart == 'sugar':
		atomtypes = ["C5'","C4'","C3'","C2'","O2'","C1'","O4'","O5'","O3'"]
	elif RNApart == 'base':
		atomtypes = ['P','OP1','OP2'] + ["C5'","C4'","C3'","C2'","O2'","C1'","O4'","O5'","O3'"]

	for i in range(len(pdbmulti[0].maxdensity)):
		pdbmulti.sort(key=lambda x: (x.maxdensity[i]))
		G1_dens,A2_dens,G3_dens,U4_dens = [],[],[],[]

		for atom in pdbmulti:
			if atom.chaintype == 'W': 
				if (atom.basetype == 'G' and atom.residuenum%5 == 1):
					if RNApart in ('phosphate','sugar'):
						if atom.atomtype in atomtypes:
							G1_dens.append(atom.maxdensity[i])
					elif RNApart in ('base'):
						if atom.atomtype not in atomtypes:
							G1_dens.append(atom.maxdensity[i])
				elif (atom.basetype == 'A' and atom.residuenum%5 == 2):
					if RNApart in ('phosphate','sugar'):
						if atom.atomtype in atomtypes:
							A2_dens.append(atom.maxdensity[i])
					elif RNApart in ('base'):
						if atom.atomtype not in atomtypes:
							A2_dens.append(atom.maxdensity[i])
				elif (atom.basetype == 'G' and atom.residuenum%5 == 3):
					if RNApart in ('phosphate','sugar'):
						if atom.atomtype in atomtypes:
							G3_dens.append(atom.maxdensity[i])
					elif RNApart in ('base'):
						if atom.atomtype not in atomtypes:
							G3_dens.append(atom.maxdensity[i])
				elif (atom.basetype == 'U' and atom.residuenum%5 == 4):
					if RNApart in ('phosphate','sugar'):
						if atom.atomtype in atomtypes:
							U4_dens.append(atom.maxdensity[i])
					elif RNApart in ('base'):
						if atom.atomtype not in atomtypes:
							U4_dens.append(atom.maxdensity[i])
				else:
					print 'error determining RNA base type --> terminating script'
					sys.exit()

		# calculate means and standard deviations for each nucleotide type, for current dataset
		G1_densMeans.append(np.mean(G1_dens))
		A2_densMeans.append(np.mean(A2_dens))
		G3_densMeans.append(np.mean(G3_dens))
		U4_densMeans.append(np.mean(U4_dens))
		G1_densStd.append(np.std(G1_dens))
		A2_densStd.append(np.std(A2_dens))
		G3_densStd.append(np.std(G3_dens))
		U4_densStd.append(np.std(U4_dens))

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	plt.errorbar(range(2,len(pdbmulti[0].mindensity)+2), G1_densMeans,yerr=G1_densStd, fmt='-o',capthick=2,label='G1')
	plt.errorbar(range(2,len(pdbmulti[0].mindensity)+2), A2_densMeans,yerr=A2_densStd, fmt='-o',capthick=2,label='A2')
	plt.errorbar(range(2,len(pdbmulti[0].mindensity)+2), G3_densMeans,yerr=G3_densStd, fmt='-o',capthick=2,label='G3')
	plt.errorbar(range(2,len(pdbmulti[0].mindensity)+2), U4_densMeans,yerr=U4_densStd, fmt='-o',capthick=2,label='U4')

	# make a legend	for line plot
	plt.legend(loc='best',fontsize=24)
	plt.xlabel('Dataset #',fontsize=24)
	plt.ylabel('Max density gain',fontsize=24)
	f.suptitle('Max Density Gain vs Dataset',fontsize=24)
	plt.setp(f.axes)
	f.savefig('maxdensgainVdatasetRNA_%s.png'%(str(RNApart)))


def ResTypeDamageErrorBara(pdbmulti,res1,atomtypes1,res2,atomtypes2,normalise):
	# this function ranks the atoms of a particular residue in terms of density metric
	# 'normalise' input is boolian and determines whether plot metrics should be normalised
	# to background Calpha density loss
	
	res1_densMeans,res2_densMeans = [],[]
	res1_densStd,res2_densStd = [],[]

	# list of doses
	doses = [1.31,3.88,6.45,9.02,11.58,14.15,16.72,19.29,21.86,24.98,33.01,40.86,48.89,56.75,64.78,72.63]
	
	for i in range(len(pdbmulti[0].mindensity)):
		
		# use Calpha atoms as measure of background deteriation of density map
		decay_normaliser = [] 
		for otheratom in pdbmulti:
			if otheratom.atomtype == 'CA':
				decay_normaliser.append(otheratom.mindensity[i])
		decay_factor = np.mean(decay_normaliser)

		res1_dens,res2_dens = [],[]
		for atom in pdbmulti:
			if (atom.basetype == res1 and atom.atomtype == atomtypes1):
				if normalise == True:
					res1_dens.append((atom.mindensity[i]-decay_factor)/decay_factor)
				else:
					res1_dens.append(atom.mindensity[i])
			elif (atom.basetype == res2 and atom.atomtype == atomtypes2):
				if normalise == True:
					res2_dens.append((atom.mindensity[i]-decay_factor)/decay_factor)
				else:
					res2_dens.append(atom.mindensity[i])

		res1_densMeans.append(np.mean(res1_dens))
		res1_densStd.append(np.std(res1_dens))
		res2_densMeans.append(np.mean(res2_dens))
		res2_densStd.append(np.std(res2_dens))

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (8, 8)})
	f = plt.figure()

	plt.errorbar(doses[0:len(pdbmulti[0].mindensity)], res1_densMeans,yerr=res1_densStd, fmt='-o',capthick=2,label=str(res1)+': '+str(atomtypes1)+ ' atoms')
	plt.errorbar(doses[0:len(pdbmulti[0].mindensity)], res2_densMeans,yerr=res2_densStd, fmt='-o',capthick=2,label=str(res2)+': '+str(atomtypes2)+ ' atoms')

	# make a legend	for line plot
	plt.legend(loc='best',fontsize=24)
	plt.xlabel('Dose (MGy)',fontsize=24)
	plt.ylabel('Max density loss',fontsize=24)
	f.suptitle('Max Density loss vs Dose',fontsize=24)
	plt.setp(f.axes)

	if normalise == True:
		f.savefig('NORMALISEDmaxdenslossVdwd_%svs%s.png' %(str(res1)+str(atomtypes1),str(res2)+str(atomtypes2)))
	else:
		f.savefig('maxdenslossVdwd_%svs%s.png' %(str(res1)+str(atomtypes1),str(res2)+str(atomtypes2)))

def gluAspDamageRateComparison(pdbmulti,formatted,split):
	# this function determines the density loss rates for each glu and asp atom of interest (carboxyl side chain)
	# and determines whether there is a significant difference of that of glu and asp in the TRAP population
	# keyword 'formatted' used to distinguish between the pdbmulti list of objects and atomsbychain_present list of 
	# objects, which is grouped into lists of equivalent atoms (different chains)
	# 'split' keyword used to keep unbound and bound TRAP rings separate (true) or treat as one entity (false) - only
	# available for case of formatted = False

	p = 10 # number of damage datasets to include in analysis

	# Create a figure instance
	sns.set(style="white", context="talk")
	sns.set_context(rc={"figure.figsize": (16, 16), "ytick.labelsize":32})
	f = plt.figure()

	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	keptAtoms = {}
	if split == True:
		keys = ['Bound Glu','Unbound Glu','Bound Asp','Unbound Asp']
	elif split == False:
		keys = ['Glu','Asp']
	for key in keys:
		keptAtoms[key] = {}

	for key in keptAtoms.keys():
		keptAtoms[key]['slope'] = []	
		keptAtoms[key]['r_squared'] = []	
		keptAtoms[key]['p_value'] = []	
		keptAtoms[key]['intercept'] = []	
		keptAtoms[key]['std_err'] = []	

	if formatted == False:
		for atom in pdbmulti:

			if split == True: # for case where unbound and bound rings treated separately
				if atom.basetype in ('GLU') and atom.atomtype == 'CD':
					if atom.chaintype in nonRNAbound:
						key = 'Unbound Glu'
					else:
						key = 'Bound Glu'
				elif atom.basetype in ('ASP') and atom.atomtype == 'CG':
					if atom.chaintype in nonRNAbound:
						key = 'Unbound Asp'
					else:
						key = 'Bound Asp'
				else:
					continue
			elif split == False: # for case where unbound and bound rings treated as one
				if atom.basetype in ('GLU') and atom.atomtype == 'CD':
					key = 'Glu'
				elif atom.basetype in ('ASP') and atom.atomtype == 'CG':
					key = 'Asp' 
				else:
					continue

			# for each set of atoms glu/asp bound/unbound calculate the rate of density loss for each atom in set
			y = np.array(atom.mindensity[0:p-1])
			x = np.array(range(2,p+1))
			slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
			keptAtoms[key]['slope'].append(abs(slope))
			keptAtoms[key]['intercept'].append(intercept)
			keptAtoms[key]['r_squared'].append(r_value**2)
			keptAtoms[key]['p_value'].append(p_value)
			keptAtoms[key]['std_err'].append(std_err)
	else:
		for atomGroup in pdbmulti:
			if (atomGroup[0].atomtype == 'CD' and atomGroup[0].basetype == 'GLU'):
				x = np.array(range(2,p+1))

				key = 'Unbound Glu'
				y_unbound = np.array(np.mean([atom.mindensity[0:p-1] for atom in atomGroup[0:11]],0))
				slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_unbound)
				keptAtoms[key]['slope'].append(slope)
				keptAtoms[key]['intercept'].append(intercept)
				keptAtoms[key]['r_squared'].append(r_value**2)
				keptAtoms[key]['p_value'].append(p_value)
				keptAtoms[key]['std_err'].append(std_err)

				key = 'Bound Glu'
				y_bound = np.array(np.mean([atom.mindensity[0:p-1] for atom in atomGroup[11:22]],0))
				slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_bound)
				keptAtoms[key]['slope'].append(slope)
				keptAtoms[key]['intercept'].append(intercept)
				keptAtoms[key]['r_squared'].append(r_value**2)
				keptAtoms[key]['p_value'].append(p_value)
				keptAtoms[key]['std_err'].append(std_err)

			elif (atomGroup[0].atomtype == 'CG' and atomGroup[0].basetype == 'ASP'):
				x = np.array(range(2,p+1))

				key = 'Unbound Asp'
				y_unbound = np.array(np.mean([atom.mindensity[0:p-1] for atom in atomGroup[0:11]],0))
				slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_unbound)
				keptAtoms[key]['slope'].append(slope)
				keptAtoms[key]['intercept'].append(intercept)
				keptAtoms[key]['r_squared'].append(r_value**2)
				keptAtoms[key]['p_value'].append(p_value)
				keptAtoms[key]['std_err'].append(std_err)

				key = 'Bound Asp'
				y_bound = np.array(np.mean([atom.mindensity[0:p-1] for atom in atomGroup[11:22]],0))
				slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_bound)
				keptAtoms[key]['slope'].append(slope)
				keptAtoms[key]['intercept'].append(intercept)
				keptAtoms[key]['r_squared'].append(r_value**2)
				keptAtoms[key]['p_value'].append(p_value)
				keptAtoms[key]['std_err'].append(std_err)

	for key in keptAtoms.keys():
		print 'Slope:'
		print 'mean: {}: {}'.format(key,np.mean(keptAtoms[key]['slope']))
		print 'std: {}: {}'.format(key,np.std(keptAtoms[key]['slope']))

		print 'p-value:'
		print 'mean: {}: {}'.format(key,np.mean(keptAtoms[key]['p_value']))
		print 'std: {}: {}'.format(key,np.std(keptAtoms[key]['p_value']))

		print 'r-squared:'
		print 'mean: {}: {}'.format(key,np.mean(keptAtoms[key]['r_squared']))
		print 'std: {}: {}'.format(key,np.std(keptAtoms[key]['r_squared']))

		print 'intercept:'
		print 'mean: {}: {}'.format(key,np.mean(keptAtoms[key]['intercept']))
		print 'std: {}: {}'.format(key,np.std(keptAtoms[key]['intercept']))

	categories_means = [np.mean(keptAtoms[key]['slope']) for key in keptAtoms.keys()]
	categories_stds = [np.std(keptAtoms[key]['slope']) for key in keptAtoms.keys()]

	# x axis is discrete number of different equivalent atom groups
	ax = plt.subplot(1, 1, 1)
	width=0.8
	x = np.arange(len(keptAtoms))
	ax.bar(x,categories_means,width,
						 yerr=categories_stds,color='#6897bb',
						 error_kw=dict(elinewidth=4,ecolor='k'))
	xTickNames = [key for key in keptAtoms.keys()]

	# Custom x-label,y-label
	sns.set(style="ticks")
	sns.despine()
	plt.setp(ax.patches,linewidth=5)  
	ax.set_xlim(-width,len(x)+width)           
	plt.xlabel('Residue type', fontsize=40)
	plt.ylabel('|Dloss| (e/A^3)', fontsize=30)
	ax.set_xticks(x+(width/2))
	xtickNames = ax.set_xticklabels(xTickNames)
	plt.setp(xtickNames, rotation=0, fontsize=32)
	f.suptitle('Max density loss slope for glu/asp',fontsize=40)
	f.savefig('maxdenslossSlopeBarplot.png')

def RNAphosphateBaseCorrelation(pdbmulti,datasetnum):
	# this function determines whether there exists a correlation between phosphate group density loss
	# and base density gain
	i = datasetnum

	# get list of RNA nucleotide numbers in pdb file
	RNAindices = []
	for atom in pdbmulti:
		if atom.chaintype == 'W':
			if atom.residuenum not in RNAindices:
				RNAindices.append(atom.residuenum) 

	phosphate_changeList,sugar_changeList,base_changeList = [],[],[]
	for index in RNAindices:
		phosphate,sugar,base = [],[],[]
		for atom in pdbmulti:
			if atom.chaintype == 'W' and atom.residuenum == index:
				if atom.atomtype in ['P','OP1','OP2']:
					phosphate.append(atom)
				elif atom.atomtype in ["C5'","C4'","C3'","C2'","O2'","C1'","O4'","O5'","O3'"]:
					sugar.append(atom)
				else:
					base.append(atom)

		phosphate_change = np.min([atom.mindensity[i] for atom in phosphate])
		sugar_change = np.max([atom.maxdensity[i] for atom in sugar])
		base_change = np.max([atom.maxdensity[i] for atom in base])
		phosphate_changeList.append(phosphate_change)
		sugar_changeList.append(sugar_change)
		base_changeList.append(base_change)

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()
	ax = plt.subplot(1,1,1)
	plt.scatter(phosphate_changeList, sugar_changeList, s=100, c='#d64d4d')
	
	f.text(0.5, 0.04, 'Phosphate Max density loss (e/A^3)', ha='center',fontsize=24)
	f.text(0.04, 0.5, 'Sugar Max density Gain (e/A^3)', va='center', rotation='vertical',fontsize=24)

	f.suptitle('Phosphate Max Density Loss Vs Sugar Max Density Gain: dataset '+str(i+2),fontsize=30)
	plt.setp(f.axes)

	slope, intercept, r_value, p_value, std_err = stats.linregress(phosphate_changeList, sugar_changeList)
	ax.text(0.05, 0.85, '\nR2: '+str("{0:.2f}".format(r_value**2)) +'\np-value: '+str("{0:.2f}\n".format(p_value))
		,color='black',transform=ax.transAxes,fontsize=24)

	f.savefig('phosphateDensLossVsSugarDensGain'+str(i+2)+'.png')


def findSolventAccessibilities(atomList):
	# script to determine correlation between solvent accessiblity and
	# damage susceptibility within the RNA binding interface. Uses new
	# processedAtom class only!

	solvAccessDict = {}
	for boundType in ('unbound','bound'):
		for resNum in (36,42):
			key = '{} Glu {}'.format(boundType,resNum)
			solvAccessDict[key] = {}
			for atomType in ('CD','OE1','OE2'):
				solvAccessDict[key][atomType] = {}
				solvAccessDict[key][atomType]['values'] = []
				for atom in atomList:
					if (atom.atomtype == atomType and atom.residuenum == resNum
						and atom.boundOrUnbound() == '{} protein'.format(boundType)):
						solvAccess = atom.findSolventAccessibility('TRAPRNA1_areaimol1.pdb')
						solvAccessDict[key][atomType]['values'].append(float(solvAccess))

	for boundType in ('unbound','bound'):
		key = '{} Asp {}'.format(boundType,39)
		solvAccessDict[key] = {}
		for atomType in ('CG','OD1','OD2'):
			solvAccessDict[key][atomType] = {}
			solvAccessDict[key][atomType]['values'] = []		
			for atom in atomList:
				if (atom.atomtype == atomType and atom.residuenum == 39
					and atom.boundOrUnbound() == '{} protein'.format(boundType)):
					solvAccess = atom.findSolventAccessibility('TRAPRNA1_areaimol1.pdb')
					solvAccessDict[key][atomType]['values'].append(float(solvAccess))

	for key in solvAccessDict.keys():
		for subkey in solvAccessDict[key].keys():
			solvAccessDict[key][subkey]['mean'] = np.mean(solvAccessDict[key][subkey]['values'])
			solvAccessDict[key][subkey]['std'] = np.std(solvAccessDict[key][subkey]['values'])
			print '{} {} mean: {} std: {}'.format(key,subkey,solvAccessDict[key][subkey]['mean'],solvAccessDict[key][subkey]['std'])
			print solvAccessDict[key][subkey]['values']
	return solvAccessDict

def findSolvAccessDamageCorrelation(atomList,densMet):
	# this function determines whether there is a correlation between damage and solvent accessibility
	# to acidic residues within the RNA binding interfaces around the TRAP ring

	# get full list of atoms to include
	atomIdentifiers = [(resNum,boundType) for resNum in (36,42) for boundType in ('unbound','bound')]
	atomIdentifiers += [(resNum,boundType) for resNum in [39] for boundType in ('unbound','bound')]

	# Create a figure instance
	sns.set_style("white")
	sns.set_context(rc={"figure.figsize": (12, 12)})
	f = plt.figure()
	ax = plt.subplot(1,1,1)

	# colours for plot defined here
	i = 0
	for ind in atomIdentifiers:
		for atom in atomList:
			if atom.residuenum == ind[0] and atom.atomtype == ind[1]:
				i+=1
				break
	jet = cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=i)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

	atomValues = {}
	counter = -1
	scatterPlots = []
	atomDict = {}
	atomValues['Identity'] = []
	for ind in atomIdentifiers:
		print '{} values as follows:'.format(ind)
		atomValues['Identity'].append(ind)
		atomValues['Solvent Access'] = []
		atomValues[densMet] = []
		for atom in atomList:
			if atom.residuenum == ind[0] and atom.atomtype in ('OE1','OE2','OD1','OD2') and atom.boundOrUnbound() =="{} protein".format(ind[1]):
				atomValues['Solvent Access'].append(float(atom.findSolventAccessibility('TRAPRNA1_areaimol1.pdb')))
				atomValues[densMet].append(atom.densMetric[densMet]['Standard']['lin reg']['slope'])
				print atomValues[densMet][-1]
		if len(atomValues['Solvent Access']) != 0:
			counter += 1
			if ind[1] == 'unbound':
				scatterPlot = plt.scatter(atomValues['Solvent Access'], atomValues[densMet], s=100, c=scalarMap.to_rgba(counter),label=ind,marker=str('o'))
			elif ind[1] == 'bound':
				scatterPlot = plt.scatter(atomValues['Solvent Access'], atomValues[densMet], s=100, c=scalarMap.to_rgba(counter),label=ind,marker=str('^'))
			scatterPlots.append(scatterPlot)

			atomDict[ind] = {}
			atomDict[ind]['Solvent Access'] = np.median(atomValues['Solvent Access'])
			atomDict[ind][densMet] = np.median(atomValues[densMet])

	plt.xlabel('Solvent Accessibility', fontsize=18)
	plt.ylabel('D{} slope'.format(densMet), fontsize=18)
	f.suptitle('Solvent Access vs D{}'.format(densMet),fontsize=30)

	# Shrink current axis by 20%
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])

	# Put a legend to the right of the current axis
	ax.legend(handles=scatterPlots,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=18)

	plt.setp(f.axes)
	f.savefig('SolventAccessVsD{}.png'.format(densMet))

	return atomDict,atomValues

def findSolvAccessDamageCorrelationChange(atomList,densMet):

	atomChangeDict = {}
	for keyType in ('Solvent Access',densMet):
		atomChangeDict['{} Change'.format(keyType)] = []
	atomDict,atomValues = findSolvAccessDamageCorrelation(atomList,densMet)
	for key in atomDict.keys():
		for otherkey in atomDict.keys():
			if (key[0:1] == otherkey[0:1] and key[1] == 'unbound' and otherkey[1] == 'bound'):
				atomChangeDict['{} Change'.format('Solvent Access')].append(np.abs(atomDict[key]['Solvent Access'] - atomDict[otherkey]['Solvent Access'])/np.mean([atomDict[key]['Solvent Access'],atomDict[otherkey]['Solvent Access']])) 
				atomChangeDict['{} Change'.format((densMet))].append(np.abs(atomDict[key][densMet] - atomDict[otherkey][densMet])) 

				print '{} unbound solvAcc: {}, bound solvAcc: {}... weighted solvAcc change: {}'.format(key[0:1],atomDict[key]['Solvent Access'],atomDict[otherkey]['Solvent Access'],atomChangeDict['Solvent Access Change'][-1])
				print '..... unbound dens slope: {}, bound dens slope: {}... weighted densmet slope: {}'.format(atomDict[key][densMet],atomDict[otherkey][densMet],atomChangeDict['{} Change'.format(densMet)][-1])
	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (12, 12)})
	f = plt.figure()
	ax = plt.subplot(1,1,1)
	scatterPlot = plt.scatter(atomChangeDict['Solvent Access Change'], atomChangeDict['{} Change'.format(densMet)], s=100)
	plt.xlabel('Solvent Accessibility Change', fontsize=18)
	plt.ylabel('D{} Slope Change'.format(densMet), fontsize=18)
	f.suptitle('Solvent Access Change vs D{} Change'.format(densMet),fontsize=30)
	plt.setp(f.axes)
	f.savefig('SolventAccessChangeVsD{}Change.png'.format(densMet))

def findSolvAccessDamageCorrelationChange2(atomList,densMet):

	# get full list of atoms to include
	# atomIdentifiers = [(resNum,atomType) for resNum in [36] for atomType in ('OE1','OE2')]
	atomIdentifiers = [(resNum,atomType) for resNum in [39] for atomType in ('OD1','OD2')]

	solvChangeList = []
	densChangeList = []
	for ind in atomIdentifiers:
		atomList.atomType 	= ind[1]
		atomList.residueNum = ind[0]
		if atomList.residueNum in (36,42):
			atomList.baseType 	= 'GLU'
		elif atomList.residueNum in [39]:
			atomList.baseType 	= 'ASP'
		atomList.getEquivalentAtoms(True)
		counter = 0
		for atom in atomList.equivAtoms:
			counter + 1
			for otheratom in atomList.equivAtoms[counter:]:
				if atom != otheratom:
					if atom.boundOrUnbound() == 'unbound protein' and otheratom.boundOrUnbound() == 'bound protein':
						print '{} vs {}'.format(atom.chaintype,otheratom.chaintype)
						solvChange = (float(atom.findSolventAccessibility('TRAPRNA1_areaimol1.pdb')) - float(otheratom.findSolventAccessibility('TRAPRNA1_areaimol1.pdb')))/float(atom.findSolventAccessibility('TRAPRNA1_areaimol1.pdb'))
						densChange = atom.densMetric[densMet]['Standard']['lin reg']['slope'] - otheratom.densMetric[densMet]['Standard']['lin reg']['slope']
						solvChangeList.append(solvChange)
						densChangeList.append(densChange)

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (12, 12)})
	f = plt.figure()
	ax = plt.subplot(1,1,1)
	scatterPlot = plt.scatter(solvChangeList,densChangeList, s=100)
	plt.xlabel('Solvent Accessibility Change', fontsize=18)
	plt.ylabel('D{} Slope Change'.format(densMet), fontsize=18)
	f.suptitle('Solvent Access Change vs D{} Change'.format(densMet),fontsize=30)
	plt.setp(f.axes)
	f.savefig('SolventAccessChangeVsD{}Change.png'.format(densMet))

def GLUASPDnetCorrel(PDBmulti,i,xWhat,yWhat,fitType,xDensMet,xNormType,yDensMet,yNormType):
	# this function is designed to determine a correlation between the extent of damage for
	# both carboxyl group side chain oxygens and carbon. The function plots a scatter plot of 
	# xWhat damage vs yWhat damage for each dataset damage level. xWhat and yWhat specify what
	# x and y scatter axes should be - takes values O1,O2,C (corresponding to GLUOE1/ASPOD1,
	# GLUOE2/ASPOE2, GLUCD/ASPCG). The density metric here is specified by xDensMet and yDensMet. 
	# Function can plot 2 separate types of correlation plots ([fitType == 'values', i = dataset
	# number] to plot for a distinct dataset number, and [fitType == 'lin reg', i = 'slope').
	# xNormType and yNormType specify whether the x and y axis density metrics should be normalised
	# or not (takes values 'Standard' or 'Calpha normalised')

	# first find all GLU/ASP carboxyl side chain oxy atoms
	keptAtoms = []
	for atom in PDBmulti:
		if atom.basetype in ('GLU','ASP'):
			if atom.atomtype in ('OE1','OE2','OD1','OD2'):
				keptAtoms.append(atom)

	# next find all GLU/ASP carboxyl side chain end carbon atoms
	for atom in PDBmulti:
		if atom.basetype in ('GLU') and atom.atomtype == 'CD':
			keptAtoms.append(atom)
		elif atom.basetype in ('ASP') and atom.atomtype == 'CG':
			keptAtoms.append(atom)

	# next create a gluOE1, gluOE2 and gluCD list such that GLU OE1 atoms in gluOE1 list
	# match to GLU OE2 atoms in gluOE2 list and GLU CD atoms in gluCD list
	gluOE1,gluOE2,gluCD=[],[],[]
	for atom in keptAtoms:
		if atom.basetype in 'GLU':
			if atom.atomtype in ('OE1'):
				gluOE1.append(atom)
				for otheratom in keptAtoms:
					if (otheratom.chaintype == atom.chaintype and
						otheratom.residuenum == atom.residuenum): 
						if otheratom.atomtype == 'OE2':
							gluOE2.append(otheratom)
						elif otheratom.atomtype == 'CD':
							gluCD.append(otheratom)
	# check that three lists are same length 
	if (len(gluOE2) != len(gluOE1)) or (len(gluOE2) != len(gluCD)):
		print 'something went wrong!'
		sys.exit()

	# same again for the ASP OD1, OD2 and CG atoms
	aspOD1,aspOD2,aspCG=[],[],[]
	for atom in keptAtoms:
		if atom.basetype in 'ASP':
			if atom.atomtype in ('OD1'):
				aspOD1.append(atom)
				for otheratom in keptAtoms:
					if (otheratom.chaintype == atom.chaintype and
						otheratom.residuenum == atom.residuenum): 
						if otheratom.atomtype == 'OD2':
							aspOD2.append(otheratom)
						elif otheratom.atomtype == 'CG':
							aspCG.append(otheratom)
	# check that three lists are same length 
	if (len(aspOD2) != len(aspOD1)) or (len(aspOD2) != len(aspCG)):
		print 'something went wrong!'
		sys.exit()

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()

	# want to create x and y lists to plot in scatter plot, GLU first
	if xWhat == 'O1':
		x = [atom.densMetric[xDensMet][xNormType][fitType][i] for atom in gluOE1]
	elif xWhat == 'O2':
		x = [atom.densMetric[xDensMet][xNormType][fitType][i] for atom in gluOE2]
	elif xWhat == 'C':
		x = [atom.densMetric[xDensMet][xNormType][fitType][i] for atom in gluCD]
	if yWhat == 'O1':
		y = [atom.densMetric[yDensMet][yNormType][fitType][i] for atom in gluOE1]
	elif yWhat == 'O2':
		y = [atom.densMetric[yDensMet][yNormType][fitType][i] for atom in gluOE2]
	elif yWhat == 'C':
		y = [atom.densMetric[yDensMet][yNormType][fitType][i] for atom in gluCD]

	gluScatter = plt.scatter(x, y, s=100, c='#d64525', alpha=0.5)

	# calculate linear regression for GLU points here
	slope, intercept, r_valueGLU, p_value, std_err = stats.linregress(x, y)
	# x2 = range(np.floor(min(x)),np.ceil(max(x))+1)
	x2 = [min(x), min(x)+ (max(x) - min(x))/2, max(x)]
	y2 = slope*np.array(x2) + np.array(intercept)
	plt.plot(x2, y2, '-',color='#d64525',linewidth=3)

	# want to create x and y lists to plot in scatter plot, ASP next
	if xWhat == 'O1':
		x = [atom.densMetric[xDensMet][xNormType][fitType][i] for atom in aspOD1]
	elif xWhat == 'O2':
		x = [atom.densMetric[xDensMet][xNormType][fitType][i] for atom in aspOD2]
	elif xWhat == 'C':
		x = [atom.densMetric[xDensMet][xNormType][fitType][i] for atom in aspCG]
	if yWhat == 'O1':
		y = [atom.densMetric[yDensMet][yNormType][fitType][i] for atom in aspOD1]
	elif yWhat == 'O2':
		y = [atom.densMetric[yDensMet][yNormType][fitType][i] for atom in aspOD2]
	elif yWhat == 'C':
		y = [atom.densMetric[yDensMet][yNormType][fitType][i] for atom in aspCG]

	aspScatter = plt.scatter(x, y, s=100, c='#35b9ac', alpha=0.5,marker='d')

	# calculate linear regression for ASP points here
	slope, intercept, r_valueASP, p_value, std_err = stats.linregress(x, y)
	x2 = [min(x), min(x)+ (max(x) - min(x))/2, max(x)]
	y2 = slope*np.array(x2) + np.array(intercept)
	plt.plot(x2, y2, '-',color='#35b9ac',linewidth=3)

	plt.legend([gluScatter,aspScatter],['GLU, R^2: '+str("{0:.2f}".format(r_valueGLU**2)),
										'ASP, R^2: '+str("{0:.2f}".format(r_valueASP**2))],
			   loc=4, borderaxespad=0.,fontsize=24)

	plt.xlabel('{} {} D{}'.format(xWhat,xNormType,xDensMet),fontsize=24)
	plt.ylabel('{} {} D{}'.format(yWhat,yNormType,yDensMet),fontsize=24)
	f.suptitle('GLU/ASP {} {} D{} vs {} {} D{} {}'.format(xWhat,xNormType,xDensMet,yWhat,yNormType,yDensMet,str(i)),fontsize=24)
	plt.setp(f.axes)
	f.savefig('GLUASP_{}_{}D{}_vs_{}_{}D{}_{}.png'.format(xWhat,xNormType.replace(" ",""),xDensMet,yWhat,yNormType.replace(" ",""),yDensMet,str(i)))

def glu50tyr62Correlation(atomList,densMet,metType,valType):
	# this function is designed to determine whether there exists a correlation between glu 50 (bound) damage
	# and Tyr 62 (unbound) damage - since these two residues are hydrogen bonded.
	damageDict = {}
	damageDict['Glu 50'] = []
	damageDict['Tyr 62'] = []
	tyr62ChainOrder = ['A','B','C','D','E','F','G','H','I','J','K']
	glu50ChainOrder = ['U','T','S','R','Q','P','O','N','M','L','V']
	for atom in atomList:
		if atom.atomtype in ['OE2'] and atom.residuenum == 50 and atom.boundOrUnbound() == 'bound protein':
			for otheratom in atomList:
				if otheratom.atomtype in ['OH'] and otheratom.residuenum == 62 and otheratom.boundOrUnbound() == 'unbound protein':
					if glu50ChainOrder.index(atom.chaintype) == tyr62ChainOrder.index(otheratom.chaintype):
						damageDict['Glu 50'].append(atom.densMetric[densMet]['Standard'][metType][valType])
						damageDict['Tyr 62'].append(otheratom.densMetric[densMet]['Standard'][metType][valType])
						break

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()
	plt.scatter(damageDict['Glu 50'], damageDict['Tyr 62'], s=100, c='#35b9ac')
	plt.xlabel('Glu 50 OE2 D{} {} {}'.format(densMet,metType,valType),fontsize=24)
	plt.ylabel('Tyr 62 OH D{} {} {}'.format(densMet,metType,valType),fontsize=24)
	f.suptitle('Glu 50 OE2 vs Tyr 62 OH D{} {} {}'.format(densMet,metType,valType),fontsize=24)
	plt.setp(f.axes)
	f.savefig('Glu50OE2vsTyr62OH_D{}_{}_{}.png'.format(densMet,metType,valType))


def asp39G1distCorrelation(atomList,densMet,metType,valType,O_C_difference):
	# this function is designed to determine whether there exists a correlation between G1 N2 density change
	# and the distance between asp39 side-chain oxygens. If 'O_C_difference' is True then Dloss difference between
	# O and C atoms in side-chain carboxyl groups plotted, otherwise O Dloss plotted only

	# Create a figure instance
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (16, 16)})
	f = plt.figure()
	ax = plt.subplot(1,1,1)

	atomresinum = [1,1,3,3]
	atomtype = ['N1','N2','N1','N2']
	otheratomtype = ['OD1','OD2','OE2','OE1']
	otheratomresinum = [39,39,36,36]
	otherotheratomtype = ['CG','CG','CD','CD']
	plotcolor = ['red','orange','green','blue']
	combinedlists = {'x':[],'y':[]}
	for i in range(0,4):
		plotDict = {'CO2 damage':[],'G:N1 damage':[],'distance':[],'O minus C damage':[]}
		for atom in atomList:
			if atom.atomtype == atomtype[i] and atom.residuenum%5 == atomresinum[i] and atom.boundOrUnbound() == 'rna':
				atomPosition = np.array([atom.X_coord,atom.Y_coord,atom.Z_coord])
				plotDict['G:N1 damage'].append(atom.densMetric[densMet]['Standard'][metType][valType])
				distList = []
				searchedAtoms = []
				for otheratom in atomList:
					if otheratom.atomtype == otheratomtype[i] and otheratom.residuenum == otheratomresinum[i] and otheratom.boundOrUnbound() == 'bound protein':
						otheratomPosition = np.array([otheratom.X_coord,otheratom.Y_coord,otheratom.Z_coord])
						distList.append(np.linalg.norm(atomPosition-otheratomPosition))
						searchedAtoms.append(otheratom)
				plotDict['distance'].append(min(distList))
				plotDict['CO2 damage'].append(searchedAtoms[distList.index(min(distList))].densMetric[densMet]['Standard'][metType][valType]) 
				
				# calculate the difference in Dloss between the oxygen and carbon in located Asp39 CO2 group
				if O_C_difference == True:
					for otherotheratom in atomList:
						if otherotheratom.atomtype == otherotheratomtype[i] and otherotheratom.residuenum == otheratomresinum[i] and otherotheratom.chaintype == searchedAtoms[distList.index(min(distList))].chaintype:
							densMetDifference = searchedAtoms[distList.index(min(distList))].densMetric[densMet]['Standard'][metType][valType] - otherotheratom.densMetric[densMet]['Standard'][metType][valType]
							plotDict['O minus C damage'].append(densMetDifference)
		if O_C_difference == True: 
			yKey = 'O minus C damage'
			yLabel = 'Closest oxygen D{} difference {} {}'.format(densMet,metType,valType)
		else:
			yKey = 'CO2 damage'
			yLabel = 'Closest oxygen D{} {} {}'.format(densMet,metType,valType)
		plt.scatter(plotDict['distance'], plotDict[yKey], s=200, c=plotcolor[i])
		combinedlists['x'].append(plotDict['distance'])
		combinedlists['y'].append(plotDict[yKey])

	# calculate linear regression values here
	xlist_1d = [x for sublist in combinedlists['x'] for x in sublist]
	ylist_1d = [x for sublist in combinedlists['y'] for x in sublist]
	slope, intercept, r_value, p_value, std_err = stats.linregress(xlist_1d, ylist_1d)
	x2 = [min(xlist_1d), min(xlist_1d)+ (max(xlist_1d) - min(xlist_1d))/2, max(xlist_1d)]
	y2 = slope*np.array(x2) + np.array(intercept)
	plt.plot(x2, y2, '-',color='grey',linewidth=3)

	plt.xlabel('G1:N1 distance',fontsize=24)
	plt.ylabel(yLabel,fontsize=24)
	ax.tick_params(axis='x', labelsize=30)
	ax.tick_params(axis='y', labelsize=30)
	f.suptitle('CO2 damage vs carboxyl proximity: D{} {} {}'.format(densMet,metType,valType),fontsize=24)
	plt.setp(f.axes)
	f.savefig('CO2DamagevCarboxylProximityD{}{}{}_{}.png'.format(densMet,metType,valType,O_C_difference))

	# output linear regression summary
	print '\n***\nLin regression summary:\nSlope: {}\nIntercept: {}\n R-squared: {}\n P-value: {}\nStd Error: {}'.format(slope,intercept,r_value**2,p_value,std_err)
	return slope,r_value

def batchRun_asp39G1distCorrelation(atomList):

	summaryDict = {}
	for densMet in ['loss','bfactor']:
		summaryDict[densMet] = {'slope':[],'r_squared':[]}
		for i in range(0,9):
			slope,r_value = asp39G1distCorrelation(atomList,densMet,'values',i,False)
			summaryDict[densMet]['slope'].append(slope)
			summaryDict[densMet]['r_squared'].append(r_value**2)

	for key in summaryDict.keys():
		print '\n***\nD{} metric: '.format(key)	
		print 'Slope: {}'.format(','.join([str(val) for val in summaryDict[key]['slope']]))
		print 'R-squared: {}'.format(','.join([str(val) for val in summaryDict[key]['r_squared']]))

def gluAspBfactorChangeOnRNABinding(atomList,densMet,percent):
	# determine the change in Bfactor upon RNA binding - which residues are most effected?

	# Create a figure instance
	sns.set(style="white", context="talk")
	sns.set_context(rc={"figure.figsize": (12, 12)})
	sns.despine()
	f = plt.figure()
	ax = plt.subplot(1,1,1)

	gluList = [(aType,'GLU',rType) for aType in ['CD','OE1','OE2'] for rType in [16,42,50,71,73]]
	aspList = [(aType,'ASP',rType) for aType in ['CG','OD1','OD2'] for rType in [8,17,29,39]]
	# find all glu and asp side-chain carboxyl atoms

	counter = -1
	colourList =['#00aedb','#d11141']
	for idSet in (gluList+aspList,[(aType,'GLU',36) for aType in ['CD','OE1','OE2']]):
		atomDict = {}
		atomDict['bfactor'],atomDict['slope'] = [],[]
		atomDict['percent bfactor'],atomDict['percent slope'] = [],[]
		counter += 1
		for id in idSet:
			atomList.atomType = id[0]
			atomList.baseType = id[1]
			atomList.residueNum = id[2]
			atomList.getEquivalentAtoms(True)
			foundAtoms = atomList.equivAtoms
			unboundList_bfac,boundList_bfac=[],[]
			unboundList_dens,boundList_dens=[],[]

			# sort atoms into bound and unbound groupings
			for atom in foundAtoms:
				if atom.boundOrUnbound() == 'unbound protein':
					unboundList_bfac.append(float(atom.Bfactor[0]))
					unboundList_dens.append(atom.densMetric[densMet]['Standard']['lin reg']['slope'])
				elif atom.boundOrUnbound() == 'bound protein':
					boundList_bfac.append(float(atom.Bfactor[0]))
					boundList_dens.append(atom.densMetric[densMet]['Standard']['lin reg']['slope'])

			atomDict['bfactor'].append(np.mean(unboundList_bfac) - np.mean(boundList_bfac))
			atomDict['percent bfactor'].append(atomDict['bfactor'][-1]/np.mean(unboundList_bfac))
			atomDict['slope'].append(np.mean(unboundList_dens) - np.mean(boundList_dens))
			atomDict['percent slope'].append(atomDict['slope'][-1]/np.mean(unboundList_dens))

		if percent == True:
			scatterPlot = plt.scatter(atomDict['percent bfactor'], atomDict['percent slope'], s=200,c=colourList[counter])

			# linearly fit to the non-glu36 points
			if counter == 0:
				slope, intercept, r_value, p_value, std_err = stats.linregress(atomDict['percent bfactor'],atomDict['percent slope'])
				x = [min(atomDict['percent bfactor'])*1.1,max(atomDict['percent bfactor'])*1.1]
				y = slope*np.array(x) + np.array(intercept)
				plt.plot(x, y, '-',color='#4b86b4',linewidth=4)
				print 'R^2: {} pval: {} stdErr: {}'.format(r_value**2,p_value,std_err)

		else:
			scatterPlot = plt.scatter(atomDict['bfactor'], atomDict['slope'], s=200,c=colourList[counter])

			# linearly fit to the non-glu36 points
			if counter == 0:
				slope, intercept, r_value, p_value, std_err = stats.linregress(atomDict['bfactor'],atomDict['slope'])
				x = [min(atomDict['bfactor'])*1.1,max(atomDict['bfactor'])*1.1]
				y = slope*np.array(x) + np.array(intercept)
				plt.plot(x, y, '-',color='#4b86b4',linewidth=4)
				print 'R^2: {} pval: {} stdErr: {}'.format(r_value**2,p_value,std_err)

	if percent == True:
		plt.xlabel('Percent Bfactor Change', fontsize=24)
		plt.ylabel('Percent D{} Slope Change'.format(densMet), fontsize=24)
		f.suptitle('Bfactor change vs D{} Change'.format(densMet),fontsize=30)
		plt.setp(f.axes)
		f.savefig('PERCENTbfacChangeVsD{}Change.png'.format(densMet))
	else:
		plt.xlabel('Bfactor Change', fontsize=24)
		plt.ylabel('D{} Slope Change'.format(densMet), fontsize=24)
		f.suptitle('Bfactor change vs D{} Change'.format(densMet),fontsize=30)
		plt.setp(f.axes)
		f.savefig('bfacChangeVsD{}Change.png'.format(densMet))


def meanOfDensitiesPerComponent():
	# small function to calculate the mean of the Dloss metric for all protein
	# RNA and solvent atoms separately and observe if they are different for each dose 
	atoms = retrieveAtomList()
	densDict = {'protein':[],'RNA':[],'solvent':[]}
	for i in range(0,10):
		print '***\nFor damage set {}:'.format(i)
		for atom in atoms.processedAtomList:
			density = atom.densMetric['loss']['Standard']['values'][i]
			if atom.chaintype not in ['W','Y']:
				densDict['protein'].append(density)
			elif atom.chaintype in ['W']:
				densDict['RNA'].append(density)
			elif atom.chaintype in ['Y']:
				densDict['solvent'].append(density)
		# calculate mean of densities, as well as std
		for key in densDict.keys():
			meanOfDensities = np.mean(densDict[key])
			stdOfDensities = np.std(densDict[key])
			print '{} ---> mean: {} std: {}'.format(key,meanOfDensities,stdOfDensities)


def asp39G1_glu36G3dist(atomList):
	# this function is designed to determine calculate the min atom distance between glu36 and G3 and also
	# asp39 and G1

	distPairs = [[39,'OD1',1,'N1'],[39,'OD2',1,'N2'],[39,'OD1',1,'N2'],[39,'OD2',1,'N1'],
			     [36,'OE1',3,'N1'],[36,'OE2',3,'N2'],[36,'OE1',3,'N2'],[36,'OE2',3,'N1']]

	for distPair in distPairs:
		distPairList = []
		for atom in atomList:
			if atom.atomtype == distPair[3] and atom.residuenum%5 == distPair[2] and atom.boundOrUnbound() == 'rna':
				atomPosition = np.array([atom.X_coord,atom.Y_coord,atom.Z_coord])
				distList = []
				searchedAtoms = []
				for otheratom in atomList:
					if otheratom.atomtype == distPair[1] and otheratom.residuenum == distPair[0] and otheratom.boundOrUnbound() == 'bound protein':
						otheratomPosition = np.array([otheratom.X_coord,otheratom.Y_coord,otheratom.Z_coord])
						distList.append(np.linalg.norm(atomPosition-otheratomPosition))
						searchedAtoms.append(otheratom)
				if distPair[0] == 39:
					print 'asp39-G1 {}-{} distance: {}'.format(distPair[1],distPair[3],min(distList))
				elif distPair[0] == 36:
					print 'glu36-G3 {}-{} distance: {}'.format(distPair[1],distPair[3],min(distList))
				distPairList.append(min(distList))
		print 'Average distance: {}\n***'.format(np.mean(distPairList))



def nonboundVboundDensityStats(PDBmulti):
	# calculate basic density summary stats for bound, non-bound and RNA

	densityDict = {'unbound protein':[],'bound protein':[],'rna':[]}
	for atom in PDBmulti:
		boundtype = atom.boundOrUnbound()
		if boundtype != 'something else':
			densityDict[boundtype].append(atom.densMetric['loss']['Standard']['values'])

	# calculate stats over each component with each dose
	csvOut = open('boundnonboundstats.csv','w')
	for i in range(0,10):
		csvOut.write('DATASET {}\n'.format(i))
		print 'For dataset {}'.format(i)
		for key in densityDict.keys():
			densVals = [dens[i] for dens in densityDict[key]]
			meanDloss = np.mean(densVals)
			stdDloss = np.std(densVals)
			ntile_10 = np.percentile(densVals,10)
			ntile_90 = np.percentile(densVals,90)
			confint_95 = mean_confidence_interval(densVals)

			print '\tFor {} component:'.format(key)
			print '\t\tMean for Dloss: {}'.format(meanDloss)
			print '\t\tSD for Dloss: {}'.format(stdDloss)
			print '\t\t10-tile for Dloss: {}'.format(ntile_10)
			print '\t\t90-tile for Dloss: {}'.format(ntile_90)
			print '\t\t95CI for Dloss: {}'.format(confint_95)

			csvOut.write('{},{},{}\n'.format(key,meanDloss,confint_95))
	csvOut.close()

def DlossHistogram(atoms,Dmetric):
	# plot a histogram of the distribution of Dloss values within a structure
	# for a given dataset number. Dmetric specifies density metric

	testStatistics = {'unbound protein':[],'bound protein':[],'rna':[]}
	skewnessDict = {'unbound protein':[],'bound protein':[],'rna':[]}
	for i in range(0,9):
		sns.set_palette("deep", desat=.6)
		sns.set_context(rc={"figure.figsize": (10, 6)})
		fig = plt.figure()

		counter = -1 
		colors = ['red','blue','green']
		for boundType in ['unbound protein','bound protein','rna']:
			counter += 1
			colorType = colors[counter]
			datax = [atm.densMetric[Dmetric]['Standard']['values'][i] for atm in atoms.processedAtomList if atm.boundOrUnbound() == boundType]
			plt.hist(datax, 300, histtype="stepfilled", alpha=.7,color=colorType)
			print '-------Dataset {}---------'.format(i+2)
			print "Component '{}' length: {}".format(boundType,len(datax))

			# test for normality in data
			print 'Testing for normality in dataset: {}'.format(boundType)
			W,critVal,SigLevel = stats.anderson(datax,dist='norm')
			print 'Test statistic: {}'.format(W)
			print 'critical value: {}'.format(critVal)
			print 'SigLevel: {}'.format(SigLevel)
			testStatistics[boundType].append(W)

			# calculate dataset skewness
			skewness = stats.skew(datax)
			print 'Calculating skewness in dataset: {}'.format(boundType)
			print 'Skewness: {}'.format(skewness)
			skewnessDict[boundType].append(skewness)

		plt.xlabel('D{} value'.format(Dmetric))
		plt.ylabel('Frequency of atoms')
		plt.title('Histrogram of D{} value per atom'.format(Dmetric))
		fig.savefig('D{}PerAtom_{}.png'.format(Dmetric,i))

	# plot line graph of dataset number against test statistics and skewness
	for i in range(0,2):
		if i == 0:
			test = testStatistics
			yLabel = 'A-D test statistic'
			title = 'Anderson-Darling test stat for D{}'.format(Dmetric)
			fileName = 'AndDarlStat_D{}.png'.format(Dmetric)
		else:
			test = skewnessDict
			yLabel = 'Skewness'
			title = 'Skewness for D{}'.format(Dmetric)
			fileName = 'Skewness_D{}.png'.format(Dmetric)

		sns.set_palette("deep", desat=.6)
		sns.set_context(rc={"figure.figsize": (10, 6)})
		fig = plt.figure()
		ax = plt.subplot(1,1,1)
		for key in test.keys():
			if key == 'unbound protein':
				label = 'Non-bound TRAP'
			elif key == 'bound protein':
				label = 'Bound TRAP'
			elif key == 'rna':
				label = 'RNA'
			plt.plot(range(2,11),test[key],label=label)
			ax.legend(loc='best')
			plt.xlabel('Dataset', fontsize=18)
			plt.ylabel(yLabel, fontsize=18)
			plt.title(title)
			fig.savefig(fileName)


def DlossHistogram_RNAcomponents(atoms,Dmetric):
	# plot a histogram of the distribution of Dloss values within a structure
	# for a given dataset number. Dmetric specifies density metric

	testStatistics = {'base':[],'sugar':[],'phosphate':[]}
	for i in range(0,9):
		sns.set_palette("deep", desat=.6)
		sns.set_context(rc={"figure.figsize": (10, 6)})
		fig = plt.figure()

		counter = -1 
		colors = ['red','blue','green']
		datax = {}
		for rnaType in ['base','sugar','phosphate']:

			if rnaType == 'phosphate':
				atomtypes = ['P','OP1','OP2']
			elif rnaType == 'sugar':
				atomtypes = ["C5'","C4'","C3'","C2'","O2'","C1'","O4'","O5'","O3'"]
			elif rnaType == 'base':
				atomtypes = ['P','OP1','OP2'] + ["C5'","C4'","C3'","C2'","O2'","C1'","O4'","O5'","O3'"]

			counter += 1
			colorType = colors[counter]
			datax[rnaType] = []
			if rnaType in ('sugar','phosphate'):
				for atm in atoms.processedAtomList:
					if atm.atomtype in atomtypes and atm.chaintype == 'W':
						datax[rnaType].append(atm.densMetric[Dmetric]['Standard']['values'][i])
			elif rnaType in ('base'):
				for atm in atoms.processedAtomList:
					if atm.atomtype not in atomtypes and atm.chaintype == 'W':
						datax[rnaType].append(atm.densMetric[Dmetric]['Standard']['values'][i])
			plt.hist(datax[rnaType], 300, histtype="stepfilled", alpha=.7,color=colorType)

			# test for normality in data
			print 'Testing for normality in dataset: {}'.format(rnaType)
			W,critVal,SigLevel = stats.anderson(datax[rnaType],dist='norm')
			print 'Test statistic: {}'.format(W)
			print 'critical value: {}'.format(critVal)
			print 'SigLevel: {}'.format(SigLevel)
			testStatistics[rnaType].append(W)

		plt.xlabel('D{} value'.format(Dmetric))
		plt.ylabel('Frequency of atoms')
		plt.title('Histrogram of D{} value per atom'.format(Dmetric))
		fig.savefig('RNAatoms_D{}PerAtom_{}.png'.format(Dmetric,i))

	# plot line graph of dataset number against test statistics
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (10, 6)})
	fig = plt.figure()
	ax = plt.subplot(1,1,1)
	for key in testStatistics.keys():
		plt.plot(range(2,11),testStatistics[key],label=key)
		ax.legend(loc='best')
		plt.xlabel('Dataset', fontsize=18)
		plt.ylabel('A-D test statistic'.format(Dmetric), fontsize=18)
		plt.title('Anderson-Darling test stat for D{}'.format(Dmetric))
		fig.savefig('RNAatoms_AndDarlStat_D{}.png'.format(Dmetric))


def DlossHistogram_RNAcomponents2(atoms,Dmetric,rnaType):
	# plot a histogram of the distribution of Dloss values within a structure
	# for a given dataset number. Dmetric specifies density metric

	testStatistics = {'G1':[],'A2':[],'G3':[],'U4':[]}
	for i in range(0,9):
		sns.set_palette("deep", desat=.6)
		sns.set_context(rc={"figure.figsize": (10, 6)})
		fig = plt.figure()

		counter = -1 
		colors = ['red','blue','green','yellow']
		datax = {}
		for nucleotideType in ['G1','A2','G3','U4']:

			if rnaType == 'phosphate':
				atomtypes = ['P','OP1','OP2']
			elif rnaType == 'sugar':
				atomtypes = ["C5'","C4'","C3'","C2'","O2'","C1'","O4'","O5'","O3'"]
			elif rnaType == 'base':
				atomtypes = ['P','OP1','OP2'] + ["C5'","C4'","C3'","C2'","O2'","C1'","O4'","O5'","O3'"]

			counter += 1
			colorType = colors[counter]
			datax[nucleotideType] = []
			if rnaType in ('sugar','phosphate'):
				for atm in atoms.processedAtomList:
					if atm.atomtype in atomtypes and atm.chaintype == 'W' and atm.residuenum%5 == int(nucleotideType[1]):
						datax[nucleotideType].append(atm.densMetric[Dmetric]['Standard']['values'][i])
			elif rnaType in ('base'):
				for atm in atoms.processedAtomList:
					if atm.atomtype not in atomtypes and atm.chaintype == 'W' and atm.residuenum%5 == int(nucleotideType[1]):
						datax[nucleotideType].append(atm.densMetric[Dmetric]['Standard']['values'][i])
			plt.hist(datax[nucleotideType], 300, histtype="stepfilled", alpha=.7,color=colorType)

			# test for normality in data
			print 'Testing for normality in dataset: {}'.format(nucleotideType)
			W,critVal,SigLevel = stats.anderson(datax[nucleotideType],dist='norm')
			print 'Test statistic: {}'.format(W)
			print 'critical value: {}'.format(critVal)
			print 'SigLevel: {}'.format(SigLevel)
			testStatistics[nucleotideType].append(W)

		plt.xlabel('D{} value'.format(Dmetric))
		plt.ylabel('Frequency of atoms')
		plt.title('Histrogram of D{} value per atom'.format(Dmetric))
		fig.savefig('RNAatoms_D{}PerAtom_{}.png'.format(Dmetric,i))

	# plot line graph of dataset number against test statistics
	sns.set_palette("deep", desat=.6)
	sns.set_context(rc={"figure.figsize": (10, 6)})
	fig = plt.figure()
	ax = plt.subplot(1,1,1)
	for key in testStatistics.keys():
		plt.plot(range(2,11),testStatistics[key],label=key)
		ax.legend(loc='best')
		plt.xlabel('Dataset', fontsize=18)
		plt.ylabel('A-D test statistic'.format(Dmetric), fontsize=18)
		plt.title('Anderson-Darling test stat for D{}'.format(Dmetric))
		fig.savefig('RNAatoms_AndDarlStat_D{}.png'.format(Dmetric))


def glu36OE2toG3O6Dist(PDBmulti):
	# function to calculate distances from Glu36 OE2 atom to O6 oxygen of 
	# G3 nucleotide - is it less than 3 Angstrom?
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']

	# calculate Glu36OE2-G3O6 and ASP39OD1-G1O6 distances for each protein chain of bound TRAP
	for choice in ([36,'OE2'],[39,'OD1']):
		distDict = {}
		for atom in PDBmulti:
			distList = []
			if atom.chaintype in RNAbound and atom.residuenum == choice[0] and atom.atomtype == choice[1]:
				for otheratom in PDBmulti:
					if otheratom.chaintype == 'W' and otheratom.atomtype == 'O6':
						atomPosition = np.array([atom.X_coord,atom.Y_coord,atom.Z_coord])
						otheratomPosition = np.array([otheratom.X_coord,otheratom.Y_coord,otheratom.Z_coord])
						distList.append(np.linalg.norm(atomPosition-otheratomPosition))
				distDict[atom.chaintype] = min(distList)
		
		# provide summary info here
		print 'Residue {}, atom {}'.format(choice[0],choice[1])
		print 'Min distance: {}'.format(min(distDict.values()))
		print 'Max distance: {}'.format(max(distDict.values()))
		print 'Mean distance: {}'.format(np.mean(distDict.values()))
		print 'Std for distances: {}\n'.format(np.std(distDict.values()))




