# -*- coding: utf-8 -*-
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def bdamchange_v_atomnum(where,PDBmulti):
	# funtion to plot Bdamage change as function of atom number
	# with each dataset number on a separate axis

	sns.set(style="white", context="talk")
	f, axes = plt.subplots(len(PDBmulti[0].meandensity), 1, figsize=(16, 20), sharex=True)

	PDBmulti.sort(key=lambda x: x.atomnum)

	for i in range(0,len(PDBmulti[0].meandensity)):
		y = [atom.bdamchange[i] for atom in PDBmulti]
		x = [atom.atomnum for atom in PDBmulti]

		# add line plot for this dataset 
		axes[i].plot(x, y)
		axes[i].set_ylabel('data '+str(i+1))

	plt.xlabel('atom number', fontsize=18)
	f.suptitle('Bdamage change vs atom number',fontsize=20)
	sns.despine(bottom=True)
	plt.setp(f.axes)
	plt.tight_layout(h_pad=3)
	plt.subplots_adjust(top=0.95)
	f.savefig(where+'Bdamage_change_vs_atom_number.png')

def bdamBfac_change_v_atomnum(where,PDBmulti):
	# function to plot Bdamage change and Bfactor change as function of 
	# atom number on same plot

	sns.set(style="white", context="talk")
	f, axes = plt.subplots(2, 1, figsize=(16, 20), sharex=True)

	PDBmulti.sort(key=lambda x: x.atomnum)

	y = [atom.bdamchange for atom in PDBmulti]
	x = [atom.atomnum for atom in PDBmulti]

	# add line plot for this dataset 
	axes[0].plot(x, y)
	axes[0].set_ylabel('Bdamage change')

	y = [atom.Bfactorchange for atom in PDBmulti]

	# add line plot for this dataset 
	axes[1].plot(x, y)
	axes[1].set_ylabel('Bfactor change')


	plt.xlabel('atom number', fontsize=18)
	f.suptitle('Bdamage/Bfactor change vs atom number',fontsize=20)
	sns.despine(bottom=True)
	plt.setp(f.axes)
	plt.tight_layout(h_pad=3)
	plt.subplots_adjust(top=0.95)
	f.savefig(where+'BdamageBfactor_change_vs_atom_number.png')

