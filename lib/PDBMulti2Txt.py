# -*- coding: utf-8 -*-

# function here is designed to read in a list of objects merged over multiple datasets
# (PDBmulti) and write a txt file containing the density metrics for each atom in the 
# structure, over the different datasets
def objlist2txt(PDBmulti,where,densmet):

	# make a new txt file in location where
	txtfile = open(where+'PDBmulti_'+str(densmet)+'.txt','w')
	txtfile.write('atomnum\t')
	txtfile.write('atominfo\t')
	for i in range(0,len(PDBmulti[0].mindensity)):
		txtfile.write('%s_%s\t' %(str(densmet),str(i+1)))
	txtfile.write('\n')

	# sort by atomnumber
	PDBmulti.sort(key=lambda x: x.atomnum)

	# for each atom in PDBmulti
	for atom in PDBmulti:
		txtfile.write('%s\t' %(str(atom.atomnum)))
		txtfile.write('%s%s%s%s\t' %(str(atom.atomtype),str(atom.residuenum),str(atom.basetype),str(atom.chaintype)))

		if densmet in ('mean','Mean'):
			for i in range(0,len(PDBmulti[0].meandensity)):
				txtfile.write('%s\t' %(str(atom.meandensity[i])))

		elif densmet in ('median','Median'):
			for i in range(0,len(PDBmulti[0].mediandensity)):
				txtfile.write('%s\t' %(str(atom.mediandensity[i])))

		elif densmet in ('min','Min'):
			for i in range(0,len(PDBmulti[0].mindensity)):
				txtfile.write('%s\t' %(str(atom.mindensity[i])))

		elif densmet in ('max','Max'):
			for i in range(0,len(PDBmulti[0].maxdensity)):
				txtfile.write('%s\t' %(str(atom.maxdensity[i])))
				
		elif densmet in ('std','Std'):
			for i in range(0,len(PDBmulti[0].stddensity)):
				txtfile.write('%s\t' %(str(atom.stddensity[i])))

		elif densmet in ('min90tile'):
			for i in range(0,len(PDBmulti[0].min90tile)):
				txtfile.write('%s\t' %(str(atom.min90tile[i])))

		elif densmet in ('max90tile'):
			for i in range(0,len(PDBmulti[0].max90tile)):
				txtfile.write('%s\t' %(str(atom.max90tile[i])))

		elif densmet in ('min95tile'):
			for i in range(0,len(PDBmulti[0].min95tile)):
				txtfile.write('%s\t' %(str(atom.min95tile[i])))

		elif densmet in ('max95tile'):
			for i in range(0,len(PDBmulti[0].max95tile)):
				txtfile.write('%s\t' %(str(atom.max95tile[i])))

		elif densmet in ('mode','Mode'):
			for i in range(0,len(PDBmulti[0].modedensity)):
				txtfile.write('%s\t' %(str(atom.modedensity[i])))

		elif densmet in ('rsd','Rsd'):
			for i in range(0,len(PDBmulti[0].rsddensity)):
				txtfile.write('%s\t' %(str(atom.rsddensity[i])))

		elif densmet in ('dipstat','Dipstat'):
			for i in range(0,len(PDBmulti[0].dipstat)):
				txtfile.write('%s\t' %(str(atom.dipstat[i])))

		txtfile.write('\n')

	txtfile.close()


