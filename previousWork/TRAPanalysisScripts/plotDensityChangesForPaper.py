
def plot(atoms):
	# plot Glu atoms of interest first
	atomType = ['CD','OE1','OE2']
	resiNum = [16,36,42,50,71,73]
	plotAtomsOfResType(atoms,atomType,resiNum,'GLU')

	# plot ASP atoms of interest next
	atomType = ['CG','OD1','OD2']
	resiNum = [8,17,29,39]
	plotAtomsOfResType(atoms,atomType,resiNum,'ASP')

	# plot Lys 37 atoms next
	atomType = ['CD','CE','NZ']
	resiNum = [37]
	plotAtomsOfResType(atoms,atomType,resiNum,'LYS')

	# plot Phe 32 atoms next
	atomType = ['CZ','CE2','CE1','CD1','CD2','CG']
	resiNum = [32]
	plotAtomsOfResType(atoms,atomType,resiNum,'PHE')

def plotAtomsOfResType(atoms,atomTypes,resiNums,residue):
	for t1 in atomTypes:
		for t2 in resiNums:
			atoms.atomType = t1
			atoms.baseType = residue
			atoms.residueNum = t2
			print '***\nPlotting for {}{} {} atoms:'.format(residue,t2,t1)
			atoms.densMetricErrorbarGraphs(True,'95confIntervalTestPlots',3,True)
