# idea is to calculate finite grid sampling of unit cell and partition dose amongst grid elements
# dependent on X-ray absorption cross-section /
# !!! make sure pdb file pre-calculated over full unit cell !!!
import numpy as np

class dosePartition(object):
	def __init__(self, pdbFileName = "",gridSamplingRate=0):
		self.pdbFileName 		= pdbFileName
		self.gridSamplingRate 	= gridSamplingRate
		self.fullBreakdown 		= False

	def calculateGrid(self):
		# read in pdb file
		pdbLines = open(self.pdbFileName,'r').readlines()
		self.symMatrix = []
		for line in pdbLines:
			if 'CRYST1' in line.split()[0]:
				self.crystDims = {}
				for i, abc_ind in enumerate(['a dim','b dim','c dim']):
					self.crystDims[abc_ind] = float(line.split()[i+1])
			for ind in range(1,4):
				if 'SCALE{}'.format(ind) in line.split()[0]:
					self.symMatrix.append(map(float,line.split()[1:4]))
			if 'ATOM' in line.split()[0]:
				break # no more action needed

		#provide some feedback
		print 'xyz dimensions as follows:\n{}'.format(self.crystDims)
		print 'pdb xyz --> fract matrix as follows:\n{}'.format(self.symMatrix)

		# create list of grid element objects
		self.gridElements = []
		for i in range(self.gridSamplingRate):
			for j in range(self.gridSamplingRate):
				for k in range(self.gridSamplingRate):
					grdElmnt = gridElement([i,j,k],self.gridSamplingRate)
					self.gridElements.append(grdElmnt)

	def findAtomGridElement(self):
		error = False # error check to determine that all atoms assigned to grid element in unit cell
		errorcount = 0
		pdbLines = open(self.pdbFileName,'r').readlines()
		print 'calculating number atoms per grid point - could take some time..'
		for line in pdbLines:
			if 'ATOM' in line.split()[0]:
				xyz = []
				for line_part in [line[30:38],line[38:46],line[46:54]]:
					xyz.append(float(line_part.strip()))
				# calculate fractional coords for atom  
				xyz_frac = np.dot(self.symMatrix,xyz)  
				xyz_frac = [x%1 for x in xyz_frac] # ensure fractional coordinate within 0-1 range

				grdElFound = False
				for grdElmnt in self.gridElements:
					if grdElmnt.a_start <= xyz_frac[0] and xyz_frac[0] <= grdElmnt.a_end:
						if grdElmnt.b_start <= xyz_frac[1] and xyz_frac[1] <= grdElmnt.b_end:
							if grdElmnt.c_start <= xyz_frac[2] and xyz_frac[2] <= grdElmnt.c_end:
								grdElmnt.atoms.append(str(line[76:78].strip()))
								grdElFound = True
				if grdElFound == False:
					error = True
					errorcount += 1

		# check whether all atoms assigned to grid element
		if error == True:
			print 'Some atoms ({} of them) not assigned to grid element --> this is bad!'.format(errorcount)

	def atomsPerGridElement(self):
		print '--------------------------------------------'
		print 'Summary of number of atoms per grid element:'
		
		if self.fullBreakdown == True:
			atomNumList = []
			for grdElmnt in self.gridElements:
				print 'x:{:.2f}-{:.2f} y:{:.2f}-{:.2f} z:{:.2f}-{:.2f}'.format(grdElmnt.a_start,grdElmnt.a_end,
																			   grdElmnt.b_start,grdElmnt.b_end,
		 																	   grdElmnt.c_start,grdElmnt.c_end)
				print '# atoms: {}'.format(len(grdElmnt.atoms))
				atomNumList.append(len(grdElmnt.atoms))

		else: # don't waste time with above loop if not required
			atomNumList = [len(grdElmnt.atoms) for grdElmnt in self.gridElements]															     	  

		print 'max atoms: {}'.format(max(atomNumList))
		print 'min atoms: {}'.format(min(atomNumList))
		print 'mean atoms: {}'.format(np.mean(atomNumList))
		print 'std atoms: {}'.format(np.std(atomNumList))

	def runToEnd(self):
		# pipe the above methods together
		self.calculateGrid()
		self.findAtomGridElement()
		self.atomsPerGridElement()

class gridElement(object):
	# a simple class for grid elements
	def __init__(self,index,gridSamplingRate):
		self.gridSamplingRate 	= gridSamplingRate
		self.index 				= index
		self.a_start 			= self.gridDivide(0,True)
		self.b_start 			= self.gridDivide(1,True)
		self.c_start 			= self.gridDivide(2,True)
		self.a_end 				= self.gridDivide(0,False)
		self.b_end 				= self.gridDivide(1,False)
		self.c_end 				= self.gridDivide(2,False)

		# set an initial empty list of atoms belonging to grid element
		self.atoms = []

	def gridDivide(self,ind,start):
		if start == True:
			return (float(self.index[ind])/self.gridSamplingRate)
		else:
			return (float(self.index[ind]+1)/self.gridSamplingRate)


