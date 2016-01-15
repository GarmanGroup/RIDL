import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve,curve_fit

class halfDoseApproximator():
	def __init__(self,atom,atoms,equivAtoms,plot,globalScaling,doseFraction,stdErrors,densityMetric,plotDir,shiftedHalfDose):
		self.atom = atom
		self.atoms = atoms
		self.equivAtoms = equivAtoms # true or false
		self.plot = plot
		self.globalScaling = globalScaling # true or false (decide whether to crudely weight for global density decay)
		self.doseFraction = doseFraction # typically 0.5 for half dose
		self.stdErrors = stdErrors # fit with raw data or mean & average per dose
		self.densMet = densityMetric # density metric to use for calculation ('loss','gain',etc..)
		self.doses = np.array([1.31,3.88,6.45,9.02,11.58,14.15,16.72,19.29,
				 			   21.86,24.98,33.01,40.86,48.89,56.75,64.78,72.63])
		self.avDloss = np.mean([atm.densMetric[self.densMet]['Standard']['values'] for atm in atoms],0)
		self.plotDir = plotDir # where saved plots will go
		self.shiftedHalfDose = shiftedHalfDose # see method for half-dose calculation below

		if self.equivAtoms == True:
			self.findEquivAtoms()
		else:
			self.equiv_atoms = [atom]

		self.printAtomSummary()
		# if self.plot == True:
		# 	self.plotData()
		self.runFitter2()
		# self.half_life_list = []
		# for i in range(2,7):
		# 	for j in range(i+1,7):		
		# 		# run solver now for current atom object
		# 		self.getParams(i,j)
		# 		half_life,success = self.runSolver()
		# 		if success == 1:
		# 			self.half_life_list.append(half_life)

		# plt.xlabel("k")
		# plt.ylabel("expression value")
		# plt.grid()
		# plt.show()

		# # provide overall summary:
		# if len(self.half_life_list) != 0:
		# 	self.halfLifeSummary()
		# else:
		# 	print 'No converged half life approximations'

	def plotExampleDecay(self,k):
		q = 1
		Y = q*(np.exp(k*self.doses[1:10])-np.exp(k*self.doses[0]))
		plt.plot(self.doses[1:10],Y)
		plt.xlabel("dose")
		plt.ylabel("Dloss")
		plt.grid()
		plt.show()

	def findEquivAtoms(self):
		nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
		RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']
		self.equiv_atoms = []
		for atm in self.atoms:
			if (atm.atomtype == self.atom.atomtype and
				atm.residuenum == self.atom.residuenum):
				if self.atom.chaintype in nonRNAbound and atm.chaintype in nonRNAbound:
					self.equiv_atoms.append(atm)
				elif self.atom.chaintype in RNAbound and atm.chaintype in RNAbound:
					self.equiv_atoms.append(atm)
				elif self.atom.chaintype in ['W','Y']:
					self.equiv_atoms.append(atm)

	def halfLifeSummary(self):
		mean_half_life = np.mean(self.half_life_list)
		std_half_life = np.mean(self.half_life_list)
		max_half_life = max(self.half_life_list)
		min_half_life = min(self.half_life_list)
		len_list = len(self.half_life_list)
		print '\nHalf life summary:'
		print '{} tests converged to a value'.format(len_list)
		print 'mean: {}'.format(mean_half_life)
		print 'std: {}'.format(mean_half_life)
		print 'range: {}--{}'.format(min_half_life,max_half_life)

	def getParams(self,i,j):
		self.d1 = self.doses[0]
		self.di = self.doses[i-1]
		self.dj = self.doses[j-1]

		if self.equivAtoms is False:
			self.diff_i = self.atom.mindensity[i-2]
			self.diff_j = self.atom.mindensity[j-2]
		else:
			self.diff_i = np.mean([atm.mindensity[i-2] for atm in self.equiv_atoms])
			self.diff_j = np.mean([atm.mindensity[j-2] for atm in self.equiv_atoms])
		if self.globalScaling == True:
			self.diff_i = self.diff_i - self.avDloss[i-2]
			self.diff_j = self.diff_j - self.avDloss[j-2]

	def runSolver(self):
		func = lambda k  : self.diff_i/self.diff_j - (np.exp(k*(self.di-self.d1))-1)/(np.exp(k*(self.dj-self.d1))-1)

		# plot it to check
		k = np.linspace(-10,10,1000)
		plt.plot(k,func(k))

		# use numerical solver to find solutions
		k_initial_guess = -1
		(k_solution,infodict,ier,mesg) = fsolve(func,k_initial_guess,full_output=True)
		print ier
		print mesg

		print "Solution is k = {}".format(k_solution)
		print "at which, the expression value is {}".format(func(k_solution))

		half_life = self.decayConstant(k_solution)
		print "Half life predicted to be {}".format(half_life)
		return half_life,ier

	def calculateHalfLife(self,fitConstants,shifted):
		# if shifted is False then half dose is dose for density to reach half of initial density
		if shifted is False:
			return -np.log((1.0/self.doseFraction)*(1 - (fitConstants[2]/fitConstants[0])))/(fitConstants[1])
		# if shifted is True then half dose is dose for density to reach av(initial density,end density limit)
		else:
			 return -np.log(1.0/self.doseFraction)/fitConstants[1]

	def runFitter(self):
		xData = self.doses[0:9]
		if self.globalScaling == True:
			yData = np.mean([(atm.densMetric[self.densMet]['Standard']['values'] - self.avDloss) for atm in self.equiv_atoms],0)
		else:
			yData = np.mean([(atm.densMetric[self.densMet]['Standard']['values']) for atm in self.equiv_atoms],0)

#		yData = [0]+list(yData[0:9])
		try:
			popt, pcov = curve_fit(DlossFunc2, xData, yData)
		except RuntimeError:
			print 'Max interations reached without convergence...'
			self.halfDose = 'NaN'
			self.halfDoseCertainty = 0
			return
		print 'Best a: {}, k: {}'.format(popt[0],popt[1])
		perr = np.sqrt(np.diag(pcov))
		print 'val/Std error: {} and {}'.format(popt[0]/perr[0],popt[1]/perr[1])
		if self.plot == True:
			plt.plot(xData,yData)
			plt.plot(xData,DlossFunc(xData,popt[0],popt[1]))
			plt.xlabel("dose")
			plt.ylabel("Dloss")
			plt.grid()
			plt.show()

		half_life = self.calculateHalfLife(popt,self.shiftedHalfDose)
		print "Half life predicted to be {}".format(half_life)
		print "certainty value: {}".format((popt[0]/perr[0])*(popt[1]/perr[1]))
		self.halfDose = half_life
		self.halfDoseCertainty = (popt[0]/perr[0])*(popt[1]/perr[1])

	def runFitter2(self):
		xData = self.doses[0:10]
		if self.globalScaling == True:
			yDataRaw = np.array([(atm.densMetric[self.densMet]['Standard']['values']/self.avDloss) for atm in self.equiv_atoms])
		else:
			yDataRaw = np.array([atm.densMetric[self.densMet]['Standard']['values'] for atm in self.equiv_atoms])
		yData = np.mean(yDataRaw,0)[0:10]
		if self.equivAtoms == True:
			weights = np.std(yDataRaw,0)[0:10] # sigma fitting weights
		else:
			weights = np.array([1]*10)

		# yData = [0]+list(yData[0:9])
		# weights = [weights[0]]+list(weights)

		# initial parameter guesses
		a = yData[0]
		k = (yData[1]-yData[0])/(a*(xData[0]-xData[1]))
		c = 0
		initGuess = [a,k,c]

		try:
			popt, pcov = curve_fit(DlossFunc3, xData, yData,sigma=weights,absolute_sigma=self.stdErrors)
			self.fitVals = popt
		except RuntimeError:
			print 'Max interations reached without convergence...'
			self.halfDose = 'NaN'
			self.halfDoseCertainty = 10**4 # big number --> no certainty 
			self.halfDoseStd = 'NaN'
			return

		print 'Best a: {}, k: {}, c: {}'.format(popt[0],popt[1],popt[2])
		perr = np.sqrt(np.diag(pcov))

		sumOfSquares = np.sum([val**2 for val in (DlossFunc3(xData,*popt)-yData)/weights])

		print 'val/Std error: {}, {}'.format(popt[0]/perr[0],popt[1]/perr[1])
		print 'Residuals: {}'.format(sumOfSquares)

		half_life = self.calculateHalfLife(popt,self.shiftedHalfDose)
		print "Half life predicted to be {}".format(half_life)
		certaintyValue = (popt[0]/perr[0])*(popt[1]/perr[1])
		print "certainty value: {}".format(sumOfSquares)
		self.halfDose = half_life
		self.halfDoseCertainty = sumOfSquares
		self.halfDoseStd = np.array([self.calculateHalfLife([popt[0],popt[1]+perr[1],popt[2]],self.shiftedHalfDose),
									 self.calculateHalfLife([popt[0],popt[1]-perr[1],popt[2]],self.shiftedHalfDose)])

		if self.plot == True:
			f = plt.figure()
			if self.equivAtoms == True:
				plt.errorbar(xData,yData,yerr=weights,fmt='--o')
			else:
				plt.plot(xData,yData,marker='.',linestyle='',markersize=10)
				
			plt.plot(xData,DlossFunc3(xData,*popt))

			# plot horizontal end density point 
			xvals = range(0,int(xData[-1])+1)
			yvals = [popt[2]]*len(xvals)
			plt.plot(xvals,yvals,':')

			# plot start density point
			plt.plot([0],[popt[0]+popt[2]],'o')

			plt.xlabel("dose")
			plt.ylabel("Density metric")
			plt.grid()
			identifier = self.atom.basetype + str(self.atom.residuenum)+'-'+ self.atom.atomtype+'-'+self.atom.chaintype
			f.suptitle('Half-dose: {} MGy'.format(self.halfDose))
			f.savefig('{}/{}.png'.format(self.plotDir,identifier))


	def plotData(self):
		xData = self.doses[1:10]
		if self.globalScaling == True:
			yDataRaw = [(atm.densMetric[self.densMet]['Standard']['values']/self.avDloss) for atm in self.equiv_atoms]
		else:
			yDataRaw = [(atm.densMetric[self.densMet]['Standard']['values']) for atm in self.equiv_atoms]
		# yData = yData[0:10]
		for i in range(len(yDataRaw)):
			plt.plot(xData,yDataRaw[i][0:9])
		plt.xlabel("dose")
		plt.ylabel("Dloss")
		plt.grid()
		plt.show()

	def printAtomSummary(self):
		print '----------------------------'
		residue = self.atom.basetype + str(self.atom.residuenum)
		atomtype = self.atom.atomtype
		if self.equivAtoms is False:
			chain = self.atom.chaintype
			print 'chain:{} residue:{} atom:{}'.format(chain,residue,atomtype)
		else:
			print 'Equiv atoms: residue:{} atom:{}'.format(residue,atomtype)

def DlossFunc(x, a, k):
	return a*(np.exp(-k*x) - np.exp(-k*1.31))

def DlossFunc2(x,a, k):
	return a*(np.exp(-k*x))

def DlossFunc3(x,a, k,c):
	return a*np.exp(-k*x) + c

def batchRun(atoms,N,globalScaling,doseFraction,stdErrors,densMetric,equivAtoms,plotDir,shiftedHalfDose):
	nonRNAbound = ['A','B','C','D','E','F','G','H','I','J','K']
	RNAbound =    ['L','M','N','O','P','Q','R','S','T','U','V']

	for resType in (['GLU','CD'],['GLU','OE1'],['GLU','OE2'],['ASP','CG'],['ASP','OD1'],['ASP','OD2']):
		csvFile = open('halfDoses_{}.csv'.format(resType[0]),'w')
		csvFile.write('Dose fraction = {}\n'.format(doseFraction))
		csvFile.write('Atom,Half-Dose,Residuals,Std-lower, Std-upper,initial-density,end-density\n')
		readAtoms = {}
		for atom in atoms:
			if (atom.basetype != resType[0] or atom.atomtype != resType[1]):
				continue
			if atom.chaintype in nonRNAbound:
				boundType = 'Non-Bound'
			elif atom.chaintype in RNAbound:
				boundType = 'Bound'
			elif atom.chaintype == 'W':
				boundType = 'RNA'
			elif atom.chaintype == 'Y':
				continue

			identifier = atom.basetype + str(atom.residuenum)+'-'+ atom.atomtype+'-'+boundType
			if equivAtoms == False:
				identifier += atom.chaintype
			if identifier in readAtoms.keys():
				continue
			else:
				half = halfDoseApproximator(atom,atoms,equivAtoms,False,globalScaling,doseFraction,stdErrors,densMetric,plotDir,shiftedHalfDose)
				readAtoms[identifier] = {'half dose':half.halfDose,'certainty':half.halfDoseCertainty,'std':half.halfDoseStd}
				if half.halfDoseCertainty < N:
					half = halfDoseApproximator(atom,atoms,equivAtoms,True,globalScaling,doseFraction,stdErrors,densMetric,plotDir,shiftedHalfDose)
					csvFile.write('{},{},{},{},{},{},{}\n'.format(identifier,half.halfDose,half.halfDoseCertainty,half.halfDoseStd[0],half.halfDoseStd[1],half.fitVals[0],half.fitVals[2]))
		csvFile.close()		
		for key in readAtoms.keys():
			if readAtoms[key]['certainty'] < N:
				halfDose = readAtoms[key]['half dose']
				certainty = readAtoms[key]['certainty']
				halfDoseStd = readAtoms[key]['std']
				print '{} --> {} (certainty:{})'.format(key,halfDose,certainty)





