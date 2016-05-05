import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from checkSeabornPresent import checkSeabornPresent as checkSns
seabornFound = checkSns()
if seabornFound is True:
	import seaborn as sns

class halfDoseApprox():
	def __init__(self,
				 atom            = [],
				 atoms           = [],
				 plot            = True,
				 doses           = [],
				 globalScaling   = False,
				 doseFraction    = 0.5,
				 densityMetric   = 'loss',
				 normType        = 'Standard',
				 plotDir         = './',
				 shiftedHalfDose = True,
				 zeroOffset      = False):
		
		self.atom 				= atom # the current atom object 
		self.atoms 				= atoms # the full list of atom objects
		self.plot 				= plot # (Boolian) plot fitted data 
		self.globalScaling 		= globalScaling # (Boolian) decide whether to crudely weight for global density decay
		self.doseFraction 		= doseFraction # typically 0.5 for half dose but other values between 0 & 1 suitable
		self.densMet 			= densityMetric # density metric to use for calculation ('loss','gain',etc..)
		self.normType 			= normType # the corresponding metric normalisation type ('Standard', etc..)
		self.doses 				= doses # list of doses 
		self.avDmetric 			= np.mean([atm.densMetric[self.densMet][normType]['values'] for atm in atoms],0)
		self.plotDir 			= plotDir # where saved plots will go
		self.shiftedHalfDose 	= shiftedHalfDose # see method for half-dose calculation below
		self.zeroOffset 		= zeroOffset # takes values 'y' (fit end density) or float n for fixed end density n 

		if seabornFound is False:
			self.plot = False

		# define fitting functions here
		self.decayFuncOffset = lambda x, *p: p[0]*np.exp(-p[1]*x) + p[2]
		if self.zeroOffset != 'y':
			self.decayFunc   = lambda x, *p: p[0]*np.exp(-p[1]*x) + self.zeroOffset

		# run the following
		success = self.checkInputs()
		if success is False:
			print 'Invalid inputs supplied'
			return

		self.printAtomSummary()
		success = self.runFitter()
		if success is True:
			self.printCalcSummary()
			self.saveAtomHalfDose()

	def checkInputs(self):
		# check the class inputs are suitable
		for attr in [self.plot,self.globalScaling,self.shiftedHalfDose]:
			if attr not in (True,False):
				return False
		if not isinstance(self.doseFraction, float):
			return False
		if self.zeroOffset != 'y':
			if not (isinstance(self.zeroOffset,float) or isinstance(self.zeroOffset,int)):
					return False
		return True

	def runFitter(self):
		# fit the data to exponentially decaying model function
		xData = np.array(self.doses)

		yData = np.array(self.atom.densMetric[self.densMet][self.normType]['values'])
		if self.globalScaling is True: yData /= self.avDmetric

		# initial parameter guesses
		a = yData[0]
		k = (yData[1]-yData[0])/(a*(xData[0]-xData[1]))
		c = 0

		initGuess = [a,k]
		if self.zeroOffset == 'y': initGuess += [c] # extra offset param if required


		# run fitting dependent on whether decay function zero-offset permitted
		try:
			if self.zeroOffset == 'y':
				popt, pcov = curve_fit(self.decayFuncOffset, xData, yData,initGuess)
			else: 
				popt, pcov = curve_fit(self.decayFunc, xData, yData,initGuess)
				popt = list(popt) + [self.zeroOffset]
		except RuntimeError:
			print 'Max interations reached without convergence...'
			self.halfDose 			= 'NaN'
			self.halfDoseCertainty 	= 10**4 # big number --> no certainty 
			self.halfDoseStd 		= 'NaN'
			return False

		self.fitParams 		= np.array(popt)
		self.fitErrors 		= np.sqrt(np.diag(pcov))
		self.fitResiduals 	= np.sum([val**2 for val in (self.decayFuncOffset(xData,*self.fitParams)-yData)])
		self.halfDose 		= round(self.calculateHalfDose(self.fitParams,self.shiftedHalfDose),2)

		lwrParams,uprParams = self.fitParams,self.fitParams
		lwrParams[1]		= self.fitParams[1]+self.fitErrors[1]
		uprParams[1]		= self.fitParams[1]-self.fitErrors[1]
		lwrHalfDoseStd 		= self.calculateHalfDose(lwrParams,self.shiftedHalfDose)
		upprHalfDoseStd 	= self.calculateHalfDose(uprParams,self.shiftedHalfDose)
		self.halfDoseStd 	= np.array([lwrHalfDoseStd,upprHalfDoseStd])
		self.calculateInitDecayRate(self.fitParams)

		if self.zeroOffset == 'y':
			self.certaintyValue = np.product(self.fitParams/self.fitErrors)
		else:
			self.certaintyValue = np.product(self.fitParams[:-1]/self.fitErrors)

		if self.plot == True:
			self.plotDecayPlot(xData,yData)
		return True

	def saveAtomHalfDose(self):
		# for current atom object save the half dose statistics as attributes
		self.atom.densMetric[self.densMet][self.normType]['Half-dose'] = {'Half-dose'       : self.halfDose,
																	   	  'Residuals'       : self.fitResiduals,
																	      'Certainty'       : self.certaintyValue,
																	      'Initial density' : self.fitParams[0],
																	      'End density'     : self.fitParams[2],
																	      'init decay rate' : self.initDecayRate}

	def plotDecayPlot(self,xData,yData):
		sns.set(context='talk',style='dark')
		f = plt.figure()
		plt.plot(xData,yData,marker='.',linestyle='',markersize=14)	
		x = np.linspace(min(xData),max(xData),100)
		plt.plot(x,self.decayFuncOffset(x,*self.fitParams))

		# plot horizontal end density point 
		xvals = range(0,int(xData[-1]+1)+1)
		yvals = [self.fitParams[2]]*len(xvals)
		plt.plot(xvals,yvals,':')

		# plot start density point
		plt.plot([0],[self.fitParams[0]+self.fitParams[2]],'o')

		plt.xlabel("Dose (MGy)",fontsize=20)
		plt.ylabel('D{} (e/A^3)'.format(self.densMet),fontsize=20)
		plt.grid()
		identifier = self.atom.getAtomID() 
		sns.despine()
		f.suptitle('{} {} D{} Half-dose: {} MGy, initial decay: {}'.format(self.atom.getAtomID(),self.normType,
																		   self.densMet,self.halfDose,self.initDecayRate))
		if self.plotDir != '':
			f.savefig('{}/{}-D{}.png'.format(self.plotDir,identifier,self.densMet))
		else:
			f.savefig('{}-D{}.png'.format(identifier,self.densMet))

	def printCalcSummary(self):
		print '----------------------------'
		print 'Half dose calculation summary:'
		print 'Best a: {}, k: {}, c: {}'.format(*self.fitParams)
		if self.zeroOffset == 'y':
			print 'val/Std error: {}, {}, {}'.format(*(self.fitParams/self.fitErrors))
		else:
			print 'val/Std error: {}, {}'.format(*(self.fitParams[:-1]/self.fitErrors))
		print 'Residuals: {}'.format(self.fitResiduals)
		print "Half life predicted to be {}".format(self.halfDose)
		print "initial decay rate: {}".format(self.initDecayRate)
		print "Certainty value: {}".format(self.certaintyValue)
		print "Half dose Std interval: {}-{}".format(*self.halfDoseStd)

	def printAtomSummary(self):
		residue 	= self.atom.basetype + str(self.atom.residuenum)
		atomtype 	= self.atom.atomtype
		chain 		= self.atom.chaintype
		print '----------------------------'
		print 'chain:{} residue:{} atom:{}'.format(chain,residue,atomtype)

	def calculateHalfDose(self,fitConstants,shift):
		# if shift is False then half dose is dose for density to reach half of initial density
		if shift is False:
			return -np.log((self.doseFraction)*(1 - (fitConstants[2]/fitConstants[0])))/(fitConstants[1])
		# if shift is True then half dose is dose for density to reach av(initial density,end density limit)
		else:
			 return -np.log(self.doseFraction)/fitConstants[1]

	def calculateInitDecayRate(self,fitConstants):
		# calculate the initial decay rate for the current exponential decay model
		# It is the derivative at dose=0
		self.initDecayRate = round(fitConstants[0]*fitConstants[1],2)



