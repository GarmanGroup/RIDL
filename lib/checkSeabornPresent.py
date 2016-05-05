import imp

def checkSeabornPresent(printText = False):

	# check whether seaborn is present and flag if not

	if printText is True:
		print 'Checking whether seaborn plotting library present...'
	try:
		imp.find_module('seaborn')
	except ImportError:
		if printText is True:
			print 'Plotting library seaborn not found..'+\
			 	  'Will not create any plots for current run.'+\
			      'Use "pip install seaborn" to install for future use'
		return False
	if printText is True:	
		print '---> success!'
	return True