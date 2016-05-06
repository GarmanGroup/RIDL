import imp

def checkSeabornPresent(printText = False):

	# check whether seaborn is present and flag if not

	if printText is True:
		print 'Checking whether seaborn plotting library present...'
	try:
		imp.find_module('seaborn')
	except ImportError:
		if printText is True:
			print 'Plotting library seaborn not found..\n'+\
			 	  'Will not create any plots for current run.\n'+\
			      'Use "pip install seaborn" to install for future use.\n'
		return False
	if printText is True:	
		print '---> success!\n'
	return True