import sys

class error():

	# a class used to print error information to command line

	def __init__(self,
				 text  = '',
				 type  = 'error'):

		# errorStr is the message to be printed to the command line
		# type specifies whether the message is an 'error' or 'warning'

		if type == 'error':
			self.errorMessage(message = text)

		elif type == 'warning':
			self.warningMessage(message = text)

		else:
			self.errorMessage(message = 'Uknown error type specified!')
			sys.exit()

	def errorMessage(self,
				     message = ''):

		# template for error message printed to screen

		print '\n****ERROR!****'
		print '{}'.format(message)
		print '\nPlease contact charles.bury@dtc.ox.ac.uk for queries regarding this failure'
		print '**************'


	def warningMessage(self,
				       message = ''):

		# template for warning message printed to screen

		print 'Warning! >>> {}'.format(message)





