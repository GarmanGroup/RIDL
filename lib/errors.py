import sys

class error():

	# a class used to print error information to command line

	def __init__(self,
				 text = '',
				 type = 'error',
				 log  = ''):

		# errorStr is the message to be printed to the command line
		# type specifies whether the message is an 'error' or 'warning'.
		# If a log file object is specified (see logFile.py), then 
		# the error message is printed into the log file

		self.log =log

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

		err = '\n****ERROR!****'+\
			  '\n{}'.format(message)+\
		      '\n\nPlease contact charles.bury@dtc.ox.ac.uk for queries regarding this failure'+\
		      '\n**************'

		self.printError(message = err)

	def warningMessage(self,
				       message = ''):

		# template for warning message printed to screen

		err = 'Warning! >>> {}'.format(message)

	def printError(self,
				   message = ''):

		# return error to command line or log file (if exists)

		if self.log != '':
			self.log.writeToLog(str   = message,
								strip = False)
		else:
			print message






