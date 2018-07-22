import sys

CONTACT_EMAIL = 'csbury@me.com'

CONTACT_EMAIL = 'csbury@me.com'

class error():

    # a class used to print error information to command line

    def __init__(self,
                 text='', type='error', log=''):

        # errorStr is the message to be printed to the command line
        # type specifies whether the message is an 'error' or 'warning'.
        # If a log file object is specified (see logFile.py), then
        # the error message is printed into the log file.
        # An 'error' is considered fatal and the script will stop, whereas
        # a 'warning' will not be considered fatal.

        self.log = log

        if type == 'error':
            self.errorMessage(message=text)

        elif type == 'warning':
            self.warningMessage(message=text)

        else:
            self.errorMessage(message='Unknown error type specified!')
            sys.exit()

    def errorMessage(self,
                     message=''):

        # template for error message printed to screen

        self.printError(
            message='\n****ERROR!****\n{}\n\nPlease contact '.format(message) +
                    '{} for queries regarding this '.format(CONTACT_EMAIL) +
                    'failure\n--->Terminating script...\n**************\n',
            killRun=True)

    def warningMessage(self,
                       message=''):

        # template for warning message printed to screen

        self.printError(message='\n***Warning!***: {}\n'.format(message),
                        killRun=False)

    def printError(self,
                   message='', killRun=True):

        # return error to command line or log file (if exists)

        if self.log != '':
            self.log.writeToLog(str=message, strip=False)
        else:
            print(message)

        if killRun:
            sys.exit()
