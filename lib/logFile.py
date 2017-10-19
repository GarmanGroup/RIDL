from time import gmtime, strftime


class logFile():

    def __init__(self,
                 fileName='untitled-log', fileDir='', printToScreen=True):

        self.logFile = fileName
        self.createLogFile()

        # fileDir is where majority of files come from
        self.fileDir = fileDir
        self.allocateDir()
        self.printToScreen = printToScreen

    def createLogFile(self):

        # create a log file

        log = open(self.logFile, 'w')
        log.write('Created at: {}\n'.format(self.getTime()))
        log.close()

    def allocateDir(self):

        # this is the specified directory where most
        # files have come from. This is written at the
        # top of the log file to avoid repetition.

        log = open(self.logFile, 'a')
        log.write('All files come from the following directory ' +
                  'unless otherwise given:\n"{}"\n'.format(self.fileDir))
        log.close()

    def writeToLog(self,
                   str='', strip=True, forcePrint=False, timeStamp=False):

        # write string to current log file

        # strip away the common directory name
        logstring = str
        if strip:
            if 'Working directory set to' not in str:
                logstring = str.replace(self.fileDir, '')

        if timeStamp is True:
            logstring = '{}\t{}'.format(self.getTime(), logstring)

        with open(self.logFile, "a") as logfile:
            logfile.write('{}\n'.format(logstring))

            if self.printToScreen or forcePrint:
                if logstring.startswith('Running'):
                    print '\n'+logstring
                else:
                    print logstring

    def getTime(self):

        # get the current date and time

        return strftime("%Y-%m-%d %H:%M:%S", gmtime())
