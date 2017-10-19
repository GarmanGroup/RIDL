import struct


class mapTools():

    # a small class to parse map files and read header information

    def __init__(self,
                 mapName='untitled.map', logFile=''):

        self.mapName = mapName
        self.log = logFile
        self.readHeader()

    def readAllHeader(self):

        # parse all header information from a
        # map file and dump to command line

        binaryMapFile = open(self.mapName, 'rb')
        for i in range(0, 30):
            print self.get4Bytes(binaryMapFile)

    def readHeader(self):

        # parse specific information from map header

        binaryMapFile = open(self.mapName, 'rb')
        self.numCols = self.get4Bytes(binaryMapFile)
        self.numRows = self.get4Bytes(binaryMapFile)
        self.numSecs = self.get4Bytes(binaryMapFile)
        binaryMapFile.seek(7*4, 0)
        self.gridsamp1 = self.get4Bytes(binaryMapFile)
        self.gridsamp2 = self.get4Bytes(binaryMapFile)
        self.gridsamp3 = self.get4Bytes(binaryMapFile)
        binaryMapFile.seek(16*4, 0)
        self.fastaxis = self.get4Bytes(binaryMapFile)
        self.medaxis = self.get4Bytes(binaryMapFile)
        self.slowaxis = self.get4Bytes(binaryMapFile)
        binaryMapFile.close()

    def getMapSize(self):

        # get map size in number of voxels

        numVoxels = self.numCols*self.numRows*self.numSecs
        return numVoxels

    def get4Bytes(self, binaryMapFile):

        # get next 4 bytes within a map file

        info = struct.unpack('=l', binaryMapFile.read(4))[0]
        return info

    def printMapInfo(self):

        # print summary map info to log file or
        # command line if log file not found

        numVoxels = self.getMapSize()
        info = '\tNumber of columns, rows, sections ......  {} {} {}\n'.format(
            self.numCols, self.numRows, self.numSecs) +\
            '\tGrid sampling on x, y, z ...............  {} {} {}\n'.format(
            self.gridsamp1, self.gridsamp2, self.gridsamp3) +\
            '\tFast, medium, slow axes ................  {} {} {}\n'.format(
            self.fastaxis, self.medaxis, self.slowaxis) +\
            '\tNumber of map voxels ...................  {}'.format(numVoxels)

        if self.log == '':
            print info
        else:
            self.log.writeToLog(str=info)
