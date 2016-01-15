import struct

class mapTools():

	def __init__(self,mapName):
		self.mapName = mapName
		self.readHeader()

	def readAllHeader(self):
		binaryMapFile = open(self.mapName,'rb')
		for i in range(0,30):
			info = self.get4Bytes(binaryMapFile)
			print info

	def readHeader(self):
		# open electron density .map file here 
		binaryMapFile = open(self.mapName,'rb')
		self.numCols = self.get4Bytes(binaryMapFile)
		self.numRows = self.get4Bytes(binaryMapFile)
		self.numSecs = self.get4Bytes(binaryMapFile)
		binaryMapFile.seek(7*4,0)
		self.gridsamp1 = self.get4Bytes(binaryMapFile)
		self.gridsamp2 = self.get4Bytes(binaryMapFile)
		self.gridsamp3 = self.get4Bytes(binaryMapFile)
		binaryMapFile.seek(16*4,0)
		self.fastaxis = self.get4Bytes(binaryMapFile)
		self.medaxis = self.get4Bytes(binaryMapFile)
		self.slowaxis = self.get4Bytes(binaryMapFile)
		binaryMapFile.close()

	def getMapSize(self):
		numVoxels = 4*(self.numCols*self.numRows*self.numSecs)
		return numVoxels

	def get4Bytes(self,binaryMapFile):
		info = struct.unpack('=l',binaryMapFile.read(4))[0]
		return info

	def printMapInfo(self):
		numVoxels = self.getMapSize()
		string 	=  	'Number of columns, rows, sections ......  {} {} {}\n'.format(self.numCols,self.numRows,self.numSecs)+\
					'Grid sampling on x, y, z ...............  {} {} {}\n'.format(self.gridsamp1,self.gridsamp2,self.gridsamp3)+\
					'Fast, medium, slow axes ................  {} {} {}\n'.format(self.fastaxis,self.medaxis,self.slowaxis)+\
					'Number of map voxels ...................  {}'.format(numVoxels)
		print string
