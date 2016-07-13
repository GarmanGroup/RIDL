class bioInfo():

	# contains basic information on amino acids etc

	def __init__(self):
		pass

	def getAminoAcids(self):

		# return ordered list of amino acids
		lst = ['ALA',
		  	   'ARG',
		  	   'ASN',
		  	   'ASP',
		  	   'CYS',
		 	   'GLN',
		 	   'GLU',
		 	   'GLY',
		 	   'HIS',
		 	   'ILE',
		  	   'LEU',
		  	   'LYS',
		  	   'MET',
		  	   'PHE',
		  	   'PRO',
		  	   'SER',
		  	   'THR',
		  	   'TRP',
		  	   'TYR',
		  	   'VAL']

		return lst

	def getNucleotides(self):

		# return ordered list of nucleic acid identifiers

		lst = ['A',
			   'C',
			   'G',
			   'T',
			   'U']

		return lst

	def getDNA(self):

		# return ordered list of DNA bases

		lst = ['DA',
			   'DC',
			   'DG',
			   'DT']

		return lst

	def getRNA(self):

		# return ordered list of RNA bases

		lst = ['A',
			   'C',
			   'G',
			   'U']
			   
		return lst