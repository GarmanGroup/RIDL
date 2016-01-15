# small helpful script to retrieve the processedAtom class list of objects from
# scratch after opening python
# NO LONGER USED SINCE PROCESSEDATOMLIST CLASS REPLACED

from eTrack_RUN import eTrack
from processedAtomList import processedAtomList

def retrieveAtomList():
	et = eTrack()
	et.runPipeline('n','n','y','n','e_Track_inputfile_TRAP.txt')
	doses = [1.31,3.88,6.45,9.02,11.58,14.15,16.72,19.29,21.86,24.98,33.01,40.86,48.89,56.75,64.78,72.63]
	atomList = processedAtomList(et.PDBmulti,10,doses[1:10])
	atomList.processAtomList()
	return atomList
