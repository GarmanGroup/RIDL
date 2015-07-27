# small helpful script to retrieve the processedAtom class list of objects from
# scratch after opening python

from eTrack_RUN import eTrack
from processedAtomList import processedAtomList

def retrieveAtomList():
	et = eTrack()
	et.runPipeline('n','n','y','n','e_Track_inputfile_TRAP.txt')
	atomList = processedAtomList(et.PDBmulti,10)
	atomList.processAtomList()
	return atomList.processedAtomList
