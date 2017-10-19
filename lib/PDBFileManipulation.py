import sys
from classHolder import StructurePDB, singlePDB


def PDBtoList(pdbFileName='', printText=False):

    # this function inputs a pdb file name and returns an list
    # of pdb objects, organised following the StructurePDB class

    PDBarray = []
    pdbin = open(pdbFileName, "r")
    lines = pdbin.readlines()
    if printText:
        print 'Reading PDB file and converting to list of objects'
    for line in lines:
        if ('ATOM' in str(line[0:6])) or ('HETATM' in str(line[0:6])):
            y = singlePDB(StructurePDB)
            y.atomnum = int(line[6:11].strip())
            y.atomtype = str(line[12:16].strip())
            y.basetype = str(line[17:20].strip())
            y.chaintype = str(line[21])
            y.residuenum = str(line[22:27].strip())
            y.X_coord = float(line[30:38].strip())
            y.Y_coord = float(line[38:46].strip())
            y.Z_coord = float(line[46:54].strip())
            y.Occupancy = str(line[54:60].strip())
            y.Bfactor = float(line[60:66].strip())
            y.atomID = str(line[76:78].strip())
            y.atomOrHetatm = str(line[0:6].strip())
            PDBarray.append(y)
        else:
            pass
    pdbin.close()
    return PDBarray


def writePDBline_DamSite(atom='', damValue=0, index=0, chain='A'):

    # script to convert atom information (in class format) to 'ATOM'
    # line format for pdb output files. Here the line is written in
    # pdb format as a 'DAM' atom, with same xyz coordinates as the
    # original atom input into the function. The Bfactor column is
    # written to be the damage metric value for that atom.

    # for 'atom' of PDBarray (ie a specific atom of the structure)
    FIELD1 = "HETATM"
    FIELD2 = " "*(5-len(str(index))) + str(index)
    BREAK1 = " "*1
    FIELD3 = str('O') + " "*(4-len('O'))
    FIELD4 = " "
    FIELD5 = str('DAM').rjust(3)
    BREAK2 = " "
    FIELD6 = str(chain)
    FIELD7 = " "*(4-len(str(index))) + str(index)
    FIELD8 = " "
    BREAK3 = " "*3
    FIELD9 = "{0:8.3f}".format(atom.X_coord)
    FIELD10 = "{0:8.3f}".format(atom.Y_coord)
    FIELD11 = "{0:8.3f}".format(atom.Z_coord)
    FIELD12 = "{0:6.2f}".format(float(1))
    FIELD13 = "{0:6.3f}".format(float(damValue))
    BREAK4 = " "*6
    FIELD14 = " "*4
    FIELD15 = str(atom.atomID).rjust(2)
    FIELD16 = " "*2

    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 +\
        BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 +\
        FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 +\
        FIELD15 + FIELD16

    # as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        sys.exit('Error: PDB ATOM line has inconsistent length (should be 80' +
                 ' characters long)\n--> aborting line writing process, ' +
                 'check atom object attributes for errors in formatting')
    return line


def writePDBline(element, Bfactor):
    # function to convert atom information (in class format) to 'ATOM'
    # line format for pdb output files.
    # The Bfactor of each atom can be specified by another metric if desired

    # for 'atom' of PDBarray (ie a specific atom of the structure):
    FIELD1 = "ATOM  "
    FIELD2 = " "*(5-len(str(element.atomnum))) + str(element.atomnum)
    BREAK1 = " "*1
    FIELD3 = str(element.atomtype) + " "*(4-len(str(element.atomtype)))
    FIELD4 = " "
    FIELD5 = str(element.basetype).rjust(3)
    BREAK2 = " "
    FIELD6 = str(element.chaintype)
    FIELD7 = " "*(4-len(str(element.residuenum))) + str(element.residuenum)
    FIELD8 = " "
    BREAK3 = " "*3
    FIELD9 = "{0:8.3f}".format(element.X_coord)
    FIELD10 = "{0:8.3f}".format(element.Y_coord)
    FIELD11 = "{0:8.3f}".format(element.Z_coord)
    FIELD12 = "{0:6.2f}".format(float(1.00))
    FIELD13 = "{0:6.2f}".format(float(Bfactor))
    BREAK4 = " "*6
    FIELD14 = " "*4
    FIELD15 = str(element.atomID).rjust(2)
    FIELD16 = " "*2

    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 +\
        BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 +\
        FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 +\
        FIELD15 + FIELD16

    # as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        sys.exit('Error: PDB ATOM line has inconsistent length (should be 80' +
                 ' characters long)\n--> aborting line writing process, ' +
                 'check atom object attributes for errors in formatting')
    return line
