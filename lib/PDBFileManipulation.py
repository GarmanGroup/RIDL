import sys
from classHolder import StructurePDB,singlePDB,multiPDB
import numpy as np
from operator import sub,truediv
from random import randint
from progbar import progress

def PDBtoList(pdbFileName = '',
              printText   = False):

    # this function inputs a pdb file name and returns an list 
    # of pdb objects, organised following the StructurePDB class

    PDBarray = []
    pdbin = open(pdbFileName, "r")
    lines = pdbin.readlines()
    if printText is True:
        print 'Reading PDB file and converting to list of objects'
    for line in lines:
        if ('ATOM' in str(line[0:6])) or ('HETATM' in str(line[0:6])):
            y               = singlePDB(StructurePDB)
            y.atomnum       = int(line[6:11].strip())
            y.atomtype      = str(line[12:16].strip())
            y.basetype      = str(line[17:20].strip())                       
            y.chaintype     = str(line[21])                     
            y.residuenum    = int(line[22:26].strip())             
            y.X_coord       = float(line[30:38].strip())                       
            y.Y_coord       = float(line[38:46].strip())                                             
            y.Z_coord       = float(line[46:54].strip())                                                 
            y.Occupancy     = str(line[54:60].strip())                                                    
            y.Bfactor       = float(line[60:66].strip())
            y.atomID        = str(line[76:78].strip())  
            y.atomOrHetatm  = str(line[0:6].strip())
            PDBarray.append(y)
        else: 
            pass
    pdbin.close()
    return PDBarray

def writePDBline_DamSite(atom     = '',
                         damValue = 0,
                         index    = 0,
                         chain    = 'A'):

    # script to convert atom information (in class format) to 'ATOM' line format for
    # pdb output files. Here the line is written in pdb format as a 'DAM' atom, 
    # with same xyz coordinates as the original atom input into the function.
    # The Bfactor column is written to be the damage metric value for that atom

    #for 'atom' of PDBarray (ie a specific atom of the structure):
    FIELD1 = "HETATM" #has length 6
    FIELD2 = " "*(5-len(str(index))) + str(index) #has length 5
    BREAK1 = " "*1 #has length 1
    FIELD3 = str('O') + " "*(4-len('O')) #has length 4
    FIELD4 = " " #has length 1
    FIELD5 = str('DAM').rjust(3) #has length 3
    BREAK2 = " " # has length 1
    FIELD6 = str(chain)
    FIELD7 = " "*(4-len(str(index))) + str(index) #has length 4
    FIELD8 = " " #has length 1
    BREAK3 = " "*3 #has length 3
    FIELD9 = "{0:8.3f}".format(atom.X_coord) #has length 8
    FIELD10 = "{0:8.3f}".format(atom.Y_coord) #has length 8
    FIELD11 = "{0:8.3f}".format(atom.Z_coord) #has length 8
    FIELD12 = "{0:6.2f}".format(float(1)) #has length 6, fix occupancy here as 1
    FIELD13 = "{0:6.3f}".format(float(damValue)) #has length 6
    BREAK4 = " "*6 # has length 6
    FIELD14 = " "*4 #has length 4
    FIELD15 = str(atom.atomID).rjust(2) #has length 2
    FIELD16 = " "*2 #has length 2
    
    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 + BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 + FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 + FIELD15 + FIELD16

    #as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        sys.exit('Error: PDB ATOM line written with inconsistent length (should be 80 characters long)\n'+\
                 '-> aborting line writing process, check atom object attributes for errors in formatting')        
    return line


def writePDBline_sudoWater(atom,Bfactor):
    # script to convert atom information (in class format) to 'ATOM' line format for
    # pdb output files. Here the line is written in pdb format as a sudo-water atom, 
    # but with same xyz coordinates as the original atom input into the function.
    # Here the Bfactor is written to be the value input by the user
    
    #for 'atom' of PDBarray (ie a specific atom of the structure):
    FIELD1 = "HETATM" #has length 6
    FIELD2 = " "*(5-len(str(atom.atomnum))) + str(atom.atomnum) #has length 5
    BREAK1 = " "*1 #has length 1
    FIELD3 = str('O') + " "*(4-len('O')) #has length 4
    FIELD4 = " " #has length 1
    FIELD5 = str('HOH').rjust(3) #has length 3
    BREAK2 = " " # has length 1
    FIELD6 = str('A')
    FIELD7 = " "*(4-len(str(atom.residuenum))) + str(atom.residuenum) #has length 4
    FIELD8 = " " #has length 1
    BREAK3 = " "*3 #has length 3
    FIELD9 = "{0:8.3f}".format(atom.X_coord) #has length 8
    FIELD10 = "{0:8.3f}".format(atom.Y_coord) #has length 8
    FIELD11 = "{0:8.3f}".format(atom.Z_coord) #has length 8
    FIELD12 = "{0:6.2f}".format(float(1)) #has length 6, fix occupancy here as 1
    FIELD13 = "{0:6.3f}".format(float(Bfactor)) #has length 6
    BREAK4 = " "*6 # has length 6
    FIELD14 = " "*4 #has length 4
    FIELD15 = str(atom.atomID).rjust(2) #has length 2
    FIELD16 = " "*2 #has length 2
    
    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 + BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 + FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 + FIELD15 + FIELD16

    #as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        print 'Error: PDB ATOM line written with inconsistent length (should be 80 characters long)'
        print '-> aborting line writing process, check atom object attributes for errors in formatting'
        sys.exit()
        
    return line





def convertPDBobject_toPDBline_fn(element,Occupancy):
    #script to convert atom information (in class format) to 'ATOM' line format for
    #pdb output files. Note that the following PDB ATOM line convention is used:
    
    #FIELD      #COLUMNS        DATA TYPE       CONTENTS                            
    #--------------------------------------------------------------------------------
    #FIELD1      1 -  6        Record name     "ATOM  "                                            
    #FIELD2      7 - 11        Integer         Atom serial number.                   
    #FIELD3     13 - 16        Atom            Atom name.                            
    #FIELD4     17             Character       Alternate location indicator.         
    #FIELD5     18 - 20        Residue name    Residue name.                         
    #FIELD6     22             Character       Chain identifier.                     
    #FIELD7     23 - 26        Integer         Residue sequence number.              
    #FIELD8     27             AChar           Code for insertion of residues.       
    #FIELD9     31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
    #FIELD10    39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
    #FIELD11    47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
    #FIELD12    55 - 60        Real(6.2)       Occupancy.                            
    #FIELD13    61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
    #FIELD14    73 - 76        LString(4)      Segment identifier, left-justified.   
    #FIELD15    77 - 78        LString(2)      Element symbol, right-justified.      
    #FIELD16    79 - 80        LString(2)      Charge on the atom.   
    
    #NOTE THAT OCCUPANCY IS NOT SPECIFIC TO THE 'ELEMENT' CHOSEN - HENCE IT CAN 
    #BE ASSIGNED TO ANOTHER METRIC IF DESIRED
    
    #for 'atom' of PDBarray (ie a specific atom of the structure):
    FIELD1 = "ATOM  " #has length 6
    FIELD2 = " "*(5-len(str(element.atomnum))) + str(element.atomnum) #has length 5
    BREAK1 = " "*1 #has length 2
    FIELD3 = str(element.atomtype) + " "*(4-len(str(element.atomtype))) #has length 4
    FIELD4 = " " #has length 1
    FIELD5 = str(element.basetype).rjust(3) #has length 3
    BREAK2 = " " # has length 1
    FIELD6 = str(element.chaintype)
    FIELD7 = " "*(4-len(str(element.residuenum))) + str(element.residuenum) #has length 4
    FIELD8 = " " #has length 1
    BREAK3 = " "*3 #has length 3
    FIELD9 = "{0:8.3f}".format(element.X_coord) #has length 8
    FIELD10 = "{0:8.3f}".format(element.Y_coord) #has length 8
    FIELD11 = "{0:8.3f}".format(element.Z_coord) #has length 8
    FIELD12 = "{0:6.2f}".format(float(Occupancy)) #has length 6
    FIELD13 = "{0:6.2f}".format(float(element.Bfactor)) #has length 6
    BREAK4 = " "*6 # has length 6
    FIELD14 = " "*4 #has length 4
    FIELD15 = str(element.atomID).rjust(2) #has length 2
    FIELD16 = " "*2 #has length 2
    
    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 + BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 + FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 + FIELD15 + FIELD16
    
    #as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        print 'Error: PDB ATOM line written with inconsistent length (should be 80 characters long)'
        print '-> aborting line writing process, check atom object attributes for errors in formatting'
        sys.exit()
        
    return line

    
def writePDBline(element,Bfactor):
    #script to convert atom information (in class format) to 'ATOM' line format for
    #pdb output files. Note that the following PDB ATOM line convention is used:
        
    #NOTE THAT Bfactor IS NOT SPECIFIC TO THE 'ELEMENT' CHOSEN - HENCE IT CAN 
    #BE ASSIGNED TO ANOTHER METRIC IF DESIRED
    
    #for 'atom' of PDBarray (ie a specific atom of the structure):
    FIELD1 = "ATOM  " #has length 6
    FIELD2 = " "*(5-len(str(element.atomnum))) + str(element.atomnum) #has length 5
    BREAK1 = " "*1 #has length 1
    FIELD3 = str(element.atomtype) + " "*(4-len(str(element.atomtype))) #has length 4
    FIELD4 = " " #has length 1
    FIELD5 = str(element.basetype).rjust(3) #has length 3
    BREAK2 = " " # has length 1
    FIELD6 = str(element.chaintype)
    FIELD7 = " "*(4-len(str(element.residuenum))) + str(element.residuenum) #has length 4
    FIELD8 = " " #has length 1
    BREAK3 = " "*3 #has length 3
    FIELD9 = "{0:8.3f}".format(element.X_coord) #has length 8
    FIELD10 = "{0:8.3f}".format(element.Y_coord) #has length 8
    FIELD11 = "{0:8.3f}".format(element.Z_coord) #has length 8
    FIELD12 = "{0:6.2f}".format(float(1.00)) #has length 6 (occupancy set to 1.00 here by default)
    FIELD13 = "{0:6.2f}".format(float(Bfactor)) #has length 6
    BREAK4 = " "*6 # has length 6
    FIELD14 = " "*4 #has length 4
    FIELD15 = str(element.atomID).rjust(2) #has length 2
    FIELD16 = " "*2 #has length 2
    
    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 + BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 + FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 + FIELD15 + FIELD16
    
    #as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        print 'Error: PDB ATOM line written with inconsistent length (should be 80 characters long)'
        print '-> aborting line writing process, check atom object attributes for errors in formatting'
        sys.exit()
        
    return line
