# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 23:39:24 2014
@author: charlie
"""
import sys
from classHolder import StructurePDB,singlePDB,multiPDB
import numpy as np
from operator import sub,truediv
from random import randint
from progbar import progress


def PDBtoList(pdbfilename,PDBarray):
    # this function inputs a pdb file name and returns an array of pdb objects, 
    # organised following the StructurePDB class
    
    pdbin = open(str(pdbfilename), "r")
    lines = pdbin.readlines()
    print 'Reading PDB file and converting to array of objects'
    for line in lines:
        if ('ATOM' in str(line[0:6])) or ('HETATM' in str(line[0:6])):
            y = singlePDB(StructurePDB)
            y.atomnum = int(line[6:11].strip())
            y.atomtype = str(line[12:16].strip())
            y.basetype = str(line[17:20].strip())                       
            y.chaintype = str(line[21])                     
            y.residuenum = int(line[22:26].strip())             
            y.X_coord = float(line[30:38].strip())                       
            y.Y_coord = float(line[38:46].strip())                                             
            y.Z_coord = float(line[46:54].strip())                                                 
            y.Occupancy = str(line[54:60].strip())                                                    
            y.Bfactor = str(line[60:66].strip())
            y.atomidentifier = str(line[76:78].strip())  
            PDBarray.append(y)
        else: 
            pass
    pdbin.close()
    return PDBarray


def getMultiDoseAtomList(data_list):
    # this function inputs a list of lists of PDB atom objects (see StructurePDB class)
    # and formats as an object of the class 'multiPDB'. It is a variant of the function 
    # above which can also cope with structures containing different numbers of atoms 
    # (say if solvent molecules/ligands are included in a subset of the structures). 
    # In this case, the smallest common substructure between all structures will be used
     
    #first check that each PDBarray contains the same number of atoms (consistency check)
    if len(data_list) > 1:
        print 'Multiple datasets detected...'
        for dataset in data_list:
            if len(dataset) != len(data_list[0]):
                print 'Not all PDB structures have same number of atoms!'\
                ' Will only include atoms common to ALL structures...'
    elif len(data_list) == 1:
        print 'Single dataset detected...'

    PDBdoses = []
    notincludedatmcounter = 0
    
    print '------------------------------------------------'
    print 'Locating common atoms to ALL datasets...'
    print 'SUMMARY:'

    i = 0
    num_atoms = len(data_list[0])

    for atom in data_list[0]:

        # unessential loading bar add-in
        i += 1
        progress(i, num_atoms, suffix='')

        atm_counter = 1
        
        Bfactor_comb = [atom.Bfactor]
        Occupancy_comb = [atom.Occupancy]
        meandensity_comb = [atom.meandensity]
        maxdensity_comb = [atom.maxdensity]
        mindensity_comb = [atom.mindensity]
        mediandensity_comb = [atom.mediandensity]
        stddensity_comb = [atom.stddensity]
        min90tile_comb = [atom.min90tile]
        max90tile_comb = [atom.max90tile]
        min95tile_comb = [atom.min95tile]
        max95tile_comb = [atom.max95tile]
        rsddensity_comb = [atom.rsddensity]
        range_comb = [atom.rangedensity]
        bdam_comb = [atom.bdam]
        numvoxels_comb = [atom.numvoxels]
        
        # list of index of same atom in each later dataset
        indexindataset = []

        # check whether atom in all datasets:
        for dataset in data_list[1:]:
            k = -1        
            for otheratom in dataset: 
                k += 1  
                if (atom.residuenum == otheratom.residuenum and 
                    atom.atomtype == otheratom.atomtype and 
                    atom.basetype == otheratom.basetype and 
                    atom.chaintype == otheratom.chaintype):
                    
                    atm_counter += 1       
                    Bfactor_comb.append(otheratom.Bfactor)
                    Occupancy_comb.append(otheratom.Occupancy)
                    meandensity_comb.append(otheratom.meandensity)
                    maxdensity_comb.append(otheratom.maxdensity)
                    mindensity_comb.append(otheratom.mindensity)
                    mediandensity_comb.append(otheratom.mediandensity)
                    stddensity_comb.append(otheratom.stddensity)
                    min90tile_comb.append(otheratom.min90tile)
                    max90tile_comb.append(otheratom.max90tile)
                    min95tile_comb.append(otheratom.min95tile)
                    max95tile_comb.append(otheratom.max95tile)
                    rsddensity_comb.append(otheratom.rsddensity)
                    range_comb.append(otheratom.rangedensity)
                    bdam_comb.append(otheratom.bdam)
                    numvoxels_comb.append(otheratom.numvoxels)
                    indexindataset.append(k)
                    break

        # remove this located atom from the later dataset now that it
        # has been located --> to make loop faster 
        for j in range(1,len(indexindataset)+1):
            if indexindataset[j-1] != -1:
                data_list[j].pop(indexindataset[j-1])
                        
        if atm_counter != len(data_list):
            print 'Atom not found in all datasets!'
            print 'res num: %s, atm type: %s, res type: %s, chain: %s '\
            %(atom.residuenum,atom.atomtype,atom.basetype,atom.chaintype)
            print '---> not including atom...'
            notincludedatmcounter += 1
            continue
      
        else:                 
            y = multiPDB()
            y.atomnum = atom.atomnum
            y.residuenum = atom.residuenum
            y.atomtype = atom.atomtype
            y.basetype = atom.basetype 
            y.chaintype = atom.chaintype
            y.X_coord = atom.X_coord
            y.Y_coord = atom.Y_coord
            y.Z_coord = atom.Z_coord
            y.atomidentifier = atom.atomidentifier
            y.numsurroundatoms = atom.numsurroundatoms
            y.numsurroundprotons = atom.numsurroundprotons

            y.Bfactor = Bfactor_comb 
            y.Occupancy = Occupancy_comb
            y.meandensity = meandensity_comb
            y.maxdensity = maxdensity_comb
            y.mindensity = mindensity_comb
            y.mediandensity = mediandensity_comb
            y.stddensity = stddensity_comb
            y.min90tile = min90tile_comb
            y.max90tile = max90tile_comb
            y.min95tile = min95tile_comb
            y.max95tile = max95tile_comb
            y.rsddensity = rsddensity_comb
            y.rangedensity = range_comb
            y.bdam = bdam_comb
            y.numvoxels = numvoxels_comb
               
            PDBdoses.append(y)
    print '\n------------------------------------------------'    
    print 'Number of atoms removed since not in all datasets: %s' %str(notincludedatmcounter)
    print '---> Finished!'
    return PDBdoses
###############################################################################



###############################################################################
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
    FIELD15 = str(element.atomidentifier).rjust(2) #has length 2
    FIELD16 = " "*2 #has length 2
    
    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 + BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 + FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 + FIELD15 + FIELD16
    
    #as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        print 'Error: PDB ATOM line written with inconsistent length (should be 80 characters long)'
        print '-> aborting line writing process, check atom object attributes for errors in formatting'
        sys.exit()
        
    return line
###############################################################################
    

###############################################################################
def convertPDBobject_toPDBline_fn_Bfactorsetter(element,Bfactor):
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
    FIELD12 = "{0:6.2f}".format(float(element.Occupancy)) #has length 6
    FIELD13 = "{0:6.2f}".format(float(Bfactor)) #has length 6
    BREAK4 = " "*6 # has length 6
    FIELD14 = " "*4 #has length 4
    FIELD15 = str(element.atomidentifier).rjust(2) #has length 2
    FIELD16 = " "*2 #has length 2
    
    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 + BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 + FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 + FIELD15 + FIELD16
    
    #as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        print 'Error: PDB ATOM line written with inconsistent length (should be 80 characters long)'
        print '-> aborting line writing process, check atom object attributes for errors in formatting'
        sys.exit()
        
    return line
###############################################################################



###############################################################################
def convertPDBobject_toPDBline_sudoWater(atom,Bfactor):
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
    FIELD15 = str(atom.atomidentifier).rjust(2) #has length 2
    FIELD16 = " "*2 #has length 2
    
    line = FIELD1 + FIELD2 + BREAK1 + FIELD3 + FIELD4 + FIELD5 + BREAK2 + FIELD6 + FIELD7 + FIELD8 + BREAK3 + FIELD9 + FIELD10 + FIELD11 + FIELD12 + FIELD13 + BREAK4 + FIELD14 + FIELD15 + FIELD16

    #as a check the length of the line should be 80 characters long:
    if len(line) != 80:
        print 'Error: PDB ATOM line written with inconsistent length (should be 80 characters long)'
        print '-> aborting line writing process, check atom object attributes for errors in formatting'
        sys.exit()
        
    return line
###############################################################################
