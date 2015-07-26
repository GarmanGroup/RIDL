# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 17:08:27 2014

@author: charlie
"""
import numpy as np
from classHolder import StructurePDB
import sys
import os
from colRowSec2Cartesian import findCartesianTranslationVectors
from PDBFileManipulation import PDBtoCLASSARRAY_v2 as pdb2list
from progbar import progress
from scipy import spatial as spa


def pdbCUR_symgen(inputpdbfile,outputpdbfile,space_group):
    # function to run CCP4 program pdbCUR to generate all symmetry related atoms of 
    # original structure 
    print '•••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••'
    print 'Determining symmetrically related atoms to original structure, using pdbCUR'

    input1 = "/Applications/ccp4-6.4.0/bin/pdbcur "+\
             "XYZIN %s " %(str(inputpdbfile))+\
             "XYZOUT %s " %(str(outputpdbfile))+\
             "SYMINFO /Applications/ccp4-6.4.0/lib/data/syminfo.lib "

    # specify to remove hydrogen atoms and pick on the most probably conformation
    # (for two conformations with 0.5 occupancy, the first - A is chosen and occupancy
    # set to 1.00). Also remove all anisou info from file - since it is not needed for 
    # current analysis
    input2 = "delhydrogen\n"+\
             "mostprob\n"+\
             "noanisou\n"+\
             "symm %s\n" %(str(space_group))+\
             "GENUNIT\n"+\
             "mkchainIDs\n"+\
             "END"

    # run input2 in pdbcur
    textinput = open('inputfile.txt','w')
    textinput.write(input2)
    textinput.close()
    os.system(input1 +' < inputfile.txt')
    os.remove('inputfile.txt')

    print '\n---> pdbcur run complete...'
    print '--------------------------'
    print 'Summary:'
    print 'Input pdb file: %s' %(str(inputpdbfile))
    print 'Output pdb file: %s' %(str(outputpdbfile))

    # determine initial number of atoms in input pdb file
    pdbin = open(str(inputpdbfile),'r')
    pdbinlines = pdbin.readlines()
    counter = 0
    for line in pdbinlines:
        if 'ATOM' in line[0:5]:
            counter += 1
    pdbin.close()
    print 'number of atoms in input pdb file: %s' %(str(counter))

    # determine number of atoms in output pdb file
    pdbout = open(str(outputpdbfile),'r')
    pdboutlines = pdbout.readlines()
    counter = 0
    for line in pdboutlines:
        if 'ATOM' in line[0:5]:
            counter += 1
    pdbout.close()
    print 'number of atoms in output pdb file: %s' %(str(counter))


def translate26cells(inputpdbfile,outpdbfile):
    # function to create an extended pdb file containing all atoms translated from
    # the original structure into the adjacent 26 unit cells surrounding the unit
    # cell the structure is situated in
    print '\n••••••••••••••••••••••••••••••••••••••••••••••••'
    print 'Determining related atoms to original structure,\n'+\
          'translated into adjacent 26 unit cells'

    pdbin = open(str(inputpdbfile),'r')
    pdbout = open(str(outpdbfile),'w')

    for line in pdbin.readlines():
        if 'CRYST1' in line[0:6]:
            a = float(line.split()[1])
            b = float(line.split()[2])
            c = float(line.split()[3])
            alpha = float(line.split()[4])
            beta = float(line.split()[5])
            gamma = float(line.split()[6])
            pdbout.write(line)

            avec,bvec,cvec = findCartesianTranslationVectors(a,b,c,alpha,beta,gamma)

        if ('ATOM' in line[0:5] or 'HETATM' in line[0:6]):

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            for i in ([-1,0,1]):
                for j in ([-1,0,1]):
                    for k in ([-1,0,1]):
                        xyz_new = np.array([x,y,z]) + i*np.array(avec) + j*np.array(bvec) + k*np.array(cvec)
                        x_new = xyz_new[0]
                        y_new = xyz_new[1]
                        z_new = xyz_new[2]

                        pdbout.write(line[0:30])
                        pdbout.write("{0:8.3f}".format(x_new))
                        pdbout.write("{0:8.3f}".format(y_new))
                        pdbout.write("{0:8.3f}".format(z_new))
                        pdbout.write(line[54:80])
                        pdbout.write('\n')
    pdbin.close()
    pdbout.close()
    print '---> success!'



def restrict14A(PDBarrayin,expanded_inputfile,outputpdbfile):
    # function to restrict the pdb file containing neighbouring 
    # atoms of original pdb structure to 14 Angstroms around
    # original structure
    print '\n••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••'
    print 'Restricting atoms to within 14 Angstroms of original structure\n'

    # find max and min x values in original PDB file
    PDBarrayin.sort(key=lambda x: (x.X_coord))
    minX = PDBarrayin[0].X_coord
    maxX = PDBarrayin[len(PDBarrayin)-1].X_coord

    # find max and min y values in original PDB file
    PDBarrayin.sort(key=lambda x: (x.Y_coord))
    minY = PDBarrayin[0].Y_coord
    maxY = PDBarrayin[len(PDBarrayin)-1].Y_coord

    # find max and min z values in original PDB file
    PDBarrayin.sort(key=lambda x: (x.Z_coord))
    minZ = PDBarrayin[0].Z_coord
    maxZ = PDBarrayin[len(PDBarrayin)-1].Z_coord

    pdbin = open(str(expanded_inputfile),'r')
    pdbout = open(str(outputpdbfile),'w')
    for line in pdbin.readlines():
        if ('ATOM' in line[0:5] or 'HETATM' in line[0:6]):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            if x - maxX < 14 and minX - x < 14:
                if y - maxY < 14 and minY - y < 14:
                    if z - maxZ < 14 and minZ - z < 14:
                        pdbout.write(line)
        else:
            pdbout.write(line)

    pdbin.close()
    pdbout.close()
    print '---> success!'



def numsurroundatoms_calculate(initialPDBfile,PDBarray,threshold):
    # function determines for each atom in structure the number of neighbouring atoms within 
    # a threshold (defined above) for all atoms. For each atom, number of contacts added
    # as class attribute for atom
    print '••••••••••••••••••••••••••••••••••••••••••••••••••••'
    print 'Calculating contact number for atoms in structure...'

    # determine the correct extended pdb file, with atoms present up to 1 unit cell 
    # away from the original structure
    inputpdbfile1 = initialPDBfile

    # determine the space group for the input pdb file:
    pdbin = open(str(inputpdbfile1),'r')
    for line in pdbin.readlines():
        if 'CRYST1' in line[0:6]:
            space_group = line[55:66]
    pdbin.close

    # run the above functions to (a) determine the symmetrically related
    # atoms to the original structure, (b) translate to determine the 
    # location of all atoms within the adjacent 26 unit cells to the 
    # original structure, and (c) to restrict to atoms only within 14 
    # Angstroms of the original structure.
    outputpdbfile1 = initialPDBfile[:-4]+'_pdbCURsymgenOUT.pdb'
    pdbCUR_symgen(inputpdbfile1,outputpdbfile1,space_group)

    outputpdbfile2 = initialPDBfile[:-4]+'_translate26cells.pdb'
    translate26cells(outputpdbfile1,outputpdbfile2)
    
    extended_pdbfile = initialPDBfile[:-4]+'_restrict14A.pdb'
    restrict14A(PDBarray,outputpdbfile2,extended_pdbfile)

    # read through extended 14A pdb file and collect all xyz coords of atoms
    # into a list allcoords. allcoords_atmtypes contains the atom identifier 
    # name for easy reference to the atom type associated with each atom found
    pdbin = open(extended_pdbfile,'r')
    allcoords = []
    allcoords_atmtypes = []
    for line in pdbin.readlines():
        if ('ATOM' in line[0:5] or 'HETATM' in line[0:6]):
            allcoords.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
            allcoords_atmtypes.append(str(line[76:78]).strip())
    pdbin.close()


    # convert here the atom names in allcoords_atmtypes into proton numbers
    print 'Locating proton number for each atom close to structure...'
    allcoords_protons = []
    for element in allcoords_atmtypes:
        atomdetailfile = open('VDVradiusfile.txt','r')
        for line in atomdetailfile.readlines():
            if element == line.split()[1]:
                allcoords_protons.append(int(line.split()[0]))
                break
        atomdetailfile.close()
    allcoords_protons = np.array(allcoords_protons)

    # check that all atoms have been assigned proton numbers in last step
    if len(allcoords_atmtypes) != len(allcoords_protons):
        print 'Not all atoms within 14A of structure successfully assigned proton numbers'
        print '---> terminating script...'
        sys.exit()
    else:
        print '---> success!'
    del allcoords_atmtypes

    counter = 0
    num_atoms = len(PDBarray)
    for atom in PDBarray:
        counter += 1

        # unessential progress bar added here
        progress(counter, num_atoms, suffix='')

        # want to determine the number of contacts (defined as number of atoms)
        # and also number of protons (to distinguish between different atom types)
        num_contacts = 0
        num_protons = 0
        atmxyz = np.array([[atom.X_coord,atom.Y_coord,atom.Z_coord]])

        # efficient distance calculation
        dist = spa.distance.cdist(np.array(allcoords),atmxyz)

        # sorted_dist = np.sort(dist,axis=None)
        sort_order = dist.argsort(axis=None)
        sorted_dist = dist[sort_order]
        sorted_allcoords_protons = allcoords_protons[sort_order]
        del dist,sort_order,atmxyz

        num_contacts = next(x[0] for x in enumerate(sorted_dist) if x[1] > threshold)
        # for element in sorted_dist:
        #     if element < threshold:
        #         num_contacts += 1
        #     else:
        #         break
        num_protons = sum(sorted_allcoords_protons[:num_contacts+1])
           
        atom.numsurroundatoms = num_contacts
        atom.numsurroundprotons = num_protons
        del num_contacts,num_protons,sorted_allcoords_protons,sorted_dist

    print '\n---> success!'



def bdamage_calculate(PDBarray):

    # function to calculate Bdamage style metric for each atom, to save bdam 
    # attribute for each atom
    print '\n•••••••••••••••••••••••••••••••••••••••••••••••••••••••'
    print 'Calculating bdam style metric for atoms in structure...\n'

    # first order by number of surrounding atoms
    PDBarray.sort(key=lambda x: x.numsurroundatoms)
    num_atoms = len(PDBarray)

    # now loop through atoms and find number of atoms in same packing density bin
    counter = 0
    atom_indices = range(0,len(PDBarray))
    for atom in PDBarray:

        # unessential loading bar add-in
        progress(counter+1, num_atoms, suffix='')
        counter += 1

        simpacking_bfactors = []
        atom_numsurroundatoms = atom.numsurroundatoms
        k = -1
        # unwantedindices list is designed to locate any atoms which have packing density
        # below the current value and remove them from the subsequent loops
        unwantedindices = []
        for atomindex in atom_indices:

            k += 1
            otheratom = PDBarray[atomindex]

            if round(atom_numsurroundatoms/10) == round(otheratom.numsurroundatoms/10):
                simpacking_bfactors.append(float(otheratom.Bfactor))
            # since atoms ordered by number of surrounding atoms, this part breaks out of 
            # loop for current atom as soon as packing density bin is larger than that of
            # current atom    
            elif round(atom_numsurroundatoms/10) < round(otheratom.numsurroundatoms/10):
                break
            else:
                unwantedindices.append(k)

        # remove the indices from the search here if unwantedindices list is nonempty
        if len(unwantedindices) != 0:
            atom_indices = [i for j, i in enumerate(atom_indices) if j not in unwantedindices]


        bdam = float(atom.Bfactor)/(np.mean(simpacking_bfactors))
        
        atom.bdam = bdam
    print '\n---> success!'


def numsurroundatms_extract(initialPDBarray,laterPDBarray):
    # function to extract numbers for surrounding atoms from intial pdb structure
    # and extend these values to the same atoms in later pdb structures (for the 
    # same damage series)

    # loop through the later dataset and assign the corresponding num of 
    # neighbouring atoms from the same atom in the initial dataset.
    # Here the seenatoms list is filled as loop progresses to speed up loop
    # by ensuring that atom in initialPDBarray cannot be called again once it
    # has been found in the laterPDBarray list
    print '\n•••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••'
    print 'Extracting number of surrounding atoms from initial PDB file...\n'
    num_atoms = len(laterPDBarray)

    # ensure atom list ordered by number of atom in structure (atomnum)
    laterPDBarray.sort(key=lambda x: x.atomnum)
    initialPDBarray.sort(key=lambda x: x.atomnum)

    initfile_indices = range(0,len(initialPDBarray))
    counter = 0
    for atom in laterPDBarray:
        counter += 1

        # unessential loading bar add-in
        progress(counter, num_atoms, suffix='')

        k = -1
        for atomindex in initfile_indices:
            k += 1
            otheratom = initialPDBarray[atomindex]
            if (atom.atomtype == otheratom.atomtype and
               atom.basetype == otheratom.basetype and
               atom.residuenum == otheratom.residuenum and
               atom.chaintype == otheratom.chaintype):
                atom.numsurroundatoms = otheratom.numsurroundatoms
                atom.numsurroundprotons = otheratom.numsurroundprotons
                break
        initfile_indices.pop(k)  
    print '\n---> success!'
      
