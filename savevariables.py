# -*- coding: utf-8 -*-
"""
Created on Sat Jan  3 22:33:08 2015

@author: charlie
"""
import cPickle as pickle
import sys
from progbar import progress

# short functions to save PDB list of atom objects containing info on 
# each dose level to save running time consuming scripts calculating 
# electron density values to each atom in structure

def save_objectlist(PDBlist,pdbname):
    # to save current PDB list in a file:

    filename = str(len(PDBlist))+'_'+str(pdbname)+'_data.pkl'
    with open(filename, 'wb') as output:
        for atom in PDBlist:
            pickle.dump(atom, output, pickle.HIGHEST_PROTOCOL)
            
    return filename

def retrieve_objectlist(filename):
    # this function retrieves a list of objects from a file, given name
    # of form filename = str(len(PDBlist))+'_'+str(pdbname)+'_data.pkl'
    print 'Retrieving dataset from .pkl file...'
    
    #to determine number of atoms saved to file from file name:
    num_atoms = (filename.split('/')[-1]).split('_')[0]
    print '\nNumber of atoms in file: ' + str(num_atoms)
        
    #to retrieve list from file to new list:
    PDBretrieved = []    
    with open(str(filename), 'rb') as input:
        for i in range(0,int(num_atoms)):
            atom = None
            atom = pickle.load(input)
            PDBretrieved.append(atom)

            # unessential loading bar add-in
            progress(i+1, num_atoms, suffix='')

                            
    # return the list of atom objects  
    PDBretrieved.sort(key=lambda x: x.atomnum)
    print '\n---> success!'
    
    return PDBretrieved
