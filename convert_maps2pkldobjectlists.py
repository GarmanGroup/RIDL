# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 13:58:08 2015

@author: lina2532
"""

from CHECK_readinatommap_efficientdensmapreadin import maps2atmdensity
from savevariables import save_objectlist


def maps2pkldobjs(where,pdbname,mapfilname1,mapfilname2):
    # function read in two maps, create a list of atom objects (of class
    # StructurePDB) and then pickle dump this to the same folder as the script
   
   PDB_objlist = maps2atmdensity(where,pdbname,mapfilname1,
                           'atom_map',mapfilname2,'density_map')

   # pickle the big list of atom objects
   pkl_filename = save_objectlist(PDB_objlist,pdbname)
   
   return pkl_filename
    