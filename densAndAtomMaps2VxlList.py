# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 13:13:30 2015

@author: lina2532
"""
import sys 
from progbar import progress

# A class for .map file voxel
class voxel_density:
    def __init__(self,atmnum=0,density=0): 
        self.atmnum = atmnum
        self.density = density
            
def combinevxl_atmanddensvals(atmmap_list,densmap_list):
    # create list of voxel objects in class voxel_density 
    # NEW VERSION OF FUNCTION - tested on 10**6 length lists (2s vs 12s)
    vxl_list = {atm:[] for atm in atmmap_list}
    for atm,dens in zip(atmmap_list,densmap_list):
        vxl_list[atm].append(dens)
    
    return vxl_list
               
def combinevxl_atmanddensvals_gen(atmmap_list,densmap_list):
    # create list of voxel objects in class voxel_density 
        
    atmmap_len = len(atmmap_list)
    
    for i in range(0,atmmap_len):
         yield voxel_density(atmmap_list[i],densmap_list[i])

