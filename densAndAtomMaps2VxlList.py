# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 13:13:30 2015

@author: lina2532
"""
import sys 
from progbar import progress

###----------------###----------------###----------------###----------------###
###############################################################################      
# A class for .map file voxel
class voxel_density:
    def __init__(self,atmnum=0,density=0): 
        self.atmnum = atmnum
        self.density = density
###############################################################################
###----------------###----------------###----------------###----------------###

def combinevxl_atmanddensvals(atmmap_list,densmap_list):
    # create list of voxel objects in class voxel_density 
    
    vxl_list = []
    appendvxl = vxl_list.append
    
    atmmap_len = len(atmmap_list)
    
    for i in range(0,atmmap_len):
        appendvxl(voxel_density(atmmap_list[i],densmap_list[i]))

        # unessential progress bar here
        progress(i+1, atmmap_len, suffix='')
        
    return vxl_list
            
               
def combinevxl_atmanddensvals_gen(atmmap_list,densmap_list):
    # create list of voxel objects in class voxel_density 
        
    atmmap_len = len(atmmap_list)
    
    for i in range(0,atmmap_len):
         yield voxel_density(atmmap_list[i],densmap_list[i])

