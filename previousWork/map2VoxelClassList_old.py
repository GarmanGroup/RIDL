# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 00:14:01 2015

@author: charlie
"""
import struct
from classHolder import MapInfo
import os
import sys

class voxel_density:
    # A class for .map file voxel
    def __init__(self,x=0,y=0,z=0,density=0,atmnum=[],atm_dist=[]):
        self.x          = x
        self.y          = y
        self.z          = z
        self.density    = density
        self.atmnum     = atmnum
        self.atm_dist   = atm_dist

def readMap(where,pdbname,mapfilename,maptype,atom_indices):

    # append to log file for this eTrack run 
    logfileName = '{}output/{}_log.txt'.format(where,pdbname)
    logfile = open(logfileName,'a')
    
    # define 'rho' electron map object
    rho = MapInfo()

    # open electron density .map file here (bmf for binary map file)
    bmf = open(where+mapfilename,'rb')
    
    # start adding header information into MapInfo class format. 
    # Note the unpacking of a struct for each byte, read as a long 'l'
    for n in ('nx','ny','nz'):
        rho.nxyz[n]   = struct.unpack('=l',bmf.read(4))[0]

    print 'Num. Col, Row, Sec: '
    print '{} {} {}'.format(*rho.nxyz.values())
    logfile.write('Num. Col, Row, Sec: {} {} {}\n'.format(*rho.nxyz.values()))
    
    rho.type = struct.unpack('=l',bmf.read(4))[0]

    for s in ('1','2','3'):
        rho.start[s] = struct.unpack('=l',bmf.read(4))[0] 

    print 'Start positions: '
    print '{} {} {}'.format(*rho.start.values())
    logfile.write('Start positions: {} {} {}\n'.format(*rho.start.values()))

    for g in ('1','2','3'):
        rho.gridsamp[g] = struct.unpack('=l',bmf.read(4))[0] 

    print 'Grid sampling:'
    print '{} {} {}'.format(*rho.gridsamp.values())
    logfile.write('Grid sampling: {} {} {}\n'.format(*rho.gridsamp.values()))

    # for cell dimensions, stored in header file as float not long 
    # integer so must account for this
    for d in ('a','b','c','alpha','beta','gamma'):
        rho.celldims[d]       = struct.unpack('f',bmf.read(4))[0]

    print 'Cell dimensions:'
    print '{} {} {}'.format(*rho.gridsamp.values()[0:3])
    print '{} {} {}'.format(*rho.gridsamp.values()[3:6])
    logfile.write('Cell dimensions: {} {} {} {} {} {}\n'.format(*rho.gridsamp.values()))

    for a in ('fast','med','slow'):
        rho.axis[a]   = struct.unpack('=l',bmf.read(4))[0] 

    print 'Fast,med,slow axes: '
    print '{} {} {}'.format(*rho.axis.values())
    logfile.write('Fast,med,slow axes: {} {} {}\n'.format(*rho.axis.values()))

    for d in ('min','max','mean'):
        rho.density[d]  = struct.unpack('f',bmf.read(4))[0] 

    print 'Density values: '
    print '{} {} {}'.format(*rho.density.values())
    logfile.write('Density values: {} {} {}\n'.format(*rho.density.values()))

    # next find .map file size, to calculate the last nx*ny*nz bytes of 
    # file (corresponding to the position of the 3D electron density 
    # array). Note factor of 4 is included since 4-byte floats used for 
    # electron density array values.
    filesize = os.path.getsize(where+mapfilename)
    densitystart = filesize - 4*(rho.nx*rho.ny*rho.nz)
    
    # next seek start of electron density data
    bmf.seek(densitystart,0)
    
    # if electron density written in shorts (I don't think it will be 
    # but best to have this case)
    if rho.type is 1:  
        print 'Untested .map type --> 1 (type 2 expected). Values read' 
        + 'as int16 ("i2?") - consult .map header in MAPDUMP(CCP4) to check'
  
        struct_fmt = '=i2'
        struct_len = struct.calcsize(struct_fmt)
        density = []
       
        while True:
            data = bmf.read(struct_len)
            if not data: break
            s = struct.unpack(struct_fmt,data)[0]
            density.append(s)    
    
    # if electron density written in floats (which is to be expected 
    # from FFT-CCP4 outputted .map file of electron density)
    if rho.type is 2:  
        struct_fmt = '=f4'
        struct_len = struct.calcsize(struct_fmt)
        density = []
        appenddens = density.append
        
        if maptype in ('atom_map'):
            atom_indices = []
            appendindex = atom_indices.append
            counter = -1
            while True:
                data = bmf.read(struct_len)
                if not data: break
                s = struct.unpack(struct_fmt,data)[0]
                counter += 1
                if int(s) == 0:
                    continue
                else:    
                    appenddens(s)
                    appendindex(counter)
        
        # efficient way to read through density map file using indices of atoms
        # from atom map file above
        elif maptype in ('density_map'):
            for i in range(0,len(atom_indices)):
                if i != 0:
                    bmf.read(struct_len*(atom_indices[i]-atom_indices[i-1] - 1))
                else:
                    bmf.read(struct_len*(atom_indices[0]))
                    
                data = bmf.read(struct_len)
                s = struct.unpack(struct_fmt,data)[0]
                appenddens(s)
                
        else:
            print 'Unknown map type --> terminating script'
            sys.exit()
         
    logfile.close()           
    bmf.close()
    
    # as a check that file has been read correctly, check that the min 
    # and max electron density stated in .map file header correspond to 
    # calculated min and max here
    # Note for the case of atom map, the min voxel val will be 0 (no atom present)
    # --> since these voxels are removed during the filtering of the map, only 
    # the map value is tested.
    # For density map, cannot currently perform a check, since there is no 
    # guarantee that the max and min voxel values may be non atom voxels and
    # thus removed
    if maptype in ('atom_map'):        
        if max(density) == rho.maxdensity:
            print 'calculated max voxel value match value stated in file header'
        else:
            print 'calculated max voxel value does NOT match value stated in file header'
            print 'have now calculated max voxel value to be: %s'\
            %str(max(density))
            sys.exit()
    
    # if each voxel value is an atom number, then want to convert to integer
    if maptype in ('atom_map'):
        density_final = [int(dens)/100 for dens in density]
    elif maptype in ('density_map'):
        density_final = density
    else:
        print 'Unknown map type --> terminating script'
        sys.exit()

    rho.vxls_val = density_final
    
    if maptype in ('atom_map'):
        return rho,atom_indices
    else:
        return rho
        


###----------------###----------------###----------------###----------------###
###############################################################################
def densmap2class_plot(voxel_list):
    
    # the following is designed to plot a histogram plot of the distibution 
    # of density values in the collection of voxels in voxel_list:
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    sns.set_palette("deep", desat=.6)
    sns.set_context(rc={"figure.figsize": (8, 4)})
    np.random.seed(9221999)
    
    # find densities from voxels
    density = [voxel.density for voxel in voxel_list]
        
    # plot a histogram of 'density' list containing all densities above:     
    plt.hist(density, 100, histtype="stepfilled", alpha=.7)
###############################################################################
###----------------###----------------###----------------###----------------###