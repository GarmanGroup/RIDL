# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 01:48:33 2014
@author: charlie
"""
import math
import sys

###############################################################################
# A class for residue information
class residue:
    
    def __init__(
        self,name="",freq=0,res_list=[],quantity=0,av_mean_density=[],
        atom_stand_dev=[]):
            
        self.name = name
        self.freq = freq
        self.res_list = res_list
        self.quantity = quantity
        self.av_mean_density = av_mean_density
        self.atom_stand_dev = atom_stand_dev
###############################################################################
   
   

###############################################################################
# A class for structure PDB file
class StructurePDB(object):
    
    def __init__(
        self,atomnum=0,residuenum=0,atomtype="",basetype="",chaintype="",
        X_coord=0,Y_coord=0,Z_coord=0,atomidentifier="",numsurroundatoms=0,
        numsurroundprotons=0,vdw_rad=0):
            
        self.atomnum = atomnum
        self.residuenum = residuenum
        self.atomtype = atomtype
        self.basetype = basetype 
        self.chaintype = chaintype
        self.X_coord = X_coord
        self.Y_coord = Y_coord
        self.Z_coord = Z_coord
        self.atomidentifier = atomidentifier
        self.numsurroundatoms = numsurroundatoms
        self.numsurroundprotons = numsurroundprotons
        self.vdw_rad = vdw_rad

    # To determine whether the current atom is a constituent of protein or nucleic acid    
    def protein_or_nucleicacid(self):
        if self.basetype in ('DA','DC','DG','DT','A','C','G','U'):
            component_type = 'nucleic acid'
        else:
            component_type = 'protein'  
        return component_type
        
    def side_or_main(self):
        if self.basetype in ('DA','DC','DG','DT','A','C','G','U'):
            if self.atomtype not in ("P","OP1","O5'","C5'","C4'","C3'","O3'","C2'","C1'","O4'","OP2"):
                sideormain = 'sidechain'
            else:
                sideormain = 'mainchain'
        else:
            if self.atomtype not in ('N','CA','C','O'):
                sideormain = 'sidechain'
            else:
                sideormain = 'mainchain'
        return sideormain
            
    def VDW_get(self):
        # determine VDW radius for atom type
        filename = 'VDVradiusfile.txt'
        VDWradfile = open(str(filename),'r')
        vdw = 'notyetdefined'
        for line in VDWradfile.readlines():
            a = str(line.split()[1]).lower()
            b = str(self.atomidentifier).lower()
            if a == b:
                # convert VDV radius for element into Angstroms as required:
                vdw = float(line.split()[5])*(0.01)
                break
            else:
                pass
        if vdw == 'notyetdefined':
            print 'Error finding vdw radius from reference file'
            print 'Check reference file...'
            print 'Atom identity given by: %s %s %s %s %s' %(str(self.atomnum),str(self.atomtype),str(self.residuenum),str(self.basetype),str(self.chaintype))
            sys.exit()
        VDWradfile.close()
        self.vdw_rad = vdw
###############################################################################
     
      
      
###############################################################################  
# A subclass for a single pdb file structure
class singlePDB(StructurePDB):
    
    def __init__(
        self,atomnum=0,residuenum=0,atomtype="",basetype="",chaintype="",
        X_coord=0,Y_coord=0,Z_coord=0,Bfactor=0,Occupancy=0,meandensity=0,
        maxdensity=0,mindensity=0,mediandensity=0,atomidentifier="",
        numsurroundatoms=0,numsurroundprotons=0,vdw_rad=0,bdam=0,bdamchange=0,
        Bfactorchange=0,numvoxels=0,stddensity=0,min90tile=0,max90tile=0,
        min95tile=0,max95tile=0,modedensity=0,rsddensity=0,dipstat=0,rangedensity=0):
            
        super(singlePDB, self).__init__(
            atomnum,residuenum,atomtype,basetype,chaintype,X_coord,Y_coord,
            Z_coord,atomidentifier,numsurroundatoms,numsurroundprotons,vdw_rad)
            
        self.Bfactor = Bfactor 
        self.Occupancy = Occupancy
        self.meandensity = meandensity
        self.maxdensity = maxdensity
        self.mindensity = mindensity
        self.mediandensity = mediandensity
        self.bdam = bdam
        self.bdamchange = bdamchange
        self.Bfactorchange = Bfactorchange        
        self.numvoxels = numvoxels
        self.stddensity = stddensity  
        self.min90tile = min90tile
        self.max90tile = max90tile
        self.min95tile = min95tile
        self.max95tile = max95tile
        self.modedensity = modedensity
        self.rsddensity = rsddensity
        self.dipstat = dipstat
        self.rangedensity = rangedensity
    
    def vdw_bfac(self):
        return 4*(math.sqrt(float(self.Bfactor) + 25))/(2 * math.pi)
###############################################################################

      
      
###############################################################################  
# A subclass for a collection of multiple different dose pdb file structures
class multiPDB(StructurePDB):
    
    def __init__(
        self,atomnum=0,residuenum=0,atomtype="",basetype="",chaintype="",
        X_coord=0,Y_coord=0,Z_coord=0,Bfactor=[],Occupancy=[],meandensity=[],
        maxdensity=[],mindensity=[],mediandensity=[],atomidentifier="",
        numsurroundatoms=0,numsurroundprotons=0,bdam=[],bdamchange=[],Bfactorchange=[],
        meandensity_norm=[],maxdensity_norm=[],mindensity_norm=[],
        mediandensity_norm=[],numvoxels=[],stddensity=[],min90tile=[],max90tile=[],
        min95tile=[],max95tile=[],modedensity=[],rsddensity=[],dipstat=[],rangedensity=[]):
            
        super(multiPDB, self).__init__(
            atomnum,residuenum,atomtype,basetype,chaintype,X_coord,Y_coord,
            Z_coord,atomidentifier,numsurroundatoms,numsurroundprotons)   
            
        self.Bfactor = Bfactor 
        self.Occupancy = Occupancy
        self.meandensity = meandensity
        self.maxdensity = maxdensity
        self.mindensity = mindensity
        self.mediandensity = mediandensity
        self.bdam = bdam
        self.bdamchange = bdamchange
        self.Bfactorchange = Bfactorchange
        self.meandensity_norm = meandensity_norm
        self.maxdensity_norm = maxdensity_norm
        self.mindensity_norm = mindensity_norm
        self.mediandensity_norm = mediandensity_norm
        self.numvoxels = numvoxels
        self.stddensity = stddensity
        self.min90tile = min90tile
        self.max90tile = max90tile
        self.min95tile = min95tile
        self.max95tile = max95tile
        self.modedensity = modedensity
        self.rsddensity = rsddensity
        self.dipstat = dipstat
        self.rangedensity = rangedensity
###############################################################################
      



###############################################################################      
#A class for .map file header info
class electron_map_info:
    def __init__(
        self, nx=0, ny=0, nz=0, type=0, start1=0, start2=0, start3=0,
        gridsamp1=0, gridsamp2=0, gridsamp3=0, celldim_a = 0, celldim_b = 0,
        celldim_c = 0, celldim_alpha = 0, celldim_beta = 0, celldim_gamma = 0,
        fast_axis=0, med_axis=0, slow_axis=0, mindensity=0, maxdensity=0,
        meandensity=0,vxls_val=[]):
            
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.type = type
        self.start1 = start1
        self.start2 = start2
        self.start3 = start3
        self.gridsamp1 = gridsamp1
        self.gridsamp2 = gridsamp2
        self.gridsamp3 = gridsamp3
        self.celldim_a = celldim_a
        self.celldim_b= celldim_b
        self.celldim_c = celldim_c
        self.celldim_alpha = celldim_alpha
        self.celldim_beta = celldim_beta
        self.celldim_gamma = celldim_gamma
        self.fast_axis = fast_axis
        self.med_axis = med_axis
        self.slow_axis = slow_axis
        self.mindensity = mindensity
        self.maxdensity = maxdensity
        self.meandensity = meandensity
        self.vxls_val = vxls_val
###############################################################################
