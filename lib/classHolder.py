# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 01:48:33 2014
@author: charlie
"""
import math
import sys
import numpy as np

class residue:
    # A class for residue information
    def __init__(self,name = "",freq = 0,res_list = [],quantity = 0,av_mean_density = [],
                 atom_stand_dev = []):
            
        self.name               = name
        self.freq               = freq
        self.res_list           = res_list
        self.quantity           = quantity
        self.av_mean_density    = av_mean_density
        self.atom_stand_dev     = atom_stand_dev


class StructurePDB(object):
    # A class for coordinate PDB file atom
    def __init__(self,atomnum = 0,residuenum = 0,atomtype = "",basetype = "",chaintype = "",
                 X_coord = 0,Y_coord = 0,Z_coord = 0,atomID = "",numsurroundatoms = 0,
                 numsurroundprotons = 0,vdw_rad = 0,atomOrHetatm = "" ):
            
        self.atomnum            = atomnum
        self.residuenum         = residuenum
        self.atomtype           = atomtype
        self.basetype           = basetype 
        self.chaintype          = chaintype
        self.X_coord            = X_coord
        self.Y_coord            = Y_coord
        self.Z_coord            = Z_coord
        self.atomID             = atomID
        self.numsurroundatoms   = numsurroundatoms
        self.numsurroundprotons = numsurroundprotons
        self.vdw_rad            = vdw_rad
        self.atomOrHetatm       = atomOrHetatm

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
        vdw      = 'notyetdefined'

        try: # only proceed if file actually found
            VDWradfile  = open(str(filename),'r')
        except IOError:
            self.vdw_rad = 'Unable to locate {} file'.format(filename)
            return

        for line in VDWradfile.readlines():
            a = str(line.split()[1]).lower()
            b = str(self.atomID).lower()
            if a == b:
                # convert VDV radius for element into Angstroms as required:
                vdw = float(line.split()[5])*(0.01)
                break
            else:
                pass
        if vdw == 'notyetdefined':
            print 'Error finding vdw radius from reference file'
            print 'Check reference file...'
            print 'Atom identity given by: {} {}'.format(self.atomnum,self.getAtomID())
            sys.exit()
        VDWradfile.close()
        self.vdw_rad = vdw

    def getAtomID(self):
        # get a unique identifier for atom within structure, not dependent on atom number
        # which may differ between different datasets of same structure
        ID = '{}-{}-{}-{}'.format(self.chaintype,self.residuenum,self.basetype,self.atomtype)
        return ID

  
class singlePDB(StructurePDB):
    # StructurePDB subclass for a single pdb file structure
    def __init__(self,atomnum = 0,residuenum = 0,atomtype = "",basetype = "",chaintype = "",
                 X_coord = 0,Y_coord = 0,Z_coord = 0,Bfactor = 0,Occupancy = 0,meandensity = 0,
                 maxdensity = 0,mindensity = 0,mediandensity = 0,atomID = "",
                 numsurroundatoms = 0,numsurroundprotons = 0,vdw_rad = 0,atomOrHetatm = "",bdam = 0,
                 bdamchange = 0,Bfactorchange = 0,numvoxels = 0,stddensity = 0,min90tile = 0,
                 max90tile = 0,min95tile = 0,max95tile = 0):
                    
        super(singlePDB, self).__init__(atomnum,residuenum,atomtype,basetype,chaintype,X_coord,Y_coord,
                                        Z_coord,atomID,numsurroundatoms,numsurroundprotons,vdw_rad,atomOrHetatm)
            
        self.Bfactor        = Bfactor 
        self.Occupancy      = Occupancy
        self.meandensity    = meandensity
        self.maxdensity     = maxdensity
        self.mindensity     = mindensity
        self.mediandensity  = mediandensity
        self.bdam           = bdam
        self.bdamchange     = bdamchange
        self.Bfactorchange  = Bfactorchange        
        self.numvoxels      = numvoxels
        self.stddensity     = stddensity  
        self.min90tile      = min90tile
        self.max90tile      = max90tile
        self.min95tile      = min95tile
        self.max95tile      = max95tile
    
    def vdw_bfac(self):
        return 4*(math.sqrt(float(self.Bfactor) + 25))/(2 * math.pi)

    def getAdditionalMetrics(self):
        self.rsddensity     = float(self.stddensity)/self.meandensity
        self.rangedensity   = np.linalg.norm(self.maxdensity - self.mindensity)
      
  
# A subclass for a collection of multiple different dose pdb file structures
class multiPDB(StructurePDB):
    
    def __init__(self,atomnum = 0,residuenum = 0,atomtype = "",basetype = "",chaintype = "",
                 X_coord = 0,Y_coord = 0,Z_coord = 0,Bfactor = [],Occupancy = [],meandensity = [],
                 maxdensity = [],mindensity = [],mediandensity = [],atomID = "",
                 numsurroundatoms = 0,numsurroundprotons = 0,bdam = [],bdamchange = [],Bfactorchange = [],
                 meandensity_norm = [],maxdensity_norm = [],mindensity_norm = [],
                 mediandensity_norm = [],numvoxels = [],stddensity = [],min90tile = [],max90tile = [],
                 min95tile = [],max95tile = [],rsddensity = [],rangedensity = []):
            
        super(multiPDB, self).__init__(atomnum,residuenum,atomtype,basetype,chaintype,X_coord,
                                       Y_coord,Z_coord,atomID,numsurroundatoms,numsurroundprotons)   
            
        self.Bfactor            = Bfactor 
        self.Occupancy          = Occupancy
        self.meandensity        = meandensity
        self.maxdensity         = maxdensity
        self.mindensity         = mindensity
        self.mediandensity      = mediandensity
        self.bdam               = bdam
        self.bdamchange         = bdamchange
        self.Bfactorchange      = Bfactorchange
        self.meandensity_norm   = meandensity_norm
        self.maxdensity_norm    = maxdensity_norm
        self.mindensity_norm    = mindensity_norm
        self.mediandensity_norm = mediandensity_norm
        self.numvoxels          = numvoxels
        self.stddensity         = stddensity
        self.min90tile          = min90tile
        self.max90tile          = max90tile
        self.min95tile          = min95tile
        self.max95tile          = max95tile
        self.rsddensity         = rsddensity
        self.rangedensity       = rangedensity


class MapInfo:
    # A class for .map file header info
    def __init__(self, nx = 0, ny = 0, nz = 0, mapType = 0, start1 = 0, start2 = 0, start3 = 0,
                 gridsamp1 = 0, gridsamp2 = 0, gridsamp3 = 0, celldim_a = 0, celldim_b = 0,
                 celldim_c = 0, celldim_alpha = 0, celldim_beta = 0, celldim_gamma = 0,
                 fast_axis = 0, med_axis = 0, slow_axis = 0, mindensity = 0, maxdensity = 0,
                 meandensity = 0,vxls_val = []):
            
        self.nxyz       = {'nx':nx,'ny':ny,'nz':nz}
        self.type       = mapType
        self.start      = {'1':start1,'2':start2,'3':start3}
        self.gridsamp   = {'1':gridsamp1,'2':gridsamp2,'3':gridsamp3}
        self.celldims   = {'a':celldim_a,'b':celldim_b,'c':celldim_c,
                          'alpha':celldim_alpha,'beta':celldim_beta,'gamma':celldim_gamma}
        self.axis       = {'fast':fast_axis,'med':med_axis,'slow':slow_axis}
        self.density    = {'min':mindensity,'max':maxdensity,'mean':meandensity}
        self.vxls_val   = vxls_val

    def getnxyzString(self):
        return self.getParamString('Num. Col, Row, Sec',
                                   ['nx','ny','nz'],self.nxyz)

    def getStartString(self):
        return self.getParamString('Start positions',
                                   ['1','2','3'],self.start)

    def getGridsampString(self):
        return self.getParamString('Grid sampling',
                                   ['1','2','3'],self.gridsamp)

    def getCelldimsString(self):
        title = 'Cell dimensions'
        s = self.getParamString(title,['a','b','c'],self.celldims) +'\n'
        s += self.getParamString(' '*len(title),['alpha','beta','gamma'],self.celldims)
        return s

    def getAxisString(self):
        return self.getParamString('Fast,med,slow axes',
                                   ['fast','med','slow'],self.axis)

    def getDensityString(self):
        return self.getParamString('Values (min, max, mean)',
                                   ['min','max','mean'],self.density)

    def getParamString(self,title,orderedKeys,mapAttr):
        s = title +(40-len(title))*'.'
        for val in orderedKeys: s += str(mapAttr[val])+' '
        return s

    def getHeaderInfo(self):
        s =  self.getnxyzString() +'\n'
        s += self.getStartString() +'\n'
        s += self.getGridsampString() +'\n'
        s += self.getCelldimsString() +'\n'
        s += self.getAxisString() +'\n'
        s += self.getDensityString()
        return s



