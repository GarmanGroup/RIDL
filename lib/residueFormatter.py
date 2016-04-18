from boxPlotter import residue_violinplotter
import sys

class residuetype:
    # for a specific residue (eg. 'MET') this class defines objects
    # for the function below

    def __init__(self,
                 name        = "",
                 atm_names   = [],
                 atms_bytype = [],
                 frequency   = 0):
            
        self.name        = name
        self.atm_names   = atm_names
        self.atms_bytype = atms_bytype
        self.frequency   = frequency


def res2atomsbytype(PDBarray,res_type):
    # for a given residue type, read through PDBarray and output an object
    # (see residue type class above)
    
    # create an object for this residue type
    y = residuetype()
    y.name = res_type
    
    # list to collect all the distinct 'residue number' + 'chain type' strings.
    spec_res_indicator = []
    
    # determine atom types present in residue (distinguished by position: CA,CB,..)
    atms_inres = []
    for atom in PDBarray:
        if atom.basetype == res_type:
            if atom.atomtype not in atms_inres:
                atms_inres.append(atom.atomtype)
            if str(atom.residuenum)+':'+str(atom.chaintype) not in spec_res_indicator:
                spec_res_indicator.append(str(atom.residuenum)+':'+str(atom.chaintype))
                
    y.atm_names =  atms_inres        
    
    # create a 2d list with all atoms of the residue split into sublists
    # by atom type
    atms_ofres = []
    for atmtype in atms_inres:
        atms_ofres.append([atom for atom in PDBarray if 
                                        (atom.atomtype == atmtype and 
                                        atom.basetype == res_type)])
   
    # add above 2d list to residue object                                   
    y.atms_bytype = atms_ofres
    
    # determines the frequency of residue in structure
    res_freq = len(spec_res_indicator)
    y.frequency = res_freq
    
    return y    
    
def densper_resatom_NOresidueclass(where    = '',
                                   PDBarray = [],
                                   plot     = False,
                                   densMet  = 'mean',
                                   pdbName  = ''):
    # this function plots violin plots of atom density for each residue 
    # type, and also outputs a list of residue objects (see residuetype
    # class)

    # KEY:
    # 'PDBarray' is list of atoms of structure
    # 'plot' (boolian) determines whether individual boxplots plotted
    # 'densmet' takes values 'mean','median','max','min' to determine 
    # metric of electron density
        
    # first ensure PDB list ordered by atom number
    PDBarray.sort(key=lambda x: (x.basetype,x.atomnum))
      
    # find list of protein residues and nucleic acids
    # present in structure (with repeated names removed)
    resi_list = [atom.basetype for atom in PDBarray]
    restype_list = []
    appendres = restype_list.append                      
    for res in resi_list:
        if res not in restype_list:
            appendres(res)
    restype_list.sort()

    residueArray = []
    for res in restype_list:
        #create an object for the residue type (in residuetype class)
        res_obj = res2atomsbytype(PDBarray,res)
                 
        print '------------SUMMARY---------------'
        print 'Total number of residues of type {}: {}'.format(res_obj.name,str(res_obj.frequency))
        print 'Maximum size of residue calculated to be: {}'.format(res_obj.atm_names) 
      
        # determine the density metric to use and create list of densities
        # from 'atms_bytype' attribute for the residue object
        dens_list = []
        
        if densmet in ('mean'):
            for atmtype in res_obj.atms_bytype:
                dens_list.append([atm.meandensity for atm in atmtype])
                
        elif densmet in ('median'):
            for atmtype in res_obj.atms_bytype:
                dens_list.append([atm.mediandensity for atm in atmtype]) 
                
        elif densmet in ('min'):
            for atmtype in res_obj.atms_bytype:
                dens_list.append([atm.mindensity for atm in atmtype])
                
        elif densmet in ('max'):
            for atmtype in res_obj.atms_bytype:
                dens_list.append([atm.maxdensity for atm in atmtype])
                
        else:
            print 'Unknown electron density metric specified'
            sys.exit()         
                
        # call the boxplotting function to plot boxplot of atom type
        # against electron density change for given residue type
        if plot is True:
            print 'Plotting now...'
            residue_violinplotter(where,dens_list,
                                  res_obj.atm_names,
                                  res_obj.name,
                                  res_obj.frequency,
                                  pdbname,
                                  densmet)
        
        residueArray.append(res_obj)
    return residueArray
    
    
def densper_res(where      = '',
                resArray   = [],
                minResNum  = 0,
                sideormain = ['sidechain','mainchain'],
                densmet    = 'min',
                pdbName    = ''):  
    # plots a boxplot for each residue detailing the electron density 
    # distn.'minResNum' is the threshold for the min number of residues
    # of a given type that need to be present in structure to be included in 
    # the end plot. 'sideormain' specifies whether 'sidechain' only, 
    # 'mainchain' only, or ['sidechain','mainchain'] are selected 'densmet'
    # specifies the density metric to be used ('mean','median','max','min')

    residue_label = []

    # The next section of code plots a boxplot for each residue/base 
    # type in structure   
    densitylist = [] 
    for res in residueArray:
        # only residue types with more than a specified frequency will be
        # included in the plot
        if res.frequency > minResNum:
            
            print str(res.name) + ' -----> '\
            + str(res.frequency)

            # need to merge 2d list res_obj of atoms by atom type into 1d list
            merged_atms = sum(res.atms_bytype,[])
        
            # determine the density metric to use and create list of densities
            # from merged_atms
            if densmet in ('mean'):
                merged_dens = [atom.meandensity for atom in merged_atms 
                                if (atom.side_or_main() in sideormain)]
                                    
            elif densmet in ('median'):
                merged_dens = [atom.mediandensity for atom in merged_atms 
                                if (atom.side_or_main() in sideormain)]
                                    
            elif densmet in ('min'):
                merged_dens = [atom.mindensity for atom in merged_atms 
                                if (atom.side_or_main() in sideormain)]
                                    
            elif densmet in ('max'):
                merged_dens = [atom.maxdensity for atom in merged_atms 
                                if (atom.side_or_main() in sideormain)]
                                 
            else:
                print 'Unknown electron density metric specified'
                sys.exit()  
                        
            densitylist.append(merged_dens)
            
            # get list of residue names for boxplot x-axis label
            residue_label.append(res.name)
            
    print 'Plotting all present residues now......'
    residue_violinplotter(where,
                          densitylist,
                          residue_label,
                          densmet,
                          'ALL'+str(sideormain),
                          pdbName,
                          densmet)
