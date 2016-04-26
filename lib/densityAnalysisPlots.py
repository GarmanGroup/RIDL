from matplotlib import pyplot as plt
from scipy import stats  
import sys

def edens_scatter(outputDir = './',
                  metrics   = ['meandensity','mindensity'],
                  PDBarray  = [],
                  pdbName   = 'untitled',
                  fileType  = '.png',
                  printText = False,
                  savefig   = True,
                  titleFont = 20,
                  axesFont  = 18,
                  edgeColor = '#FFFFFF',
                  spotColor = '#a93b1c'):

    # plot scatter plot of two selected metrics. 'metrics' of 
    # form ['meandensity','mediandensity'] to choose two metrics 
    # of electron density per atom to plot against each other 
    # in a scatter plot.
    
    valsPerParam = []
    for metric in metrics:   
        valPerAtom = [getattr(atom,metric) for atom in PDBarray]
        valsPerParam.append(valPerAtom)
        
    # check two generated lists of same length:
    if len(valsPerParam[0]) != len(valsPerParam[1]):
        print 'Error: lists not same length for scatter plot'
        sys.exit()   
            
    scatter1 = plt.figure()
    plt.scatter(x          = valsPerParam[0],
                y          = valsPerParam[1],
                color      = spotColor,
                edgecolors = edgeColor,
                marker     = 'o',
                s          = 100)

    # calculate linear regression and R-squared value
    slope, intercept, r_value, p_value, std_err = stats.linregress(*valsPerParam)

    infoStr = '------------------------------------------------------\n'+\
              'Scatter plot: {} density vs {} density\n'.format(*metrics)+\
              'r-squared: {}\n'.format(r_value**2)+\
              'p-value: {}\n'.format(p_value)

    if printText is True:
        print infoStr

    scatter1.suptitle('{} vs {} density'.format(*metrics), fontsize = titleFont)      
    plt.xlabel('{} density'.format(metrics[0]),fontsize = axesFont)   
    plt.ylabel('{} density'.format(metrics[1]),fontsize = axesFont)

    figName = '{}plots/{}_{}_vs_{}{}'.format(outputDir,
                                             pdbName,
                                             metrics[0],
                                             metrics[1],
                                             fileType)
    if savefig is True:
        scatter1.savefig(figName)

    return infoStr





###############################################################################
def numneighbours_scatter(where,PDBarray,pdbname):
    # function plots a scatter plot of number of neighbouring atoms against
    # number of neighbouring protons to determine correlation between two metrics
    
    numatms_peratom = [atom.numsurroundatoms for atom in PDBarray]
    numprotons_peratom = [atom.numsurroundprotons for atom in PDBarray]

    #check two generated lists of same length:
    if len(numatms_peratom) != len(numprotons_peratom):
        print 'Error: lists not same length for scatter plot'
        sys.exit()   
            
    scatter1 = plt.figure()
    plt.scatter(numatms_peratom,numprotons_peratom)
    
    scatter1.suptitle('# neighbouring atoms vs protons',fontsize=20)      
    plt.xlabel('# atoms',fontsize=18)  
    plt.ylabel('# protons',fontsize=16)

    figname = where+pdbname+ '_'+'numneighbouratomsVprotons.png'
    scatter1.savefig(figname)
###############################################################################


###############################################################################
def Npercentmostdamaged(var,PDBarray,percent):
    densityfile = open('atom_edensity.txt','r')
    counter = 0
    num_atoms = len(PDBarray)
    res_list = []
    for line in densityfile.readlines():
        if str(var) == line.split()[0]:
            counter += 1
            if counter > num_atoms*(float(100-percent)/100):
                res_list.append(str(line.split()[4])[0:len(line.split()[4])-1])
    densityfile.close()
    
    aminoacids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    DNAbases = ['DA','DC','DG','DT']   
    restype_list  = aminoacids + DNAbases
    
    res_list_counter = [0]*len(restype_list)
    counter = 0
    for res in res_list:
        index = restype_list.index(res)
        newlist = []
        for i in range(0,index):
            newlist.append(res_list_counter[i])
        newlist.append(res_list_counter[index]+1)
        for i in range(index+1,len(restype_list)):
            newlist.append(res_list_counter[i])
        res_list_counter = newlist

    #want fractions of 100 for pie chart: 
    res_list_counter_normalised = [(float(x)*float(100)) / float(sum(res_list_counter)) for x in res_list_counter]
    explode = [0.2]*len(restype_list)
    pie1 = plt.figure()
    plt.pie(res_list_counter_normalised, explode = explode, labels=restype_list,
                autopct='%1.1f%%', shadow=True, startangle=90)


    pie1.suptitle(str(percent)+'%-tile ' +str(var) + ' density' ,fontsize=20)
    figname = str(percent)+'%' +str(var) + ' density_piechart'+'.png'
    pie1.savefig(figname)
###############################################################################
    

###############################################################################
def bfac_scatter(pdbfilename,PDBarray,mainorside,restypes,atomtypes,densmet,bmet):
    ##simply plots scatter plot of two numerical attributes of atom objects of StructurePDB class
    #'PDBarray' is list of atoms
    #'mainorside' specifies 'mainchain', 'sidechain'
    #'restypes' specifies 'GLU','MET','HIS','DA',... etc
    #'atomtypes' specifies atom types eg: 'C','O','S',...etc
    #'densmet' specifies electron density metric per atom ('mean','min',...etc)
    #'bmet' specifies the Bfactor associated metric per atom ('Bfactor','Bdamage',...etc)
    
    X = []
    Y = []
    
    for atom in PDBarray:
        if atom.side_or_main() in (mainorside) and atom.basetype in (restypes) and atom.atomidentifier in (atomtypes):
            if densmet in ('mean','Mean'):
                X.append(atom.meandensity) 
            elif densmet in ('median','Median'):
                X.append(atom.mediandensity) 
            elif densmet in ('min','Min'):
                X.append(atom.mindensity) 
            elif densmet in ('max','Max'):
                X.append(atom.maxdensity) 
            else:
                print 'Unrecognised variable name, stopping program!'
                sys.exit()
            
            if bmet in ('Bfactor'):
                Y.append(atom.Bfactor) 
            elif bmet in ('Bfactorchange'):
                Y.append(atom.Bfactorchange) 
            elif bmet in ('Bdamage'):
                Y.append(atom.bdam) 
            elif bmet in ('Bdamagechange'):
                Y.append(atom.bdamchange) 
            else:
                print 'Unrecognised variable name, stopping program!'
                sys.exit()
               
    scatter1 = plt.figure()
    plt.scatter(X,Y)    
    scatter1.suptitle(str(densmet) + ' electron density' + ' vs ' + str(bmet),fontsize=20)
    plt.xlabel(str(densmet) + ' electron density',fontsize=18)
    plt.ylabel(str(bmet),fontsize=16)
    
    figname = str(pdbfilename) + '_' + str(densmet) + ' electron density' + ' vs ' + str(bmet) + '.jpg'
    scatter1.savefig(figname)
############################################################################### 
