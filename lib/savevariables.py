from future import standard_library
standard_library.install_aliases()
from progbar import progress

import sys
if sys.version_info[0] < 3:
    import cPickle as pickle
else:
    import pickle

# small functions to save PDB list of atom objects containing info on
# each dose level to save running time consuming scripts calculating
# electron density values to each atom in structure


def save_objectlist(PDBlist, pdbName):
    # to save current PDB list in a file:
    filename = str(len(PDBlist))+'_'+pdbName+'_data.pkl'
    with open(filename, 'wb') as output:
        for atom in PDBlist:
            pickle.dump(atom, output, pickle.HIGHEST_PROTOCOL)
    return filename


def retrieve_objectlist(fileName='untitled.pkl', loadBar=False, logFile=''):

    # this function retrieves a list of objects
    # from a file, given name of form filename =
    # str(len(PDBlist))+'_'+str(pdbName)+'_data.pkl'

    ln = 'Retrieving dataset from .pkl file...'
    if logFile != '':
        logFile.writeToLog(str=ln)
    else:
        print(ln)

    checkFileFormat(fileName)

    # to determine number of atoms saved to file from file name:
    num_atoms = (fileName.split('/')[-1]).split('_')[0]
    ln = 'Number of atoms in file: ' + str(num_atoms)

    if logFile != '':
        logFile.writeToLog(str=ln)
    else:
        print(ln)

    # to retrieve list from file to new list:
    PDBretrieved = []
    with open(fileName, 'rb') as input:
        for i in range(0, int(num_atoms)):
            atom = None
            atom = pickle.load(input)
            PDBretrieved.append(atom)

            if loadBar:
                progress(i+1, num_atoms, suffix='')

    # return the list of atom objects
    PDBretrieved.sort(key=lambda x: x.atomnum)

    return PDBretrieved


def saveGenericObject(obj=[], fileName='untitled'):

    # save a generic object to file 'fileName'_data.pkl

    filename = fileName+'_data.pkl'
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
    return filename


def checkFileFormat(fileName):
    if fileName.split('_')[-1] != 'data.pkl':
        sys.exit('.pkl file of wrong format to retrieve. ' +
                 'File {} supplied'.format(fileName))


def retrieveGenericObject(
                          fileName='untitled.pkl'):

    # retrieve a generic object from
    # file name '{fileName}_data.pkl'

    checkFileFormat(fileName)
    with open(fileName, 'rb') as input:
        obj = pickle.load(input)
    return obj
