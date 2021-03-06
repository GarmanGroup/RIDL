from __future__ import division
from struct import unpack, calcsize
from classHolder import MapInfo
from errors import error
import os
import mmap
import numpy as np
from functools import reduce


def readMap(dirIn='./', dirOut='./', mapName='untitled.map',
            mapType='atom_map', atomInds=[], log='',
            standardise=False, fixMaxMapVal=''):

    # a function to read in a .map file of either density or atom-tagged type

    # define 'rho' electron map object
    rho = MapInfo()

    mapName = dirIn + mapName
    filesize = os.path.getsize(mapName)
    log.writeToLog(str='Map file of size {} bytes to be read'.format(filesize))

    # open electron density .map file here (bmf for binary map file)
    if os.name != 'nt':
        with open(mapName) as f:
            bmf = mmap.mmap(
                f.fileno(), 0, prot=mmap.PROT_READ, flags=mmap.MAP_PRIVATE)
    else:
        with open(mapName, "r+b") as f:
            bmf = mmap.mmap(f.fileno(), 0)

    # start adding header information into MapInfo class format.
    # Note the unpacking of a struct for each byte, read as a long 'l'
    for n in ('nx', 'ny', 'nz'):
        rho.nxyz[n] = unpack('=l', bmf.read(4))[0]

    rho.type = unpack('=l', bmf.read(4))[0]

    for s in ('fast', 'med', 'slow'):
        rho.start[s] = unpack('=l', bmf.read(4))[0]

    for g in ('1', '2', '3'):
        rho.gridsamp[g] = unpack('=l', bmf.read(4))[0]

    for d in ('a', 'b', 'c', 'alpha', 'beta', 'gamma'):
        # cell dims stored as float not long int
        rho.celldims[d] = unpack('f', bmf.read(4))[0]

    for a in ('fast', 'med', 'slow'):
        rho.axis[a] = unpack('=l', bmf.read(4))[0]

    for d in ('min', 'max', 'mean'):
        rho.density[d] = unpack('f', bmf.read(4))[0]

    s = rho.getHeaderInfo(tab=True)

    # write to log file if specified
    if log != '':
        if os.path.exists(log.logFile):
            for l in s.split('\n'):
                log.writeToLog(l, priority='minor')
    else:
        print(s)

    # calculate the last nx*ny*nz bytes of file (corresponding to
    # the position of the 3D electron density array). Note factor
    # of 4 is included since 4-byte floats used for electron
    # density array values.

    densitystart = filesize - 4*(reduce(lambda x, y: x*y, list(rho.nxyz.values())))

    # if sys.version_info[0] >= 3:
    #     import functools
    #     densitystart = filesize - 4*(functools.reduce(lambda x, y: x*y, list(rho.nxyz.values())))
    # else:
    #     densitystart = filesize - 4*(reduce(lambda x, y: x*y, list(rho.nxyz.values())))

    # get symmetry operations from map header
    for j in range(23, 58):
        if j != 53 and j < 57:

            if j != 55:
                val = unpack('=l', bmf.read(4))[0]
            else:
                val = unpack('=f', bmf.read(4))[0]

            if j == 24:
                numSymBytes = val

            if j == 55:
                mapStdev = round(val, 5)

        elif j < 57:
            for i in range(4):
                val = unpack('c', bmf.read(1))[0]
        else:
            for i in range(10):
                for k in range(80):
                    char = unpack('c', bmf.read(1))[0]
    symOps = []

    for j in range(numSymBytes // 80):
        line = ''
        for k in range(80):
            char = unpack('=c', bmf.read(1))[0].decode("utf-8")
            line += char
        symOps.append(line)
    rho.curateSymOps(symOps)

    # next seek start of electron density data
    bmf.seek(densitystart, 0)

    # if electron density written in shorts this is not
    # currently expected, so program will halt
    if rho.type == 1:
        error(text='Untested .map type --> 1 (type 2 expected). ' +
                   'Values read as int16 ("i2?") - consult .map header in ' +
                   'MAPDUMP(CCP4) to check', log=log, type='error')

    # if electron density written in floats (which is to be expected
    # from FFT-CCP4 outputted .map file of electron density)
    if rho.type == 2:

        # the length of the each data chunk to read map voxels below
        struct_fmt = '=f'
        struct_len = calcsize(struct_fmt)
        # should get struct_len = 4
        if struct_len != 4:
            error(text='Bad data chunk length assigned when reading map',
                  log=log, type='error')

        density = []
        appenddens = density.append

        if mapType in ('atom_map'):
            atomInds = []
            appendindex = atomInds.append
            counter = -1
            while True:
                data = bmf.read(struct_len)
                if not data:
                    break
                s = unpack(struct_fmt, data)[0]
                counter += 1
                if int(s) == 0:
                    continue
                else:
                    appenddens(s)
                    appendindex(counter)
            log.writeToLog(str='# voxels in total : {}'.format(counter + 1))

        # efficient way to read through density map file
        # using indices of atoms from atom map file above
        elif mapType in ('density_map'):

            for i in range(0, len(atomInds)):
                if i != 0:
                    bmf.read(struct_len*(atomInds[i]-atomInds[i-1] - 1))
                else:
                    bmf.read(struct_len*(atomInds[0]))

                data = bmf.read(struct_len)
                s = unpack(struct_fmt, data)[0]
                appenddens(s)

            # check that resulting list of same length as atomInds
            if len(density) != len(atomInds):
                error(text='Failure to process the density map ' +
                           'using atom-tagged map', log=log, type='error')

            # UNCOMMENT TO INSTEAD READ IN ALL VOXELS
            # while True:
            #     data = bmf.read(struct_len)
            #     if not data: break
            #     s = unpack(struct_fmt,data)[0]
            #     appenddens(s)

        else:
            error(text='Unknown map type!', log=log, type='error')
    bmf.close()

    # as a check that file has been read correctly, check that the min
    # and max electron density stated in .map file header correspond to
    # calculated min and max here.
    # NOTE for the case of atom map, the min voxel val will be 0 (no
    # atom present)--> since these voxels are removed during the filtering
    # of the map, only the map value is tested.
    # For density map, cannot currently perform a check, since there is no
    # guarantee that the max and min voxel values may be non atom voxels and
    # thus removed
    if mapType in ('atom_map'):
        if max(density) == rho.density['max']:
            log.writeToLog(
                str='calculated max voxel value match value ' +
                    'stated in file header')
        else:
            error(
                text='Calculated max voxel value:{} does NOT '.format(
                    max(density)) +
                'match value stated in file header:{}'.format(
                    rho.density['max']),
                log=log, type='error')

    # if each voxel value is an atom number, then want to convert to integer
    if mapType in ('atom_map'):
        density_final = [int(dens // 100) for dens in density]
    elif mapType in ('density_map'):
        density_final = density
    else:
        error(text='Unknown map type!', log=log, type='error')

    # provide option to standardise density to map mean and standard deviation
    if standardise and mapType in ('density_map'):
        density_final = (np.array(density_final)-rho.density['mean'])/mapStdev

        if fixMaxMapVal != '':
            density_final /= np.max(np.absolute(density_final))
            if fixMaxMapVal == 'rand':
                fixMaxMapVal = np.random.uniform(6, 10)
                log.writeToLog(
                    str='Random max map value chosen to be {}'.format(
                        fixMaxMapVal))
            density_final *= fixMaxMapVal

        density_final = density_final.tolist()

    rho.vxls_val = density_final

    if mapType in ('atom_map'):
        return rho, atomInds
    else:
        return rho
