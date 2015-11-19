# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 20:55:21 2014
@author: charlie
"""
import numpy as np
import math

###############################################################################
def map2CartCoords(rho,edensity_index,aCartesianVector,bCartesianVector,cCartesianVector):
      
    #currently, start positions ordered by columns, rows, sections as follows (fast,medium,slow axes): (bac)
    startposition_unordered = [rho.start1,rho.start2,rho.start3]

    #to change order to correspond to a,b,c directions of unit cell (in that order!) do the following: (abc)
    currentindexorder = [rho.fast_axis,rho.med_axis,rho.slow_axis]
    index_a = currentindexorder.index(min(currentindexorder))
    index_c = currentindexorder.index(max(currentindexorder))
    index_b = list(set([0,1,2])-set([index_a,index_c]))[0]

    STARTPOSITION_ordered = [startposition_unordered[index_a],startposition_unordered[index_b],startposition_unordered[index_c]]

    #similarly to above, the grid sampling is ordered by columns, rows, and sections as follows:
    gridsampling = [rho.gridsamp1,rho.gridsamp2,rho.gridsamp3]
    #want to convert this to a,b,c unit cell directions (by switching order!). 
    #Eg for P65 3CLC case, columns -> b, rows -> x, and sections -> z, so want
    #order of grid sampling to be rows,columns,sections for use below:
    #GRIDSAMPLING_ordered = [gridsampling[index_a],gridsampling[index_b],gridsampling[index_c]]
    GRIDSAMPLING_ordered = gridsampling

    #next need to convert the unit cell start a,b,c position into a x,y,z
    #cartesian position
    START = np.array(aCartesianVector)*(float(STARTPOSITION_ordered[0])/float(GRIDSAMPLING_ordered[0])) + np.array(bCartesianVector)*(float(STARTPOSITION_ordered[1])/float(GRIDSAMPLING_ordered[1])) + np.array(cCartesianVector)*(float(STARTPOSITION_ordered[2])/float(GRIDSAMPLING_ordered[2]))

    #find the cartesian position using the cartisian basis vectors. Note that 
    #index vector is written as [fastindex,medindex,slowindex] so corresponds 
    #to columns, rows, sections or equivalently b ,a, c (for P65 3CLC case). To 
    #use below, index order must be converted to a,b,c
    #relativeCartesianposition = numpy.array(aCartesianVector)*(float(edensity_index[index_a])/float(GRIDSAMPLING_ordered[0])) + numpy.array(bCartesianVector)*(float(edensity_index[index_b])/float(GRIDSAMPLING_ordered[1])) + numpy.array(cCartesianVector)*(float(edensity_index[index_c])/float(GRIDSAMPLING_ordered[2]))
    #relativeCartesianposition = numpy.array(aCartesianVector)*float(edensity_index[index_a]) + numpy.array(bCartesianVector)*float(edensity_index[index_b]) + numpy.array(cCartesianVector)*float(edensity_index[index_c])
    x = float(edensity_index[index_a])/(rho.gridsamp1)
    y = float(edensity_index[index_b])/(rho.gridsamp2)
    z = float(edensity_index[index_c])/(rho.gridsamp3)

    rrr = np.array([aCartesianVector,bCartesianVector,cCartesianVector])
    sss = np.transpose(rrr)

    xyz = np.dot(sss,np.array([x,y,z]))

    #add these coordinates to the 'position' matrix. Note, must be
    #translated to START to remove relative scale
    # position = relativeCartesianposition + START
    return xyz
    ###############################################################################
    
    
    
    
###############################################################################
def findCartesianTranslationVectors(a,b,c,alpha,beta,gamma):
    #Function to convert from unit cell parameters to cartesian coordinates

    #convert the angles from degrees to radians
    alpha = math.radians(alpha)
    beta = math.radians(beta)
    gamma = math.radians(gamma)

    #Define parameter v which is the final input of the conversion matrix below
    v = math.sqrt(1 - math.pow((math.cos(alpha)),2) - math.pow((math.cos(beta)),2) - math.pow((math.cos(gamma)),2) + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma))

    #Define the elements of the conversion matrix
    a11 = a
    a12 = b * math.cos(gamma)
    a13 = c * math.cos(beta)

    a21 = 0
    a22 = b * math.sin(gamma)
    a23 = c * ((math.cos(alpha) - (math.cos(beta) * math.cos(gamma)))/(math.sin(gamma)))

    a31 = 0
    a32 = 0
    a33 = c * (v / math.sin(gamma))

    #Convert the lattice basis vector in each direction to cartesian coordinates
    aCartesianVector = [a11,a21,a31]
    bCartesianVector = [a12,a22,a32]
    cCartesianVector = [a13,a23,a33]

    return aCartesianVector, bCartesianVector, cCartesianVector
###############################################################################
