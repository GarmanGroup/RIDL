# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 10:14:29 2015

@author: charlie
"""
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import math

def vxlsperatm_hist(pdbname,where,vxlsperatom):
    # histogram plot of number of voxels per atom
    print '\n------------------------------------------------'
    print 'Plotting histogram of number of voxels per atom...' 
    sns.set_palette("deep", desat=.6)
    sns.set_context(rc={"figure.figsize": (10, 6)})
    fig = plt.figure()

    datax = [len(atm.vxls) for atm in vxlsperatom]

    plt.hist(datax, 300, histtype="stepfilled", alpha=.7);
    plt.xlabel('Voxels per atom')
    plt.ylabel('Frequency')
    plt.title('Histrogram of voxels per atom')
    fig.savefig(str(where)+'output/plots/'+str(pdbname)+'vxlsperatm_hist.png')
    

def vxlsperatm_kde(pdbname,where,vxlsperatom):
    # kde plot of number of voxels per atom
    print '\n------------------------------------------------'
    print 'Plotting kde plot of number of voxels per atom...' 
    sns.set_palette("deep", desat=.6)
    sns.set_context(rc={"figure.figsize": (10, 6)})
    fig = plt.figure()

    datax = [len(atm.vxls) for atm in vxlsperatom]

    sns.kdeplot(np.array(datax), shade=True);
    plt.xlabel('Voxels per atom')
    plt.ylabel('Frequency')
    plt.title('kde plot of voxels per atom')
    fig.savefig(str(where)+'output/plots/'+str(pdbname)+'vxlsperatm_kde.png')
