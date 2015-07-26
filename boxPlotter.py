# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 23:56:28 2014

@author: charlie
"""
import numpy as np 
import matplotlib as mpl 

## agg backend is used to create plot as a .png file
mpl.use('agg')

import matplotlib.pyplot as plt 

def residue_boxplotter(densitylist_atomorder,atoms_labels,residue_name,quantity):

    # Create a figure instance
    fig = plt.figure()
    
    # Create an axes instance
    ax = fig.add_subplot(111)
    
    ## add patch_artist=True option to ax.boxplot() 
    ## to get fill color
    bp = ax.boxplot(densitylist_atomorder, patch_artist=True)
    
    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = '#1b9e77' )
    
    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)
    
    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)
    
    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)
    
    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)
    
    ## Custom x-axis labels
    ax.set_xticklabels(atoms_labels)
    
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    ## Custom title,x-label,y-label    
    fig.suptitle(str(residue_name)+' sampled over '+str(quantity)+' residues',
                 fontsize=20)
    plt.xlabel('Atoms', fontsize=18)
    plt.ylabel('Electron Density Change', fontsize=16)
    
    ## Save the figure
    fig.savefig(str(residue_name)+'_boxplot.png', bbox_inches='tight')
    
    

import seaborn as sns
import matplotlib.pyplot as plt

def residue_violinplotter(where,densitylist_atomorder,atoms_labels,
                          residue_name,quantity,pdbname,densmet):

    # Create a figure instance
    fig = plt.figure()
    
    # Create an axes instance
    ax = fig.add_subplot(111)    
    sns.violinplot(densitylist_atomorder, color="coolwarm_r", lw=2);
    
    ## Custom x-axis labels
    ax.set_xticklabels(atoms_labels)
    
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    ## Custom title,x-label,y-label    
    fig.suptitle(str(residue_name)+' sampled over '+str(quantity)+' residues',
                 fontsize=20)
                 
    plt.xlabel('Atoms', fontsize=18)
    plt.ylabel('Electron Density Change', fontsize=16)
    
    ## Save the figure
    fig.savefig(str(where)+'output/plots/' + str(pdbname) + '_'+str(residue_name)+'_'+str(densmet)+'_boxplot.png',
                bbox_inches='tight')
    
    
