# ETRACK

A collection of scripts to calculate per-atom density change metrics within a specific damage MX investigation.
Suitable for any MX experiment in which datasets are collected on the **same crystal** over **multiple doses**.
**NOTE: These scripts are currently under development and updated regularly..**

## Usage

During MX data collection, when a protein or nucleic acid crystal is exposed to increasing doses of radiation, localised radiation-induced chemical changes can occur within the crystalline macromolecules, even with doses as low as a few MGy. These *specific damage* manifestations can ultimately leading to false biological interpretations within structures during subsequent model building if not accounted for. 

Localised chemical changes within a macromolecule correlate with shifts in the electron density attributed to particular atoms within the crystal with increasing dose. In fact, radiation-induced decarboxylation of glutamate and aspartate residues can be detected by considering the deterioration in density local to the carboxylate group at different dose states, as is frequently performed using F(dn)-F(d1) Fourier difference maps between different accumulated dose states d1 and dn within a single crystal.

To quantitatively investigate dose-dependent density changes for any individual atoms within a structure, these **ETRACK** scripts can be utilised. This pipeline calculates metrics to quantify the damage susceptibility of each refined atom within a structure, including max density loss *Dloss*.

The scripts require the following to run:

- A series of input *.pdb* and *.mtz* files corresponding to a damage series collected from an individual crystal
- Python 2.7 
- The CCP4 suite downloaded (version non-specific, but tested on CCP4-6.4.0 and CCP4-6.5)
- The *seaborn* python plotting library
- A list of calculated doses for the series for full quantitative analysis (can be set to increasing integers if unknown)

First, for each dataset within the selected damage series, the *.pdb* and *.mtz* files are used to create compatible SFALL-output atom-tagged *.map* files and a Fourier difference map Fn-F1 *.map*. The *runProcessFiles.py script performs this task, using the following input file:

`dir /Users/charlie/DPhil/YEAR2/JAN/Weik2000_ETRACK`
`INITIALDATASET`
`name1 1qidinit`
`mtz1 /Users/charlie/DPhil/PDBredo_damageSeries/Weik2000/1qid/1qid.mtz`
`mtzlabels1 P_1qid`
`pdb1 /Users/charlie/DPhil/PDBredo_damageSeries/Weik2000/1qid/1qid.pdb`
`RfreeFlag1 FreeR_flag`
`LATERDATASET`
`name2 1qid`
`mtz2 /Users/charlie/DPhil/PDBredo_damageSeries/Weik2000/1qid/1qid.mtz`
`mtzlabels2 P_1qid`
`pdb2 /Users/charlie/DPhil/PDBredo_damageSeries/Weik2000/1qid/1qid.pdb`
`PHASEDATASET`
`name3 1qidinit`
`mtz3 /Users/charlie/DPhil/PDBredo_damageSeries/Weik2000/1qid/1qid.mtz`
`mtzlabels3 C_1qid`
`MAPINFO`
`sfall_VDWR 1`
`densMapType DIFF`
`FFTmapWeight True`
