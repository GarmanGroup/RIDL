# ETRACK

A collection of scripts to calculate per-atom density change metrics within a specific damage MX investigation.
Suitable for any MX experiment in which datasets are collected on the **same crystal** over **multiple doses**.
**NOTE: These scripts are currently under development and updated regularly..**


## Brief background

During MX data collection, when a protein or nucleic acid crystal is exposed to increasing doses of radiation, localised radiation-induced chemical changes can occur within the crystalline macromolecules, even with doses as low as a few MGy. These *specific damage* manifestations can ultimately leading to false biological interpretations within structures during subsequent model building if not accounted for. 

Localised chemical changes within a macromolecule correlate with shifts in the electron density attributed to particular atoms within the crystal with increasing dose. In fact, radiation-induced decarboxylation of glutamate and aspartate residues can be detected by considering the deterioration in density local to the carboxylate group at different dose states, as is frequently performed using F(dn)-F(d1) Fourier difference maps between different accumulated dose states d1 and dn within a single crystal.

To quantitatively investigate dose-dependent density changes for any individual atoms within a structure, these **ETRACK** scripts can be utilised. This pipeline calculates metrics to quantify the damage susceptibility of each refined atom within a structure, including max density loss *Dloss*.

## Usage

The scripts require the following to run:

- A series of input *.pdb* and *.mtz* files corresponding to a damage series collected from an individual crystal
- Python 2.7 (main testing performed on 2.7.10)
- The CCP4 suite downloaded (version non-specific, but tested on 6.4.0 and 6.5)
- The *seaborn* python plotting library
- A list of calculated doses for the series is ideal but not essential

To show how the scripts can be run, use the *getPDBseries.py* script to retrieve a pdb series from pdb_redo. Here the first 3 entries within the *Weik M, et al. (2000) PNAS 97(2):623–8* have been selected as a test damage series. In python:

```python
import os
os.system('mkdir testOutput') # make a directory for this test run
from getPDBSeries import getSeries

for i in ['1qid','1qie','1qif']:
	p = getSeries(i,'initial','./testOutput/')
```
The *initial* parameter indicates that files have been downloaded from pdb_redo **before** re-refinement.

For each dataset within this selected damage series, the *.pdb* and *.mtz* files are used to generated compatible SFALL-output atom-tagged *.map* files and a Fourier difference map Fn-F1 *.map* calculated between the initial dataset (*1qid* here) and higher-dose dataset (*1qie* or *1qif* here). The *runProcessFiles.py script performs this task, but it requires an input file. Luckily for this test damage series, the input file can be written for each dataset:

```python
from runProcessFiles import process
p = process()
p.writeTestInputFile(1) # generate the 1st set of maps: a 1qid - 1qid difference map
p.writeTestInputFile(2) # generate the 2nd set of maps: a 1qie - 1qid difference map
p.writeTestInputFile(3) # generate the first set of maps: a 1qif - 1qid difference map
```

For the general case this input file must be written manually, specifying correct mtz column label information as required. However there are two functions `p.writeTemplateInputFile()` and `p.writePDBredoInputFile('1qid','1qie','./testOutput',./testOutput/ETRACK)` that can speed up this process.

After writing an input file, the script assigns this input file name as the *current* input file. After `p.writeTestInputFile(3)` is run above, the current input file is then *testInput3.txt*. The contents of this input file can be printed using `p.printInputFile()`. The process the 3 test input files here in the correct order, use:

```python
p.setInputFile('testInput1.txt') 
p.run()
p.setInputFile('testInput2.txt')
p.run()
p.setInputFile('testInput3.txt')
p.run()
```

Within the directory `./testOutput/ETRACK/` the results of this run should be detailed. For each of the 3 runs, an atom-tagged map file e.g. `1qie_atoms.map` and corresponding Fourier difference map file `1qie_density.map` have been generated. A compressed directory `1qie-1qidinit_additionalFiles.tar` contains CCP4 log files and intermediate run files, and additionally 2 run logs `1qie-1qidinit_runLog_1.txt` and `1qie-1qidinit_runLog_2.txt` are generated to indicate the success or failure of the run.

**The rest of this example is to come!...**





