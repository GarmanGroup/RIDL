# -*- coding: utf-8 -*-

# Script to run ccp4 programs pdbCUR, SFALL, FFT, MAPMASK
# pdbCUR --> to remove alternate conformations from .pdb file (pick most occupied)
# SFALL --> to create map file with voxel values given as atom numbers (of most contributing atom)
# FFT --> to create difference density map file over unit cell with same grid sampling as SFALL output atom map
# MAPMASK --> to crop/extend the FFT density map to same extent as atom map file (over asymmetric unit)

import os
import struct
import sys

# A class for pdbcur output
class pdbcur_output:
    def __init__(self,inputpdbfile="",outputpdbfile=""):

		self.inputpdbfile = inputpdbfile
		self.outputpdbfile = outputpdbfile

# A class for SFALL output
class SFALL_output:
    def __init__(self,inputpdbfile="",outputmapfile="",
    			 symmetrygroup=0,VDWR=0,gridsamp1="",
    			 gridsamp2="",gridsamp3="",fastaxis="",
    			 medaxis="",slowaxis=""):

		self.inputpdbfile = inputpdbfile
		self.outputmapfile = outputmapfile
		self.symmetrygroup = symmetrygroup
		self.VDWR = VDWR
		self.gridsamp1 = gridsamp1
		self.gridsamp2 = gridsamp2
		self.gridsamp3 = gridsamp3
		self.fastaxis = fastaxis
		self.medaxis = medaxis
		self.slowaxis = slowaxis

# A class for FFT output
class FFT_output:
    def __init__(self,inputmergedmtzfile="",outputmapfile="",
    			 gridsamp1="",gridsamp2="",gridsamp3="",
    			 fastaxis="",medaxis="",slowaxis=""):

		self.inputmergedmtzfile = inputmergedmtzfile
		self.outputmapfile = outputmapfile
		self.gridsamp1 = gridsamp1
		self.gridsamp2 = gridsamp2
		self.gridsamp3 = gridsamp3
		self.fastaxis = fastaxis
		self.medaxis = medaxis
		self.slowaxis = slowaxis

def pdbCUR_run(inputpdbfile,outputpdbfile):

	input1 = "/Applications/ccp4-6.4.0/bin/pdbcur "+\
			 "XYZIN %s " %(str(inputpdbfile))+\
			 "XYZOUT %s " %(str(outputpdbfile))

	# specify to remove hydrogen atoms and pick on the most probably conformation
	# (for two conformations with 0.5 occupancy, the first - A is chosen and occupancy
	# set to 1.00). Also remove all anisou info from file - since it is not needed for 
	# current analysis
	input2 = "delhydrogen\n"+\
			 "mostprob\n"+\
			 "noanisou\n"+\
			 "END"

	# run input2 in pdbcur
	textinput = open('inputfile.txt','w')
	textinput.write(input2)
	textinput.close()
	os.system(input1 +' < inputfile.txt')
	os.remove('inputfile.txt')


	print '\n---> pdbcur run complete...'
	print '--------------------------'
	print 'Summary:'
	print 'Input pdb file: %s' %(str(inputpdbfile))
	print 'Output pdb file: %s' %(str(outputpdbfile))

	# determine initial number of atoms in input pdb file
	pdbin = open(str(inputpdbfile),'r')
	pdbinlines = pdbin.readlines()
	counter = 0
	for line in pdbinlines:
		if 'ATOM' in line[0:5]:
			counter += 1
	pdbin.close()
	print 'number of atoms in input pdb file: %s' %(str(counter))

	# determine number of atoms in output pdb file
	pdbout = open(str(outputpdbfile),'r')
	pdboutlines = pdbout.readlines()
	counter = 0
	for line in pdboutlines:
		if 'ATOM' in line[0:5]:
			counter += 1
	pdbout.close()
	print 'number of atoms in output pdb file: %s' %(str(counter))

	print '--------------------------'
	# save details as object of class pdbcur_output
	pdbcur_runoutput = pdbcur_output(inputpdbfile,outputpdbfile)

	return pdbcur_runoutput

def SFALL_run(inputpdbfile,outputmapfile,symmetrygroup,sfall_GRID):

	VDWR = 1
	title = 'run of sfall'

	input1 = "/Applications/ccp4-6.4.0/bin/sfall "+\
			 "XYZIN %s " %(str(inputpdbfile))+\
			 "ATOMSF /Users/charlie/Desktop/atomsf_addNplus1.lib "+\
		 	 "MAPOUT %s " %(str(outputmapfile))+\
			 "SYMINFO /Applications/ccp4-6.4.0/lib/data/syminfo.lib "

	if len(sfall_GRID) == 0:
		input2 = "MODE ATMMAP ATMMOD\n"+\
			 	 "SYMMETRY %s\n" %(str(symmetrygroup))+\
			 	 "VDWR %s\n" %(str(VDWR))+\
			 	 "title %s\n" %(str(title))+\
   			 	 "END"
	else:
		input2 = "MODE ATMMAP ATMMOD\n"+\
			 	 "SYMMETRY %s\n" %(str(symmetrygroup))+\
			 	 "VDWR %s\n" %(str(VDWR))+\
			 	 "title %s\n" %(str(title))+\
			 	 "GRID %s %s %s\n" %(str(sfall_GRID[0]),str(sfall_GRID[1]),str(sfall_GRID[2]))+\
			 	 "END"

	# run input2 in sfall
	textinput = open('inputfile.txt','w')
	textinput.write(input2)
	textinput.close()
	os.system(input1 +' < inputfile.txt')
	os.remove('inputfile.txt')

	# open electron density .map file here 
	binarymapfile = open(str(outputmapfile),'rb')
	binarymapfile.seek(7*4,0)
	gridsamp1 = struct.unpack('=l',binarymapfile.read(4))[0]
	gridsamp2 = struct.unpack('=l',binarymapfile.read(4))[0]
	gridsamp3 = struct.unpack('=l',binarymapfile.read(4))[0]
	binarymapfile.seek(16*4,0)
	fastaxis = struct.unpack('=l',binarymapfile.read(4))[0]
	medaxis = struct.unpack('=l',binarymapfile.read(4))[0]
	slowaxis = struct.unpack('=l',binarymapfile.read(4))[0]
	binarymapfile.close()

	print '\n---> SFALL run complete...'
	print '--------------------------'
	print 'Summary:'
	print 'Input pdb file: %s' %(str(inputpdbfile))
	print 'Output map file: %s' %(str(outputmapfile))
	print 'Grid sampling 1: %s' %(str(gridsamp1))
	print 'Grid sampling 2: %s' %(str(gridsamp2))
	print 'Grid sampling 3: %s' %(str(gridsamp3))

	xyz = ['X','Y','Z']
	print 'Fast axis: %s' %(str(xyz[int(fastaxis)-1]))
	print 'Medium axis: %s' %(str(xyz[int(medaxis)-1]))
	print 'Slow axis: %s' %(str(xyz[int(slowaxis)-1]))
	print '--------------------------'


	SFALL_runoutput = SFALL_output(inputpdbfile,outputmapfile,
    			 				   symmetrygroup,VDWR,gridsamp1,
    			                   gridsamp2,gridsamp3,str(xyz[int(fastaxis)-1]),
    			                   str(xyz[int(medaxis)-1]),str(xyz[int(slowaxis)-1]))

	return SFALL_runoutput

def FFT_run(SFALL_runoutput,inputmergedmtzfile,outputmapfile,mtzlabels):

	Fobs_dam = mtzlabels[0]
	SIGobs_dam = mtzlabels[1]
	Fobs_init = mtzlabels[2]
	SIGobs_init = mtzlabels[3]
	PHIC_dam = mtzlabels[4]
	FOM_dam = mtzlabels[5]
	FOM_init = mtzlabels[6]

	input1 = "/Applications/ccp4-6.4.0/bin/fft "+\
			 "HKLIN %s " %(str(inputmergedmtzfile))+\
		 	 "MAPOUT %s " %(str(outputmapfile))+\
		 	 "SYMINFO /Applications/ccp4-6.4.0/lib/data/syminfo.lib "\

	# if FOM damaged exists, sub this line in below:
	# "LABIN F1=%s SIG1=%s F2=%s SIG2=%s " %(str(Fobs_dam),str(SIGobs_dam),str(Fobs_init),str(SIGobs_init)) +\
	#       "PHI=%s W=%s W2=%s\n" %(str(PHIC_dam),str(FOM_dam),str(FOM_init)) +\

	input2 = "AXIS %s %s %s\n" %(str(SFALL_runoutput.fastaxis),str(SFALL_runoutput.medaxis),str(SFALL_runoutput.slowaxis))+\
			 "title FTT DENSITY MAP RUN\n"+\
			 "grid %s %s %s\n" %(str(SFALL_runoutput.gridsamp1),str(SFALL_runoutput.gridsamp2),str(SFALL_runoutput.gridsamp3))+\
			 "xyzlim 0 1 0 1 0 1\n"+\
			  "LABIN F1=%s SIG1=%s F2=%s SIG2=%s " %(str(Fobs_dam),str(SIGobs_dam),str(Fobs_init),str(SIGobs_init)) +\
			       "PHI=%s W=%s\n" %(str(PHIC_dam),str(FOM_init)) +\
			 "SCALE F1 1.0 0.0 F2 1.0 0.0\n"+\
			 "END"

	# run input2 in sfall
	textinput = open('inputfile.txt','w')
	textinput.write(input2)
	textinput.close()
	os.system(input1 +' < inputfile.txt')
	os.remove('inputfile.txt')

	# open electron density .map file here 
	binarymapfile = open(str(outputmapfile),'rb')
	binarymapfile.seek(7*4,0)
	gridsamp1 = struct.unpack('=l',binarymapfile.read(4))[0]
	gridsamp2 = struct.unpack('=l',binarymapfile.read(4))[0]
	gridsamp3 = struct.unpack('=l',binarymapfile.read(4))[0]
	binarymapfile.seek(16*4,0)
	fastaxis = struct.unpack('=l',binarymapfile.read(4))[0]
	medaxis = struct.unpack('=l',binarymapfile.read(4))[0]
	slowaxis = struct.unpack('=l',binarymapfile.read(4))[0]
	binarymapfile.close()

	# save run details for map in object here
	xyz = ['X','Y','Z']
	FFT_runoutput = FFT_output(inputmergedmtzfile,outputmapfile,gridsamp1,gridsamp2,gridsamp3,
								str(xyz[int(fastaxis)-1]),str(xyz[int(medaxis)-1]),str(xyz[int(slowaxis)-1]))

	# summary of run to command line
	print '\n---> FFT run complete...'
	print '--------------------------'
	print 'Summary:'
	print 'Input merged mtz file: %s' %(str(inputmergedmtzfile))
	print 'Output map file: %s' %(str(outputmapfile))
	print 'Grid sampling 1: %s' %(str(gridsamp1))
	print 'Grid sampling 2: %s' %(str(gridsamp2))
	print 'Grid sampling 3: %s' %(str(gridsamp3))

	xyz = ['X','Y','Z']
	print 'Fast axis: %s' %(str(xyz[int(fastaxis)-1]))
	print 'Medium axis: %s' %(str(xyz[int(medaxis)-1]))
	print 'Slow axis: %s' %(str(xyz[int(slowaxis)-1]))
	print '--------------------------'

	return FFT_runoutput

def mapmask_run(inputmapfile,inputmapfile2,outputmapfile):

	input1 = "/Applications/ccp4-6.4.0/bin/mapmask "+\
			 "MAPIN %s " %(str(inputmapfile))+\
			 "MAPLIM %s " %(str(inputmapfile2))+\
		 	 "MAPOUT %s " %(str(outputmapfile))+\
		 	 "SYMINFO /Applications/ccp4-6.4.0/lib/data/syminfo.lib "

	input2 = "EXTEND\n"+\
			 "XYZLIM MATCH\n"+\
			 "END"

	# run input2 in sfall
	textinput = open('inputfile.txt','w')
	textinput.write(input2)
	textinput.close()
	os.system(input1 +' < inputfile.txt')
	os.remove('inputfile.txt')

	# summary of run to command line
	print '\n---> MAPMASK run complete...'
	print '--------------------------'

def mapmask_run_confine2AS(inputmapfile,outputmapfile):

	input1 = "/Applications/ccp4-6.4.0/bin/mapmask "+\
			 "MAPIN %s " %(str(inputmapfile))+\
		 	 "MAPOUT %s " %(str(outputmapfile))+\
		 	 "SYMINFO /Applications/ccp4-6.4.0/lib/data/syminfo.lib "

	input2 = "EXTEND\n"+\
			 "XYZLIM ASU\n"+\
			 "END"

	# run input2 in sfall
	textinput = open('inputfile.txt','w')
	textinput.write(input2)
	textinput.close()
	os.system(input1 +' < inputfile.txt')
	os.remove('inputfile.txt')

	# summary of run to command line
	print '\n---> MAPMASK run complete...'
	print '--------------------------'

def atmmap_densmap_check(SFALL_runoutput,FFT_runoutput):
	# this function determines whether the atom map and density map calculated using SFALL and FFT
	# are compatible - meaning the grid dimensions/filtering are the same and the ordering of the
	# fast, medium, and slow axes are identical.

	print 'Checking that atom map (SFALL) and density map (FFT) are compatible...'

	if (SFALL_runoutput.gridsamp1 != FFT_runoutput.gridsamp1 or
		SFALL_runoutput.gridsamp2 != FFT_runoutput.gridsamp2 or
		SFALL_runoutput.gridsamp3 != FFT_runoutput.gridsamp3):
		print 'Incompatible grid sampling found...'
		sys.exit()
	if (SFALL_runoutput.fastaxis != FFT_runoutput.fastaxis or
		SFALL_runoutput.medaxis != FFT_runoutput.medaxis or
		SFALL_runoutput.slowaxis != FFT_runoutput.slowaxis):
		print 'Incompatible fast,med,slow axes ordering found...'
		sys.exit()
	if (SFALL_runoutput.gridsamp1 != FFT_runoutput.gridsamp1 or
		SFALL_runoutput.gridsamp2 != FFT_runoutput.gridsamp2 or
		SFALL_runoutput.gridsamp3 != FFT_runoutput.gridsamp3):
		print 'Incompatible grid sampling found...'
		sys.exit()

	# find number of bytes per file used in map (removing header which
	# may vary in size...)
	binarymapfile = open(str(SFALL_runoutput.outputmapfile),'rb')
	nx = struct.unpack('=l',binarymapfile.read(4))[0]
	ny = struct.unpack('=l',binarymapfile.read(4))[0]
	nz = struct.unpack('=l',binarymapfile.read(4))[0]
	densitystart_atm = 4*(nx*ny*nz)

	binarymapfile = open(str(FFT_runoutput.outputmapfile),'rb')
	nx = struct.unpack('=l',binarymapfile.read(4))[0]
	ny = struct.unpack('=l',binarymapfile.read(4))[0]
	nz = struct.unpack('=l',binarymapfile.read(4))[0]
	densitystart_dens = 4*(nx*ny*nz)

	if densitystart_atm != densitystart_dens:
		print 'Incompatible map file sizes'
		sys.exit()

	print '---> success!'

def renumber_pdbfile(pdbfilein,pdbfileout):
	# reorder atoms in pdb file since some may be missing now after
	# pdbcur has been used

	pdbin = open(pdbfilein,'r')
	pdbout = open(pdbfileout,'w')

	counter = 0
	for line in pdbin.readlines():
		if ('ATOM' not in line[0:4] and 'HETATM' not in line[0:6]):
			pdbout.write(line)
		else:
			counter += 1
			pdbout.write(line[0:6])
			new_atomnum = " "*(5-len(str(counter))) + str(counter) #has length 5
			pdbout.write(new_atomnum)
			pdbout.write(line[11:80]+'\n')

	return pdbfileout

def get_atomanddensmaps():
	# read in inputfile
	inputfile = open('inputfile_ccp4jobs.txt','r')
	inputlines = inputfile.readlines()

	sfall_GRID = []
	for line in inputlines:
		if 'pdbIN' == line.split()[0]:
			pdbname = line.split()[1]
		if 'runname' == line.split()[0]:
			runname = line.split()[1]
		if 'dataset' == line.split()[0]:
			N = line.split()[1]
		if 'sfall_symm' == line.split()[0]:
			sfall_symmetrygroup = line.split()[1]
		if 'sfall_GRID' == line.split()[0]:
			sfall_GRIDnx = line.split()[1]
			sfall_GRIDny = line.split()[2]
			sfall_GRIDnz = line.split()[3]
			sfall_GRID = [sfall_GRIDnx,sfall_GRIDny,sfall_GRIDnz]
		if 'mtzIN' == line.split()[0]:
			fft_inputmergedmtzfile = line.split()[1]
		if 'fft_Fobs_dam' == line.split()[0]:
			Fobs_dam = line.split()[1]
		if 'fft_SIGobs_dam' == line.split()[0]:
			SIGobs_dam = line.split()[1]
		if 'fft_Fobs_init' == line.split()[0]:
			Fobs_init = line.split()[1]
		if 'fft_SIGobs_init' == line.split()[0]:
			SIGobs_init = line.split()[1]
		if 'fft_PHIC_dam' == line.split()[0]:
			PHIC_dam = line.split()[1]
		if 'fft_FOM_dam' == line.split()[0]:
			FOM_dam = line.split()[1]
		if 'fft_FOM_init' == line.split()[0]:
			FOM_init = line.split()[1]
		if 'foldername' == line.split()[0]:
			where = line.split()[1]
		if 'END' == line.split()[0]:
			break

	pdbcur_inputpdbfile = pdbname
	pdbcur_outputpdbfile = where+'pdbcur_'+runname+N+'.pdb'
	pdbreorder_outpdbfile = where+'pdbreorder_'+runname+N+'.pdb'
	sfall_outputmapfile = where+runname+N+'_atoms.map'
	FFT_outputmapfile = where+runname+N+'_density.map'

	# sfall_outputmapfileCROPPED = where+runname+N+'_atoms_CROPPED.map'

	# run pdbcur job
	pdbcur_runoutput = pdbCUR_run(pdbcur_inputpdbfile,pdbcur_outputpdbfile)

	# run my own reordering script
	inputpdbfile = pdbcur_runoutput.outputpdbfile
	reorderedpdb = renumber_pdbfile(inputpdbfile,pdbreorder_outpdbfile)

	# run sfall job
	SFALL_runoutput = SFALL_run(reorderedpdb,sfall_outputmapfile,sfall_symmetrygroup,sfall_GRID)

	# run fft job
	mtzlabels = [Fobs_dam,SIGobs_dam,Fobs_init,SIGobs_init,PHIC_dam,FOM_dam,FOM_init]
	FFT_runoutput = FFT_run(SFALL_runoutput,fft_inputmergedmtzfile,FFT_outputmapfile,mtzlabels)

	# crop to asymmetric unit:
	mapmask_run_confine2AS(SFALL_runoutput.outputmapfile,SFALL_runoutput.outputmapfile)
	# mapmask_run_confine2AS(SFALL_runoutput.outputmapfile,sfall_outputmapfileCROPPED)

	mapmask_run_confine2AS(FFT_runoutput.outputmapfile,FFT_runoutput.outputmapfile)
	
	# run MAPMASK job to crop FFT density map to same grid 
	# sampling dimensions as SFALL atom map
	mapmask_run(FFT_runoutput.outputmapfile,SFALL_runoutput.outputmapfile,FFT_runoutput.outputmapfile)
	# mapmask_run(FFT_runoutput.outputmapfile,sfall_outputmapfileCROPPED,FFT_runoutput.outputmapfile)

	# check consistency between atom and density maps generated
	atmmap_densmap_check(SFALL_runoutput,FFT_runoutput)











	

