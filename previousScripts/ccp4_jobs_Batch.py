from ccp4_jobs import get_atomanddensmaps as run_ccp4Pipeline

def run():
	# state number of datasets within damage data series here
	numDatasets = 16

	for i in range(2,numDatasets+1):
		inputfile = open('inputfile_ccp4jobs.txt','w')
		inputfile.write('pdbIN /Users/charlie/DPhil/YEAR1/MAY/TRAPwork/TRAPd_Reprocess/phenixrefine_isotropicBfactors_run/trapd'+str(i)+'/TRAPd_May2015_reprocess_refine_'+str(i+11)+'.pdb\n')
		inputfile.write('runname TRAPRNAdamage\n')
		inputfile.write('foldername ../TRAPwork/eTrack_TRAPd/\n')
		inputfile.write('dataset '+str(i)+'\n')
		inputfile.write('sfall_symm 5\n')
		inputfile.write('mtzIN /Users/charlie/DPhil/YEAR1/MAY/TRAPwork/TRAPd_Reprocess/scaleit_run/trapd'+str(i)+'/trapd'+str(i)+'_cad_scaleit1.mtz\n')
		inputfile.write('fft_Fobs_dam F-obs-filtered_'+str(i)+'\n')
		inputfile.write('fft_SIGobs_dam SIGF-obs-filtered_'+str(i)+'\n')
		inputfile.write('fft_Fobs_init F_trapd_dose1\n')
		inputfile.write('fft_SIGobs_init SIGF_trapd_dose1\n')
		inputfile.write('fft_PHIC_dam PHIF-model_'+str(i)+'\n')
		inputfile.write('fft_FOM_init FOM\n')
		inputfile.write('fft_FOM_dam FOM_'+str(i)+'\n')
		inputfile.write('END')
		inputfile.close()

		# now run the ccp4 pipeline
		run_ccp4Pipeline()
