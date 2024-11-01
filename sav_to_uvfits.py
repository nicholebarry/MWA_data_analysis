from pyuvdata import UVData
UV = UVData()
fhd_prefix = '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/'

# Construct the list of files
filepath='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23_longrunstyle.txt'
with open(filepath) as fp: 
	for obsid in fp:
		obsid = obsid.rstrip('\n')
		print 'Working on coverting ' + obsid
		fhd_files = [fhd_prefix + f for f in ['vis_data/'+obsid+'_flags.sav', 'vis_res/'+obsid+'_vis_XX.sav', 'metadata/'+obsid+'_params.sav',
                                      'vis_res/'+obsid+'_vis_YY.sav', 'metadata/'+obsid+'_settings.txt']]
		print fhd_files
		UV.read_fhd(fhd_files)
		UV.write_uvfits(fhd_prefix+'uvfits/'+obsid+'.uvfits', spoof_nonessential=True)
