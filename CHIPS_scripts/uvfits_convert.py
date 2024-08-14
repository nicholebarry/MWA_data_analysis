from pyuvdata import UVData
import numpy
import os
import copy
import sys
import subprocess
from optparse import OptionParser



#********************************
def main():

	#Parse the command line inputs. 
	parser=OptionParser()
	parser.add_option("-c", "--chips_format", dest="chips_format", \
		help="Convert to CHIPS format: uvfits per coarse band")
	parser.add_option("-e", "--fhd_epp_format", dest="fhd_epp_format", \
		help="Convert to FHD/eppsilon format: one uvfits")
	parser.add_option("-s", "--specific_dir", dest="specific_dir", \
		help="Specific fhd directory")
	parser.add_option("-o", "--obsfile_name", dest="obsfile_name", \
		help="Name of the file containing the observation list")
	parser.add_option("-d", "--make_dirty", dest="make_dirty", \
		help="Make dirty uvfits.",default='')
	parser.add_option("-m", "--make_model", dest="make_model", \
		help="Make model uvfits",default='')
	parser.add_option("-r", "--make_residual", dest="make_residual", \
		help="Make residul uvfits. Waiting on a pyuvdata update for this option.",default='')
	(options, args)=parser.parse_args()
	chips_format=options.chips_format
	fhd_epp_format=options.fhd_epp_format
	specific_dir=options.specific_dir
	obsfile_name=options.obsfile_name
	make_dirty=options.make_dirty
	make_model=options.make_model
	make_residual=options.make_residual

	try: 
		chips_format
	except NameError:
		chips_format = 0
	try:
		fhd_epp_format
	except NameError:
		fhd_epp_format = 1
		
        common_dir = '/fred/oz048/MWA/CODE/FHD/'

	try:
		specific_dir
	except NameError:
		specific_dir = 'fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z'
	try:
		obsfile_name
	except NameError:
		obsfile_name = '/home/nbarry/MWA/FHD/obs_list/zenith.txt'
	try:
		make_dirty
	except NameError:
		make_dirty = 1
	try:
		make_model
	except NameError:
		make_model = 1
	
	if str(chips_format) == '1':
		out = fhd_vis_to_chips_format(common_dir,specific_dir,obsfile_name,make_dirty,make_model,make_residual)
	
	if str(fhd_epp_format) == '1':
		out = rts_vis_to_fhd_format(common_dir,specific_dir,obsfile_name)

#********************************

#********************************
#Turn FHD sav directories into CHIPS compatible uvfits per coarse band
#Example call:
#python uvfits_convert.py -c 1 -e 0 -s fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z -o ~/MWA/FHD/obs_list/Aug23_longrunstyle.txt -d 1 -m 1
def fhd_vis_to_chips_format(common_dir,specific_dir,obsfile_name,make_dirty,make_model,make_residual):

	obsinfo_dir = '/fred/oz048/MWA/CODE/CHIPS/chips/obsinfo_files/'
	script_flag = 0

	#Get obsids to download
	obsfile = open(obsfile_name, "r")
	obsids = [line.split( ) for line in obsfile.readlines()]
	obsids = [obs[0] for obs in obsids]
	obsfile.close()

	if not os.path.exists(common_dir + specific_dir + '/CHIPS_input'):
   		os.mkdir(common_dir + specific_dir + '/CHIPS_input')

	for obs_id in obsids:
    		#Make the dir in the FHD output folder
    		if not os.path.exists(common_dir + specific_dir + '/CHIPS_input/' + obs_id):
	        	os.mkdir(common_dir + specific_dir + '/CHIPS_input/' + obs_id)

    		#Set CHIPS shell script flag to make obs info if not present already
    		if (not os.path.exists(obsinfo_dir + 'obsid' + obs_id + '.txt')) or \
        		(not os.path.exists(obsinfo_dir + 'obsidra' + obs_id + '.txt')) or \
        		(not os.path.exists(obsinfo_dir + 'obsidtsys' + obs_id + '.txt')):
            		script_flag = 1

                if (os.path.exists(common_dir + specific_dir + '/CHIPS_input/' + obs_id + '/uv_dirty_' +str(1).zfill(2) + '.uvfits')) and \
			(os.path.exists(common_dir + specific_dir + '/CHIPS_input/' + obs_id + '/uv_dirty_' +str(24).zfill(2) + '.uvfits')):
                        print('Uv_dirty already exists for ' + obs_id)
                        continue

    		# Construct the list of files
    		fhd_prefix_vis = common_dir + specific_dir + '/vis_data/' + obs_id + '_'
    		fhd_prefix_set = common_dir + specific_dir + '/metadata/' + obs_id + '_'

		if (not os.path.exists(fhd_prefix_vis + 'vis_XX.sav')) or \
			(not os.path.exists(fhd_prefix_vis + 'vis_YY.sav')):
			print(obs_id + ' is missing')
			continue

    		fhd_files = [fhd_prefix_vis + 'flags.sav', fhd_prefix_vis + 'vis_XX.sav', fhd_prefix_set + 'params.sav',
                                       fhd_prefix_vis + 'vis_YY.sav', fhd_prefix_vis + 'vis_model_XX.sav',
                                       fhd_prefix_vis + 'vis_model_YY.sav', fhd_prefix_set + 'settings.txt']

    		# Read in all fhd files into UV class
    		if make_dirty or make_residual:
        		UV_dirty = UVData()
        		UV_dirty.read_fhd(fhd_files)
    		if make_model or make_residual:
        		UV_model = UVData()
        		UV_model.read_fhd(fhd_files,use_model=1)

		# Find the upper and lower bound for each coarse band given the input is FHD-typical (384 freq channels)
    		mod_list = numpy.array([(x % 16) for x in range(384)])
    		top_coarse_ind = numpy.where(mod_list == 15) 
    		bottom_coarse_ind = numpy.where(mod_list == 0) 

		# Split up UV class into coarse bands and save uvfits
    		for coarse_i in range(len(bottom_coarse_ind[0])):
        		if make_dirty:
            			UV1 = copy.deepcopy(UV_dirty)
            			UV1.select(frequencies=UV_dirty.freq_array[0][bottom_coarse_ind[0][coarse_i]:top_coarse_ind[0][coarse_i]+1])
            			UV1.write_uvfits(common_dir + specific_dir + '/CHIPS_input/' + obs_id +
                		'/uv_dirty_'+str(coarse_i+1).zfill(2)+'.uvfits', spoof_nonessential=True)
        		if make_model:
            			UV1 = copy.deepcopy(UV_model)
            			UV1.select(frequencies=UV_model.freq_array[0][bottom_coarse_ind[0][coarse_i]:top_coarse_ind[0][coarse_i]+1])
            			UV1.write_uvfits(common_dir + specific_dir + '/CHIPS_input/' + obs_id +
                		'/uv_model_'+str(coarse_i+1).zfill(2)+'.uvfits', spoof_nonessential=True)
        		if make_residual:
            			UV1 = copy.deepcopy(UV_dirty - UV_model)
            			UV1.select(frequencies=UV_dirty.freq_array[0][bottom_coarse_ind[0][coarse_i]:top_coarse_ind[0][coarse_i]+1])
            			UV1.write_uvfits(common_dir + specific_dir + '/CHIPS_input/' + obs_id +
                		'/uv_res_'+str(coarse_i+1).zfill(2)+'.uvfits', spoof_nonessential=True)

	# Run obs info generator for CHIPS input if it doesn't already exist
	if script_flag:
    		subprocess.call(['/fred/oz048/MWA/CODE/CHIPS/chips/scripts/make_obs_info_ozstar.sh',obsfile_name])

#********************************

#********************************
#Turn RTS uvfits per coarse band outputs into FHD compatible full uvfits
#Example call:
#python uvfits_convert.py -c 0 -e 1 -s RTS_compare -o ~/MWA/FHD/obs_list/Aug23_longrunstyle.txt
def rts_vis_to_fhd_format(common_dir,specific_dir,obsfile_name):

        #Get obsids to download
        obsfile = open(obsfile_name, "r")
        obsids = [line.split( ) for line in obsfile.readlines()]
        obsids = [obs[0] for obs in obsids]
        obsfile.close()

	# Remove uvdump files afterwards
	cleanup = 1

        for obs_id in obsids:
        	if not os.path.exists(common_dir + specific_dir + '/' + obs_id):
			print('Download uvdump files for ' + str(obs_id) + ' first')
			continue
		else: 
			dir_obs = common_dir + specific_dir + '/' + obs_id + '/'
	
		if os.path.exists(dir_obs + obs_id + '.uvfits'):
			print('Uvfits already exists for ' + obs_id)
			continue

		if not os.path.exists(dir_obs + 'uvdump_01.uvfits'):
			print('Download uvdump files for ' + str(obs_id) + ' first using uvdump naming convention.')
			continue		
	
		UV_total = UVData()

		for coarse_i in range(24):
			UV = UVData()
			UV.read(dir_obs + 'uvdump_' + str(coarse_i+1).zfill(2)+'.uvfits')
			if coarse_i == 0:
				UV_total = UV
			else:
				UV_total = UV_total + UV

		print("Writing " + dir_obs + obs_id + ".uvfits")
		UV_total.write_uvfits(dir_obs + obs_id + '.uvfits', spoof_nonessential=True)
		
		if cleanup:
			files = os.listdir(dir_obs)
			for file_i in files:
				if file_i.startswith("uvdump"):
					os.remove(dir_obs + file_i)

#********************************

#********************************

if __name__ == '__main__':
	main()
