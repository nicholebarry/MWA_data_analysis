import numpy as np
import os
import glob
from pyuvdata import UVData
import argparse


#********************************
def main():


	# Parse the command line inputs. 
	parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=
		"Apply van vleck corrections to gpubox files iteratively and then save a singular uvfits file \n" +
		"Example call: \n"
		"python van_vleck_corrector.py --obsfile_name=zenith.txt --data_path='/fred/oz048/MWA/data/2013/van_vleck_corrected/' ")
	parser.add_argument("-o", "--obs_id", required=True,
		help="Name of the observation")
	parser.add_argument("-d", "--data_path", required=True,
		help="Path to the data (gpubox files, metadata, and mwaf files)")
	parser.add_argument("-i", "--integrate_only", required=False, default=False,
		help="Combine coarse channel uvfits into one uvfits")
	args = parser.parse_args()

	obs_id = args.obs_id
	# Get obsids to download
	#obsfile = open(args.obsfile_name, "r")
	#obsids = [line.split( ) for line in obsfile.readlines()]
	#obsids = [obs[0] for obs in obsids]
	#obsfile.close()

	# List of gpu files to cycle through
	#coarse_list = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19',
	#	'20','21','22','23','24']
	coarse_list = [['01','02','03','04','05','06','07','08','09','10','11','12'],['13','14','15','16','17','18','19',
		'20','21','22','23','24']]

	iter_flag=0
	if not args.integrate_only:
		for coarse_id in coarse_list:
			UV = van_vleck_corrections(obs_id, coarse_id, args.data_path)
			UV = integration(UV)
			if iter_flag == 0:
				print("Writing " + args.data_path + obs_id + "_firsthalf.uvfits")
				UV.write_uvfits(args.data_path + obs_id + '_firsthalf.uvfits', spoof_nonessential=True)
			else:
				print("Writing " + args.data_path + obs_id + "_secondhalf.uvfits")
				UV.write_uvfits(args.data_path + obs_id + '_secondhalf.uvfits', spoof_nonessential=True)
			iter_flag = iter_flag + 1
			#print("Writing " + args.data_path + obs_id + ".uvfits")
			#UV.write_uvfits(args.data_path + obs_id + '_' + coarse_id + '.uvfits', spoof_nonessential=True)
			#if coarse_id == '01':
			#	UV_total = UV
			#else:
			#	UV_total = UV_total + UV
	else:
		#for coarse_id in coarse_list:
		#	UV = UVData()
		#	UV.read_uvfits(args.data_path + obs_id + '_' + coarse_id + '.uvfits')
		#	if coarse_id == '01':
		#		UV_total = UV
		#		#UV_add = UV_total
		#	else:
		#		UV_total = UV_total + UV
		#		if coarse_id == '12':
		#			UV_total = UV_total + UV
		#			UV_total = time_integration(UV_total)
		#			UV_total.write_uvfits(args.data_path + obs_id + '_firsthalf.uvfits', spoof_nonessential=True)
		#		else:
		#			if coarse_id == '13':
		#				UV_total = UV
		#			else:
		#				UV_total = UV_total + UV
		UV = UVData()
		UV_total = UVData()
		#UV_total = time_integration(UV_total)
		#UV_total.write_uvfits(args.data_path + obs_id + '_secondhalf.uvfits', spoof_nonessential=True)
		UV.read_uvfits(args.data_path + obs_id + '_firsthalf.uvfits')
		UV_total.read_uvfits(args.data_path + obs_id + '_secondhalf.uvfits')
		UV_total = UV_total + UV
		UV_total = time_integration(UV_total)
		UV_total.write_uvfits(args.data_path + obs_id + '.uvfits', spoof_nonessential=True)

	#print("Writing " + args.data_path + obs_id + ".uvfits")
	#UV_total.write_uvfits(data_obs + obs_id + '.uvfits', spoof_nonessential=True)

#********************************
def van_vleck_corrections(obs_id, coarse_id, data_path):

	UV = UVData()

	#gpubox = glob.glob(data_path + obs_id + '*_gpubox'+coarse_id+'*')
	gpubox = [glob.glob(data_path + obs_id + '*_gpubox'+coarse_channel+'*') for coarse_channel in coarse_id]
	gpubox = [ item for elem in gpubox for item in elem]
	mwaf = [glob.glob(data_path + obs_id + '_' + coarse_channel+'.mwaf') for coarse_channel in coarse_id]
	mwaf = [ item for elem in mwaf for item in elem]
	filelist = [data_path + obs_id + '.metafits', data_path + obs_id + '_metafits_ppds.fits', *mwaf, *gpubox]
	#filelist = [data_path + obs_id + '.metafits', data_path + obs_id + '_metafits_ppds.fits', data_path + obs_id + '_' + coarse_id + '.mwaf', *gpubox]
	#gpubox = glob.glob(data_path + obs_id + '*_gpubox*')
	#mwaf = glob.glob(data_path + obs_id + '_*.mwaf')
	#filelist = [data_path + obs_id + '.metafits', data_path + obs_id + '_metafits_ppds.fits', *mwaf, *gpubox]

	UV.read_mwa_corr_fits(filelist,use_aoflagger_flags=True,correct_cable_len=True,correct_van_vleck=True,phase_to_pointing_center=True)

	return UV

#********************************
def integration(UV):

	# Should be a nice, round number for the MWA
	#time_int = np.mean(UV.integration_time)
	#if time_int != 2.0:
	#	average_factor = int(2.0 / time_int)
	#	UV.downsample_in_time(n_times_to_avg=average_factor,keep_ragged=False)

	freq_int = UV.channel_width
	if freq_int != 80000.0:
		average_factor = int(80000.0 / freq_int)
		UV.frequency_average(average_factor)

	return UV

#********************************

#********************************
def time_integration(UV):

        # Should be a nice, round number for the MWA
        time_int = np.mean(UV.integration_time)
        if time_int != 2.0:
               average_factor = int(2.0 / time_int)
               UV.downsample_in_time(n_times_to_avg=average_factor,keep_ragged=False)

        #freq_int = UV.channel_width
        #if freq_int != 80000.0:
        #        average_factor = int(80000.0 / freq_int)
        #        UV.frequency_average(average_factor)

        return UV

#********************************

#********************************

if __name__ == '__main__':
    main()
