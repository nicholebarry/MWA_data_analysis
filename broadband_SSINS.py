#!/usr/bin/python

from SSINS import SS
from SSINS import MF
from SSINS import INS
import os
import numpy as np
from pyuvdata import UVData, UVFlag
from pyuvdata import utils as uvutils



#********************************
#Script that moves updates a few necessary stats

#Script written by Nichole Barry, May 2017.

def main():

	#test!!
	#obsids = ['1064768840']

	#directory of original uvfits
	orig_dir = '/fred/oz048/MWA/data/2013/van_vleck_corrected/'
	prefix = 'SSINS/data/'

	#Threshold for flagging
	sig_thresh = 5

	#Get obsids to download
	obsfile_name='/home/nbarry/MWA/pipeline_scripts/bash_scripts/ozstar/obs_list/1_plustwo.txt'
	obsfile = open(obsfile_name, "r")
	obsids = [line.split( ) for line in obsfile.readlines()]
	obsids = [obs[0] for obs in obsids]
	obsfile.close()
	nonredundant_obsids = list(set(obsids))
	if len(obsids) != len(nonredundant_obsids):
		print("WARNING: Obs list contains redundant entries.")
		obsids = nonredundant_obsids



	for obsid in obsids:

		obs_index = obsids.index(obsid)

		ss = SS()
		ss.read(orig_dir+obsid +'.uvfits', diff=True)
		ss.apply_flags(flag_choice='original')
		ins = INS(ss)
		mf = MF(ins.freq_array, sig_thresh, streak=True, narrow=False, shape_dict={})
		mf.apply_match_test(ins)

		uvd = UVData()
		uvd.read(orig_dir+obsid+'.uvfits')
		uvf = UVFlag(uvd, waterfall=True, mode='flag')
		flags = ins.mask_to_flags()
		flags = ins.flag_uvf(uvf)
		uvutils.apply_uvflag(uvd, flags)
		uvd.write_uvfits(orig_dir + prefix + obsid + '.uvfits')




if __name__ == '__main__':
	main()
