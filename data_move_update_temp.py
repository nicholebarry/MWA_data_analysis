#!/usr/bin/python

import glob,os
from astropy.time import Time
from astropy.io import fits
import psycopg2
import sys
import socket
from optparse import OptionParser
import subprocess
import datetime
import time
import zipfile
import numpy as np
import shutil
import ntpath

#********************************
#Script that moves old uvfits, metafits, and qs files from locations on the MIT cluster
#to places that people are not utilizing while updating a few necessary stats along 
#the way.

#Script written by Nichole Barry, May 2017.

def main():

	#Assume cotter version is the same for files in eor-14 batch folder, and version, subversion
	cotter_version = 'version 2.4 (2014-08-07)'
	version=4
	subversion=1
	old_path = '/nfs/mwa-14/r1/EoRuvfits/jd2456528v4_1/'
	old_path2 = '/nfs/eor-14/r1/EoRuvfits/jd2456528v4_1/'
	#folder_name = 'jd2456528v4_1'

	all_nodes = ["eor-03", "eor-04", "eor-05","eor-06", "eor-07", "eor-08", "eor-10", "eor-11"]#, "eor-12"]#, "eor-13", "eor-14"]
	all_nodes = ["/nfs/" + nodename + "/r1/" for nodename in all_nodes]


	try:
		folder_name
	except NameError:
		#Get obsids to download
		obsfile_name='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/obs_list/datafile_update.txt'
		obsfile = open(obsfile_name, "r")
		obsids = [line.split( ) for line in obsfile.readlines()]
		obsids = [obs[0] for obs in obsids]
		obsfile.close()
		nonredundant_obsids = list(set(obsids))
		if len(obsids) != len(nonredundant_obsids):
			print "WARNING: Obs list contains redundant entries."
			obsids = nonredundant_obsids
	else:	
		#Get obsids to download directly from the folder
		old_path = old_path + folder_name + '/'
		old_path2 = old_path2 + folder_name + '/'
		temp_path = '/nfs/mwa-11/r1/EoRuvfits/jd2456528v4_1'
		obsids_path = glob.glob(temp_path + '/*')
		obsids = [ os.path.splitext(os.path.basename(obs))[0] for obs in obsids_path]

	#Find the obsids' save directory names
	t = Time([int(obsid) for obsid in obsids], format="gps", scale="utc")
	jds = t.jd
	jds = [int(jd) for jd in jds]
	save_directories = ["EoRuvfits/jd" + str(jd) + "v"+ str(version) + "_" + str(subversion) + "/" for jd in jds]

	#Connect to the database on eor-00 using the mwa username.
	try:
		conn = psycopg2.connect(database='mwa_qc',user='mwa',password='BowTie',host='eor-00.mit.edu')
	except:
		print 'Could not connect to mwa database.'
		sys.exit(1)

	#Module only works on the MIT cluster.
	if not 'mit.edu' in socket.gethostname():
		print 'Sorry, this script is currently only supported on eor-xx.mit.edu machines.'
		sys.exit(1)  
		
	cur = conn.cursor()

	for obsid in obsids:

		obs_index = obsids.index(obsid)

		save_path = all_nodes[7] + save_directories[obs_index] + obsid #put approximately 500 obsids on specified node if following test turns up false
		for loc_node in all_nodes:
			loc_path = loc_node + save_directories[obs_index]
			if os.path.isdir(loc_path): #checks to see if the directory exists
				save_path = loc_path + obsid
		#test!!
		#print save_path

		#Open up the metafits file that was made with the uvfits file (assumes they are in the same location)
		metafits_file = save_path + '/' + obsid + ".metafits"
		if not os.path.exists(metafits_file):
			print metafits_file + ' does not exist, skipping'
			continue
		hdu_list_metafits = fits.open(metafits_file)
		header_metafits = hdu_list_metafits[0].header

		#Get the frequency center and bandwidth (assumes contiguous frequency channels)
		freq_cent = header_metafits['FREQCENT']
		bandwidth = header_metafits['BANDWDTH']
		top_freq_mhz = "{0:.2f}".format(freq_cent + bandwidth/2)
		bottom_freq_mhz = "{0:.2f}".format(freq_cent - bandwidth/2)

		#Close the metafits file
		hdu_list_metafits.close()

    		if not os.path.exists(save_path):
        		os.makedirs(save_path)
		new_path = save_path + "/" + obsid + ".uvfits"
		#shutil.move(old_path + obsid + ".uvfits", new_path)
		#shutil.move(old_path + obsid + ".metafits", save_path + "/" + obsid + ".metafits")
		#shutil.move(old_path + obsid + ".qs", save_path + "/" + obsid + ".qs")

		#print new_path + ', ' + cotter_version + ', ' + bottom_freq_mhz + ', ' + top_freq_mhz
		#print "/nfs/mwa-14/r1/EoRuvfits/batch/" + obsid + ".uvfits" + ', ' + obsid
		print old_path + obsid + ".uvfits",obsid
		#Create the database row, and fill it with the inputs. 
		cur.execute("UPDATE uvfits SET (path,cotter_version,bottom_freq_mhz,top_freq_mhz)=(%s,%s,%s,%s) WHERE (path,obsid)=(%s,%s);", \
			(new_path,cotter_version,bottom_freq_mhz,top_freq_mhz,old_path + obsid + ".uvfits",obsid))
		conn.commit()

	#Commit all the cur.execute, and close the connection.
	conn.commit()
	cur.close()
	conn.close()

if __name__ == '__main__':
	main()
