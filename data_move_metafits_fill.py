#!/usr/bin/python

import os
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

#********************************
#Script that moves uvfits, makes metafits at MIT cluster
#and fills the database with necessary information

#Script written by Nichole Barry, May 2017.

def main():

	#test!!
	#obsids = ['1064768840']

	#Assume cotter version is the same for files in eor-14 batch folder, and version, subversion
	version=5
	subversion=1

	all_nodes = ["eor-10"]#, "eor-13", "eor-14"]
	all_nodes = ["/nfs/" + nodename + "/r1/" for nodename in all_nodes]

	#Get obsids to download
	obsfile_name='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/obs_list/season2_zenith_rest.txt'
	obsfile = open(obsfile_name, "r")
	obsids = [line.split( ) for line in obsfile.readlines()]
	obsids = [obs[0] for obs in obsids]
	obsfile.close()
	nonredundant_obsids = list(set(obsids))
	if len(obsids) != len(nonredundant_obsids):
		print "WARNING: Obs list contains redundant entries."
		obsids = nonredundant_obsids

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

		os.system('python /nfs/mwa-00/h1/nbarry/MWA_Tools/scripts/make_metafits.py -f /nfs/mwa-10/r1/EoRuvfits/EoR2014/' + obsid + '.uvfits -g ' + obsid)

		save_path = all_nodes[0] + save_directories[obs_index] + obsid #put approximately 500 obsids on specified node if following test turns up false
		for loc_node in all_nodes:
			loc_path = loc_node + save_directories[obs_index]
			if os.path.isdir(loc_path): #checks to see if the directory exists
				save_path = loc_path + obsid
		#test!!
		print save_path

		#Open up the metafits file that was made with the uvfits file (assumes they are in the same location)
		metafits_file = "/nfs/mwa-10/r1/EoRuvfits/EoR2014/" + obsid + ".metafits"
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

		uvfits_file = "/nfs/mwa-10/r1/EoRuvfits/EoR2014/" + obsid + ".uvfits"
		hdu_list_uvfits = fits.open(uvfits_file)
		header_uvfits = hdu_list_uvfits[0].header

		#Get the frequency center and bandwidth (assumes contiguous frequency channels)
		cotter_version = header_uvfits['COTVER']
		hdu_list_uvfits.close()

    		if not os.path.exists(save_path):
        		os.makedirs(save_path)
		new_path = save_path + "/" + obsid + ".uvfits"
		shutil.move("/nfs/eor-10/r1/EoRuvfits/EoR2014/" + obsid + ".uvfits", new_path)
		shutil.move("/nfs/eor-10/r1/EoRuvfits/EoR2014/" + obsid + ".metafits", save_path + "/" + obsid + ".metafits")

		#Get the time for the timestamp in UTC
		timestamp = str(datetime.datetime.utcnow())

		db_comment = 'Season 2 (2014) EoR0 zenith pointing'

		#Check to make sure that a similar uvfits file does not already exist, throw warning if it does.
		cur.execute("SELECT uvfits.obsid FROM uvfits WHERE (obsid,version,subversion,cotter_version,bottom_freq_mhz,top_freq_mhz)=(%s,%s,%s,%s,%s,%s);", \
			(obsid,version,subversion,cotter_version,bottom_freq_mhz,top_freq_mhz))
		if cur.fetchall():
			print "WARNING: A uvfits file for obsid " + obsid + ", version " + version + ", subversion " + subversion + \
				", cotter " + cotter_version + ", and frequency range " + bottom_freq_mhz + "-" + top_freq_mhz + " already exists."

		#Create the database row, and fill it with the inputs. 
		cur.execute("INSERT INTO uvfits(obsid,version,subversion,path,cotter_version,timestamp,comment,bottom_freq_mhz,top_freq_mhz) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s);", \
			(obsid,version,subversion,new_path,cotter_version,timestamp,db_comment,bottom_freq_mhz,top_freq_mhz))
		conn.commit()

	#Commit all the cur.execute, and close the connection.
	conn.commit()
	cur.close()
	conn.close()

if __name__ == '__main__':
	main()
