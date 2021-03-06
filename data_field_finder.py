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
#Script to find field of database query

#Script written by Nichole Barry, May 2017.

def main():

	version=4

	#Get obsids to download
	obsfile_name='/nfs/mwa-00/h1/nbarry/database_query.txt'
	obsfile = open(obsfile_name, "r")
	obsids = [line.split( ) for line in obsfile.readlines()]
	obsids = [obs[0] for obs in obsids]
	obsfile.close()
	nonredundant_obsids = list(set(obsids))
	if len(obsids) != len(nonredundant_obsids):
		print "WARNING: Obs list contains redundant entries."
		obsids = nonredundant_obsids

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

		cur.execute("SELECT uvfits.path, uvfits.obsid FROM uvfits WHERE (obsid,version)=(%s,%s);", \
			(obsid,version))
		#path = cur.fetchall()
		for path1, obsid1 in cur.fetchall():
			path = path1
			#print path1

		#Open up the metafits file that was made with the uvfits file (assumes they are in the same location)
		#metafits_file = "/nfs/mwa-14/r1/EoRuvfits/batch/" + obsid + ".metafits"
		uvfits_file = path
		if not os.path.exists(uvfits_file):
			print metafits_file + ' does not exist, skipping'
			continue
		hdu_list_uvfits = fits.open(uvfits_file)
		header_uvfits = hdu_list_uvfits[0].header

		#Get the frequency center and bandwidth (assumes contiguous frequency channels)
		gridname = header_uvfits['GRIDNAME']
		print obsid, gridname

		#Close the metafits file
		hdu_list_uvfits.close()

	#Commit all the cur.execute, and close the connection.
	conn.commit()
	cur.close()
	conn.close()

if __name__ == '__main__':
	main()
