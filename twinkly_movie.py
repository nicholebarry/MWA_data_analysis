from astropy.io import fits
import numpy as np

for time_i in range(2,51):
	if time_i == 47:
		continue
	time_step = str(time_i).zfill(3) 
	image_file = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_channel_7_'+time_step+'/output_data/1061313496_uniform_Dirty_XX.fits'
	hdu_list = fits.open(image_file)
	image_data1 = hdu_list[0].data
	if time_i == 2:
		image_full = image_data1
	else:
		image_full = image_full + image_data1
	hdu_list.close()

image_full = image_full / 52.

for time_i in range(2,51):
	if time_i == 47:
		continue
	time_step = str(time_i).zfill(3) 
	image_file = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_channel_7_'+time_step+'/output_data/1061313496_uniform_Dirty_XX.fits'
	hdu_list = fits.open(image_file)
	image_data1 = hdu_list[0].data
	final_image = abs(image_data1 - image_full)

	outfile = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_channel_7_'+time_step+'/dirty_xx_avesub'+time_step+'.fits'
	hdu = fits.PrimaryHDU(final_image)
	hdu.writeto(outfile, clobber=True)

