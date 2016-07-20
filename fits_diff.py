from astropy.io import fits

image_file = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_nb_sim_perfect_farextent1_maxbaseline/output_data/1061316176_uniform_Dirty_XX.fits'
hdu_list = fits.open(image_file)
hdu_list.info()
image_data1 = hdu_list[0].data
hdu_list.close()

image_file = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_nb_sim_perfect_farextent2_maxbaseline/output_data/1061316176_uniform_Dirty_XX.fits'
hdu_list = fits.open(image_file)
hdu_list.info()
image_data2 = hdu_list[0].data
hdu_list.close()

final_image= image_data1 - image_data2

outfile = '/nfs/mwa-00/h1/nbarry/fits_diff.fits'
hdu = fits.PrimaryHDU(final_image)
hdu.writeto(outfile, clobber=True)
