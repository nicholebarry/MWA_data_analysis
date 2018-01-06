from astropy.io import fits
import numpy as np

image_file = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_channel_7_midobs_-3/output_data/1061313496_uniform_Dirty_XX.fits'
hdu_list = fits.open(image_file)
image_data1 = hdu_list[0].data
hdu_list.close()


image_file = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_channel_7_midobs_-2/output_data/1061313496_uniform_Dirty_XX.fits'
hdu_list = fits.open(image_file)
image_data2 = hdu_list[0].data
hdu_list.close()

#final_image= abs(image_data1 - image_data2) / image_data2
#zero_inds = np.where( abs(image_data2) < .05 )
#final_image[zero_inds]=0
#final_image_crop = final_image[799:1247,799:1247]

final_image = abs(image_data1 - image_data2)

outfile = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_channel_7_midobs/diff_data/dirty_xx_-2.fits'
hdu = fits.PrimaryHDU(final_image)
hdu.writeto(outfile, clobber=True)

#optional
import matplotlib.pyplot as plt
plt.hist(final_image_crop, bins='auto')
plt.title("Percent difference between GLEAM and KGS in primary")
plt.ylabel("counts")
plt.xlabel("Percent difference between restored images")
plt.axis([0, 1, 0, 50])
plt.show
