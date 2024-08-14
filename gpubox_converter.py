from pyuvdata import UVData
import numpy as np
UV = UVData()

# Construct the list of files
data_path = '/fred/oz048/MWA/data/1061316296/rts_thermal_noise/'
filelist = [data_path + i for i in ['1061316296.metafits',
     '1061316296_thermal01_00.fits','1061316296_thermal01_01.fits',
     '1061316296_thermal02_00.fits','1061316296_thermal02_01.fits',
     '1061316296_thermal03_00.fits','1061316296_thermal03_01.fits',
     '1061316296_thermal04_00.fits','1061316296_thermal04_01.fits',
     '1061316296_thermal05_00.fits','1061316296_thermal05_01.fits',
     '1061316296_thermal06_00.fits','1061316296_thermal06_01.fits']]
#     '1061316296_thermal07_00.fits','1061316296_thermal07_01.fits',
#     '1061316296_thermal08_00.fits','1061316296_thermal08_01.fits',
#     '1061316296_thermal09_00.fits','1061316296_thermal09_01.fits',
#     '1061316296_thermal10_00.fits','1061316296_thermal10_01.fits',
#     '1061316296_thermal11_00.fits','1061316296_thermal11_01.fits',
#     '1061316296_thermal12_00.fits','1061316296_thermal12_01.fits']]
#     '1061316296_thermal13_00.fits','1061316296_thermal13_01.fits',
#     '1061316296_thermal14_00.fits','1061316296_thermal14_01.fits',
#     '1061316296_thermal15_00.fits','1061316296_thermal15_01.fits',
#     '1061316296_thermal16_00.fits','1061316296_thermal16_01.fits',
#     '1061316296_thermal17_00.fits','1061316296_thermal17_01.fits',
#     '1061316296_thermal18_00.fits','1061316296_thermal18_01.fits']]
#     '1061316296_thermal19_00.fits','1061316296_thermal19_01.fits',
#     '1061316296_thermal20_00.fits','1061316296_thermal20_01.fits',
#     '1061316296_thermal21_00.fits','1061316296_thermal21_01.fits',
#     '1061316296_thermal22_00.fits','1061316296_thermal22_01.fits',
#     '1061316296_thermal23_00.fits','1061316296_thermal23_01.fits',
#     '1061316296_thermal24_00.fits','1061316296_thermal24_01.fits']]

# Use the `read` method, optionally specify the file type. Can also use the
# file type specific `read_mwa_corr_fits` method, but only if reading files
# from a single observation.
# Apply cable corrections and phase data before writing to uvfits
# Skip routine time/frequency flagging - see flag_init and associated keywords in documentation
UV.read(filelist, correct_cable_len=True, phase_to_pointing_center=True, flag_init=False,data_array_dtype=np.complex64)
#UV.read(filelist, file_type='mwa_corr_fits', correct_cable_len=True, phase_to_pointing_center=True, flag_init=False, data_array_dtype=np.complex64)
#UV.read_mwa_corr_fits(filelist, correct_cable_len=True, phase_to_pointing_center=True, flag_init=False)

# Write out uvfits file
write_file = '/fred/oz048/MWA/data/1061316296/rts_thermal_noise/1061316296_1.uvfits'
UV.write_uvfits(write_file, spoof_nonessential=True)
