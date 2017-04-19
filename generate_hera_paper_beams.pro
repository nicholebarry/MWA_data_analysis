pro generate_HERA_PAPER_beams

	antenna_paper = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_paper_imag_test2/beams/zen.2457458.16694.xx.uvAU_antenna.sav','antenna')
	antenna_hera = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_hera_imag_test/beams/zen.2457458.16694.xx.uvAU_antenna.sav','antenna')
	
	hera_inds = [80,104,96,64,53,31,65,88,9,20,89,43,105,22,81,10,72,112,97]
	paper_inds = [1,3,4,13,15,16,23,26,37,38,41,42,46,47,49,50,56,57,58,59,61,63,66,67,70,71,73,74,82,83,87,90,98,99,103,106,114,115,116,117,118,119,120,121,122,123,124,125,126,127]
	
	paper_hex = [2,21,45,17,68,62,0,113,84,100,85,54,69,40,101,102,44,14,86]
	paper_pol = [25,19,48,29,24,28,55,34,27,51,35,75,18,76,5,77,32,78,30,79,33,91,6,92,52,93,7,94,12,95,8,107,11,108,36,109,60,110,39,111]
	
	antenna = antenna_paper
	antenna[hera_inds] = antenna_hera[hera_inds-1]
	
	save, antenna, filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_hera_imag_test/heraxpaper_antenna.sav'
	
end