pro paper_hera_cal_plots

	filename = STRARR(2)
	filename[0] ='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_paper_imag_test7/calibration/zen.2457458.55666.xx.uvU_cal.sav'
	filename[1] ='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_hera_x_paper_test/calibration/zen.2457458.55666.xx.uvU_cal.sav'
	
	hera_inds = [80,104,96,64,53,31,65,88,9,20,89,43,105,22,81,10,72,112,97]
	paper_inds = [1,3,4,13,15,16,23,26,37,38,41,42,46,47,49,50,56,57,58,59,61,63,66,67,70,71,73,74,82,83,87,90,98,99,103,106,114,115,116,117,118,119,120,121,122,123,124,125,126,127]
	paper_hex = [2,21,45,17,68,62,0,113,84,100,85,54,69,40,101,102,44,14,86]
	paper_pol = [25,19,48,29,24,28,55,34,27,51,35,75,18,76,5,77,32,78,30,79,33,91,6,92,52,93,7,94,12,95,8,107,11,108,36,109,60,110,39,111]
	
	for file_i=0,1 do begin
		restore, filename[file_i]
		gain = *cal.gain[0]
		n_tile = N_elements(gain[0,*])
		mean_gain = FLTARR(1024,n_tile)
		
		for tile_i=0, n_tile-1 do begin
			resistant_mean, abs(gain[*,tile_i]),2,res_mean
			mean_gain[*,tile_i] = res_mean
		endfor
		
		if file_i EQ 0 then begin
			abs_gain = abs(gain)/mean_gain
			abs_gain[*,[hera_inds,paper_hex,paper_pol]] = 0
		endif else begin
			abs_gain2 = abs(gain)/mean_gain
			abs_gain2[*,[paper_hex,paper_pol]] = 0
		endelse
	endfor
	
	stop
	cgPS_Open,'/nfs/eor-00/h1/nbarry/hera_tile88typical_littlespikes.pdf',/quiet,/nomatch
	cgplot, freq, abs_gain2[*,hera_inds[7]], yrange=[0,10], xtitle='Frequency (MHz)', ytitle = 'Normalized gain amplitude', title= 'HERA hex, tile 88'
	cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
	
	cgPS_Open,'/nfs/eor-00/h1/nbarry/hera_tile31great.pdf',/quiet,/nomatch
	cgplot, freq, abs_gain2[*,hera_inds[5]], yrange=[0,2], xtitle='Frequency (MHz)', ytitle = 'Normalized gain amplitude', title= 'HERA hex, tile 31'
	cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
	
	cgPS_Open,'/nfs/eor-00/h1/nbarry/hera_tile80typical.pdf',/quiet,/nomatch
	cgplot, freq, abs_gain2[*,hera_inds[0]], yrange=[0,10], xtitle='Frequency (MHz)', ytitle = 'Normalized gain amplitude', title= 'HERA hex, tile 80'
	cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
	
	cgPS_Open,'/nfs/eor-00/h1/nbarry/paper_tile4great.pdf',/quiet,/nomatch
	cgplot, freq, abs_gain[*,paper_inds[2]], yrange=[0,2], xtitle='Frequency (MHz)', ytitle = 'Normalized gain amplitude', title= 'PAPER imaging, tile 4'
	cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
	
	cgPS_Open,'/nfs/eor-00/h1/nbarry/paper_tile123badrfi.pdf',/quiet,/nomatch
	cgplot, freq, abs_gain[*,paper_inds[45]], yrange=[0,6], xtitle='Frequency (MHz)', ytitle = 'Normalized gain amplitude', title= 'PAPER imaging, tile 123'
	cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
	
	cgPS_Open,'/nfs/eor-00/h1/nbarry/paper_tile71refl_large.pdf',/quiet,/nomatch
	cgplot, freq, abs_gain[*,paper_inds[25]], yrange=[0,6], xtitle='Frequency (MHz)', ytitle = 'Normalized gain amplitude', title= 'PAPER imaging, tile 71'
	cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
	
	cgPS_Open,'/nfs/eor-00/h1/nbarry/paper_tile26_refl.pdf',/quiet,/nomatch
	cgplot, freq, abs_gain[*,paper_inds[7]], yrange=[0,1.5], xtitle='Frequency (MHz)', ytitle = 'Normalized gain amplitude', title= 'PAPER imaging, tile 26'
	cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
end