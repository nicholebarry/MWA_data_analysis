pro beam_image_maker

	if keyword_set(create) then begin
		restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e3/beams/1061316296_beams.sav'
		print, "beam 1"
		restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e3/metadata/1061316296_obs.sav'
		print, 'obs 1"
		beam_arr_thresh1e3 = beam_image_cube(obs,psf,n_freq_bin=384)
		save, beam_arr_thresh1e3, filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e3/1061316296_beam_image.sav'
		
		restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e2_interpolate/beams/1061316296_beams.sav'
		print, "beam 1"
		restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e2_interpolate/metadata/1061316296_obs.sav'
		print, 'obs 1"
		beam_arr_thresh1e2_int = beam_image_cube(obs,psf,n_freq_bin=384)
		save, beam_arr_thresh1e2_int, filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e2_interpolate/1061316296_beam_image.sav'
		
		restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh5e2/beams/1061316296_beams.sav'
		print, "beam 1"
		restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh5e2/metadata/1061316296_obs.sav'
		print, 'obs 1"
		beam_arr_thresh5e2 = beam_image_cube(obs,psf,n_freq_bin=384)
		save, beam_arr_thresh5e2, filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh5e2/1061316296_beam_image.sav'
		
		restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e4/beams/1061316296_beams.sav'
		print, "beam 1"
		restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e4/metadata/1061316296_obs.sav'
		print, 'obs 1"
		beam_arr_thresh1e4 = beam_image_cube(obs,psf,n_freq_bin=384)
		save, beam_arr_thresh1e4, filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e4/1061316296_beam_image.sav'
		
		restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e1/beams/1061316296_beams.sav'
		print, "beam 1"
		restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e1/metadata/1061316296_obs.sav'
		print, 'obs 1"
		beam_arr_thresh1e1 = beam_image_cube(obs,psf,n_freq_bin=384)
		save, beam_arr_thresh1e1, filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e1/1061316296_beam_image.sav'
	endif
	
	restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e3/1061316296_beam_image.sav'
	restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e2/1061316296_beam_image.sav'
	restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e2_interpolate/1061316296_beam_image.sav'
	restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh5e2/1061316296_beam_image.sav'
	restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e4/1061316296_beam_image.sav'
	restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_beam2b_thresh1e4/metadata/1061316296_obs.sav'
	
	thresh1e3 = FLTARR(384)
	thresh1e2 = FLTARR(384)
	thresh1e2_int = FLTARR(384)
	thresh5e2 = FLTARR(384)
	thresh1e4 = FLTARR(384)

	
	for freq_i=0, 383 do thresh1e3[freq_i] = (*beam_arr_thresh1e3[0,freq_i])[1023,1023]
	for freq_i=0, 383 do thresh1e2[freq_i] = (*beam_arr_thresh1e2[0,freq_i])[1023,1023]
	for freq_i=0, 383 do thresh1e2_int[freq_i] = (*beam_arr_thresh1e2_int[0,freq_i])[1023,1023]
	for freq_i=0, 383 do thresh5e2[freq_i] = (*beam_arr_thresh5e2[0,freq_i])[1023,1023]
	for freq_i=0, 383 do thresh1e4[freq_i] = (*beam_arr_thresh1e4[0,freq_i])[1023,1023]

	
	freq_arr = (*obs.baseline_info).freq / 1e6
	
	cgPS_Open,'/nfs/eor-00/h1/nbarry/beamvalues_zenith.png',/quiet,/nomatch
	cgplot, freq_arr,thresh1e3, yrange = [.95,1.05], color='black', title = 'Beam value in center pixel at zenith', ytitle='Beam value', xtitle='Frequency (MHz)' , charsize=1
	cgoplot, freq_arr,thresh1e2, color='blue'
	cgoplot, freq_arr,thresh1e2_int, color='green'
	cgoplot, freq_arr,thresh5e2, color='purple'
	cgoplot, freq_arr,thresh1e4, color='teal'

	cglegend, title=['1e3','1e2','1e2 interpolated','5e2','1e4'], color=['black','blue','green','purple','teal'],location=[.65,.75],charsize=1
	cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
	
	
end