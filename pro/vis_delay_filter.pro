pro vis_delay_filter, vis_model_ptr, vis_weight_ptr, params, obs
	; This is a script to generate delay spectra from visibilities
	; TODO: Contruct file_path_fhd from dir and obsid, and then use standard fhd_save_io
	; TODO: Handle backward compatibility (see read flags lines)
	; TODO: Filter out flagged baselines better (ie, bl_use array)
	; TODO: Put data into physical units
	;save, vis_model_ptr, filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/vis_model_ptr.sav'
	i_comp = Complex(0,1)
	
	;;;;; prepare data
	
	;u,v,w are in light travel time in seconds
	freq_arr = (*obs.baseline_info).freq
	freq_res = obs.freq_res
	n_pol = obs.n_pol
	nfreq = n_elements(freq_arr)
	nbl = n_elements(params.uu)
	flags = fltarr(nfreq,nbl,n_pol)
	
	;for poli=0,1 do flags[*,*,poli] = *flag_arr[poli] ; TODO: handle backward compatibility
	;undefine_fhd,flag_arr
	;vis_weight_ptr_input = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_GLEAM_beam2b_filtered/vis_data/1061316296_flags.sav','vis_weights')
	for pol_i=0,n_pol-1 do flags[*,*,pol_i] = *vis_weight_ptr[pol_i]
	
	
	data = Complex(fltarr(nfreq,nbl,n_pol)) ; Stack pols
	for pol_i=0,n_pol-1 do data[*,*,pol_i] = *vis_model_ptr[pol_i]
	
	; Only keep unflagged baselines
	flag_test = Total(Total(flags>0,1),n_pol)
	bi_use=where(flag_test eq 2*nfreq)
	
	;test with removing zeroed visibilities instead
	;total_data = total(total(data,1),2)
	;bi_use = where(total_data ne 0)
	
	data = data[*,bi_use,*]
	uu = params.uu[bi_use]
	vv = params.vv[bi_use]
	ww = params.ww[bi_use]
	bb = sqrt(uu^2.+vv^2.+ww^2.)
	nbl = n_elements(bi_use)
	undefine_fhd,flags
	; Phase to zenith (see Danny for explanation)
	w_mat = freq_arr#ww ; This should now be in wavelengths
	for pol_i=0,n_pol-1 do data[*,*,pol_i] *= exp(i_comp * 2. * !pi * w_mat)

	; FFT
	spectra = shift(fft(data,dim=1),nfreq/2,0,0) ; Shift only in fft direction.
	
	;Cut at the horizon
	tau_cut=1.
	;Calculate upper and lower delay limits for each baseline
	lower_limit = nfreq*(.5 - bb * tau_cut * freq_res)
	upper_limit = nfreq*(.5 + bb * tau_cut * freq_res)
	
	for freq_i=0, nfreq-1 do begin
		mask_high = freq_i/upper_limit
		mask_high_inds = where(mask_high GT 1., n_count)
		if n_count GT 0 then spectra[freq_i,mask_high_inds,*] = 0
		;print, freq_i
		;if n_count GT 0 then print, mask_high_inds[0],  mask_high_inds[N_elements(mask_high_inds)-1]
		mask_low = freq_i/lower_limit
		mask_low_inds = where(mask_low LT 1., n_count)
		if n_count GT 0 then spectra[freq_i,mask_low_inds,*] = 0
	;if n_count GT 0 then print, mask_low_inds[0],  mask_low_inds[N_elements(mask_low_inds)-1]
	endfor
	
	masked_data = fft(spectra,dim=1,/inverse)
	
	; UnPhase from zenith?
	w_mat = freq_arr#ww ; This should now be in wavelengths
	for pol_i=0,n_pol-1 do masked_data[*,*,pol_i] *= exp(-i_comp * 2. * !pi * w_mat)
	
	for pol_i=0, n_pol-1 do (*vis_model_ptr[pol_i])[*,bi_use] = masked_data[*,*,pol_i]
	
	;save, vis_model_ptr, filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/vis_model_ptr_filtered.sav'
	
	return

end