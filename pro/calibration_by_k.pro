pro calibration_by_k, vis_ptr, vis_model_ptr, vis_weight_ptr, obs, params, cal_k_coverage = cal_k_coverage, cal_k_wedge=cal_k_wedge

	n_pol=2
	n_freq=384.
	c = 299792458.
	
	uu = params.uu * c ; in meters
	vv = params.vv * c ; in meters
	
	z0_freq = 1420.40E6 ;; Hz
	freq = (*obs.baseline_info).freq
	redshifts = z0_freq/freq - 1 ;; frequencies will be identical if kx, ky, kz match
	;mean_redshift = mean(redshifts)
	
	cosmology_measures, redshifts, wedge_factor = wedge_factor, comoving_dist_los=comov_dist_los, Ez=Ez
	
	mean_comoving_dist_los = mean(comov_dist_los)
	
	kx = 2. * !PI * uu / mean_comoving_dist_los
	ky = 2. * !PI * vv / mean_comoving_dist_los
	kperp = sqrt(abs(kx)^2. + abs(ky)^2.)
	
	;*****calculate wedge
	;; assume 20 degrees from pointing center to first null
	source_dist = 20d * !dpi / 180d
	fov_amp = wedge_factor * source_dist
	max_theta= 10d ;random guess!
	
	;; calculate angular distance to horizon
	horizon_amp = wedge_factor * ((max_theta+90d) * !dpi / 180d)
	
	wedge_amp = [fov_amp, horizon_amp] ; to be multiplied by k perp
	wedge_names = ['fov', 'horizon']
	;*****
	
	;*****FT along frequency axis
	for pol_i=0, n_pol-1 do begin
		*vis_ptr[pol_i] = FFT(*vis_ptr[pol_i],DIMENSION=1);,/center)
		*vis_model_ptr[pol_i] = FFT(*vis_model_ptr[pol_i],DIMENSION=1);,/center)
	endfor
	;*****
	
	;***Failed eta conversions
	;deltaeta = 1. / (n_freq * obs.freq_res)
	;eta = FINDGEN(n_freq) * deltaeta
	
	;hubble_param = .71
	;kpar = (eta * 2. * !PI * .71 * z0_freq * Ez) / (c * (1 + mean_redshift)^2.)
	;****
	
	comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
	comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
	
	z_mpc_delta = float(mean(comov_los_diff))
	z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
	kz_mpc_range =  (2.*!pi) / (z_mpc_delta)
	kz_mpc_delta = (2.*!pi) / z_mpc_length
	kz = findgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - kz_mpc_range/2.
	folded_kz = abs(kz)
	
	;Filtering
	kperp_cut = kperp * wedge_amp[0]
	kperp_mask = FLTARR(n_freq, N_elements(kperp))
	kperp_mask[*,*] = 1.
	
	if keyword_set(cal_k_wedge) then begin
		if cal_k_wedge GT 0 then begin
			for freq_i=0, n_freq-1 do begin
				cut_freq_i = where(kperp_cut[freq_i,*] LT folded_kz[freq_i],n_count)
				if n_count GT 0 then kperp_mask[freq_i,cut_freq_i] = 0
			endfor
		endif else begin
			for freq_i=0, n_freq-1 do begin
				cut_freq_i = where(kperp_cut[freq_i,*] GT folded_kz[freq_i],n_count)
				if n_count GT 0 then kperp_mask[freq_i,cut_freq_i] = 0
			endfor
		endelse
	endif
	
	if keyword_set(cal_k_coverage) then begin
		coverage_cut = where(kperp GT .1,n_count)
		if n_count GT 0 then kperp_mask[*,coverage_cut] = 0
	endif
	
	for pol_i=0, n_pol-1 do *vis_weight_ptr[pol_i] = (*vis_weight_ptr[pol_i]) * (kperp_mask)
	
	return
	
end