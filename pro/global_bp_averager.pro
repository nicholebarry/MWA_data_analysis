pro global_bp_averager

  ; Options
  reconstruct_bp = 0
  flag_long_tiles = 0

  ; Dir, obs/cal to be averaged
  ;dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_cal_stop_novanvleck/'
  ;dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_stop/'
  dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_2014_data_pointing_beam_coarse_corr_no_ao_nodiffuse_tenthperc_noFoV_cal_stop/'
  cal_dir = 'calibration/'
  obs_dir = 'metadata/'

  ; Parameters for the observation id files
  dir_obsids = '/home/nbarry/MWA/pipeline_scripts/bash_scripts/ozstar/obs_list/2014/'
  day = '2014_1'
;  parsednames=[day+'_minustwo',day+'_minusone',day+'_zenith',day+'_plusone',day+'_plustwo']
  parsednames=[day+'_plusone']
  n_pointings = N_elements(parsednames)
  obs_ptr=PTRARR(n_pointings,/allocate_heap)

  ; Get observation ids and put them in a pointing pointer
  FOR poi_i=0,n_pointings-1 DO BEGIN
    readcol, dir_obsids + parsednames[poi_i]+'.txt', obs_temp, format='A', /silent
    *obs_ptr[poi_i]=obs_temp
  ENDFOR

  ; Get a sample obs file so that the polyfit params can be extracted
  obs = getvar_savefile(dir + obs_dir + (*obs_ptr[0])[0] + '_obs.sav','obs')
  freq_use=(*obs.baseline_info).freq_use
  freq_use=where(freq_use,nf_use)
  pre_dig_inds = where((*obs.baseline_info).freq[freq_use] LT 187.515E6,n_count)
  f_d = max(pre_dig_inds)
  f_end = N_elements(freq_use)-1
  amp_degree = 2
  phase_degree = 1

  ; Typical observational parameters
  n_freq = 384
  n_tile = 128
  n_pol = 2

  ; Pointing loop
  for poi_i=0, n_pointings - 1 do begin

    obs_ids = *obs_ptr[poi_i]
    n_obs = N_elements(obs_ids)
 
    ; Observation pre-loop to check that the files exist and curate a more representative list
    n_obs_checked=0
    for obs_i=0,n_obs-1 do begin
      if ~file_test(dir + cal_dir + (*obs_ptr[poi_i])[obs_i] + '_cal.sav') then begin
        print, dir + cal_dir + (*obs_ptr[poi_i])[obs_i] + "_cal.sav not found, skipping"
      endif else begin
        if n_obs_checked EQ 0 then obs_ids_checked = (*obs_ptr[poi_i])[obs_i] $
          else obs_ids_checked = [obs_ids_checked,(*obs_ptr[poi_i])[obs_i]]
        n_obs_checked += 1
      endelse
    endfor
    n_obs = n_obs_checked
    obs_ids = obs_ids_checked 

    ; Init the various gain arrays that will be used
    gain_arr = complex(fltarr(n_pol,n_freq,n_tile,n_obs))
    gain_orig = complex(fltarr(n_pol,n_freq,n_tile,n_obs))
    gain_bandpass = complex(fltarr(n_pol,n_freq,n_tile,n_obs))
    polyfit = complex(fltarr(n_pol,n_freq,n_tile,n_obs))
    polyfit_global = fltarr(n_pol,n_freq,n_tile,n_obs)

    ; Observation loop
    for obs_i=0,n_obs-1 do begin

      ; Restore the cal file
;      if ~file_test(dir + cal_dir + (*obs_ptr[poi_i])[obs_i]) then continue
      restore, dir + cal_dir + obs_ids[obs_i] + '_cal.sav'

      ; Polarization loop
      for pol_i=0,1 do begin
        gain_arr[pol_i,*,*,obs_i] = *cal.gain[pol_i]
        gain_orig[pol_i,*,*,obs_i] = *cal.gain[pol_i] + *cal.gain_residual[pol_i]

        ;Reconstruct per-tile polyfit (digital gain jump style)
        for tile_i=0, n_tile-1 do begin
          gain_fit=fltarr(n_freq)
          fit_params1= (*cal.amp_params[pol_i,tile_i])[0,*]
          fit_params2= (*cal.amp_params[pol_i,tile_i])[1,*]
          FOR di=0L,amp_degree-1 DO gain_fit[freq_use[0]:freq_use[f_d]] += fit_params1[di]*findgen(freq_use[f_d])^di
          FOR di=0L,amp_degree-1 DO gain_fit[freq_use[f_d+1]:freq_use[f_end]] += $
            fit_params2[di]*(findgen(freq_use[f_end] - freq_use[f_d+1]+1) + freq_use[f_d+1])^di

          phase_fit=fltarr(n_freq)
          FOR di=0L,phase_degree DO phase_fit+=(*cal.phase_params[pol_i,tile_i])[0,di]*findgen(n_freq)^di
        
          ; Per-tile polyfit
          polyfit[pol_i,*,tile_i,obs_i]=gain_fit*Exp(Complex(0,1)*phase_fit)
        endfor

        ; if it exists easily, grab the global bandpass
        ; otherwise, reset the gains for reconstruction
        if ~keyword_set(reconstruct_bp) then begin
          gain_bandpass[pol_i,*,*,obs_i] = *cal.gain_bandpass[pol_i]
        endif else *cal.gain[pol_i] = *cal.gain[pol_i] + *cal.gain_residual[pol_i]

      endfor ; end pol loop

      ; if it doesn't exist easily, reconstruct the global bandpass
      ; takes time!
      if keyword_set(reconstruct_bp) then begin
        uvfits_read,hdr,params,layout,vis_arr,vis_weights,file_path_vis='/fred/oz048/MWA/data/2013/van_vleck_corrected/' + obs_ids[obs_i] +'.uvfits',$
          silent=silent,error=error,_Extra=extra
        restore, dir + obs_dir + obs_ids[obs_i] + '_obs.sav'
        vis_auto=vis_extract_autocorr(obs,vis_arr = vis_arr,/time_average,auto_tile_i=auto_tile_i)
        cal_autos_init=cal_auto_ratio(obs,cal,auto_ratio=auto_ratio,vis_auto=vis_auto,auto_tile_i=auto_tile_i,/divide)

        if keyword_set(flag_long_tiles) then begin
          long_tiles = [78,79,87,88,95,96,97,98,104,112,113,122,123,124]
          obs_flag_long = obs
          (*obs_flag_long.baseline_info).tile_use[long_tiles] = 0
          cal_bandpass=vis_cal_bandpass(cal_autos_init,obs_flag_long,params,cal_remainder=cal_remainder,auto_ratio_calibration=1)
        endif else cal_bandpass=vis_cal_bandpass(cal_autos_init,obs,params,cal_remainder=cal_remainder,auto_ratio_calibration=1)
        
        for pol_i=0,n_pol-1 do begin
          gain_bandpass[pol_i,*,*,obs_i] = *cal_bandpass.gain[pol_i]
        endfor
      endif

      ; Once the bandpass is for-sure reconstructed, fit a polyfit to the global bandpass
      for pol_i=0, n_pol-1 do begin
        gain = abs(reform(gain_bandpass[pol_i,freq_use,0,obs_i])) ;all tiles are the same
        ;Only fit for amplitude if amp_degree is set and greater than zero
        ;Fit pre- and post-digital gain jump separately in highband MWA data
        gain_fit=fltarr(n_freq)
        fit_params1=poly_fit(freq_use[0:f_d],gain[0:f_d],amp_degree-1)
        fit_params2=poly_fit(freq_use[f_d+1:f_end],gain[f_d+1:f_end],amp_degree-1)
        FOR di=0L,amp_degree-1 DO gain_fit[freq_use[0]:freq_use[f_d]] += fit_params1[di]*findgen(freq_use[f_d])^di
        FOR di=0L,amp_degree-1 DO gain_fit[freq_use[f_d+1]:freq_use[f_end]] += $
          fit_params2[di]*(findgen(freq_use[f_end] - freq_use[f_d+1]+1) + freq_use[f_d+1])^di

        polyfit_global[pol_i, *,*,obs_i] =  rebin(gain_fit,n_freq,n_tile)

      endfor

    endfor ; end observation loop


    gain_auto_ratio = abs(gain_arr / gain_bandpass / polyfit)
    
    res_bandpass = abs(gain_bandpass) / polyfit_global
    zeros = where(res_bandpass EQ 0, n_count)
    if n_count GT 0 then res_bandpass[zeros] = !Values.F_NAN ;to exclude from the mean

    ; This is a little redundant given that each tile has the same value, but it won't affect the mean, only the runtime
    ; Rebin (i.e. repeat values) so that division is easy later
    ave_res_bandpass = rebin( mean(res_bandpass,dimension=4,/NAN), n_pol,n_freq,n_tile,n_obs)

    ; Reconstruct the gains with the averaged bandpass instead
    new_gains = (ave_res_bandpass * polyfit_global) * abs(polyfit) * gain_auto_ratio * exp(Complex(0,1)*atan(gain_arr,/phase))

    ; Save the reconstructing gains in per obs cal files for ease of read-in
    for obs_i=0, n_obs-1 do begin
      restore, dir + cal_dir + obs_ids[obs_i] + '_cal.sav'
      for pol_i=0, n_pol-1 do *cal.gain[pol_i] = reform(new_gains[pol_i,*,*,obs_i])
      for pol_i=0, n_pol-1 do *cal.gain_residual[pol_i] = reform(gain_orig[pol_i,*,*,obs_i] - new_gains[pol_i,*,*,obs_i])
      if keyword_set(reconstruct_bp) then cal = structure_update(cal, gain_bandpass=cal.gain)
      for pol_i=0, n_pol-1 do *cal.gain_bandpass[pol_i] = reform((ave_res_bandpass * polyfit_global)[pol_i,*,*,obs_i])
      save, cal, filename = dir + 'cal_transfer/' + obs_ids[obs_i] + '_cal.sav'
    endfor


  endfor ; end pointing loop





stop

end
