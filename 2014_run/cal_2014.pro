pro cal_2014

  ;Generate statistics based on the number of converged frequencies and the number of iterations
  ;Perform high-level observation flagging if the calibration is poor
  ;Create pointing-based calibration solutions from remaining observations

  dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_calstop/'
  cal_dir = 'calibration/'
  obs_dir = 'metadata/'

  ; Parameters for the observation id files
  dir_obsids = '/home/nbarry/MWA/pipeline_scripts/bash_scripts/ozstar/obs_list/2014/'
  day = ['10','11','13','14','15','16','18','19','21','22','23']
  n_days = N_elements(day)
  pointingnames=['_minustwo_ssins','_minusone_ssins',$
    '_zenith_ssins','_plusone_ssins','_plustwo_ssins']
  n_pointings = N_elements(pointingnames)
  obs_ptr=PTRARR(n_days,n_pointings,/allocate_heap)

  max_n_obs=0
  ; Get observation ids and put them in a pointing pointer
  FOR poi_i=0,n_pointings-1 DO BEGIN
    max_n=0
    FOR day_i=0,n_days-1 DO BEGIN
      readcol, dir_obsids + day[day_i] + pointingnames[poi_i]+'.txt', obs_temp, format='A', /silent
      *obs_ptr[day_i,poi_i]=obs_temp
      max_n_obs = max([N_elements(obs_temp),max_n_obs])
    ENDFOR
  ENDFOR
  

  ; Get a sample obs file so that the polyfit params can be extracted
  obs = getvar_savefile(dir + obs_dir + (*obs_ptr[0,0])[0] + '_obs.sav','obs')
  freq_use=(*obs.baseline_info).freq_use
  freq_use=where(freq_use,nf_use)
  pre_dig_inds = where((*obs.baseline_info).freq[freq_use] LT 187.515E6,n_count)
  f_d = max(pre_dig_inds)
  f_end = N_elements(freq_use)-1
  amp_degree = 2
  phase_degree = 1

  ; Typical observational parameters
  n_freq = obs.n_freq
  n_tile = obs.n_tile
  n_pol = obs.n_pol

  n_conv = lonarr(n_days,n_pointings,n_pol,max_n_obs)
  conv_iter = fltarr(n_days,n_pointings,n_pol,n_freq,n_tile,max_n_obs)
  n_obs_arr = fltarr(n_days,n_pointings)

  gain_arr = complex(fltarr(n_days,n_pointings,n_pol,n_freq,n_tile,max_n_obs))
  final_gain = complex(fltarr(n_days,n_pointings,n_pol,n_freq,n_tile,max_n_obs))
  gain_res = complex(fltarr(n_days,n_pointings,n_pol,n_freq,n_tile,max_n_obs))
  stddev_res = fltarr(n_days,n_pointings,n_pol,n_tile,max_n_obs)
  gain_bandpass = complex(fltarr(n_days,n_pointings,n_pol,n_freq,n_tile,max_n_obs))
  polyfit = complex(fltarr(n_days,n_pointings,n_pol,n_freq,n_tile,max_n_obs))
  polyfit_global = fltarr(n_days,n_pointings,n_pol,n_freq,n_tile,max_n_obs) ;polyfit of the global bandpass to remove extra large-band structure
  first_defined_tile = fltarr(n_days,n_pointings,n_pol,max_n_obs)
  first_defined_tile[*]=-1 ;set this flag in order to query whether some obs toward the end of the array exist 

  for day_i=0, n_days-1 do begin
    for poi_i=0, n_pointings-1 do begin

      obs_ids = *obs_ptr[day_i,poi_i]
      n_obs = N_elements(obs_ids)
      n_obs_arr[day_i,poi_i] = n_obs

      for obs_i=0,n_obs-1 do begin

        restore, dir + cal_dir + obs_ids[obs_i] + '_cal.sav'
        ;Gather statistics and gains
        for pol_i=0,n_pol-1 do begin
          gain_arr[day_i,poi_i,pol_i,*,*,obs_i] = *cal.gain[pol_i]
          gain_res[day_i,poi_i,pol_i,*,*,obs_i] = *cal.gain_residual[pol_i]
          gain_bandpass[day_i,poi_i,pol_i,*,*,obs_i] = *cal.gain_bandpass[pol_i]
          stddev_res[day_i,poi_i,pol_i,*,obs_i] = stddev(abs(*cal.gain_residual[pol_i]), /NAN, dimension=1)
          n_conv[day_i,poi_i,pol_i,obs_i] = cal.n_converged[pol_i]
          conv_iter[day_i,poi_i,pol_i,*,*,obs_i] = *cal.conv_iter[pol_i]

          ;Reconstruct polyfit that was fit
          for tile_i=0, n_tile-1 do begin
            if (first_defined_tile[day_i,poi_i,pol_i,obs_i] EQ -1) AND (mean(gain_bandpass[day_i,poi_i,pol_i,*,tile_i,obs_i]) GT 0) $
             then first_defined_tile[day_i,poi_i,pol_i,obs_i]=tile_i

            gain_fit=fltarr(n_freq)
            fit_params1= (*cal.amp_params[pol_i,tile_i])[0,*];poly_fit(freq_use[0:f_d],gain[0:f_d],amp_degree-1)
            fit_params2= (*cal.amp_params[pol_i,tile_i])[1,*];poly_fit(freq_use[f_d+1:f_end],gain[f_d+1:f_end],amp_degree-1)
            FOR di=0L,amp_degree-1 DO gain_fit[freq_use[0]:freq_use[f_d]] += fit_params1[di]*findgen(freq_use[f_d])^di
            FOR di=0L,amp_degree-1 DO gain_fit[freq_use[f_d+1]:freq_use[f_end]] += $
            fit_params2[di]*(findgen(freq_use[f_end] - freq_use[f_d+1]+1) + freq_use[f_d+1])^di

            phase_fit=fltarr(n_freq)
            FOR di=0L,phase_degree DO phase_fit+=(*cal.phase_params[pol_i,tile_i])[0,di]*findgen(n_freq)^di

            polyfit[day_i,poi_i,pol_i,*,tile_i,obs_i]=gain_fit*Exp(Complex(0,1)*phase_fit)
          endfor ;end tile loop

          zero_mean_individual_res = abs((gain_arr[day_i,poi_i,pol_i,*,*,obs_i] + gain_res[day_i,poi_i,pol_i,*,*,obs_i]) / $
            polyfit[day_i,poi_i,pol_i,*,*,obs_i] / gain_bandpass[day_i,poi_i,pol_i,*,*,obs_i]) - 1
          zeros = where(zero_mean_individual_res EQ 0)
          zero_mean_individual_res[zeros]=!Values.F_NAN
          deviation =  stddev(zero_mean_individual_res,/nan)

          inds = where(zero_mean_individual_res GT deviation*4,n_count)
          if n_count GT 0 then begin
            bad_tile = inds / n_freq
            gain_arr[day_i,poi_i,pol_i,*,bad_tile,obs_i] = 0
          endif

        endfor ;end pol loop
      endfor ;end obs loop
    endfor ;end pointing loop
  endfor ;end day loop

zeros = where(gain_arr EQ 0)
gain_arr[zeros]=!Values.F_NAN
;Final amplitude fit with the individual tile polyfit (without reflections) and the global obs bandpass removed, and zero mean
; Will have antenna-dependent parameters (i.e. cable reflections) picked up from the autos
zero_mean_individual_amps = abs(gain_arr / polyfit / gain_bandpass) - 1

  ;Set zeros in the bandapss, including coarse band edges and flagged tiles, to 0 to exclude them from averaging?
  zeros = where(~finite(gain_bandpass))
  gain_bandpass[zeros]=0.

;Polyfit is usually fit per tile, but we want to average together the global bandpass
;find a polyfit per global bandpass (per observation / pol)
for day_i=0, n_days-1 do begin
  for poi_i=0, n_pointings-1 do begin
    n_obs = n_obs_arr[day_i,poi_i]
    for obs_i=0,n_obs-1 do begin
      for pol_i=0,n_pol-1 do begin
        first_tile = first_defined_tile[day_i,poi_i,pol_i,obs_i]
        gain = abs(reform(gain_bandpass[day_i,poi_i,pol_i,freq_use,first_tile,obs_i])) ;all tiles are the same
        ;Only fit for amplitude if amp_degree is set and greater than zero
        ;Fit pre- and post-digital gain jump separately in highband MWA data
        gain_fit=fltarr(n_freq)
        fit_params1=poly_fit(freq_use[0:f_d],gain[0:f_d],amp_degree-1)
        fit_params2=poly_fit(freq_use[f_d+1:f_end],gain[f_d+1:f_end],amp_degree-1)
        FOR di=0L,amp_degree-1 DO gain_fit[freq_use[0]:freq_use[f_d]] += fit_params1[di]*findgen(freq_use[f_d])^di
        FOR di=0L,amp_degree-1 DO gain_fit[freq_use[f_d+1]:freq_use[f_end]] += $
          fit_params2[di]*(findgen(freq_use[f_end] - freq_use[f_d+1]+1) + freq_use[f_d+1])^di

        polyfit_global[day_i,poi_i,pol_i,*,*,obs_i] =  rebin(gain_fit,n_freq,n_tile)
      endfor
    endfor
  endfor
endfor

;Remove the polyfit fit to the bandpass to get just the small-scale global variations in amplitude across the array
gains = abs(gain_bandpass) / polyfit_global
;ave_gains = fltarr(n_pointings,n_pol,n_freq,n_tile,max_n_obs)
ave_gains = fltarr(n_days,n_pointings,n_pol,n_freq,n_tile,max_n_obs)

;Find the average small-scale fluctuations that exist across the entire array 
for poi_i=0, n_pointings-1 do begin
  for pol_i=0,n_pol-1 do begin
    ;Find the first defined tile that is common across the days/obs 
    result = histogram(first_defined_tile[*,poi_i,pol_i,*],binsize=1,min=0)
    temp = max(result,first_tile)
    ;ave across obs
    if n_days GT 1 then ave_gains[*,poi_i,pol_i,*,*,*] = rebin(mean(reform(gains[*,poi_i,pol_i,*,first_tile,*]),dimension=3,/NAN),n_days,n_freq,n_tile,max_n_obs) 
    if n_days EQ 1 then ave_gains[*,poi_i,pol_i,*,*,*] = rebin(mean(reform(gains[*,poi_i,pol_i,*,first_tile,*]),dimension=2,/NAN),n_days,n_freq,n_tile,max_n_obs) 
    ;ave_gains[*,poi_i,pol_i,*,*,*] = rebin(mean(reform(gains[*,poi_i,pol_i,*,first_tile,*]),dimension=3,/NAN),n_days,n_freq,n_tile,max_n_obs) 
    ;full average across days/obs
    ;ave_gains[poi_i,pol_i,*,*,*] = rebin(mean(mean(reform(gains[*,poi_i,pol_i,*,first_tile,*]),dimension=3,/NAN),dimension=1,/NAN),n_freq,n_tile,max_n_obs) 
  endfor
endfor

new_gain_bandpass = polyfit_global * ave_gains
;for day_i=0,n_days-1 do new_gain_bandpass[day_i,*,*,*,*,*] *= ave_gains

;new_gain_bandpass = polyfit_global
;for day_i=0,n_days-1 do new_gain_bandpass[day_i,*,*,*,*,*] *= ave_gains

;Reconstruct the gains using the auto per-tile fluctuations, the previously calculated phases, and the now averaged 
;small scale fluctuations that exist across the entire array.
final_gain = ((zero_mean_individual_amps+1) * abs(polyfit) * new_gain_bandpass)* exp(Complex(0,1)*atan(gain_arr,/phase))

;don't do anything new, but flagged tiles
final_gain = ((zero_mean_individual_amps+1) * abs(polyfit) * gain_bandpass)* exp(Complex(0,1)*atan(gain_arr,/phase))


filepath_transfer = dir + 'cal_transfer_with23/'

for day_i=0,n_days-1 do begin
  for poi_i=0,n_pointings-1 do begin
    obs_ids = *obs_ptr[day_i,poi_i]
    n_obs = N_elements(obs_ids)
    for obs_i=0,n_obs-1 do begin
      restore, dir + cal_dir + obs_ids[obs_i] + '_cal.sav'
      for pol_i=0,n_pol-1 do begin
        *cal.gain[pol_i] = reform(final_gain[day_i,poi_i,pol_i,*,*,obs_i])
      endfor
      save, cal, filename = filepath_transfer + obs_ids[obs_i] + '_cal.sav'
    endfor
  endfor
endfor


end

