pro cal_metadata


  ;dir = ['/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_Jack_hbeam_LoBES_extraremoved_ondata_cal_stop/','/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_Jack_hbeam_LoBES_newshaplet_ondata_cal_stop/','/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_Jack_hbeam_GLEAM_newshaplet_ondata_cal_stop/']
  dir = ['/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_calstop/']
  n_filepaths = N_elements(dir)
  cal_dir = 'calibration/'
  obs_dir = 'metadata/'

  ; Parameters for the observation id files
  ;dir_obsids = '/home/nbarry/MWA/pipeline_scripts/bash_scripts/ozstar/obs_list/'
  dir_obsids = '/home/nbarry/MWA/pipeline_scripts/bash_scripts/ozstar/obs_list/2014/'
  day = '10'
  parsednames=[day+'_minustwo_ssins',day+'_minusone_ssins',day+'_zenith_ssins',day+'_plusone_ssins',day+'_plustwo_ssins']
;  parsednames=[day+'_zenith']
  n_pointings = N_elements(parsednames)
  obs_ptr=PTRARR(n_pointings,/allocate_heap)

  max_n_obs=0
  ; Get observation ids and put them in a pointing pointer
  FOR poi_i=0,n_pointings-1 DO BEGIN
    readcol, dir_obsids + parsednames[poi_i]+'.txt', obs_temp, format='A', /silent
    *obs_ptr[poi_i]=obs_temp
    max_n_obs = max([N_elements(obs_temp),max_n_obs])
  ENDFOR

  ; Get a sample obs file so that the polyfit params can be extracted
  obs = getvar_savefile(dir[0] + obs_dir + (*obs_ptr[0])[0] + '_obs.sav','obs')
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



;gain_arr = complex(fltarr(n_filepaths,n_pol,n_freq,128,n_pointings,max_n_obs))
;gain_res = complex(fltarr(n_filepaths,n_pol,n_freq,128,n_pointings,max_n_obs))
;gain_bandpass = complex(fltarr(n_filepaths,n_pol,n_freq,128,n_pointings,max_n_obs))
n_conv = lonarr(n_filepaths,n_pol,n_pointings,max_n_obs)
conv_iter = fltarr(n_filepaths,n_pol,n_freq,128,n_pointings,max_n_obs)

;polyfit = complex(fltarr(3,2,n_freq,128,n_obs))
;polyfit_global = fltarr(3,2,n_freq,128,n_obs)

;cal_save_point=1
;if ~keyword_set(cal_save_point) then begin

for filepath_i=0, n_filepaths-1 do begin
for poi_i=0, n_pointings-1 do begin

  obs_ids = *obs_ptr[poi_i]
  n_obs = N_elements(obs_ids)

  for obs_i=0,n_obs-1 do begin


    restore, dir[filepath_i] + cal_dir + obs_ids[obs_i] + '_cal.sav'
    for pol_i=0,1 do begin
;      gain_arr[filepath_i,pol_i,*,*,obs_i] = *cal.gain[pol_i]
;      gain_res[filepath_i,pol_i,*,*,obs_i] = *cal.gain_residual[pol_i]
      n_conv[filepath_i,pol_i,poi_i,obs_i] = cal.n_converged[pol_i]
      conv_iter[filepath_i,pol_i,*,*,poi_i,obs_i] = *cal.conv_iter[pol_i]

    endfor


  endfor
endfor
endfor


stop

end
