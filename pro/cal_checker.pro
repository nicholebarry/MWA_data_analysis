pro cal_checker

cal_filepath = ['/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_windowTx2_cal_stop/calibration/','/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_haslam_windowTx2_cal_stop/calibration/'];,'/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_cal_stop/calibration/']
obs_filepath = ['/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_windowTx2_cal_stop/metadata/','/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_haslam_windowTx2_cal_stop/metadata/'];,'/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_cal_stop/metadata/']
n_filepaths = N_elements(cal_filepath)

obs = getvar_savefile('/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_stop/metadata/1061317400_obs.sav','obs')
freq_use=(*obs.baseline_info).freq_use
freq_use=where(freq_use,nf_use)
pre_dig_inds = where((*obs.baseline_info).freq[freq_use] LT 187.515E6,n_count)
f_d = max(pre_dig_inds)
f_end = N_elements(freq_use)-1
amp_degree = 2
phase_degree = 1

obs_ids = ['1061317400','1061317520','1061317640','1061317888','1061318128','1061318248','1061318376','1061318496','1061318736','1061318864']
n_obs = N_elements(obs_ids)
n_freq = 384

gain_arr = complex(fltarr(3,2,n_freq,128,n_obs))
gain_res = complex(fltarr(3,2,n_freq,128,n_obs))
gain_bandpass = complex(fltarr(3,2,n_freq,128,n_obs))
n_conv = lonarr(3,2,n_obs)
conv_iter = fltarr(3,2,n_freq,128,n_obs)

polyfit = complex(fltarr(3,2,n_freq,128,n_obs))
polyfit_global = fltarr(3,2,n_freq,128,n_obs)

;cal_save_point=1
if ~keyword_set(cal_save_point) then begin

for filepath_i=0, n_filepaths-1 do begin
  for obs_i=0,n_obs-1 do begin


    restore, cal_filepath[filepath_i] + obs_ids[obs_i] + '_cal.sav'
    for pol_i=0,1 do begin
      gain_arr[filepath_i,pol_i,*,*,obs_i] = *cal.gain[pol_i]
      gain_res[filepath_i,pol_i,*,*,obs_i] = *cal.gain_residual[pol_i]
      n_conv[filepath_i,pol_i,obs_i] = cal.n_converged[pol_i]
      conv_iter[filepath_i,pol_i,*,*,obs_i] = *cal.conv_iter[pol_i]

      ;Reconstruct polyfit
      for tile_i=0, 127 do begin
        gain_fit=fltarr(n_freq)
        fit_params1= (*cal.amp_params[pol_i,tile_i])[0,*];poly_fit(freq_use[0:f_d],gain[0:f_d],amp_degree-1)
        fit_params2= (*cal.amp_params[pol_i,tile_i])[1,*];poly_fit(freq_use[f_d+1:f_end],gain[f_d+1:f_end],amp_degree-1)
        FOR di=0L,amp_degree-1 DO gain_fit[freq_use[0]:freq_use[f_d]] += fit_params1[di]*findgen(freq_use[f_d])^di
        FOR di=0L,amp_degree-1 DO gain_fit[freq_use[f_d+1]:freq_use[f_end]] += $
          fit_params2[di]*(findgen(freq_use[f_end] - freq_use[f_d+1]+1) + freq_use[f_d+1])^di

        phase_fit=fltarr(n_freq)
        FOR di=0L,phase_degree DO phase_fit+=(*cal.phase_params[pol_i,tile_i])[0,di]*findgen(n_freq)^di
        
        polyfit[filepath_i,pol_i,*,tile_i,obs_i]=gain_fit*Exp(Complex(0,1)*phase_fit)
      endfor

    *cal.gain[pol_i] = *cal.gain[pol_i] + *cal.gain_residual[pol_i]

    endfor


;    uvfits_read,hdr,params,layout,vis_arr,vis_weights,file_path_vis='/fred/oz048/MWA/data/2013/van_vleck_corrected/' + obs_ids[obs_i] +'.uvfits',$
;      n_pol=n_pol,silent=silent,error=error,_Extra=extra
;    restore, obs_filepath[filepath_i] +  obs_ids[obs_i] + '_obs.sav'
;    vis_auto=vis_extract_autocorr(obs,vis_arr = vis_arr,/time_average,auto_tile_i=auto_tile_i)
;    cal_autos_init=cal_auto_ratio(obs,cal,auto_ratio=auto_ratio,vis_auto=vis_auto,auto_tile_i=auto_tile_i,/divide)
;    cal_bandpass=vis_cal_bandpass(cal_autos_init,obs,params,cal_remainder=cal_remainder,auto_ratio_calibration=1)

;    for pol_i=0,1 do begin
;      gain_bandpass[filepath_i,pol_i,*,*,obs_i] = *cal_bandpass.gain[pol_i]
;    endfor

  endfor
endfor
stop
;gain = gain_res + gain_arr




zeros = where(gain_arr EQ 0)
gain_arr[zeros]=!Values.F_NAN
amps = abs(gain_arr / polyfit / gain_bandpass) - 1

endif else restore, 'cal_save_point.sav'

zeros = where(~finite(gain_bandpass))
gain_bandpass[zeros]=0.


;Polyfit is usually fit per tile, but we want to average together the global bandpass
;find a polyfit per global bandpass (per observation / pol)
for filepath_i=0, n_filepaths-1 do begin
  for obs_i=0,n_obs-1 do begin
    for pol_i=0,1 do begin
        gain = abs(reform(gain_bandpass[filepath_i,pol_i,freq_use,0,obs_i])) ;all tiles are the same
        ;Only fit for amplitude if amp_degree is set and greater than zero
        ;Fit pre- and post-digital gain jump separately in highband MWA data
          gain_fit=fltarr(n_freq)
          fit_params1=poly_fit(freq_use[0:f_d],gain[0:f_d],amp_degree-1)
          fit_params2=poly_fit(freq_use[f_d+1:f_end],gain[f_d+1:f_end],amp_degree-1)
          FOR di=0L,amp_degree-1 DO gain_fit[freq_use[0]:freq_use[f_d]] += fit_params1[di]*findgen(freq_use[f_d])^di
          FOR di=0L,amp_degree-1 DO gain_fit[freq_use[f_d+1]:freq_use[f_end]] += $
            fit_params2[di]*(findgen(freq_use[f_end] - freq_use[f_d+1]+1) + freq_use[f_d+1])^di

          polyfit_global[filepath_i,pol_i,*,*,obs_i] =  rebin(gain_fit,n_freq,128)

    endfor
  endfor
endfor

gains = abs(gain_bandpass) / polyfit_global
ave_gains = fltarr(3,2,n_freq,128,n_obs)

for filepath_i=0, n_filepaths-1 do begin
  for pol_i=0,1 do begin
    ave_gains[filepath_i,pol_i,*,*,*] = rebin(mean(reform(gains[filepath_i,pol_i,*,0,*]),dimension=2,/NAN),n_freq,128,n_obs)
  endfor
endfor

new_gain_bandpass = ave_gains * polyfit_global

gains = ((amps+1) * abs(polyfit) * new_gain_bandpass)* exp(Complex(0,1)*atan(gain_arr,/phase))

;only save the diffuse files
filepath_transfer = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_stop/cal_transfer/'

for obs_i=0,n_obs-1 do begin
  restore, cal_filepath[0] + obs_ids[obs_i] + '_cal.sav'
  for pol_i=0,1 do begin
    *cal.gain[pol_i] = reform(gains[0,pol_i,*,*,obs_i])
  endfor
  save, cal, filename = filepath_transfer + obs_ids[obs_i] + '_cal.sav'
endfor

stop

end
