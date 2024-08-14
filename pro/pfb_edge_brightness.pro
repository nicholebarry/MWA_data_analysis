pro pfb_edge_brightness


  n_tile=128
  n_freq=384
  n_pol=2

  dir_path = ['/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max100_int/calibration/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max150_int/calibration/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max200_int/calibration/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max250_int/calibration/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_int/calibration/']

  obs_path = ['/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max100_int/metadata/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max150_int/metadata/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max200_int/metadata/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_max250_int/metadata/',$
    '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_int/metadata/']

  n_dir = N_elements(dir_path)

  obs_ids = ['1061317400','1061317520','1061317640','1061317888','1061318128','1061318248','1061318376','1061318496','1061318736','1061318864']
  n_obs = N_Elements(obs_ids)

  rgbcolors=$
      [[228,134,161],$
;      [214,63,91],$
;      [164,79,86],$
;      [200,79,48],$
;      [217,144,103],$
      [219,133,46],$
;      [137,106,42],$
;      [198,167,62],$
;      [153,181,52],$
;      [159,178,106],$
      [82,115,55],$
;      [73,147,61],$
;      [95,196,87],$
;      [81,183,144],$
;      [70,174,206],$
      [100,138,212],$
;      [99,101,214],$
;      [101,87,160],$
;      [186,126,223],$
;      [194,150,217],$
      [164,79,201],$
;      [210,85,181],$
;      [156,79,136],$
      [209,67,131]]
  rgbcolors=reverse(rgbcolors)
  color_num = [10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56]

  for n_i=0, N_elements(rgbcolors[0,*])-1 do $
    TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]
  color_byte = [10B,12B,14B,16B,18B,20B,22B,24B,26B,28B,30B,32B,34B,36B,38B,40B,42B,44B,46B,48B,50B,52B,54B,56B]

  gain_amp = fltarr(n_pol,n_obs,n_dir)
  gain_ind_amp = fltarr(n_pol,n_obs,48,128,n_dir)
  gain_dip = fltarr(n_pol,n_obs,n_dir)
  gain_ind_dip = fltarr(n_pol,n_obs,48,128,n_dir)
  gain_renorm_stddev = fltarr(n_pol,n_obs,48,128,n_dir)
  
  for dir_i=0, n_dir-1 do begin
 
  ; Init the various gain arrays that will be used
  gain_arr = complex(fltarr(n_pol,n_freq,n_tile,n_obs))
  ;gain_bandpass = complex(fltarr(n_pol,n_freq,n_tile,n_obs))
  gain_renorm = fltarr(n_pol,n_freq,n_tile,n_obs)
  gain_renorm2 = fltarr(n_pol,n_freq,n_tile,n_obs)
  gain_mode_fit = complex(fltarr(n_pol,n_freq,n_tile,n_obs))

  ; Observation loop
  for obs_i=0,n_obs-1 do begin

    ; Restore the cal file
    restore, dir_path[dir_i] + obs_ids[obs_i] + '_cal.sav'

    for pol_i=0,1 do begin
      gain_arr[pol_i,*,*,obs_i] = *cal.gain[pol_i]
;      gain_arr[pol_i,*,*,obs_i] = *cal.gain[pol_i]+*cal.gain_residual[pol_i]
      ;gain_bandpass[pol_i,*,*,obs_i] = *cal.gain_bandpass[pol_i]

      for tile_i=0,n_tile-1 do begin
        if ptr_valid(cal.mode_params[pol_i,tile_i]) then begin
          amp_use = (*cal.mode_params[pol_i,tile_i])[1]
          mode_i = (*cal.mode_params[pol_i,tile_i])[0]
          phase_use = (*cal.mode_params[pol_i,tile_i])[2]
          gain_mode_fit[pol_i,*,tile_i,obs_i]=amp_use*exp(-Complex(0,1)*2.*!Pi*(mode_i*findgen(n_freq)/n_freq)+Complex(0,1)*phase_use) + 1
;gain_mode_fit[pol_i,*,tile_i,obs_i] = 1.
        endif else gain_mode_fit[pol_i,*,tile_i,obs_i] = 1.
      endfor

    endfor

  endfor


  ; Get a sample obs file so that the correct freq can be extracted
  obs = getvar_savefile('/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam_diffuse_cal_transfer_int/metadata/' + obs_ids[0] + '_obs.sav','obs')
  freq_use=(*obs.baseline_info).freq_use
  freq_use=where(freq_use,nf_use)
  
  subband_low =  where((findgen(384) mod 16) EQ 0) + 1
  subband_high =  where((findgen(385) mod 16) EQ 0) -2
  ;subband_mean_freq = FLTARR(14,6)

  for obs_i=0, n_obs-1 do begin
    for pol_i=0, n_pol-1 do begin
      gain_fit=fltarr(n_freq,n_tile)
       for subband_i=0, 23 do begin
         for tile_i=0, 127 do begin
           subband_gain = reform(abs(gain_arr[pol_i,subband_low[subband_i]:subband_high[subband_i+1],tile_i,obs_i] / gain_mode_fit[pol_i,subband_low[subband_i]:subband_high[subband_i+1],tile_i,obs_i] ))
           fit_params_sub=poly_fit(INDGEN(12)+1,subband_gain[1:12],1)
           gain_fit[subband_low[subband_i]:subband_high[subband_i+1],tile_i] = fit_params_sub[0]*findgen(14)^0 + fit_params_sub[1]*findgen(14)^1
       gain_ind_amp[pol_i,obs_i,subband_i*2,tile_i,dir_i] = abs(gain_fit[subband_low[subband_i]+1,tile_i])
       gain_ind_amp[pol_i,obs_i,subband_i*2+1,tile_i,dir_i] = abs(gain_fit[subband_low[subband_i]+1,tile_i])
         endfor
       endfor
       restore, obs_path[dir_i] + obs_ids[obs_i] + '_obs.sav'
       tile_use = where((*obs.baseline_info).tile_use EQ 1)
       tile_flag = where((*obs.baseline_info).tile_use EQ 0)
       gain_amp[pol_i,obs_i,dir_i] = mean(abs(gain_arr[pol_i,*,tile_use,obs_i]),/nan)
       ;gain_amp[pol_i,obs_i,dir_i] = mean(gain_fit[*,tile_use],/nan)
       gain_renorm2[pol_i,*,*,obs_i] = reform(abs(gain_arr[pol_i,*,*,obs_i] / gain_mode_fit[pol_i,*,*,obs_i])) / gain_fit 
       gain_renorm[pol_i,*,*,obs_i] = reform(abs(gain_arr[pol_i,*,*,obs_i])) / gain_fit
       gain_renorm[pol_i,*,tile_flag,*]=!Values.F_NAN

      gain_ind_dip[pol_i,obs_i,INDGEN(24)*2,*,dir_i] = reform(gain_renorm[pol_i,subband_low,*,obs_i])
      gain_ind_dip[pol_i,obs_i,INDGEN(24)*2+1,*,dir_i] = reform(gain_renorm[pol_i,subband_high[1:24],*,obs_i])

      ;gain_renorm_stddev[pol_i,obs_i,*,*,dir_i] = rebin(stddev(reform(gain_renorm[pol_i,*,*,obs_i]),dimension=1),48,128)
      gain_renorm_stddev[pol_i,obs_i,*,*,dir_i] = rebin(reform(stddev(reform(gain_renorm[pol_i,*,*,obs_i]),dimension=1,/nan),1,128),48,128)

    endfor
  endfor

  gain_dip[*,*,dir_i] = mean(mean([[gain_renorm[*,subband_low,*,*]],[gain_renorm[*,subband_high[1:24],*,*]]],/nan,dimension=2),dimension=2,/nan)
  gain_dip[*,*,dir_i] = mean(mean([[gain_renorm[*,subband_low,*,*]],[gain_renorm[*,subband_high[1:24],*,*]]],/nan,dimension=2),dimension=2,/nan)
  ;gain_dip[*,*,dir_i] = mean([[gain_renorm[*,subband_low,3,*]],[gain_renorm[*,subband_high[1:24],3,*]]],/nan,dimension=2)

  endfor


cgps_open,dir_path[0] + 'van_vleck.png',/quiet,/nomatch

  cgplot, [reform(mean(gain_amp[1,*,0],/nan))],[reform(mean(gain_dip[1,*,0],/nan))],psym=16,yrange=[.98,1.0],xrange=[1.45,2.0],xtitle='Average amplitude',ytitle='Average Edge Channel Gain (Normalized)'
  ;cgplot, [reform(mean(gain_amp[1,*,0],/nan))],[reform(mean(gain_dip[1,*,0],/nan))],psym=16,yrange=[.98,1.0],xrange=[1.45,1.7],xtitle='Average amplitude',ytitle='Average Edge Channel Gain (Normalized)'
  cgoplot, [reform(mean(gain_amp[1,*,*],/nan,dimension=2))],[reform(mean(gain_dip[1,*,*],/nan,dimension=2))]
  cgoplot, [reform(mean(gain_amp[0,*,*],/nan,dimension=2))],[reform(mean(gain_dip[0,*,*],/nan,dimension=2))]


  for dir_i=0,n_dir-1 do cgoplot, [reform((gain_amp[1,*,dir_i]))],[reform((gain_dip[1,*,dir_i]))],psym=16,color=color_byte[dir_i],sym_size=4
  for dir_i=0,n_dir-1 do cgoplot, [reform((gain_amp[0,*,dir_i]))],[reform((gain_dip[0,*,dir_i]))],psym=16,color=color_byte[dir_i],sym_size=4
  ;for dir_i=0,n_dir-1 do cgoplot, [reform(mean(gain_amp[1,*,dir_i],/nan))],[reform(mean(gain_dip[1,*,dir_i],/nan))],psym=16,color=color_byte[dir_i],sym_size=4
  ;for dir_i=0,n_dir-1 do cgoplot, [reform(mean(gain_amp[0,*,dir_i],/nan))],[reform(mean(gain_dip[0,*,dir_i],/nan))],psym=16,color=color_byte[dir_i],sym_size=4

  cglegend, Title=['<100 $\lambda$ calibration', '<150 $\lambda$ calibration', '<200 $\lambda$ calibration', '<250 $\lambda$ calibration', '<711 $\lambda$ calibration'], color=color_byte[0:4], location = [.155,.85], psym=16, charsize=1.5,Length=0.0

cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage


cgps_open,dir_path[0] + 'all_edges_with_fit_and_errors.png',/quiet,/nomatch

n_tile_use=N_elements(tile_use)
cgplot, reform(gain_ind_amp[*,*,*,tile_use,0],2*10*48.*n_tile_use),reform(gain_ind_dip[*,*,*,tile_use,0],2*10*48.*n_tile_use),psym=12,yrange=[0.5,1.5],charsize=1, ytitle='Edge Channel Gain (Normalized)', xtitle='Fit gain of next channel' 

for dir_i=0, n_dir-1 do cgoplot, reform(gain_ind_amp[*,*,*,tile_use,dir_i],2*10*48.*n_tile_use),reform(gain_ind_dip[*,*,*,tile_use,dir_i],2*10*48.*n_tile_use),psym=12,color=color_byte[dir_i],ERR_YHigh=reform(gain_renorm_stddev[*,*,*,tile_use,dir_i],2*10*48.*n_tile_use),ERR_YLow=reform(gain_renorm_stddev[*,*,*,tile_use,dir_i],2*10*48.*n_tile_use)

temp_x = reform(gain_ind_amp[*,*,*,tile_use,*],2*10*48.*n_tile_use*5)  
temp_y = reform(gain_ind_dip[*,*,*,tile_use,*],2*10*48.*n_tile_use*5)
temp_s = reform(gain_renorm_stddev[*,*,*,tile_use,*],2*10*48.*n_tile_use*5)
inds =where(finite(temp_x))
temp_x = temp_x[inds]
temp_y = temp_y[inds]
temp_s = temp_s[inds]
inds =where(temp_x NE 0)   
temp_x = temp_x[inds]   
temp_y = temp_y[inds]   
temp_s = temp_s[inds]   
;fit_params=poly_fit(temp_x,temp_y,1)    
fit_params=poly_fit(temp_x,temp_y,1,measure_errors=temp_s,sigma=sigma_fit)    
fit_line = fit_params[0]*findgen(7)^0 + fit_params[1]*findgen(7)^1
cgoplot, fit_line, thickness=3

  cglegend, Title=['<100 $\lambda$ calibration', '<150 $\lambda$ calibration', '<200 $\lambda$ calibration', '<250 $\lambda$ calibration', '<711 $\lambda$ calibration'], color=color_byte[0:4], location = [.155,.85], psym=12, charsize=1.5,Length=0.0

cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
stop
;  cgps_open,'/nfs/eor-00/h1/nbarry/pfb/gain_residuals/temp_season_digjump.pdf',/quiet,/nomatch



end
