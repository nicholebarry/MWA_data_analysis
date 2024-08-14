pro variance_noise_factor, pol=pol,deluv=deluv
  ;Procedure to make frequency-averaged and night-averaged apparent frame data and variance cubes from pre-run data.
  ;Here, the apparent frame is defined as data / weights and variance / weights^2.

  if ~keyword_set(pol) then pol='XX'
  
  if ~keyword_set(deluv) then deluv = '' else deluv = '_deluv' + deluv
  
  ;vis_sigma = sum(sum(vis_noise(freq, obs)*n_vis_freq(freq, obs))/sum(n_freq_vis(freq, obs))
  ;fhd_path='/nfs/mwa-05/r1/EoRuvfits/EoR2013/fhd_nb_2013longrun_savedbp'
  fhd_path='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun'
  ps_dir='ps'
  obs_name='beardsley_thesis_list';_forlarge'
  ;obs_name='Aug23zenith_thesis
  filter_name='_tk'
  ;filter_name=''
  ;obs_name='Oct23'
  
  metadata_struct = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_cube'+pol+'__even_odd_joint_info.idlsave','metadata_struct')
  freq = metadata_struct.frequencies
  
  z0_freq = 1420.40;E6 ;; Hz
  ;freq = (*obs.baseline_info).freq
  redshifts = z0_freq/freq - 1 ;; frequencies will be identical if kx, ky, kz match
  
  cosmology_measures, redshifts, wedge_factor = wedge_factor, comoving_dist_los=comov_dist_los, Ez=Ez
  
  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  
  z_mpc_delta = float(mean(comov_los_diff))
  z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
  kz_mpc_range =  (2.*!pi) / (z_mpc_delta)
  kz_mpc_delta = (2.*!pi) / z_mpc_length
  kz = findgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - kz_mpc_range/2.
  folded_kz = abs(kz)
  

  
  
  ;Grab the odd and even residual uvf data cubes and add them together. Sum them in frequency
  ;data_odd = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_odd_cube'+pol+deluv+filter_name + '_res_uvf.idlsave','data_cube')
  ;data_even = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_even_cube'+pol+deluv +filter_name+'_res_uvf.idlsave','data_cube')
  ;data_model_odd = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_odd_cube'+pol+deluv+filter_name + '_model_uvf.idlsave','data_cube')
  ;data_model_even = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_even_cube'+pol+deluv +filter_name+'_model_uvf.idlsave','data_cube')
  data_odd = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_odd_cube'+pol+deluv+filter_name + '_res_uvf.idlsave','data_cube')
  data_even = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_even_cube'+pol+deluv +filter_name+'_res_uvf.idlsave','data_cube')
  data_evenodd = data_odd + data_even
  ;data_evenodd = (data_dirty_odd -data_model_odd)+ (data_dirty_even -data_model_even)
  data_cube = total(data_evenodd,3)
  
  ;Grab the odd and even variance uvf data cubes and add them together. Sum them in frequency
  variance_odd = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_odd_cube'+pol+deluv +filter_name+'_weights_uvf.idlsave','variance_cube')
  variance_even = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_even_cube'+pol+deluv +filter_name+'_weights_uvf.idlsave','variance_cube')
  variance_evenodd = variance_odd + variance_even
  variance_cube = total(variance_evenodd,3)
  
  ;Grab the odd and even weights uvf data cubes and add them together. Sum them in frequency
  ;Note: there is a bit of a choice when it comes to squaring the weights for creating the apparent variances.
  ;Here, we have chosen to square then sum rather than sum then square. This is wholly dependent on your choice
  ;of variance definition and the desired result.
  weights_odd = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_odd_cube'+pol+deluv +filter_name+'_weights_uvf.idlsave','weights_cube')
  weights_even = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_even_cube'+pol+deluv +filter_name+'_weights_uvf.idlsave','weights_cube')
  weights_evenodd = weights_odd + weights_even
  weights_cube = total(weights_evenodd,3)
  weights_cube_squared = total(abs(weights_evenodd^2.),3)
  
  xarr = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_odd_cube'+pol+deluv +filter_name+'_weights_uv_plane.idlsave','xarr')
  yarr = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_odd_cube'+pol+deluv +filter_name+'_weights_uv_plane.idlsave','yarr')
  kperp_lambda_conv = getvar_savefile(fhd_path+'/'+ps_dir+'/Combined_obs_'+obs_name+'_odd_cube'+pol+deluv +filter_name+'_weights_uv_plane.idlsave','kperp_lambda_conv')
  y_arr = kperp_lambda_conv * yarr
  x_arr = kperp_lambda_conv * xarr
  
  
     mean_comoving_dist_los = mean(comov_dist_los)
    
    kx = 2. * !PI * x_arr / mean_comoving_dist_los
    kx = kx[where(kx GE 0)]
    ky = 2. * !PI * y_arr / mean_comoving_dist_los
    kperp = sqrt(abs(kx)^2. + abs(ky)^2.)
    
        ;*****calculate wedge
    ;; assume 20 degrees from pointing center to first null
    source_dist = 20d * !dpi / 180d
    fov_amp = wedge_factor * source_dist
    max_theta= 10d ;random guess!
    
    ;; calculate angular distance to horizon
    horizon_amp = wedge_factor * ((max_theta+90d) * !dpi / 180d)
    
    wedge_amp = [[fov_amp], [horizon_amp]] ; to be multiplied by k perp
    wedge_names = ['fov', 'horizon']
    
        ;Filtering
    kperp_cut = kperp#transpose(wedge_amp[*,1])
  
  
  kz_data_cube = fft(data_evenodd*weight_invert(weights_evenodd), dimension=3) * N_elements(freq) * z_mpc_delta
  
  ;; put k0 in middle of cube
  kz_data_cube = shift(kz_data_cube, [0,0,N_elements(kz)/2])
  
  n_val = round(kz / kz_mpc_delta)
  kz[where(n_val eq 0)] = 0
  
  kz_mpc = kz[where(n_val ge 0)]
  n_kz = n_elements(kz_mpc)
  
  ;; these an and bn calculations don't match the standard
  ;; convention (they ares down by a factor of 2) but they make more sense
  ;; and remove factors of 2 we'd otherwise have in the power
  ;; and variance calculations
  ;; note that the 0th mode will have higher noise because there's half as many measurements going into it
  a1_0 = kz_data_cube[*,*,where(n_val eq 0)]
  a1_n = (kz_data_cube[*,*, where(n_val gt 0)] + kz_data_cube[*,*, reverse(where(n_val lt 0))])/2.
  b1_n = complex(0,1) * (kz_data_cube[*,*, where(n_val gt 0)] - kz_data_cube[*,*, reverse(where(n_val lt 0))])/2.
  
  n_u = (size(kz_data_cube))[1]
  n_v = (size(kz_data_cube))[2]
  data_sum_1 = complex(fltarr(n_u, n_v, n_kz))
  data_sum_1[*, *, 0] = a1_0
  data_sum_1[*, *, 1:n_kz-1] = a1_n + b1_n
  
  for kz_i=0, 95 do quick_image, alog(abs(data_sum_1[*,*,kz_i])), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', charsize=1, cb_title='Log(Jy)', data_range=[-2,10],$
    title='Magnitude of longrun residual (xx), kz ' + strtrim(kz_i,2), savefile='/nfs/mwa-00/h1/nbarry/UVfromFHD/kz/res_data_longrun_'+string(strtrim(kz_i,2),FORMAT='(I03)'),/png
    
  if pol EQ 'XX' then begin
    apparent_data_XX = data_cube/weights_cube
    apparent_variance_XX = variance_cube / weights_cube_squared
    quick_image, alog10(abs(apparent_data_XX)), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted model in uv plane (xx) in log(Jy), longrun savedbp', $
      charsize=1, window=1, data_range = [0,2],cb_title='Log(Jy)', savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_model_XX_longrun_savedpb_tk_log',/png
    ;quick_image, abs(apparent_variance_XX), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted variance in uv plane (xx) in Jy, longrun', $
    ;  charsize=1, window=2;,savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_variance_XX_longrun_zenith_tk',/png
      
    mid = ((size(apparent_data_XX))[1]-1)/2+1
    quick_image, alog10(abs(apparent_data_XX[mid-50:mid+50,0:50])), x_arr[mid-50:mid+50], y_arr[0:50], xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted model in uv plane (xx) in log(Jy), longrun savedbp', $
      charsize=1, window=3, data_range = [0,2],cb_title='Log(Jy)', savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_model_XX_longrun_savedbp_zoom_tk_log',/png
    ;quick_image, abs(apparent_variance_XX[mid-50:mid+50,0:50]), x_arr[mid-50:mid+50], y_arr[0:50], xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted variance in uv plane (xx) in Jy, longrun zoom', $
    ;  charsize=1, window=4;,savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_variance_XX_longrun_zenith_tk',/png
    stop
    ;save, apparent_data_XX, x_arr, y_arr, apparent_variance_XX, filename = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_datavariance_XX_fullnight_none.sav'
    stop
  endif else begin
    apparent_data_YY = data_cube/weights_cube
    apparent_variance_YY = variance_cube / weights_cube_squared
    quick_image, abs(apparent_data_YY), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted res in uv plane (yy) in Jy, longrun zenith', $
      charsize=1, window=1,data_range = [0,1E1],savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_data_YY_longrun_zenith_tk',/png
    quick_image, abs(apparent_variance_YY), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted variance in uv plane (yy) in Jy, longrun zenith', $
      charsize=1, window=2,savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_variance_YY_longrun_zenith_tk',/png
    ;save, apparent_data_YY, x_arr, y_arr, apparent_variance_YY, filename = '/nfs/mwa-00/h1/nbarry/UVfromFHD/uv_weighted_variance_YY_bh.sav'
      
    mid = ((size(apparent_data_YY))[1]-1)/2+1
    quick_image, abs(apparent_data_YY[mid-50:mid+50,0:50]), x_arr[mid-50:mid+50], y_arr[0:50], xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted res in uv plane (yy) in Jy, longrun zoom', $
      charsize=1, window=3, data_range = [0,1E2];, savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_data_YY_longrun_zoom_tk',/png
    quick_image, abs(apparent_variance_YY[mid-50:mid+50,0:50]), x_arr[mid-50:mid+50], y_arr[0:50], xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted variance in uv plane (yy) in Jy, longrun zoom', $
      charsize=1, window=4;,savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_variance_XX_longrun_zenith_tk',/png
    stop
    stop
  endelse
  
  n_kx = N_elements(x_arr)
  n_ky = N_elements(y_arr)
  n_freq = 192
  
  ;; get sigma^2 into Jy^2
  sigma2_cube1_factor = rebin(reform(vis_sigma, 1, 1, n_freq), n_kx, n_ky, n_freq)^2. ; where sigma^2 is "true" variance
  
  stop
  
end