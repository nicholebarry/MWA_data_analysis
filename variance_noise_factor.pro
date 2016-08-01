pro variance_noise_factor, pol=pol
;Procedure to make frequency-averaged and night-averaged apparent frame data and variance cubes from pre-run data. 
;Here, the apparent frame is defined as data / weights and variance / weights^2.

  if ~keyword_set(pol) then pol='XX'
  
  ;vis_sigma = sum(sum(vis_noise(freq, obs)*n_vis_freq(freq, obs))/sum(n_freq_vis(freq, obs))
  fhd_path='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_March2016_small_through_firstpass_largeind'
  obs_name='obs_id_6296'
  
  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/'+obs_name+'.txt'
  readcol, filename, obsids, format='A', /silent
  
  vis_noise=FLTARR(N_elements(obsids),384)
  nf_vis=FLTARR(N_elements(obsids),384)
  
  ;Find the calculated noise for each obsid
  for obs_i=0, N_elements(obsids)-1 do begin
    obs = getvar_savefile(fhd_path+'/metadata/'+strtrim(string(obsids[obs_i]),2)+'_obs.sav','obs')
    vis_noise[obs_i,*] = (*obs.vis_noise)[0,*] ;grab XX for now
    nf_vis[obs_i,*] = obs.nf_vis
  endfor
  
  vis_sigma = total(vis_noise*nf_vis)/total(nf_vis)
  
  ;Grab the odd and even residual uvf data cubes and add them together. Sum them in frequency
  data_odd = getvar_savefile(fhd_path+'/ps/Combined_obs_'+obs_name+'_odd_cube'+pol+'_res_uvf.idlsave','data_cube')
  data_even = getvar_savefile(fhd_path+'/ps/Combined_obs_'+obs_name+'_even_cube'+pol+'_res_uvf.idlsave','data_cube')
  data_evenodd = data_odd + data_even
  data_cube = total(data_evenodd,3)
  
  ;Grab the odd and even variance uvf data cubes and add them together. Sum them in frequency
  variance_odd = getvar_savefile(fhd_path+'/ps/Combined_obs_'+obs_name+'_odd_cube'+pol+'_weights_uvf.idlsave','variance_cube')
  variance_even = getvar_savefile(fhd_path+'/ps/Combined_obs_'+obs_name+'_even_cube'+pol+'_weights_uvf.idlsave','variance_cube')
  variance_evenodd = variance_odd + variance_even
  variance_cube = total(variance_evenodd,3)
  
    ;Grab the odd and even weights uvf data cubes and add them together. Sum them in frequency
    ;Note: there is a bit of a choice when it comes to squaring the weights for creating the apparent variances.
    ;Here, we have chosen to square then sum rather than sum then square. This is wholly dependent on your choice
    ;of variance definition and the desired result.
  weights_odd = getvar_savefile(fhd_path+'/ps/Combined_obs_'+obs_name+'_odd_cube'+pol+'_weights_uvf.idlsave','weights_cube')
  weights_even = getvar_savefile(fhd_path+'/ps/Combined_obs_'+obs_name+'_even_cube'+pol+'_weights_uvf.idlsave','weights_cube')
  weights_evenodd = weights_odd + weights_even
  weights_cube = total(weights_evenodd,3)
  weights_cube_squared = total(abs(weights_evenodd^2.),3)
  
  xarr = getvar_savefile(fhd_path+'/ps/Combined_obs_'+obs_name+'_odd_cube'+pol+'_weights_uv_plane.idlsave','xarr')
  yarr = getvar_savefile(fhd_path+'/ps/Combined_obs_'+obs_name+'_odd_cube'+pol+'_weights_uv_plane.idlsave','yarr')
  kperp_lambda_conv = getvar_savefile(fhd_path+'/ps/Combined_obs_'+obs_name+'_odd_cube'+pol+'_weights_uv_plane.idlsave','kperp_lambda_conv')
  y_arr = kperp_lambda_conv * yarr
  x_arr = kperp_lambda_conv * xarr
  
  if pol EQ 'XX' then begin
    apparent_data_XX = data_cube/weights_cube
    apparent_variance_XX = variance_cube / weights_cube_squared
    quick_image, abs(apparent_data_XX), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Magnitude of frequency weighted res in uv plane (xx) in Jy', $
      charsize=1, window=1, data_range = [0,1E1], savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/uv_weighted_data_XX_bh_large',/png
    quick_image, abs(apparent_variance_XX), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Magnitude of frequency weighted variance in uv plane (xx) in Jy', $
      charsize=1, window=2,savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/uv_weighted_variance_XX_bh_large',/png
    save, apparent_data_XX, x_arr, y_arr, apparent_variance_XX, filename = '/nfs/mwa-00/h1/nbarry/UVfromFHD/uv_weighted_variance_XX_bh_large.sav'
    stop
  endif else begin
    apparent_data_YY = data_cube/weights_cube
    apparent_variance_YY = variance_cube / weights_cube_squared
    quick_image, abs(apparent_data_YY), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Magnitude of frequency weighted res in uv plane (yy) in Jy', $
      charsize=1, window=1,data_range = [1E0,1E1],savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/uv_weighted_data_YY_bh',/png
    quick_image, abs(apparent_variance_YY), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Magnitude of frequency weighted variance in uv plane (yy) in Jy', $
      charsize=1, window=2,savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/uv_weighted_variance_YY_bh',/png
    save, apparent_data_YY, x_arr, y_arr, apparent_variance_YY, filename = '/nfs/mwa-00/h1/nbarry/UVfromFHD/uv_weighted_variance_YY_bh.sav'
    stop
  endelse
  
  n_kx = N_elements(x_arr)
  n_ky = N_elements(y_arr)
  n_freq = 192
  
  ;; get sigma^2 into Jy^2
  sigma2_cube1_factor = rebin(reform(vis_sigma, 1, 1, n_freq), n_kx, n_ky, n_freq)^2. ; where sigma^2 is "true" variance
  
  stop
  
end