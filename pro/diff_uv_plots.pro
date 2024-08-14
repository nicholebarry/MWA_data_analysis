pro diff_UV_plots, folder_names, obs_names_in, ps_foldernames=ps_foldernames, pol=pol,deluv=deluv
  ;Procedure to make frequency-averaged and night-averaged apparent frame data and variance cubes from pre-run data.
  ;Here, the apparent frame is defined as data / weights and variance / weights^2.

  if ~keyword_set(pol) then pol='XX'
  
  if ~keyword_set(deluv) then deluv = '' else deluv = '_deluv' + deluv
  
  if N_elements(folder_names) LT 2 then message, 'Must have two folder_names to difference'
  folder_names = mit_folder_locs(folder_names)
  
  if N_elements(ps_foldernames) EQ 0 then ps_foldernames=['ps','ps'] else if N_elements(ps_foldernames) EQ 1 then ps_foldernames=[ps_foldernames,ps_foldernames]
  if N_elements(obs_names_in) EQ 1 then obs_names_in=[obs_names_in,obs_names_in]
  
  ;filter_name='_tk'
  filter_name=''
  
  data_cube = PTRARR(2,/allocate)
  variance_cube = PTRARR(2,/allocate)
  weights_cube = PTRARR(2,/allocate)
  weights_cube_squared = PTRARR(2,/allocate)
  
  for cube_i=0,1 do begin
    data_odd = getvar_savefile(folder_names[cube_i]+'/'+ps_foldernames[cube_i]+'/Combined_obs_'+obs_names_in[cube_i]+'_odd_cube'+pol+deluv+filter_name + '_res_uvf.idlsave','data_cube')
    data_even = getvar_savefile(folder_names[cube_i]+'/'+ps_foldernames[cube_i]+'/Combined_obs_'+obs_names_in[cube_i]+'_even_cube'+pol+deluv +filter_name+'_res_uvf.idlsave','data_cube')
    data_evenodd = data_odd + data_even
    *data_cube[cube_i] = total(data_evenodd,3)
    
    ;Grab the odd and even variance uvf data cubes and add them together. Sum them in frequency
    variance_odd = getvar_savefile(folder_names[cube_i]+'/'+ps_foldernames[cube_i]+'/Combined_obs_'+obs_names_in[cube_i]+'_odd_cube'+pol+deluv +filter_name+'_weights_uvf.idlsave','variance_cube')
    variance_even = getvar_savefile(folder_names[cube_i]+'/'+ps_foldernames[cube_i]+'/Combined_obs_'+obs_names_in[cube_i]+'_even_cube'+pol+deluv +filter_name+'_weights_uvf.idlsave','variance_cube')
    variance_evenodd = variance_odd + variance_even
    *variance_cube[cube_i] = total(variance_evenodd,3)
    
    ;Grab the odd and even weights uvf data cubes and add them together. Sum them in frequency
    ;Note: there is a bit of a choice when it comes to squaring the weights for creating the apparent variances.
    ;Here, we have chosen to square then sum rather than sum then square. This is wholly dependent on your choice
    ;of variance definition and the desired result.
    weights_odd = getvar_savefile(folder_names[cube_i]+'/'+ps_foldernames[cube_i]+'/Combined_obs_'+obs_names_in[cube_i]+'_odd_cube'+pol+deluv +filter_name+'_weights_uvf.idlsave','weights_cube')
    weights_even = getvar_savefile(folder_names[cube_i]+'/'+ps_foldernames[cube_i]+'/Combined_obs_'+obs_names_in[cube_i]+'_even_cube'+pol+deluv +filter_name+'_weights_uvf.idlsave','weights_cube')
    weights_evenodd = weights_odd + weights_even
    *weights_cube[cube_i] = total(weights_evenodd,3)
    *weights_cube_squared[cube_i] = total(abs(weights_evenodd^2.),3)
    
    xarr = getvar_savefile(folder_names[cube_i]+'/'+ps_foldernames[cube_i]+'/Combined_obs_'+obs_names_in[cube_i]+'_odd_cube'+pol+deluv +filter_name+'_weights_uv_plane.idlsave','xarr')
    yarr = getvar_savefile(folder_names[cube_i]+'/'+ps_foldernames[cube_i]+'/Combined_obs_'+obs_names_in[cube_i]+'_odd_cube'+pol+deluv +filter_name+'_weights_uv_plane.idlsave','yarr')
    kperp_lambda_conv = getvar_savefile(folder_names[cube_i]+'/'+ps_foldernames[cube_i]+'/Combined_obs_'+obs_names_in[cube_i]+'_odd_cube'+pol+deluv +filter_name+'_weights_uv_plane.idlsave','kperp_lambda_conv')
    y_arr = kperp_lambda_conv * yarr
    x_arr = kperp_lambda_conv * xarr
  endfor
  
  data_diff = alog10(abs(*data_cube[0]/*weights_cube[0])) - alog10(abs(*data_cube[1]/*weights_cube[1]))
  data_ratio = alog10(abs(*data_cube[0]/*weights_cube[0])) * weight_invert(alog10(abs(*data_cube[1]/*weights_cube[1])))
  
  data_diff_rebin = rebin(data_diff,(size(data_diff))[1] * (size(data_diff))[2])
  resistant_mean, data_diff_rebin, 8, res_mean, goodvec=goodvec
  data_diff_range = minmax(data_diff_rebin[goodvec])
  
  data_ratio_rebin = rebin(data_ratio,(size(data_ratio))[1] * (size(data_ratio))[2])
  resistant_mean, data_ratio_rebin, 6, res_mean, goodvec=goodvec
  data_ratio_range = minmax(data_ratio_rebin[goodvec])
  
  note = file_basename(folder_names[0],'') + ' - ' + file_basename(folder_names[1],'') + ', ' + obs_names_in[0] + ' - ' + obs_names_in[1]
  stop
  quick_image, ((data_diff)), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag diff of weighted res in uv plane (xx) in log(Jy)', note=note, $
    charsize=1, window=1,cb_title='Log(Jy)',savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_res_XX_longrun_difftest_tk_log',/png, data_range=data_diff_range
  note = file_basename(folder_names[0],'') + ' \ ' + file_basename(folder_names[1],'') + ', ' + obs_names_in[0] + ' \ ' + obs_names_in[1]
  quick_image, ((data_ratio)), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag ratio of weighted res in uv plane (xx) in log(Jy)', note=note,$
    charsize=1, window=1,cb_title='Log(Jy)',savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_res_XX_longrun_ratiotest_tk_log',/png,data_range=data_ratio_range
    
  if keyword_set(zoom) then begin
    mid = ((size(data_diff))[1]-1)/2+1
    
    data_diff_rebin = rebin(data_diff[mid-50:mid+50,0:50],(size(data_diff[mid-50:mid+50,0:50]))[1] * (size(data_diff[mid-50:mid+50,0:50]))[2])
    resistant_mean, data_diff_rebin, 8, res_mean, goodvec=goodvec
    data_diff_range = minmax(data_diff_rebin[goodvec])
    
    data_ratio_rebin = rebin(data_ratio[mid-50:mid+50,0:50],(size(data_ratio[mid-50:mid+50,0:50]))[1] * (size(data_ratio[mid-50:mid+50,0:50]))[2])
    resistant_mean, data_ratio_rebin, 6, res_mean, goodvec=goodvec
    data_ratio_range = minmax(data_ratio_rebin[goodvec])
    
    note = file_basename(folder_names[0],'') + ' - ' + file_basename(folder_names[1],'') + ', ' + obs_names_in[0] + ' - ' + obs_names_in[1]
    quick_image, ((data_diff[mid-50:mid+50,0:50])), x_arr[mid-50:mid+50], y_arr[0:50], xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag diff of weighted res in uv plane (xx) in log(Jy)', $
      charsize=1, window=3, cb_title='Log(Jy)', savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_res_XX_longrun_difftest_zoom_tk_log',note=note,/png,data_range=data_diff_range
    note = file_basename(folder_names[0],'') + ' \ ' + file_basename(folder_names[1],'') + ', ' + obs_names_in[0] + ' \ ' + obs_names_in[1]
    quick_image, ((data_ratio[mid-50:mid+50,0:50])), x_arr[mid-50:mid+50], y_arr[0:50], xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag ratio of weighted res in uv plane (xx) in log(Jy)', $
      charsize=1, window=3, cb_title='Log(Jy)', savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_res_XX_longrun_ratiotest_zoom_tk_log',data_range=data_ratio_range,note=note,/png
  endif
  
;  if pol EQ 'XX' then begin
;    apparent_data_XX = data_cube/weights_cube
;    apparent_variance_XX = variance_cube / weights_cube_squared
;    quick_image, alog10(abs(apparent_data_XX)), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted model in uv plane (xx) in log(Jy), longrun savedbp', $
;      charsize=1, window=1, data_range = [0,2],cb_title='Log(Jy)', savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_model_XX_longrun_savedpb_tk_log',/png
;    ;quick_image, abs(apparent_variance_XX), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted variance in uv plane (xx) in Jy, longrun', $
;    ;  charsize=1, window=2;,savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_variance_XX_longrun_zenith_tk',/png
;
;    mid = ((size(apparent_data_XX))[1]-1)/2+1
;    quick_image, alog10(abs(apparent_data_XX[mid-50:mid+50,0:50])), x_arr[mid-50:mid+50], y_arr[0:50], xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted model in uv plane (xx) in log(Jy), longrun savedbp', $
;      charsize=1, window=3, data_range = [0,2],cb_title='Log(Jy)', savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_model_XX_longrun_savedbp_zoom_tk_log',/png
;    ;quick_image, abs(apparent_variance_XX[mid-50:mid+50,0:50]), x_arr[mid-50:mid+50], y_arr[0:50], xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted variance in uv plane (xx) in Jy, longrun zoom', $
;    ;  charsize=1, window=4;,savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_variance_XX_longrun_zenith_tk',/png
;    stop
;    ;save, apparent_data_XX, x_arr, y_arr, apparent_variance_XX, filename = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_datavariance_XX_fullnight_none.sav'
;    stop
;  endif else begin
;    apparent_data_YY = data_cube/weights_cube
;    apparent_variance_YY = variance_cube / weights_cube_squared
;    quick_image, abs(apparent_data_YY), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted res in uv plane (yy) in Jy, longrun zenith', $
;      charsize=1, window=1,data_range = [0,1E1],savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_data_YY_longrun_zenith_tk',/png
;    quick_image, abs(apparent_variance_YY), x_arr, y_arr, xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted variance in uv plane (yy) in Jy, longrun zenith', $
;      charsize=1, window=2,savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_variance_YY_longrun_zenith_tk',/png
;    ;save, apparent_data_YY, x_arr, y_arr, apparent_variance_YY, filename = '/nfs/mwa-00/h1/nbarry/UVfromFHD/uv_weighted_variance_YY_bh.sav'
;
;    mid = ((size(apparent_data_YY))[1]-1)/2+1
;    quick_image, abs(apparent_data_YY[mid-50:mid+50,0:50]), x_arr[mid-50:mid+50], y_arr[0:50], xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted res in uv plane (yy) in Jy, longrun zoom', $
;      charsize=1, window=3, data_range = [0,1E2];, savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_data_YY_longrun_zoom_tk',/png
;    quick_image, abs(apparent_variance_YY[mid-50:mid+50,0:50]), x_arr[mid-50:mid+50], y_arr[0:50], xtitle='U (wavelengths)', ytitle='V (wavelengths)', title='Mag of weighted variance in uv plane (yy) in Jy, longrun zoom', $
;      charsize=1, window=4;,savefile = '/nfs/mwa-00/h1/nbarry/UVfromFHD/updated_July2016/uv_weighted_variance_XX_longrun_zenith_tk',/png
;    stop
;    stop
;  endelse
;
;  stop
  
end