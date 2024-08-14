pro diffuse_normalization

  dirty_uv_arr_even = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_baseline_cut_50_2017/1061316296_even_gridded_uvf.sav','dirty_uv_arr')
  dirty_uv_arr_odd = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_baseline_cut_50_2017/1061316296_odd_gridded_uvf.sav','dirty_uv_arr')
  weights_uv_arr_even = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_baseline_cut_50_2017/1061316296_even_gridded_uvf.sav','weights_uv_arr')
  weights_uv_arr_odd = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_baseline_cut_50_2017/1061316296_odd_gridded_uvf.sav','weights_uv_arr')
  
  weighted_data = ((*dirty_uv_arr_even[0,0])+(*dirty_uv_arr_odd[0,0]))*weight_invert((*weights_uv_arr_even[0,0])+(*weights_uv_arr_odd[0,0]))
  for freq_i=1,191 do weighted_data = weighted_data + ((*dirty_uv_arr_even[0,freq_i])+(*dirty_uv_arr_odd[0,freq_i]))*weight_invert((*weights_uv_arr_even[0,freq_i])+(*weights_uv_arr_odd[0,freq_i]))
  
    weighted_data_YY = ((*dirty_uv_arr_even[1,0])+(*dirty_uv_arr_odd[1,0]))*weight_invert((*weights_uv_arr_even[1,0])+(*weights_uv_arr_odd[1,0]))
  for freq_i=1,191 do weighted_data_YY = weighted_data_YY + ((*dirty_uv_arr_even[1,freq_i])+(*dirty_uv_arr_odd[1,freq_i]))*weight_invert((*weights_uv_arr_even[1,freq_i])+(*weights_uv_arr_odd[1,freq_i]))
  
  weights_tot=(*weights_uv_arr_even[0,0])
  for freq_i=1,191 do weights_tot = weights_tot + (*weights_uv_arr_even[0,freq_i])
  ;dirty_tot=(*dirty_uv_arr_even[0,0])
  ;for freq_i=1,191 do dirty_tot = dirty_tot + (*dirty_uv_arr_even[0,freq_i])
  min_ind=where(weights_tot LT 1, n_count)
  col = min_ind mod 1200.
  row = min_ind / 1200.
  if n_count GT 0 then weights_tot[min_ind]=0
  range = [495,705]
  ;range=[500,700]
  image = (shift(fft(weighted_data[range[0]:range[1],range[0]:range[1]],/inverse),(range[1]-range[0])/2.,(range[1]-range[0])/2.))/192./(2.*!Pi)^3. ;freq averaged
    image_YY = (shift(fft(weighted_data_YY[range[0]:range[1],range[0]:range[1]],/inverse),(range[1]-range[0])/2.,(range[1]-range[0])/2.))/192./(2.*!Pi)^3. ;freq averaged

  
    N=(range[1]-range[0])
  T= .5 ;kbinsize of half a wavelength
  X = FINDGEN((N - 1)/2) + 1
  x_axis = [0.0, X, N/2, -N/2 + X]/(N*T)*180./!Pi
  
    beam_arr = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_baseline_cut_50_2017_beam/1061316296_even_gridded_imagecube.sav', 'beam_arr')
  beam = *beam_arr[0,0]
  for freq_i=1,191 do beam = beam + *beam_arr[0,freq_i]
  beam = beam/192. ; freq averaged
    beam_YY = *beam_arr[0,0]
  for freq_i=1,191 do beam_YY = beam_YY + *beam_arr[0,freq_i]
  beam_YY = beam_YY/192. ; freq averaged
  
    image_new = (abs((image)*weight_invert(beam[range[0]:range[1],range[0]:range[1]])+image_YY*weight_invert(beam_YY[range[0]:range[1],range[0]:range[1]])))/2.
  stop
  
  ;res_image_arr_even = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_baseline_cut_50_2017/1061316296_even_gridded_imagecube.sav', 'residual_arr1')
  ;weights_image_arr_even = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_baseline_cut_50_2017/1061316296_even_gridded_imagecube.sav','weights_arr1')
  ;res_image_arr_odd = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_baseline_cut_50_2017/1061316296_odd_gridded_imagecube.sav', 'residual_arr1')
  ;weights_image_arr_odd = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_baseline_cut_50_2017/1061316296_odd_gridded_imagecube.sav','weights_arr1')
  
  ;weighted_dirty_zero = (*res_image_arr_even[0,0])/(*weights_image_arr_even[0,0])+(*res_image_arr_odd[0,0])/(*weights_image_arr_odd[0,0])
  
  weights_image_tot=(*weights_image_arr_even[0,0])
  for freq_i=1,191 do weights_image_tot = weights_image_tot + (*weights_image_arr_even[0,freq_i])
  res_image_tot=(*res_image_arr_even[0,0])
  for freq_i=1,191 do res_image_tot = res_image_tot + (*res_image_arr_even[0,freq_i])
  original_image = (abs(res_image_tot*weight_invert(weights_image_tot)))[500:700,500:700]
  
  stop
  
end