pro degrade_healpix_sav, filename_in,int_list,nside_out

  evenodd = ['even','odd']
  pols = ['XX','YY']
  n_freq=192.
  
  for evenodd_i=0, 1 do begin
    for pol_i=0,1 do begin
      restore, filename_in + '/Combined_obs_'+int_list+'_'+evenodd[evenodd_i]+'_cube'+pols[pol_i]+'.sav'
      
      for freq_i=0, n_freq-1 do begin
        dirty_temp = reform(dirty_cube[*,freq_i])
        model_temp = reform(model_cube[*,freq_i])
        variance_temp = reform(variance_cube[*,freq_i])
        weights_temp = reform(weights_cube[*,freq_i])
        beam_squared_temp = reform(beam_squared_cube[*,freq_i])
        
        hpx_inds_temp = hpx_inds
        ud_grade_cut4,hpx_inds_temp, dirty_temp, nside_in=nside, nside_out=nside_out,order_in='RING'
        if freq_i EQ 0 then dirty_cube_downgrade = dirty_temp else dirty_cube_downgrade = [[dirty_cube_downgrade],[temporary(dirty_temp)]]
        hpx_inds_temp = hpx_inds
        ud_grade_cut4,hpx_inds_temp, model_temp, nside_in=nside, nside_out=nside_out,order_in='RING'
        if freq_i EQ 0 then model_cube_downgrade = model_temp else model_cube_downgrade = [[model_cube_downgrade],[temporary(model_temp)]]
        hpx_inds_temp = hpx_inds
        ud_grade_cut4,hpx_inds_temp, variance_temp, nside_in=nside, nside_out=nside_out,order_in='RING'
        if freq_i EQ 0 then variance_cube_downgrade = variance_temp else variance_cube_downgrade = [[variance_cube_downgrade],[temporary(variance_temp)]]
        hpx_inds_temp = hpx_inds
        ud_grade_cut4,hpx_inds_temp, weights_temp, nside_in=nside, nside_out=nside_out,order_in='RING'
        if freq_i EQ 0 then weights_cube_downgrade = weights_temp else weights_cube_downgrade = [[weights_cube_downgrade],[temporary(weights_temp)]]
        hpx_inds_temp = hpx_inds
        ud_grade_cut4,hpx_inds_temp, beam_squared_temp, nside_in=nside, nside_out=nside_out,order_in='RING'
        if freq_i EQ 0 then beam_squared_cube_downgrade = beam_squared_temp else beam_squared_cube_downgrade = [[beam_squared_cube_downgrade],[temporary(beam_squared_temp)]]
      endfor
      hpx_inds = hpx_inds_temp
      nside=nside_out
      beam_squared_cube=beam_squared_cube_downgrade
      dirty_cube=dirty_cube_downgrade
      model_cube=model_cube_downgrade
      variance_cube=variance_cube_downgrade
      weights_cube=weights_cube_downgrade

      save, beam_squared_cube, dirty_cube,frequencies, hpx_inds, model_cube, nside, n_avg, obs_arr, variance_cube, weights_cube,$
        filename = filename_in + '/Combined_obs_'+int_list+'_'+strtrim(nside_out,2)+'_'+evenodd[evenodd_i]+'_cube'+pols[pol_i]+'.sav'
    endfor
  endfor
  
  
end