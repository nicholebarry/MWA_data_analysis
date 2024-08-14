pro weights_corr_sim

  filebase = '/fred/oz048/MWA/CODE/FHD/'

  ;recalculate_info=1

  if keyword_set(recalculate_info) then begin
    for run_num=6,9 do begin
      ps_wrapper, filebase + 'fhd_nb_arrsim_flat_'+strtrim(run_num,2)+'_full/', 'zenith',refresh_binning = 1, $
        cube_power_info = cube_power_info,/png,wt_cutoffs=0,/refresh_ps
      save, cube_power_info, filename=filebase + 'fhd_nb_arrsim_flat_'+strtrim(run_num,2)+'_full/cube_power_info.sav'
    endfor
  endif

  nbsl_lambda2 = FLTARR(7)
  wt_ave_power = FLTARR(7)
  run_num_index = [3,4,5,6,7,8,9]

  for run_num_i=0,6 do begin
    cube_power_info = getvar_savefile(filebase + 'fhd_nb_arrsim_flat_'+strtrim(run_num_index[run_num_i],2)+'_full/cube_power_info.sav','cube_power_info')
    if keyword_set(cube_power_info) then begin
    nbsl_lambda2[run_num_i] = cube_power_info.nbsl_lambda2[0]
    wt_ave_power[run_num_i] = cube_power_info.wt_ave_power[0]
    endif
    undefine, cube_power_info
  endfor

  stop

end
