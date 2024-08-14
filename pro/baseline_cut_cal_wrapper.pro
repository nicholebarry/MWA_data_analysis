pro baseline_cut_cal_wrapper

  ; parse command line args
  compile_opt strictarr
  args = Command_Line_Args(count=nargs)
  IF keyword_set(args) then begin
    obs_id = args[0]
    output_directory = args[1]
    version = args[2]
  endif
  
  vis_XX = getvar_savefile(output_directory + 'fhd_' + version + '/vis_data/'+obs_id+'_vis_XX.sav', 'vis_ptr')
  vis_YY = getvar_savefile(output_directory + 'fhd_' + version + '/vis_data/'+obs_id+'_vis_YY.sav', 'vis_ptr')
  vis_ptr = PTRARR(2,/allocate)
  *vis_ptr[0] = *vis_XX
  *vis_ptr[1] = *vis_YY
  
  vis_model_XX = getvar_savefile(output_directory + 'fhd_' + version + '/vis_data/'+obs_id+'_vis_model_XX.sav', 'vis_model_ptr')
  vis_model_YY = getvar_savefile(output_directory + 'fhd_' + version + '/vis_data/'+obs_id+'_vis_model_YY.sav', 'vis_model_ptr')
  
  vis_weight_ptr = getvar_savefile(output_directory + 'fhd_' + version + '/vis_data/'+obs_id+'_flags.sav','vis_weights')
  
  obs = getvar_savefile(output_directory + 'fhd_' + version + '/metadata/'+obs_id+'_obs.sav', 'obs')
  params = getvar_savefile(output_directory + 'fhd_' + version + '/metadata/'+obs_id+'_params.sav', 'params')
  cal = getvar_savefile(output_directory + 'fhd_' + version + '/calibration/'+obs_id+'_cal.sav', 'cal')
  cal.time_avg = 0
  
  vis_model_ptr = PTRARR(2,/allocate)
  *vis_model_ptr[0]=*vis_model_XX
  *vis_model_ptr[1]=*vis_model_YY
  
  apply=1
  if keyword_set(apply) then begin
  
    n_freq=cal.n_freq
    tile_A_i=cal.tile_A-1 ;tile numbering starts at 1
    tile_B_i=cal.tile_B-1 ;tile numbering starts at 1
    n_baselines=Long(N_Elements(tile_A_i))
    
    inds_A=Rebin(Lindgen(n_freq),n_freq,n_baselines,/sample)+Rebin(transpose(tile_A_i)*n_freq,n_freq,n_baselines)
    inds_B=Rebin(Lindgen(n_freq),n_freq,n_baselines,/sample)+Rebin(transpose(tile_B_i)*n_freq,n_freq,n_baselines)
    
    for pol_i=0, 1 do begin
      vis_cal=*vis_ptr[pol_i]
      gain_arr = (*cal.gain[pol_i])
      vis_gain=gain_arr[inds_A]*Conj(gain_arr[inds_B])
      vis_cal*=(vis_gain)
      *vis_ptr[pol_i] = vis_cal
    endfor
  endif
  
  cal.min_cal_baseline = 20.;7.
  cal.max_cal_baseline = 150.
  
  ;(*obs.baseline_info).tile_use[76]=0
  cal2=vis_calibrate_subroutine_mod(vis_ptr, vis_model_ptr, vis_weight_ptr, obs, params, cal)
  ;cal2=vis_cal_temp2(vis_ptr, vis_model_ptr, vis_weight_ptr, obs, params, cal)
  save, cal2, filename=output_directory + 'fhd_' + version + '/selected_cal/GT7LT20_GT150/'+obs_id+'_GT7LT20_GT150_cal.sav'
  
end