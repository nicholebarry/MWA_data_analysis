pro temp

  vis_XX = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/vis_data/1061316176_vis_XX.sav', 'vis_ptr')
  vis_YY = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/vis_data/1061316176_vis_YY.sav', 'vis_ptr')
  ;vis_XX = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/vis_data/1061316296_vis_XX.sav', 'vis_ptr')
  ;vis_YY = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/vis_data/1061316296_vis_YY.sav', 'vis_ptr')
  vis_ptr = PTRARR(2,/allocate)

  zero_inds_XX = where(*vis_XX EQ 0,nx)
  zero_inds_YY = where(*vis_YY EQ 0,ny)
  
  vis_model_XX = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/vis_data/1061316176_vis_model_XX.sav', 'vis_model_ptr')
  vis_model_YY = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/vis_data/1061316176_vis_model_YY.sav', 'vis_model_ptr')
  ;vis_model_XX = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/vis_data/1061316296_vis_model_XX.sav', 'vis_model_ptr')
  ;vis_model_YY = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/vis_data/1061316296_vis_model_YY.sav', 'vis_model_ptr')
  
  vis_weight_ptr = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/vis_data/1061316176_flags.sav','flag_arr')
  ;vis_weight_ptr = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/vis_data/1061316296_flags.sav','vis_weights')
  
  if nx GT 0 then (*vis_weight_ptr[0])[zero_inds_XX]=0
  if ny GT 0 then (*vis_weight_ptr[1])[zero_inds_YY]=0
  
  obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/metadata/1061316176_obs.sav', 'obs')
  params = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/metadata/1061316176_params.sav', 'params')
  cal = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/calibration/1061316176_cal.sav', 'cal')
  cal3 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/calibration/double_stef_cal_step1.sav', 'cal2')
  ;obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/1061316296_obs.sav', 'obs')
  ;params = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/1061316296_params.sav', 'params')
  ;cal = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/calibration/1061316296_cal.sav', 'cal')
  cal.time_avg = 1
  
  ;*vis_XX = *vis_model_XX

  ;*vis_YY=*vis_model_YY
    *vis_ptr[0]=*vis_XX
  *vis_ptr[1]=*vis_YY

  
  vis_model_ptr = PTRARR(2,/allocate)
  *vis_model_ptr[0]=*vis_model_XX
  *vis_model_ptr[1]=*vis_model_YY
  
  apply=1
  if keyword_set(apply) then begin
    ;cal = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included/calibration/1061316176_cal.sav', 'cal')
    for pol_i=0,1 do begin
    ;;  (*cal.gain[pol_i]) = complex(FLTARR(384, 128)) +1.; + Complex(0,1)*1.
      (*cal.gain[pol_i])[*,0] = reform((*cal.gain[pol_i])[*,0]) + real_part(.2*exp(Complex(0,1)*30.*FINDGEN(384)/384.))
    endfor
    
    n_freq=cal.n_freq
    tile_A_i=cal.tile_A-1 ;tile numbering starts at 1
    tile_B_i=cal.tile_B-1 ;tile numbering starts at 1
    n_baselines=Long(N_Elements(tile_A_i))
    
    inds_A=Rebin(Lindgen(n_freq),n_freq,n_baselines,/sample)+Rebin(transpose(tile_A_i)*n_freq,n_freq,n_baselines)
    inds_B=Rebin(Lindgen(n_freq),n_freq,n_baselines,/sample)+Rebin(transpose(tile_B_i)*n_freq,n_freq,n_baselines)
    
    for pol_i=0, 1 do begin
      vis_cal=*vis_ptr[pol_i]
      gain_arr = (*cal.gain[pol_i])
      gain_arr2 = (*cal3.gain[pol_i]); + (*cal.gain_residual[pol_i])
      vis_gain=gain_arr[inds_A]*Conj(gain_arr[inds_B])
      vis_gain2=gain_arr2[inds_A]*Conj(gain_arr2[inds_B])
      vis_cal*=(vis_gain)
      vis_cal*=weight_invert(vis_gain2)
      *vis_ptr[pol_i] = vis_cal
    endfor
  endif
  
  ;(*obs.baseline_info).tile_use[76]=0
  ;cal2=vis_calibrate_subroutine_tau_g2(vis_ptr, vis_model_ptr, vis_weight_ptr, obs, params, cal)
  cal2=vis_calibrate_subroutine(vis_ptr, vis_model_ptr, vis_weight_ptr, obs, params, cal)
  ;cal2=vis_cal_temp2(vis_ptr, vis_model_ptr, vis_weight_ptr, obs, params, cal)
  stop
  
end