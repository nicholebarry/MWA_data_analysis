pro apparent_frame_from_uvf_cubes
  ;Must be run from gridengine since the memory load is too high

  dirty_even = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/calibration_sim/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod/1061316176_even_gridded_uvf.sav','dirty_uv_arr')
  dirty_odd = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/calibration_sim/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod/1061316176_odd_gridded_uvf.sav','dirty_uv_arr')
  model_even = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/calibration_sim/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod/1061316176_even_gridded_uvf.sav','model_uv_arr')
  model_odd = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/calibration_sim/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod/1061316176_odd_gridded_uvf.sav','model_uv_arr')
  
  res_evenodd_XX = complex(FLTARR(1200,1200,192))
  for freq_i=0, 191 do res_evenodd_XX[*,*,freq_i] = (*dirty_even[0,freq_i] - *model_even[0,freq_i]) + (*dirty_odd[0,freq_i] - *model_odd[0,freq_i])
  undefine, dirty_even, model_even
  
  res_XX = total(res_evenodd_XX,3)
  
  weights_even = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/calibration_sim/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod/1061316176_even_gridded_uvf.sav','weights_uv_arr')
  weights_odd = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/calibration_sim/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod/1061316176_odd_gridded_uvf.sav','weights_uv_arr')
  
  weights_evenodd_XX = complex(FLTARR(1200,1200,192))
  for freq_i=0, 191 do weights_evenodd_XX[*,*,freq_i] = *weights_even[0,freq_i] + *weights_odd[0,freq_i]
  
  weights_XX = total(weights_evenodd_XX,3)
  
  undefine, res_evenodd_XX, weights_evenodd_XX
  
  data_XX = res_XX/weights_XX
  save, data_XX,res_XX,weights_XX, filename='/nfs/mwa-00/h1/nbarry/uvf_cube_apparent_data.sav'
  
end