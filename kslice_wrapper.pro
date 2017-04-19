pro kslice_wrapper

;obs_name='beardsley_thesis_list'
obs_name= 'obs_id_6296'

  power_filename = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_fivegrid_zenithmissing_2048/ps/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_model_xx_averemove_bh_power.idlsave'
  restore, power_filename
  
  power_filename = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_fivegrid_zenithmissing_4096/ps/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_model_xx_averemove_bh_power.idlsave'
  power_3D_4096 = getvar_savefile(power_filename, 'power_3D')
  
 ;   power_filename = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_4000_2017/ps/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_res_xx_averemove_bh_power.idlsave'
 ; weights_quarter_power = getvar_savefile(power_filename, 'weights_3D')
  
 ;     power_filename = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_4000_2017/ps/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_dirty_xx_averemove_bh_power.idlsave'
 ; dirty_power = getvar_savefile(power_filename, 'power_3D')
  
  metadata_struct = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_fivegrid_zenithmissing_2048/ps/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_info.idlsave','metadata_struct')
  freq = metadata_struct.frequencies
  
  z0_freq = 1420.40;E6 ;; Hz
  redshifts = z0_freq/freq - 1 ;; frequencies will be identical if kx, ky, kz match
  cosmology_measures, redshifts, wedge_factor = wedge_factor, comoving_dist_los=comov_dist_los, Ez=Ez
  
  ;*****calculate wedge
  ;; assume 20 degrees from pointing center to first null
  source_dist = 20d * !dpi / 180d
  fov_amp = wedge_factor * source_dist
  max_theta= 10d ;random guess!
  
  ;; calculate angular distance to horizon
  horizon_amp = wedge_factor * ((max_theta+90d) * !dpi / 180d)
  data_range_set = [1E5,1E15]
  
  for kz_i=0, N_elements(freq)/2-1 do begin
  wedge_amp = [[fov_amp[kz_i]], [horizon_amp[kz_i]]] ; to be multiplied by k perp ;WHY was this the 10th index before?
  plotfile='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_fivegrid_zenithmissing_2048/uvkz/2048-4096/uvkz_diff__'+string(strtrim(kz_i,2),FORMAT='(I03)')
  ;noise_ratio is power=reform(noise_3D[*,*,kz_i])*sqrt(reform(weights_3D[*,*,kz_i]))
  kpower_slice_plot, power=reform(power_3D[*,*,kz_i])-reform(power_3D_4096[*,*,kz_i]), xarr=kx_mpc, yarr=ky_mpc, kz_mpc = kz_mpc[kz_i], /baseline_axis, /plot_wedge_line, $
    slice_axis=2, kperp_lambda_conv = kperp_lambda_conv, wedge_amp=wedge_amp, /linear_axes,charsize_in=1,plotfile=plotfile,/png
  endfor
  
end