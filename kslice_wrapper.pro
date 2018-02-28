pro kslice_wrapper

  obs_name='beardsley_thesis_list_minustwo'
  ;obs_name= 'Aug23_longrunstyle'
  ;obs_name= 'obs_id_6296'
  ;obs_name= '1061316296'
  ;obs_name='season2_1089577208'
  ;obs_name='season2_zenith_calcut_1-100'
  ;obs_name='beardsley_thesis_list_int_chunk1'
  
 ;power_filename = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal/ps/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_tk_res_yy_averemove_bh_power.idlsave'
  ;power_filename = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_uvf_pskspan300/ps_150/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_dirty_yy_averemove_bh_power.idlsave'
  power_filename1 = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_bh/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_maxuv100_bh_ch10-127_dirty_yy_swbh_power.idlsave'
  ;power_filename = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_Aug2017_autocal1_v2/ps/'+obs_name+'_gridded_uvf__even_odd_joint_res_yy_averemove_bh_power.idlsave'
  restore, power_filename1
  ;power_3D_1 = getvar_savefile(power_filename1, 'power_3D')
  
  ;     power_filename = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_4000_2017/ps/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_dirty_xx_avceremove_bh_power.idlsave'
  ; dirty_power = getvar_savefile(power_filename, 'power_3D')
  
  ;metadata_struct = getvar_savefile('/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_uvf_pskspan300/ps_uvf/'+obs_name+'_gridded_uvf__even_odd_joint_uvf_info.idlsave','metadata_struct')
  metadata_struct = getvar_savefile('/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_bh/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_info.idlsave','metadata_struct')
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
  
  ;power_weights = (8*(sigma2_1)^2d)
  ;power_1_num = abs(data_sum_1)^2. - abs(data_diff_1)^2.
  
  ;for kz_i=0, N_elements(freq)/2-1 do begin
  kz_i=8
  wedge_amp = [[fov_amp[kz_i]], [horizon_amp[kz_i]]]
  ;plotfile='/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_uvf_pskspan300/ps_150/kzslice_weights_'+string(strtrim(kz_i,2),FORMAT='(I03)')
  ;noise_ratio is power=reform(noise_3D[*,*,kz_i])*sqrt(reform(weights_3D[*,*,kz_i]))
  ;thermal noise is reform(weight_invert(sqrt(weights_3d[*,*,kz_i])))
  ;power is power=power_3d[*,*,kz_i]
  kpower_slice_plot, power=power_3d[*,*,kz_i], xarr=kx_mpc, yarr=ky_mpc, kz_mpc = kz_mpc, /baseline_axis, /plot_wedge_line, $
    slice_axis=2, kperp_lambda_conv = kperp_lambda_conv, wedge_amp=wedge_amp, /linear_axes,charsize_in=1,window=3,data_range=[10^(3.),10^(15.)];,plotfile=plotfile,/png;, data_range=[10^3.,10^15.]
  stop
;endfor
  
end