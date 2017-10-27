function translate_1d_power,k_edges,k_bin,hubble_param,power,k_out=k_out
  kscale = k_edges[1:*]/hubble_param
  k_mid = kscale[1:*] - k_bin/hubble_param
  ;k_out = [min(kscale), k_mid, max(kscale)]
  k_out = [min(kscale)-k_bin/hubble_param/2, k_mid, max(kscale),max(kscale)+k_bin/hubble_param]
  power_out = hubble_param^3. * [power[1],power[1:-1],power[-1]] / (2.*!Pi^2)
  
  return,power_out
end

pro limit_compare

  hubble_param = 0.71
  
  filename_1 = '/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps_jan2016/'
  filename_2 = '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/ps/'
  
  obsname_1='wedge_cut_plus_res_cut'
  obsname_2='beardsley_thesis_lis_noautocutt'
  analysis_1='ch0-95_res_yy_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-60_kpar0.15-200_1dkpower'
  analysis_2='tk_ch10-127_res_yy_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_kpar0.15-200_1dkpower'
  power_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_'+analysis_1+'.idlsave', 'power')
  k_edges_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_'+analysis_1+'.idlsave', 'k_edges')
  k_bin_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_'+analysis_1+'.idlsave', 'k_bin')
  power_2 = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_'+analysis_2+'.idlsave', 'power')
  k_edges_2 = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_'+analysis_2+'.idlsave', 'k_edges')
  k_bin_2 = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_'+analysis_2+'.idlsave', 'k_bin')
  
  delta_power_1 = k_edges_1^3. * power_1 / (2.*!Pi^2)
  delta_power_2 = k_edges_2^3. * power_2 / (2.*!Pi^2)
  
  cgplot, k_edges_1, abs(delta_power_1),/ylog,/xlog,psym=10,charsize=1, xtitle = 'k (h / Mpc)',yrange=[10^3.,10^8.],xrange=[.1,1.5],$
    ytitle='$\Delta$^2 (mK^2)',title = '1D power 10-60$\lambda$ res XX'
  cgoplot, k_edges_2, delta_power_2,/ylog,/xlog,psym=10, color='blue'
  
  ;ps_wrapper, '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/', 'beardsley_thesis_lis_noautocutt', /png, plot_1d_delta=1, image_filter_name='Tukey', freq_ch_range=[10,127], kperp_range_lambda=[10,60]
  
  stop
  
end