function translate_1d_power,k_edges,k_bin,hubble_param,power,k_out=k_out
  kscale = k_edges[1:*]/hubble_param
  k_mid = kscale[1:*] - k_bin/hubble_param
  ;k_out = [min(kscale), k_mid, max(kscale)]
  k_out = [min(kscale)-k_bin/hubble_param/2, k_mid, max(kscale),max(kscale)+k_bin/hubble_param]
  power_out = hubble_param^3. * [power[1],power[1:-1],power[-1]] / (2.*!Pi^2)
  
  return,power_out
end

pro limit_compare

  hubble_param = 0.71d
  
  filename_1 = '/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps_jan2016/'
  filename_2 = '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/ps/'
  
  obsname_1='wedge_cut_plus_res_cut'
  obsname_2='beardsley_thesis_list'
  analysis_1='ch0-95_res_xx_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-60_kpar0.15-200_1dkpower'
  analysis_2='ch0-95_res_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_kpar0.15-200_1dkpower'
  power_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_'+analysis_1+'.idlsave', 'power')
  k_edges_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_'+analysis_1+'.idlsave', 'k_edges')
  k_bin_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_'+analysis_1+'.idlsave', 'k_bin')
  power_2 = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_'+analysis_2+'.idlsave', 'power')
  k_edges_2 = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_'+analysis_2+'.idlsave', 'k_edges')
  k_bin_2 = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_'+analysis_2+'.idlsave', 'k_bin')
  
  k_mid_1 = (k_edges_1 - k_bin_1/2d)/hubble_param
  k_mid_2 = (k_edges_2 - k_bin_1/2d)/hubble_param; + k_bin_2/2d)/hubble_param
  delta_power_1 = (hubble_param*k_mid_1)^3d * power_1 / (2.*!Pi^2.)
  delta_power_2 = (hubble_param*k_mid_2)^3d * power_2 / (2.*!Pi^2.)
  
  delta_power_1_arr=[delta_power_1[0], delta_power_1, delta_power_1[N_elements(delta_power_1)-1]]
  k_mid_1=[k_edges_1[0],k_edges_1+k_bin_1/2d,max(k_edges_1)+k_bin_1/2d]/hubble_param
  
  delta_power_2_arr=[delta_power_2[0], delta_power_2, delta_power_2[N_elements(delta_power_2)-1]]
  k_mid_2=[k_edges_2[0],k_edges_2+k_bin_2/2d,max(k_edges_2)+k_bin_2/2d]/hubble_param
  
  ;power_out_1 = translate_1d_power(k_edges_1,k_bin_1,hubble_param,power_1,k_out=k_out_1)
  
  cgplot, (k_mid_1), abs(delta_power_1_arr),/ylog,/xlog,psym=10,charsize=1, xtitle = 'k (h / Mpc)',yrange=[10^3.,10^9.],xrange=[.1,1.5],$
    ytitle='$\Delta$^2 (mK^2)',title = '1D power 10-60$\lambda$ res XX'
  cgoplot, (k_mid_2), abs(delta_power_2_arr),/ylog,/xlog,psym=10, color='blue'
  
  ;ps_wrapper, '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal/', 'beardsley_thesis_list_firstthird', /png, plot_1d_delta=1,wedge_angles=120d , freq_ch_range=[0,95], kperp_range_lambda_1dave=[10,60], wedge_angles=120d
  
  stop
  
end