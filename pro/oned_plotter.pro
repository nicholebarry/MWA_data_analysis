pro oned_plotter

  plot_weights=1
  
  filename_24 = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_bh/'
  filename_1 = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_bh/'
  filename_2 = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_bh/'
  filename_3 = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_256bh/'
  ;filename_4 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_uvf_pskspan300_v2_noflags_large/ps/'
  ;obsname= 'wedge_cut_plus_res_cut'
  ;obsname_1='beardsley_thesis_list'
  obsname='beardsley_thesis_list_zenith'
  ;obsname='Aug23_longrunstyle'
  ;obsname_1='Aug23_longrunstyle'
  obsname_1='beardsley_thesis_list_minustwo'
  obsname_2='beardsley_thesis_list_plustwo'
  obsname_3='beardsley_thesis_list_256'
  ;obsname_4='obs_id_6296'
  cube_type='res'
  
  
  full_path_24=filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_maxuv100_bh_ch10-127_'+cube_type+'_yy_swbh_dencorr_no_110deg_wedge_kperplambda20-80_1dkpower.idlsave'
  power_24 = getvar_savefile(full_path_24, 'power')
  k_edges = getvar_savefile(full_path_24, 'k_edges')
  if keyword_set(plot_weights) then weights_24 = getvar_savefile(full_path_24, 'weights')
  n_k=n_elements(k_edges)
  k_24 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
  
  full_path_1=filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_maxuv100_bh_ch10-127_'+cube_type+'_yy_swbh_dencorr_no_110deg_wedge_kperplambda20-80_1dkpower.idlsave'
  power_1 = getvar_savefile(full_path_1, 'power')
  k_edges = getvar_savefile(full_path_1, 'k_edges')
  ;power_1 = getvar_savefile(filename_1 +'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_maxuv100_han_ch5-120_'+cube_type+'_yy_averemove_swbh_dencorr_no_115deg_wedge_kperplambda20-75_1dkpower.idlsave', $
  ;  'power')
  ;k_edges = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_maxuv100_han_ch5-120_'+cube_type+'_yy_averemove_swbh_dencorr_no_115deg_wedge_kperplambda20-75_1dkpower.idlsave', $
  ;  'k_edges')
  if keyword_set(plot_weights) then weights_1 = getvar_savefile(full_path_1, 'weights')
  n_k=n_elements(k_edges)
  k_1 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
  
  full_path_2=filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_maxuv100_bh_ch10-127_'+cube_type+'_yy_swbh_dencorr_no_110deg_wedge_kperplambda20-80_1dkpower.idlsave'
  power_2 = getvar_savefile(full_path_2, 'power')
  k_edges = getvar_savefile(full_path_2, 'k_edges')
  if keyword_set(plot_weights) then weights_2 = getvar_savefile(full_path_2, 'weights')
  n_k=n_elements(k_edges)
  k_2 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
  ;
  full_path_3=filename_3 + 'Combined_obs_'+obsname_3+'_cubeXX__even_odd_joint_maxuv100_bh_ch5-120_'+cube_type+'_yy_averemove_swbh_dencorr_no_115deg_wedge_kperplambda20-75_1dkpower.idlsave'
  power_3 = getvar_savefile(full_path_3, 'power')
  k_edges = getvar_savefile(full_path_3, 'k_edges')
  if keyword_set(plot_weights) then weights_3 = getvar_savefile(full_path_3, 'weights')
  n_k=n_elements(k_edges)
  k_3 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
  
  ;      power_4 = getvar_savefile(filename_4 + 'Combined_obs_'+obsname_4+'_cubeXX__even_odd_joint_bn_dirty_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
  ;    'power')
  ;  k_edges = getvar_savefile(filename_4 + 'Combined_obs_'+obsname_4+'_cubeXX__even_odd_joint_bn_dirty_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
  ;    'k_edges')
  ;  n_k=n_elements(k_edges)
  ;  k_4 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
  ;
  power_eor = getvar_savefile('/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave','power')
  k_centers_eor = getvar_savefile('/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave','k_centers')
  
  ;if keyword_set(rebin_eor) then begin
  ;  for bin_i=0, N_elements(k_centers_eor)-2 do begin
  ;  where((k_1 GT k_centers_eor[bin_i]) AND (k_1 LT k_centers_eor[bin_i]
  ;  endfor
  ;endif
  
  delta=1
  if keyword_set(delta) then begin
    if keyword_set(power_24) then power_24=power_24*(k_24^3.)/(2.*!pi^2.)
    if keyword_set(power_1) then power_1=power_1*(k_1^3.)/(2.*!pi^2.)
    if keyword_set(power_2) then power_2=power_2*(k_2^3.)/(2.*!pi^2.)
    if keyword_set(power_3) then power_3=power_3*(k_3^3.)/(2.*!pi^2.)
    if keyword_set(power_eor) then power_eor=power_eor*(k_centers_eor^3.)/(2.*!pi^2.)
    if keyword_set(plot_weights) then begin
      if keyword_set(power_24) then weights_24=(k_24^3.)/(2.*!pi^2.)/sqrt(weights_24)
      if keyword_set(power_1) then weights_1=(k_1^3.)/(2.*!pi^2.)/sqrt(weights_1)
      if keyword_set(power_2) then weights_2=(k_2^3.)/(2.*!pi^2.)/sqrt(weights_2)
      if keyword_set(power_3) then weights_3=(k_3^3.)/(2.*!pi^2.)/sqrt(weights_3)
    endif
  endif
  
  if ~keyword_set(delta) then ytitle = 'P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)' else ytitle = '$\Delta$$\up2$ (mK$\up2$)'
  
  ;cgPS_Open,'/nfs/eor-00/h1/nbarry/1d_yy.png',/quiet,/nomatch
  cgplot, k_24, (power_24),/ylog,/xlog,psym=10,xrange=[.1,1.2],charsize=1.3, xtitle = 'k (!8h!x / Mpc)',yrange=[10^2., 10^8.],$
    ytitle=ytitle,title = 'YY 1D '+cube_type+', 20-80 $\lambda$ with 110$\deg$ wedge'
  if keyword_set(plot_weights) then cgoplot, k_24, (weights_24),/ylog,/xlog,psym=10,linestyle=2
  cgoplot, k_1, (power_1),/ylog,/xlog,psym=10,$
    xrange=[.001,1.3],charsize=1, xtitle = 'k (h / Mpc)',yrange=[10^2., 10^14.],ytitle='P (mK^2 Mpc^3 / h^3)',title = '1D power, instrumental XX', color='blue'
  if keyword_set(plot_weights) then cgoplot, k_1, (weights_1),/ylog,/xlog,psym=10,linestyle=2, color='blue'
  cgoplot, k_2, power_2,/ylog,/xlog,psym=10,color='green',$
    xrange=[.01,1.2],charsize=1.3, xtitle = 'k (!8h!x / Mpc)',yrange=[10^1., 10^8.],ytitle='P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)',title = '1D res power, 10-50'
  if keyword_set(plot_weights) then cgoplot, k_2, (weights_2),/ylog,/xlog,psym=10,linestyle=2, color='green'
  ;cgoplot, k_3, power_3,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='pink'
  ;cgoplot, k_4, power_4,/ylog,/xlog,psym=10, charsize=1, color='purple'
  ;cgoplot, k_5, power_5,/ylog,/xlog,psym=10,charsize=1, color='magenta'
  ;cgoplot, k_6, power_6,/ylog,/xlog,psym=10,charsize=1, color='dark green'
  cgoplot, k_centers_eor, power_eor,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='purple'
  ;cglegend, title=['Floor, all','Floor, no sidelobes', 'Floor, no nulls','Floor, no primary', 'EoR'], color=['black','blue','green','pink','purple'], location=[.55,.8], charsize=1
  ;cglegend, title=['1','2', '4','6','8','12','16', 'EoR'], color=['black','blue','green','pink','sky blue','magenta','dark green','purple'], location=[.65,.85], charsize=1
  cglegend, title=['zenith','minustwo','plustwo','EoR'], color=['black','blue','green','purple'], location=[.15,.85], charsize=1
  ;  cglegend, title=['UV BH + image BH', 'EoR'], color=['black','purple'], location=[.65,.85], charsize=1
  
  ;cglegend, title=['input Gaussian EoR','recovered Gaussian EoR, 300 wavelength','recovered Gaussian EoR, 60 wavelength'], color=['purple','black','blue'], location=[.15,.80], charsize=1.3, thick=1
  ;cgoplot, k_24, power_1-power_24,/ylog,/xlog, psym=10, color='green'
  ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  stop
  
end