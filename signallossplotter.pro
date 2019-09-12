power_eor = getvar_savefile('/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave','power')
k_centers_eor = getvar_savefile('/Users/nabarry/MWA/FHD/catalog_data/eor_power_1d.idlsave','k_centers')

power_signalloss_nobh = getvar_savefile('/Users/nabarry/Combined_obs_zenith_cubeXX__even_odd_joint_res_xx_averemove_swbh_dencorr_no_horizon_wedge_kperplambda10-50_1dkpower_test.idlsave','power')
k_edges_signalloss_nobh = getvar_savefile('/Users/nabarry/Combined_obs_zenith_cubeXX__even_odd_joint_res_xx_averemove_swbh_dencorr_no_horizon_wedge_kperplambda10-50_1dkpower_test.idlsave','k_edges')

power_signalloss_bh = getvar_savefile('/Users/nabarry/MWA/data/fhd_nb_sim_perfect/ps/data/1d_binning/Combined_obs_zenith_cubeXX__even_odd_joint_bh_res_xx_averemove_swbh_nodensitycorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave','power')
k_edges_signalloss_bh = getvar_savefile('/Users/nabarry/MWA/data/fhd_nb_sim_perfect/ps/data/1d_binning/Combined_obs_zenith_cubeXX__even_odd_joint_bh_res_xx_averemove_swbh_nodensitycorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave','k_edges')

power_perf_bh = getvar_savefile('/Users/nabarry/MWA/data/fhd_nb_sim_perfect/ps_orig/data/1d_binning/Combined_obs_zenith_cubeXX__even_odd_joint_bh_res_xx_averemove_swbh_nodensitycorr_no_horizon_wedge_1dkpower.idlsave','power')
k_edges_perf_bh = getvar_savefile('/Users/nabarry/MWA/data/fhd_nb_sim_perfect/ps_orig/data/1d_binning/Combined_obs_zenith_cubeXX__even_odd_joint_bh_res_xx_averemove_swbh_nodensitycorr_no_horizon_wedge_1dkpower.idlsave','k_edges')
;IDL>/fred/oz048/MWA/CODE/FHD/fhd_nb_sim_perfect/ps/data/1d_binning/Combined_obs_zenith_cubeXX__even_odd_joint_bh_res_xx_averemove_swbh_dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave Combined_obs_zenith_cubeXX__even_odd_joint_bh_res_xx_averemove_swbh_dencorr_no_horizon_wedge_kperplambda10-50_1dkpower_orig.idlsave

n_k=n_elements(k_edges_perf_bh)
k_perf_bh = (k_edges_perf_bh[1:(n_k-1)]+k_edges_perf_bh[0:(n_k-2)])/2.

n_k=n_elements(k_edges_signalloss_bh)
k_signalloss_bh = (k_edges_signalloss_bh[1:(n_k-1)]+k_edges_signalloss_bh[0:(n_k-2)])/2.

n_k=n_elements(k_edges_signalloss_nobh)
k_signalloss_nobh = (k_edges_signalloss_nobh[1:(n_k-1)]+k_edges_signalloss_nobh[0:(n_k-2)])/2.

ytitle = 'P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)'

cgPS_Open,'/Users/nabarry/1d_compare.png',/quiet,/nomatch
cgplot, k_perf_bh, (power_perf_bh),/ylog,/xlog,psym=10,xrange=[.05,1.2],charsize=1.3, xtitle = 'k (!8h!x / Mpc)',yrange=[10^3., 10^7.],ytitle=ytitle
cgoplot, k_signalloss_bh, (power_signalloss_bh),/ylog,/xlog,psym=10, color='green'
cgoplot, k_signalloss_nobh, (power_signalloss_nobh),/ylog,/xlog,psym=10, color='blue'
cgoplot, k_centers_eor,power_eor,/ylog,/xlog,psym=10, color='purple'
cglegend, titles=['Foregrounds, bh','Signal loss sim, bh','Signal loss sim, postage stamp','Input EoR'], color=['black','green','blue','purple'], charsize=1, location=[.45,.85]
cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage