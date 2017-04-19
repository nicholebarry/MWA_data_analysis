pro oned_plotter

;	power_24 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_through_firstpass/ps_Tukey_2_alpha5_norm/Combined_obs_Aug23_cubeXX__even_odd_joint_deluv2_res_xx_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
;		'power')
;	power_1 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_through_firstpass2/ps_Tukey_2/Combined_obs_Aug23_cubeXX__even_odd_joint_deluv2_res_xx_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
;		'power')
;	k_1 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_through_firstpass2/ps_Tukey_2/Combined_obs_Aug23_cubeXX__even_odd_joint_deluv2_res_xx_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
;		'k_edges')
;	k_24 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_through_firstpass/ps_Tukey_2_alpha5_norm/Combined_obs_Aug23_cubeXX__even_odd_joint_deluv2_res_xx_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
;		'k_edges')
		
;		power_24 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_4000_2017/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
;    'power')
;      k_24 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_4000_2017/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
;    'k_edges')
;  power_24 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_4000_2017_gauss/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
;    'power')
;  k_24 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_4000_2017_gauss/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
;    'k_edges')
  ;    power_2 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_4000_2017_gauss_norm/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
  ;  'power')
  ;k_2 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_delay_4000_2017_gauss_norm/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
  ;  'k_edges')
;  power_2 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/calibration_sim/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod/ps_regular/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_model_xx_bh_dencorr_1dkpower.idlsave', $
;    'power')
;      k_2 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/calibration_sim/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod/ps_regular/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_model_xx_bh_dencorr_1dkpower.idlsave', $
;    'k_edges')
power_24 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_longrun_cal_amp_dig/ps/Combined_obs_Aug23_longrunstyle_cubeXX__even_odd_joint_dirty_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
    'power')
      k_24 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_longrun_cal_amp_dig/ps/Combined_obs_Aug23_longrunstyle_cubeXX__even_odd_joint_dirty_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
    'k_edges')
    power_2 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_longrun_cal_amp_dig_conv/ps/Combined_obs_Aug23_longrunstyle_cubeXX__even_odd_joint_dirty_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
    'power')
      k_2 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_longrun_cal_amp_dig_conv/ps/Combined_obs_Aug23_longrunstyle_cubeXX__even_odd_joint_dirty_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
    'k_edges')
;  power_1 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_delaycut_1source/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_model_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
;    'power')
;  k_1 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_delaycut_1source/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_model_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
;    'k_edges')

		
	;cgPS_Open,'/nfs/eor-00/h1/nbarry/1_1000_1000.png',/quiet,/nomatch	
	cgplot, k_24, power_24,/ylog,/xlog,psym=10,xrange=[.1,1.3],charsize=1, xtitle = 'k (h / Mpc)',yrange=[10^3., 10^15.],$
		ytitle='P (mK^2 Mpc^3 / h^3)',title = 'Res 1D power 10-50$\lambda$ XX'
	;cgoplot, k_1, power_1,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='blue'
	cgoplot, k_2, power_2,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='green'
	;cgoplot, k_24, power_1-power_24,/ylog,/xlog, psym=10, color='green'
	;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
	stop
	
end