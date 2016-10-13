pro oned_plotter

	power_24 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_GLEAM_beam2b_filtered/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_model_xx_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
		'power')
	power_1 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_GLEAM_beam2b/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_model_xx_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
		'power')
	k_1 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_GLEAM_beam2b/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_model_xx_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
		'k_edges')
	k_24 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_model_GLEAM_beam2b_filtered/ps/Combined_obs_obs_id_6296_cubeXX__even_odd_joint_model_xx_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
		'k_edges')
	
	
		
	cgPS_Open,'/nfs/eor-00/h1/nbarry/1Dfiltered.png',/quiet,/nomatch	
	cgplot, k_1, power_1,/ylog,/xlog,psym=10, yrange=[10^4., 10^15.],xrange=[.005,2.],charsize=1, xtitle = 'k (h / Mpc)',$
		ytitle='P (mK^2 Mpc^3 / h^3)',title = 'Model 1D power 10-50$\lambda$ with filtering'
	cgoplot, k_24, power_24,/ylog,/xlog,psym=10, yrange=[10^4., 10^15.],xrange=[.005,2.],charsize=1, color='blue'
	cgoplot, k_24, power_1-power_24,/ylog,/xlog, psym=10, color='green'
	cglegend, title=['unfiltered','horizon delay filter', 'unfiltered - filtered'],color=['black','blue','green'], location=[.5,.85], charsize=1
	cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
	
end