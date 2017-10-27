pro oned_plotter

filename_24 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_goldenset_bubbles_over_all/ps/'
filename_1 = '/nfs/mwa-04/r1/EoRuvfits/analysis/calibration_sim/fhd_nb_goldenset_bubbles/ps/'
;obsname= 'wedge_cut_plus_res_cut'
;obsname_1='beardsley_thesis_list'
obsname='obs_id_6296'
obsname_1='obs_id_6296'
  ;power_24 = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
  ;  'power')
  ;k_24 = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
  ;  'k_edges')
  power_24 = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
    'power')
  k_24 = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
    'k_edges')
  power_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
    'power')
  k_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
    'k_edges')

  power_eor = getvar_savefile('/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave','power')
  k_centers_eor = getvar_savefile('/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave','k_centers')
  
  ;cgPS_Open,'/nfs/eor-00/h1/nbarry/hash.png',/quiet,/nomatch
  cgplot, k_24, (power_24),/ylog,/xlog,psym=10,xrange=[.001,1.3],charsize=1, xtitle = 'k (h / Mpc)',yrange=[10^2., 10^7.],$
    ytitle='P (mK^2 Mpc^3 / h^3)',title = '1D power 10-60$\lambda$ res XX'
  cgoplot, k_1, (power_1),/ylog,/xlog,psym=10,xrange=[.001,1.3],charsize=1, xtitle = 'k (h / Mpc)',yrange=[10^2., 10^7.],$
    ytitle='P (mK^2 Mpc^3 / h^3)',title = '1D power 10-50$\lambda$ res XX', color='blue'
  ;cgoplot, k_2, power_2,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='green'
  ;cgoplot, k_3, power_3,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='pink'
  ;cgoplot, k_4, power_4,/ylog,/xlog,psym=10, charsize=1, color='sky blue'
  ;cgoplot, k_5, power_5,/ylog,/xlog,psym=10,charsize=1, color='magenta'
  ;cgoplot, k_6, power_6,/ylog,/xlog,psym=10,charsize=1, color='dark green'
  cgoplot, k_centers_eor, power_eor,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='purple'
  ;cglegend, title=['Floor, all','Floor, no sidelobes', 'Floor, no nulls','Floor, no primary', 'EoR'], color=['black','blue','green','pink','purple'], location=[.55,.8], charsize=1
  ;cglegend, title=['1','2', '4','6','8','12','16', 'EoR'], color=['black','blue','green','pink','sky blue','magenta','dark green','purple'], location=[.65,.85], charsize=1
  ;cglegend, title=['std','autocal'], color=['black','blue'], location=[.25,.85], charsize=1
  ;cgoplot, k_24, power_1-power_24,/ylog,/xlog, psym=10, color='green'
  ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  stop
  
end