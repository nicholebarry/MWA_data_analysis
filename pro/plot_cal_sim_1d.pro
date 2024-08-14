pro plot_cal_sim_1D

;first run commands like:
;kperp_range_lambda_1dave=[10,20] ;my goto
;ps_diff_wrapper, ['fhd_nb_sim_savedfit_cal_eor_ones_short_baselines_included_zenithpointing','fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included'], ['obs_id_6176','obs_id_6176'], data_range=[-10^10.,10^10.], data_min_abs=10.,pols='xx',cube_types='res',/invert_colorbar,plot_1d=1, kperp_range_lambda_1dave=kperp_range_lambda_1dave
;ps_diff_wrapper, ['fhd_nb_sim_cablefit_cal_eor_ones_short_baselines_included','fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included'], ['obs_id_6176','obs_id_6176'], data_range=[-10^10.,10^10.], data_min_abs=10.,pols='xx',cube_types='res',/invert_colorbar,plot_1d=1, kperp_range_lambda_1dave=kperp_range_lambda_1dave
;ps_diff_wrapper, ['fhd_nb_sim_smoothfit_cal_eor_ones_short_baselines_included','fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included'], ['obs_id_6176','obs_id_6176'], data_range=[-10^10.,10^10.], data_min_abs=10.,pols='xx',cube_types='res',/invert_colorbar,plot_1d=1, kperp_range_lambda_1dave=kperp_range_lambda_1dave
;ps_diff_wrapper, ['fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included','fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included'], ['obs_id_6176','obs_id_6176'], data_range=[-10^10.,10^10.], data_min_abs=10.,pols='xx',cube_types='res',/invert_colorbar,plot_1d=1, kperp_range_lambda_1dave=kperp_range_lambda_1dave

;ps_diff_wrapper, ['fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included_zenithpointing_1000','fhd_nb_sim_savedfit_cal_eor_ones_short_baselines_included_zenithpointing_1000'], ['obs_id_6176','obs_id_6176'], data_range=[-10^10.,10^10.], data_min_abs=10.,pols='xx',cube_types='res',plot_1d=1, kperp_range_lambda_1dave=kperp_range_lambda_1dave
;ps_diff_wrapper, ['fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included_zenithpointing_2000','fhd_nb_sim_savedfit_cal_eor_ones_short_baselines_included_zenithpointing_2000'], ['obs_id_6176','obs_id_6176'], data_range=[-10^10.,10^10.], data_min_abs=10.,pols='xx',cube_types='res',plot_1d=1, kperp_range_lambda_1dave=kperp_range_lambda_1dave
;ps_diff_wrapper, ['fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included_zenithpointing_3000','fhd_nb_sim_savedfit_cal_eor_ones_short_baselines_included_zenithpointing_3000'], ['obs_id_6176','obs_id_6176'], data_range=[-10^10.,10^10.], data_min_abs=10.,pols='xx',cube_types='res',plot_1d=1, kperp_range_lambda_1dave=kperp_range_lambda_1dave
;ps_diff_wrapper, ['fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included_zenithpointing','fhd_nb_sim_savedfit_cal_eor_ones_short_baselines_included_zenithpointing'], ['obs_id_6176','obs_id_6176'], data_range=[-10^10.,10^10.], data_min_abs=10.,pols='xx',cube_types='res',plot_1d=1, kperp_range_lambda_1dave=kperp_range_lambda_1dave
;ps_diff_wrapper, ['fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included_zenithpointing_5000','fhd_nb_sim_savedfit_cal_eor_ones_short_baselines_included_zenithpointing_5000'], ['obs_id_6176','obs_id_6176'], data_range=[-10^10.,10^10.], data_min_abs=10.,pols='xx',cube_types='res',plot_1d=1, kperp_range_lambda_1dave=kperp_range_lambda_1dave
;ps_diff_wrapper, ['fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included_zenithpointing_6000','fhd_nb_sim_savedfit_cal_eor_ones_short_baselines_included_zenithpointing_6000'], ['obs_id_6176','obs_id_6176'], data_range=[-10^10.,10^10.], data_min_abs=10.,pols='xx',cube_types='res',plot_1d=1, kperp_range_lambda_1dave=kperp_range_lambda_1dave



  ;dir = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_cal_eor_ones_short_baselines_included__'
  ;power_savefiles = [dir + 'overfit_minus_perfect/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_res_xx_dencorr_1dkpower.idlsave', $
  ;  dir + 'cablefit_minus_perfect/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_res_xx_dencorr_1dkpower.idlsave', $
  ;  dir + 'smoothfit_minus_perfect/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_res_xx_dencorr_1dkpower.idlsave', $
  ;  dir + 'savedfit_zenithpointing_minus_perfect/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_res_xx_dencorr_1dkpower.idlsave', $
  ;  '/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave']
    
      dir = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_cal_eor_ones_short_baselines_included_zenithpointing_'
  power_savefiles = [dir + '1000__perfect_minus_savedfit/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_res_xx_dencorr_1dkpower.idlsave', $
    dir + '2000__perfect_minus_savedfit/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_res_xx_dencorr_1dkpower.idlsave', $
    dir + '3000__perfect_minus_savedfit/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_res_xx_dencorr_1dkpower.idlsave', $
    dir + '_perfect_minus_savedfit/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_res_xx_dencorr_1dkpower.idlsave', $
    dir + '5000__perfect_minus_savedfit/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_res_xx_dencorr_1dkpower.idlsave', $
    dir + '6000__perfect_minus_savedfit/Combined_obs_obs_id_6176_cubeXX__even_odd_joint_res_xx_dencorr_1dkpower.idlsave', $
    '/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave']
    
 ;psyms = [intarr(5)+10];, -3]
psyms = [intarr(7)+10]
    
    ;kpower_1d_plots, power_savefiles, title = ' ', data_range = [1E0,1E8],  kpar_power=1, yaxis_type = 'sym_log', k_range=[1E-1,1.6E0],psyms=psyms, hinv=1, $
    ;  names = ['Per-frequency, per-antenna','Cable reflection fit only','Low-order polynomial fit','Time and antenna averaged', 'EoR simulation'], $
    ;  plotfile='/nfs/mwa-00/h1/nbarry/kperp_averaged_window_cal_sim1_5.pdf', pdf=1
    
   kpower_1d_plots, power_savefiles, title = ' ', data_range = [1E0,1E8],  kpar_power=1, yaxis_type = 'sym_log', k_range=[1E-1,1.6E0],psyms=psyms, hinv=1, $
      names = ['1000, averaged','2000, averaged','3000, averaged','4000, averaged','5000, averaged','6000, averaged', 'EoR simulation'], $
      plotfile='/nfs/mwa-00/h1/nbarry/kperp_averaged_window_cal_savedfits.pdf', png=1
    
    ;        kpower_1d_plots, compare_files.savefiles_1d[j,i], window_num=window_num, start_multi_params = start_multi_params, multi_pos = pos_use, $
    ;      names=compare_files.titles[j], hinv = hinv, png = png, eps = eps, pdf = pdf, plotfile = plotfiles_use, note = note_use, yaxis_type = axis_type_1d, $
    ;      title='hello world', kpar_power=1
    
end