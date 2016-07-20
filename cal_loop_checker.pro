Pro cal_loop_checker


  vis_XX_model = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_beamperchannel_unflagged/vis_data/1061316176_vis_model_XX.sav', 'vis_model_ptr') ;restore array of calibrated visibilities
  vis_YY_model = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_beamperchannel_unflagged/vis_data/1061316176_vis_model_YY.sav', 'vis_model_ptr')
  
  vis_ptr=PTRARR(2,/allocate) ;correct pol format
  *vis_ptr[0] = *vis_XX_model
  *vis_ptr[1] = *vis_YY_model
  
  vis_XX_model =GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_beamperchannel_novisflagbasic_modelnoflag_noeor_weightfix/vis_data/1061316296_vis_model_XX.sav', 'vis_model_ptr')
  vis_YY_model = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_beamperchannel_novisflagbasic_modelnoflag_noeor_weightfix/vis_data/1061316296_vis_model_YY.sav', 'vis_model_ptr')
  
  vis_model_arr=PTRARR(2,/allocate) ;correct pol format
  *vis_model_arr[0] = *vis_XX_model
  *vis_model_arr[1] = *vis_YY_model
  
  flag_ptr =GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_beamperchannel_novisflagbasic_modelnoflag_noeor_weightfix/vis_data/1061316296_flags.sav', 'flag_arr')
  
  obs =GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_beamperchannel_novisflagbasic_modelnoflag_noeor_weightfix/metadata/1061316296_obs.sav', 'obs')
  params =GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_beamperchannel_novisflagbasic_modelnoflag_noeor_weightfix/metadata/1061316296_params.sav', 'params')
  
  
  catalog_file_path=filepath('MRC_full_radio_catalog.fits',root=rootdir('FHD'),subdir='catalog_data')
  calibration_catalog_file_path=filepath('mwa_calibration_source_list.sav',root=rootdir('FHD'),subdir='catalog_data')
  
  IF Keyword_Set(calibration_catalog_file_path) THEN catalog_use=calibration_catalog_file_path
  IF ~Keyword_Set(calibration_source_list) THEN $
    calibration_source_list=generate_source_cal_list(obs,psf,catalog_path=catalog_use,_Extra=extra)
  cal=fhd_struct_init_cal(obs,params,source_list=calibration_source_list,$
    catalog_path=catalog_use,transfer_calibration=transfer_calibration,_Extra=extra)
    
    stop
    
  calibration_flag_iterate=0
  FOR iter=0,calibration_flag_iterate DO BEGIN
    t2_a=Systime(1)
    IF iter LT calibration_flag_iterate THEN preserve_flag=1 ELSE preserve_flag=preserve_visibilities
    cal=vis_calibrate_subroutine(vis_ptr,vis_model_arr,flag_ptr,obs,params,cal,$
      preserve_visibilities=0,_Extra=extra)
    t3_a=Systime(1)
    t2+=t3_a-t2_a
    
    IF Keyword_Set(flag_calibration) THEN vis_calibration_flag,obs,cal,n_tile_cut=n_tile_cut,_Extra=extra
  ;IF n_tile_cut EQ 0 THEN BREAK
  ENDFOR
  
  stop
  
End