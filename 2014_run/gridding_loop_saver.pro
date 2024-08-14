pro gridding_loop_saver

  ;When the gridding_loop is called with an obs ID passed via a bash script, 
  ; parse command line args to get the obs ID. Otherwise, run with a hardcoded
  ; obs ID for testing purposes.
  compile_opt strictarr
  args = Command_Line_Args(count=nargs)
  IF keyword_set(args) then begin
    obs_id = args[0]
  endif else begin
    obs_id= '1088544584';'1088716176' ;A minustwo pointing obs ID
  endelse

  ;Starting memory and time for usage statistics
  t0=Systime(1)
  start_mem = memory(/current)

  ;Define the output directory (with obs ID appended)
  file_path_fhd = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_pyfhd/'+obs_id
  file_path_vis = '/fred/oz048/MWA/data/2014/van_vleck_corrected/coarse_corr_no_ao/' + obs_id +'.uvfits'

  ;Define the directories to the calibrated data visibilities, model visibilities, and metadata required for gridding.
  ;Also define the directory to the hyperresolved beam (>10Gb for single gauss, >40Gb for full beam)
  ;Each pointing requires a different beam, hence the hardcoded pointing grouping. Can be coded in for future runs
  data_dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_calstop/'
  cal_dir = 'calibration/'
  meta_dir = 'metadata/'
  beam_dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam/beams/'
  model_transfer = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_calstop/woden_models/combined/'
  pointing = 'zenith';'minustwo'
  pointing_num = '0';'-2'

  ;Load the raw visibilities and generate the PARAMS structure
  n_pol = 2 ;don't need cross-pols
  uvfits_read,hdr,params,layout,vis_arr,vis_weights,file_path_vis=file_path_vis,n_pol=n_pol

  ;Reform the params based on the appropriate metadata and expected format (obs will be overwritten with the saved
  ; run, but the params will be modified accordingly)

  obs=fhd_struct_init_obs(file_path_vis,hdr,params,n_pol=n_pol,dimension=2048) 

  ;Load the PSF structure (>10Gb, structures in file: PSF, OBS). Note that this is a representative OBS,
  ; and will need to be overwritten with the correct OBS
  restore, beam_dir + 'gauss_beam_pointing'+pointing_num+'.sav'

  ;Load the OBS structure (check to make sure this saved obs struct is the same as what is req for gridding)
  restore, data_dir + meta_dir + obs_id + '_obs.sav'

  ;Run basic flagging on the raw visibilities based on MWA instrumental parameters. 
  vis_weights=vis_flag_basic(vis_weights,obs,params,n_pol=n_pol,n_freq=n_freq,freq_start=freq_start,$
    freq_end=freq_end,tile_flag_list=tile_flag_list,vis_ptr=vis_arr,_Extra=extra)
  vis_weights_update,vis_weights,obs,psf,params,_Extra=extra

  ;Load the calibration CAL structure
  restore, data_dir + cal_dir + obs_id + '_cal.sav'
  ;Apply the calibration solutions to the raw visibilities
  vis_cal=vis_calibration_apply(vis_arr,cal)

  ;Load the model visibilities
  vis_model_arr = vis_model_transfer(obs,params,model_transfer)

  ;Calculate the expected noise given the visibilities and store in the OBS struct
  vis_noise_calc,obs,vis_arr,vis_weights

  ;Save metadata, calibrated vis, and vis weights to output directory
  ;fhd_save_io,status_str,obs,var='obs',/compress,file_path_fhd=file_path_fhd,_Extra=extra
  ;fhd_save_io,status_str,params,var='params',/compress,file_path_fhd=file_path_fhd,_Extra=extra
  ;vis_export,obs,status_str,vis_arr,vis_weights,file_path_fhd=file_path_fhd,/compress

  ;Check that not all data as been flagged
  IF obs.n_vis EQ 0 THEN BEGIN
    print,"All data flagged! Returning."
    error=1
    IF Keyword_Set(!Journal) THEN Journal ;write and close log file if present
    RETURN
  ENDIF

  ;Some hardcoded variables
  rephase_weights=1
  n_avg=2
  ps_kbinsize=0.5
  kbinsize=ps_kbinsize
  FoV_use=!RaDeg/kbinsize
  ps_kspan=200
  dimension_use=ps_kspan/kbinsize
  restrict_hpx_inds='EoR0_high_healpix_inds_3x.idlsave'
  cube_name=['hpx_even','hpx_odd']
  n_freq=obs.n_freq
  ps_nfreq_avg=n_freq
  n_freq_use=Floor(n_freq/n_avg)
  nside_use=1024.
  obs_out=fhd_struct_update_obs(obs,n_pol=n_pol,beam_nfreq_avg=ps_nfreq_avg,FoV=FoV_use,dimension=dimension_use)
  beam_arr=beam_image_cube(obs_out,psf,n_freq=n_freq_use,beam_mask=beam_mask,/square,beam_threshold=beam_threshold)
  hpx_cnv=healpix_cnv_generate(obs_out,file_path_fhd=file_path_fhd,nside=nside_use,restore_last=0,/no_save,$
    mask=beam_mask,hpx_radius=hpx_radius,restrict_hpx_inds=restrict_hpx_inds,_Extra=extra)
  nside=nside_use
  hpx_inds=hpx_cnv.inds
  n_hpx=N_Elements(hpx_inds)
  vis_weights_use=Pointer_copy(vis_weights)
  vis_weights_update,vis_weights_use,obs_out,psf,params,_Extra=extra
  n_iter=2
  vis_weights_use=split_vis_weights(obs_out,vis_weights_use,bi_use=bi_use,_Extra=extra)
  vis_noise_calc,obs_out,vis_arr,vis_weights_use,bi_use=bi_use
  uvf_name = ['even','odd']
  residual_flag=0
  model_flag=1
  dirty_flag=1

    pol_names=obs.pol_names
  vis_filepath=file_path_fhd+'_vis_'
    n_freq=obs.n_freq
  n_pol=obs.n_pol
  dimension=obs.dimension
  degpix=obs.degpix
  residual_flag=obs.residual
   freq_bin_i2=Floor(lindgen(n_freq)/n_avg)
   nf=Max(freq_bin_i2)+1L
   model_return=1
   rephase_use=phase_shift_uv_image(obs_out,/to_orig_phase)
   init_arr=Fltarr(dimension,dimension)
      variance_flag=1 ;initialize
      weights_flag=1 ;initialize
      preserve_visibilities=1
stop
  ;Save metadata, calibrated vis, and vis weights to output directory
  fhd_save_io,status_str,obs_out,var='obs',/compress,file_path_fhd=file_path_fhd,_Extra=extra
  fhd_save_io,status_str,params,var='params',/compress,file_path_fhd=file_path_fhd,_Extra=extra
  vis_export,obs,status_str,vis_arr,vis_weights,file_path_fhd=file_path_fhd,/compress
  vis_export,obs,status_str,vis_model_arr,vis_weights,file_path_fhd=file_path_fhd,/compress,model=1
  fhd_save_io,status_str,vis_weights,var='vis_weights',/compress,file_path_fhd=file_path_fhd,_Extra=extra
  save, weights_flag, variance_flag, model_return, bi_use, preserve_visibilities, filename = file_path_fhd + '_variables.sav'

  ;    vis_ptr,vis_weights_use[pol_i],obs_out,0,psf_out,params,timing=t_grid0,fi_use=fi_use,bi_use=bi_use,$
  ;          polarization=pol_i,weights=weights_holo,variance=variance_holo,silent=1,mapfn_recalculate=0,$
  ;          model_ptr=model_ptr,n_vis=n_vis,/preserve_visibilities,model_return=model_return, _Extra=extra

  ;healpix_snapshot_cube_generate,obs,status_str,psf,cal,params,vis_arr,$
  ;  vis_model_arr=vis_model_arr,file_path_fhd=file_path_fhd,vis_weights=vis_weights,cmd_args=cmd_args,$
  ;  n_avg=2,ps_kbinsize=0.5,ps_kspan=200,restrict_hpx_inds='EoR0_high_healpix_inds_3x.idlsave',_Extra=extra

  ;print mem and time used
  run_report, start_mem, t0, silent=silent


end