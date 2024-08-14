pro visibility_subtract

file_path='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_noeor_ones_maxcalsources_nod_zenithpointing_notileflag_5outof10_2_psf200/'
file_path_alternate='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_noeor_ones_maxcalsources_nod_zenithpointing_notileflag_5outof10_2_psf1000/'


  ;so you certainly could use visibility_grid (*_wrap also handles saving files etc), then FFT the gridded uv to get an image
  ;also, I suggest using dirty_image_generate for your FFTs
  ;even without any normalization factors, it will handle the proper shifts and such, to keep images consistent

  ;vis_ptr=PTRARR(2,/allocate)
  ;vis_model_ptr=PTRARR(2,/allocate)
  vis_residual=PTRARR(2,/allocate)
  
  ;restore XX vis
  ;vis_model_ptr = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_noeor_ones_maxcalsources_nod_zenithpointing_notileflag_5outof10/vis_data/1061316176_vis_model_XX.sav','vis_model_ptr')
  ;vis_ptr = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_noeor_ones_maxcalsources_nod_zenithpointing_notileflag_5outof10/vis_data/1061316176_vis_XX.sav','vis_ptr')
  
  ;make residual visbilities
  ;vis = *vis_ptr
  ;vis_model = *vis_model_ptr
  ;*vis_residual[0] = vis - vis_model
  
  ;undefine, vis, vis_model, vis_ptr, vis_model_ptr
  
  ;restore YY vis
  vis_model_ptr = getvar_savefile(file_path+'vis_data/1061316176_vis_model_XX.sav','vis_model_ptr')
  vis_model_ptr_dim = getvar_savefile(file_path_alternate+'vis_data/1061316176_vis_model_XX.sav','vis_model_ptr')
  vis_ptr = getvar_savefile(file_path+'vis_data/1061316176_vis_XX.sav','vis_ptr')
  
  ;make residual visbilities
  vis = *vis_ptr
  vis_model = *vis_model_ptr
  vis_model_dim = *vis_model_ptr_dim
  ;*vis_residual[0] = vis - vis_model - vis_model_dim
    *vis_residual[0] = vis_model - vis_model_dim
  
  
  undefine, vis, vis_model, vis_ptr, vis_model_ptr, vis_model_ptr_dim, vis_model_dim
  
    vis_model_ptr = getvar_savefile(file_path+'vis_data/1061316176_vis_model_YY.sav','vis_model_ptr')
  vis_model_ptr_dim = getvar_savefile(file_path_alternate+'vis_data/1061316176_vis_model_YY.sav','vis_model_ptr')
  vis_ptr = getvar_savefile(file_path+'vis_data/1061316176_vis_YY.sav','vis_ptr')
    ;make residual visbilities
  vis = *vis_ptr
  vis_model = *vis_model_ptr
  vis_model_dim = *vis_model_ptr_dim
  ;*vis_residual[0] = vis - vis_model - vis_model_dim
    *vis_residual[1] = vis_model - vis_model_dim
      undefine, vis, vis_model, vis_ptr, vis_model_ptr, vis_model_ptr_dim, vis_model_dim
    
  
  ;restore flags
  vis_flags = getvar_savefile(file_path+'vis_data/1061316176_flags.sav','flag_arr')
  
  ;make residual visbilities
  ;vis = *vis_ptr
  ;vis_model = *vis_model_ptr
  ;*vis_residual = vis - vis_model
  
  ; *vis_residual[1] = *vis_ptr[1] - *vis_model_ptr[1]
  
  ;get obs structure
  obs = getvar_savefile(file_path+'metadata/1061316176_obs.sav','obs')
  status = getvar_savefile(file_path+'metadata/1061316176_status.sav','status_str')
  psf = getvar_savefile(file_path+'beams/1061316176_beams.sav','psf')
  params = getvar_savefile(file_path+'metadata/1061316176_params.sav','params')
  weights_XX = getvar_savefile(file_path+'grid_data/1061316176_uv_weights_XX.sav','weights_uv')
  weights_YY = getvar_savefile(file_path+'grid_data/1061316176_uv_weights_YY.sav','weights_uv')
  weights = PTRARR(2,/allocate)
  *weights[0]=weights_XX
  *weights[1]=weights_YY
  
  uv_residual_xx = visibility_grid(vis_residual[0], vis_flags[0], obs, status, psf, params, polarization=0)
    uv_residual_yy = visibility_grid(vis_residual[1], vis_flags[1], obs, status, psf, params, polarization=1)
  
  
  
  ;it would have to be constructed
  ;BUT, for a zenith observation, the value you want from it is 1
  ;so, you can just define a n_pol pointer array
  ;with a fltarr(dimension,elements)+1 for each element
  dimension=Float(obs.dimension)
  elements=Float(obs.elements)
  beam_base=PTRARR(2,/allocate)
  *beam_base[0] = fltarr(dimension,elements)+1
  *beam_base[1] = fltarr(dimension,elements)+1
  
  normalization_factor = get_image_renormalization(obs, weights=weights, degpix=obs.degpix,image_filter_fn='filter_uv_uniform',beam_base=beam_base)
  
  image_residual_xx = dirty_image_generate(uv_residual_xx, normalization=normalization_factor,image_filter_fn='filter_uv_uniform')
    image_residual_yy = dirty_image_generate(uv_residual_yy, normalization=normalization_factor,image_filter_fn='filter_uv_uniform')
  
  stop
  quick_image, abs(image_residual[712:1336,712:1336])
  stop
  
end