pro gridded_stokes_cubes, cube_file, savefile=savefile

  restore, cube_file
  jones = fhd_struct_init_jones(obs_out,beam_model_version=2)
  n_pol=4 
  n_freq=192

  stokes_dirty_uv = PTRARR(n_pol,n_freq,/allocate)
  stokes_model_uv = PTRARR(n_pol,n_freq,/allocate)
  stokes_variance_uv = PTRARR(n_pol,n_freq,/allocate)
  stokes_weights_uv = PTRARR(n_pol,n_freq,/allocate)

  for freq_i=0, n_freq-1 do begin
    inst_dirty_image = PTRARR(n_pol,/allocate)
    inst_model_image = PTRARR(n_pol,/allocate)
    inst_variance_image = PTRARR(n_pol,/allocate)
    inst_weights_image = PTRARR(n_pol,/allocate)
    for pol_i=0,n_pol-1 do begin
      *inst_dirty_image[pol_i] = fft_shift(FFT(fft_shift(*dirty_uv_arr[pol_i,freq_i])))
      *inst_model_image[pol_i] = fft_shift(FFT(fft_shift(*model_uv_arr[pol_i,freq_i])))
      *inst_variance_image[pol_i] = fft_shift(FFT(fft_shift(*variance_uv_arr[pol_i,freq_i])))
      *inst_weights_image[pol_i] = fft_shift(FFT(fft_shift(*weights_uv_arr[pol_i,freq_i])))
    endfor
   
    stokes_dirty_image=stokes_cnv(inst_dirty_image,jones,obs_out)
    stokes_model_image=stokes_cnv(inst_model_image,jones,obs_out)
    stokes_variance_image=stokes_cnv(inst_variance_image,jones,obs_out)
    stokes_weights_image=stokes_cnv(inst_weights_image,jones,obs_out)

    for pol_i=0,n_pol-1 do begin
      *stokes_dirty_uv[pol_i,freq_i]=fft_shift(FFT(fft_shift(*stokes_dirty_image[pol_i]),/inverse))
      *stokes_model_uv[pol_i,freq_i]=fft_shift(FFT(fft_shift(*stokes_model_image[pol_i]),/inverse))
      *stokes_variance_uv[pol_i,freq_i]=fft_shift(FFT(fft_shift(*stokes_variance_image[pol_i]),/inverse))
      *stokes_weights_uv[pol_i,freq_i]=fft_shift(FFT(fft_shift(*stokes_weights_image[pol_i]),/inverse))
    endfor
  endfor

  ;rename to match expected eppsilon inputs
  dirty_uv_arr = stokes_dirty_uv
  model_uv_arr = stokes_model_uv
  variance_uv_arr = stokes_variance_uv
  weights_uv_arr = stokes_weights_uv

  save, filename=savefile, dirty_uv_arr, hpx_is_defined, model_uv_arr, obs_out, variance_uv_arr, weights_uv_arr

end
