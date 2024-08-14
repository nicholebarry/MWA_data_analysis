pro kslice_wrapper, uvf_input=uvf_input

  obs_name= '10_11_13_14_15_16_18_19_21_22_23_24_25_26_28_29_30_31_32_33_34_35_36_5metric_pointing2';'1125953248'
  filedir='/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_redo_redo_CasA_N13_rescaled_smoothcal/' 
  cubedir='ps/data/uvf_cubes/' 

  ;power_filename1 = filedir+cubedir+'Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_fullimg_model_xx_averemove_swbh_kcube.idlsave'
  ;power_filename1 = filedir+cubedir+obs_name+'_gridded_uvf__even_odd_joint_fullimg_model_xx_averemove_swbh_power.idlsave'
  power_filename1 = filedir+cubedir+'Combined_obs_'+obs_name+'_even_cubeXX_noimgclip_res_uvf.idlsave'
;  weights_filename1 = filedir+cubedir+'Combined_obs_'+obs_name+'_even_cubeXX_fullimg_weights_uvf.idlsave'
  restore, power_filename1
;  restore, weights_filename1
  ;power_3D_1 = getvar_savefile(power_filename1, 'power_3D')
 

; metadata_struct = getvar_savefile(filedir+'ps/' +obs_name+'_gridded_uvf__even_odd_joint_uvf_info.idlsave','metadata_struct')
  metadata_struct = getvar_savefile(filedir+'ps/Combined_obs_'+obs_name+'_cubeXX__even_odd_joint_info.idlsave','metadata_struct')
  freq = metadata_struct.frequencies
  pix_area_rad = (!dtor*metadata_struct.degpix)^2.  

  z0_freq = 1420.40;E6 ;; Hz
  redshifts = z0_freq/freq - 1 ;; frequencies will be identical if kx, ky, kz match
  cosmology_measures, redshifts, wedge_factor = wedge_factor, comoving_dist_los=comov_dist_los, Ez=Ez
z_mpc_mean = z_mpc(freq*1e6)
kperp_lambda_conv = z_mpc_mean / (2.*!pi) 
 
  ;*****calculate wedge
  ;; assume 20 degrees from pointing center to first null
  source_dist = 20d * !dpi / 180d
  fov_amp = wedge_factor * source_dist
  max_theta= 10d ;random guess!
  
  ;; calculate angular distance to horizon
  horizon_amp = wedge_factor * ((max_theta+90d) * !dpi / 180d)
  data_range_set = [10^(3.),10^(8.)]

  if keyword_set(uvf_input) then begin
    data_range_set=[.1,100.] 
    restore, filedir + '1061316296_even_gridded_uvf.sav'
    n_k = (size((*model_uv_arr[0,0])))[2]
    kx_rad_vals = (FINDGEN(n_k)-n_k/2.)/metadata_struct.degpix
    ky_rad_vals = (FINDGEN(n_k/2.))/metadata_struct.degpix
  endif
 
  ;power_weights = (8*(sigma2_1)^2d)
  ;power_1_num = abs(data_sum_1)^2. - abs(data_diff_1)^2.
x_arr = KX_RAD_VALS / z_mpc_mean ;kx_mpc
y_arr = KY_RAD_VALS / z_mpc_mean;ky_mpc
data = data_cube ;power_3d
;weights = weights_cube
;x_arr = kx_mpc
;y_arr = ky_mpc
;data = abs(data_sum_1);noise_3D;weight_invert(sqrt(weights_3d));power_3D

;stop
  for kz_i=1, N_elements(freq)/2-1 do begin
  ;kz_i=8

  wedge_amp = [[fov_amp[kz_i]], [horizon_amp[kz_i]]]
  ;plotfile='/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_uvf_pskspan300/ps_150/kzslice_weights_'+string(strtrim(kz_i,2),FORMAT='(I03)')
  ;noise_ratio is power=reform(noise_3D[*,*,kz_i])*sqrt(reform(weights_3D[*,*,kz_i]))
  ;thermal noise is reform(weight_invert(sqrt(weights_3d[*,*,kz_i])))
  ;power is power=power_3d[*,*,kz_i]
  plotfile=filedir+cubedir+string(kz_i,format='(I02)')+'.png'
  if keyword_set(uvf_input) then begin
    data = abs((*model_uv_arr[0,kz_i])[*,n_k/2:n_k-1]/(*weights_uv_arr[0,kz_i])[*,n_k/2:n_k-1])
    inds_nan = where(~finite(data),/null,n_count)
    if n_count GT 0 then data[inds_nan]=0
    kpower_slice_plot, power=data, xarr=x_arr, yarr=y_arr, $
      kz_mpc = kz_mpc, /baseline_axis, /plot_wedge_line, slice_axis=2, kperp_lambda_conv = kperp_lambda_conv, wedge_amp=wedge_amp, $
      /linear_axes,charsize_in=1,window=3,plotfile=plotfile,/png
    ;data_range=data_range_set, data_min_abs=data_range_set[0] $
  endif else $
    kpower_slice_plot, power=abs(data[*,*,kz_i]), xarr=x_arr, yarr=y_arr, kz_mpc = kz_mpc, /baseline_axis, /plot_wedge_line, color_profile='abs',$
      slice_axis=2, kperp_lambda_conv = kperp_lambda_conv, wedge_amp=wedge_amp, /linear_axes,charsize_in=1,window=3,plotfile=plotfile,/png,data_range=data_range_set
  
endfor
  
end
