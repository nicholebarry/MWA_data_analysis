pro sim_check
;  no_conv=1
  gauss_weights=1
  dir_name='sim_20_10th_conv_gauss'
;  regenerate=1

  if keyword_set(regenerate) then begin
    ;make flat uvf plane
    restore, '/Users/nabarry/MWA/data/fhd_nb_pipeline_paper/ps/Combined_obs_Aug23_longrunstyle_cubeXX__even_odd_joint_info.idlsave'
    restore, '/Users/nabarry/MWA/data/fhd_nb_pipeline_paper/ps/data/uvf_cubes/Combined_obs_Aug23_longrunstyle_odd_cubeXX_bh_weights_uvf.idlsave'

    u_arr = (FINDGEN(215)-107)*.5
    v_arr = (FINDGEN(215)-107)*.5 ;eppsilon expects full v plane input
    freq_arr = metadata_struct.frequencies
    ;flat uvf plane
    uv_flat = eor_sim(u_arr,v_arr,freq_arr, flat_sigma=1,kx_mpc=kx_mpc,ky_mpc=ky_mpc,kz_mpc=kz_mpc,power_3d=power_3d)
    if keyword_set(gauss_weights) then uv_flat = uv_flat * 1E9
    ;uv_flat = eor_sim(u_arr,v_arr,freq_arr,kx_mpc=kx_mpc,ky_mpc=ky_mpc,kz_mpc=kz_mpc,power_3d=power_3d)
    ;power_3d is in mK^2*Mpc^3 (I just removed those factors of h), so add h to get mK^2*Mpc^3*h^-3
    ;power_3d = power_3d * (.71)^3.

    ;noise only uvf plane of sigma 10
    ;uv_noise = randomn(seed, N_elements(u_arr), N_elements(v_arr), N_elements(freq_arr))*10.+Complex(0,1)*randomn(seed, N_elements(u_arr), N_elements(v_arr), N_elements(freq_arr))*10.
    ;uv_noise2 = randomn(seed, N_elements(u_arr), N_elements(v_arr), N_elements(freq_arr))*10.+Complex(0,1)*randomn(seed, N_elements(u_arr), N_elements(v_arr), N_elements(freq_arr))*10.

    ;create correlated pixels by convolving with a gaussian of sigma 2 pixels
    uv_flat_noise = (uv_flat); + temporary(uv_noise)
    uv_flat_noise2 = temporary(uv_flat); + temporary(uv_noise2)
    uv_flat_noise_corr = uv_flat_noise
    uv_flat_noise_corr2 = uv_flat_noise
    if ~keyword_set(no_conv) then begin
      for freq_i=0,N_elements(freq_arr)-1 do begin
        uv_flat_noise_corr[*,*,freq_i] = gauss_smooth(abs(uv_flat_noise[*,*,freq_i]), 2, /edge_truncate) * exp(Complex(0,1)*atan(uv_flat_noise[*,*,freq_i],/phase))
        uv_flat_noise_corr2[*,*,freq_i] = gauss_smooth(abs(uv_flat_noise2[*,*,freq_i]), 2,  /edge_truncate) * exp(Complex(0,1)*atan(uv_flat_noise[*,*,freq_i],/phase))
      endfor
    endif

    ;****make weights
    if keyword_set(gauss_weights) then begin
      ;35 pixel sigma 2d gauss
      gauss2d = gaussian_function([35,35],/normalize) ;211,211
      gauss_int = total(gauss2d)
      gauss2d_sized = FLTARR(215,215)
      gauss2d_sized_full = FLTARR(215,215) + 1.
      gauss2d_sized[2:212, 2:212] = (gauss2d) ;* 1E9
      gauss_ratio = total(gauss2d_sized/gauss2d_sized_full,/NAN) / total(gauss2d_sized_full,/NAN)
      gauss2d_sized = gauss2d_sized * 1E9
      uv_weights_one = gauss2d_sized;[*,107:214] * exp(Complex(0,1)*atan(weights_cube[*,*,0],/phase))
      uv_variance_one = (gauss2d_sized[*,*])^2. / 1E9
    endif else begin
      uv_weights_one = FLTARR(N_elements(u_arr),N_elements(v_arr))
      uv_variance_one = uv_weights_one
      uv_weights_one[*,*] = 1.;1E9 * exp(Complex(0,1)*atan(weights_cube[*,*,0],/phase))
      uv_variance_one[*,*] = 1.;abs(uv_weights_one)
    endelse


    uv_variance = uv_flat_noise
    uv_weights = uv_flat_noise
    ;variance should be convolved with gauss squared, but whatever for now
    if ~keyword_set(no_conv) then begin
      uv_variance_one = gauss_smooth(uv_variance_one[*,*], 2,  /edge_truncate)
      uv_weights_one = gauss_smooth(uv_weights_one[*,*], 2,  /edge_truncate)
    endif

    for freq_i=0,N_elements(freq_arr)-1 do begin
      uv_variance[*,*,freq_i] = uv_variance_one[*,*]
      uv_weights[*,*,freq_i] = uv_weights_one[*,*]
    endfor

    restore, '/Users/nabarry/obs.sav'
    for pol_i=0,1 do (*obs.vis_noise)[pol_i,*] = 1.

    PRIMARY_BEAM_AREA=1
    PRIMARY_BEAM_SQ_AREA=PTRARR(2)
    ;PRIMARY_BEAM_SQ_AREA[*]=Ptr_new((.5^2. * obs.freq_center * 2. * !pi / 3E8))
    ;PRIMARY_BEAM_SQ_AREA[*]=Ptr_new(1.2816781) ;calculated from Healpix pixel number and (!dtor*degpix)^2
    ;PRIMARY_BEAM_SQ_AREA[*]=Ptr_new(4.) ;calculated from (2.*!pi)^2./(kx_mpc_delta*ky_mpc_delta) / (z_mpc_mean)^2. ;wait, this is just the area...I want area squared
    PRIMARY_BEAM_SQ_AREA[*]=Ptr_new(4.^2.)
    if keyword_set(gauss_weights) then PRIMARY_BEAM_SQ_AREA[*]=Ptr_new((4.*gauss_ratio)^2.)
    updated_values = {PRIMARY_BEAM_SQ_AREA:PRIMARY_BEAM_SQ_AREA, PRIMARY_BEAM_AREA:PRIMARY_BEAM_AREA}
    obs = structure_update(obs,_Extra=updated_values)
    res_uv_arr = uv_flat_noise_corr
    weights_uv_arr = uv_weights
    variance_uv_arr = uv_variance
    save, res_uv_arr, weights_uv_arr, variance_uv_arr,obs,filename='/Users/nabarry/MWA/data/'+dir_name+'/1061316296_even_gridded_uvf.sav'
    res_uv_arr = uv_flat_noise_corr2
    save, res_uv_arr, weights_uv_arr, variance_uv_arr,obs,filename='/Users/nabarry/MWA/data/'+dir_name+'/1061316296_odd_gridded_uvf.sav'
  endif else begin
    restore, '/Users/nabarry/MWA/data/sim_noconv_gauss/1061316296_even_gridded_uvf.sav'
    
    if ~keyword_set(no_conv) then begin
      sigma=1.0
      for freq_i=0,191 do begin
        res_uv_arr[*,*,freq_i] = gauss_smooth(abs(res_uv_arr[*,*,freq_i]), sigma, /edge_truncate) * exp(Complex(0,1)*atan(res_uv_arr[*,*,freq_i],/phase))
        weights_uv_arr[*,*,freq_i] = gauss_smooth(weights_uv_arr[*,*,freq_i], sigma, /edge_truncate)
        variance_uv_arr[*,*,freq_i] = gauss_smooth(variance_uv_arr[*,*,freq_i], sigma, /edge_truncate)
      endfor
    endif
    
    save, res_uv_arr, weights_uv_arr, variance_uv_arr,obs,filename='/Users/nabarry/MWA/data/'+dir_name+'/1061316296_even_gridded_uvf.sav'
    save, res_uv_arr, weights_uv_arr, variance_uv_arr,obs,filename='/Users/nabarry/MWA/data/'+dir_name+'/1061316296_odd_gridded_uvf.sav'
  endelse
  ;stop
  ps_wrapper, '/Users/nabarry/MWA/data/'+dir_name+'/','1061316296', uvf_input=1,wt_cutoffs=0
  ;get 1D power spectrum from flat power
  k_edges = getvar_savefile($
    '/Users/nabarry/MWA/data/fhd_nb_pipeline_paper/ps/data/1d_binning/Combined_obs_Aug23_longrunstyle_cubeXX__even_odd_joint_bh_res_xx_averemove_swbh_dencorr_kperplambda10-50_1dkpower.idlsave'$
    ,'k_edges')

  restore, '/Users/nabarry/MWA/data/sim_check/options.sav'
  restore, '/Users/nabarry/MWA/data/sim_check/kperp_lambda_conv.sav'

  oned_out = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc,k_edges,binning_1d_options=binning_1d_options, plot_options=plot_options, kperp_lambda_conv, hubble_param = .71)
  stop
  
  ;*********some plotting pro
  power_noconv = getvar_savefile('/Users/nabarry/MWA/data/sim_noconv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  k_edges = getvar_savefile('/Users/nabarry/MWA/data/sim_noconv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','k_edges')
  power_5_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_5_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_6_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_6_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_7_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_7_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_8_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_8_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_9_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_9_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_10_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_10_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_11_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_11_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_12_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_12_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_13_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_13_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_14_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_14_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_15_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_15_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_20_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_20_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  power_30_10th_conv = getvar_savefile('/Users/nabarry/MWA/data/sim_30_10th_conv/ps/data/1d_binning/1061316296_gridded_uvf__even_odd_joint_res_xx_averemove_swbh_nodensitycorr_kperplambda10-50_1dkpower.idlsave','power')
  
  cgLoadCT, 25, clip=[0,190], Ncolors=14
  Device, Decomposed=0
  
  cgPS_Open,'/Users/nabarry/MWA/data/flat_power_corr.png',/quiet,/nomatch
  cgplot, k_edges, power_noconv, yrange=[2*10^5.,2.5*10^6.], psym=10, xtitle='k (!8h!x / Mpc)',ytitle = 'P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)',charsize=1, title='Flat Power, correlation comparison'
  cgoplot, k_edges, power_5_10th_conv, color=1, psym=10
  cgoplot, k_edges, power_6_10th_conv, color=2, psym=10
  cgoplot, k_edges, power_7_10th_conv, color=3, psym=10
  cgoplot, k_edges, power_8_10th_conv, color=4, psym=10
  cgoplot, k_edges, power_9_10th_conv, color=5, psym=10
  cgoplot, k_edges, power_10_10th_conv, color=6, psym=10
  cgoplot, k_edges, power_11_10th_conv, color=7, psym=10
  cgoplot, k_edges, power_12_10th_conv, color=8, psym=10
  cgoplot, k_edges, power_13_10th_conv, color=9, psym=10
  cgoplot, k_edges, power_14_10th_conv, color=10, psym=10
  cgoplot, k_edges, power_15_10th_conv, color=11, psym=10
  cgoplot, k_edges, power_20_10th_conv, color=12, psym=10
  cgoplot, k_edges, power_30_10th_conv, color=13, psym=10
  
  cgLegend, Title=['No correlation', '2$\sigma$=0.5','2$\sigma$=0.6','2$\sigma$=0.7','2$\sigma$=0.8','2$\sigma$=0.9',$
  '2$\sigma$=1.0', '2$\sigma$=1.1', '2$\sigma$=1.2','2$\sigma$=1.3','2$\sigma$=1.4','2$\sigma$=1.5','2$\sigma$=2.0','2$\sigma$=3.0'], $
    Color=['black',1,2,3,4,5,6,7,8,9,10,11,12,13], Location=[0.4,0.62],charsize=1.0
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    Device, Decomposed=1
    
  power_mean_noconv = mean(power_noconv,/NAN) 
  power_mean_5_10th = mean(power_5_10th_conv,/NAN) / power_mean_noconv
  power_mean_6_10th = mean(power_6_10th_conv,/NAN) / power_mean_noconv
  power_mean_7_10th = mean(power_7_10th_conv,/NAN) / power_mean_noconv
  power_mean_8_10th = mean(power_8_10th_conv,/NAN) / power_mean_noconv
  power_mean_9_10th = mean(power_9_10th_conv,/NAN) / power_mean_noconv
  power_mean_10_10th = mean(power_10_10th_conv,/NAN) / power_mean_noconv
  power_mean_11_10th = mean(power_11_10th_conv,/NAN) / power_mean_noconv
  power_mean_12_10th = mean(power_12_10th_conv,/NAN) / power_mean_noconv
  power_mean_13_10th = mean(power_13_10th_conv,/NAN) / power_mean_noconv
  power_mean_14_10th = mean(power_14_10th_conv,/NAN) / power_mean_noconv
  power_mean_15_10th = mean(power_15_10th_conv,/NAN) / power_mean_noconv
  power_mean_20_10th = mean(power_20_10th_conv,/NAN) / power_mean_noconv
  power_mean_30_10th = mean(power_30_10th_conv,/NAN) / power_mean_noconv
  power_mean_noconv_norm = power_mean_noconv / power_mean_noconv 
  
  mean_power=[power_mean_noconv_norm,power_mean_5_10th,power_mean_6_10th,power_mean_7_10th,power_mean_8_10th,power_mean_9_10th,power_mean_10_10th,$
    power_mean_11_10th,power_mean_12_10th,power_mean_13_10th,power_mean_14_10th,power_mean_15_10th,power_mean_20_10th,power_mean_30_10th]
  x_corr = [0,.5,.6,.7,.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,2.0,3.0]
  cgPS_Open,'/Users/nabarry/MWA/data/flat_power_2sigma_norm.png',/quiet,/nomatch
  cgplot, x_corr, mean_power, xtitle='2$\sigma$',ytitle = 'P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)',charsize=1, title='Flat Power, correlation comparison',yrange = [0,1] ;yrange=[1.5*10^6.,2.5*10^6.]
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
end