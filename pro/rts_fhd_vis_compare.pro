pro rts_fhd_vis_compare

  rts_vis_file='/fred/oz048/MWA/CODE/FHD/RTS_compare/1061316296/1061316296.uvfits'
  uvfits_read,hdr_rts,params_rts,layout_rts,vis_rts,flag_rts,file_path_vis=rts_vis_file,n_pol=2
stop  
  n_freq=384
  fhd_vis_file='/fred/oz048/MWA/data/2013/v5_1/1061316296.uvfits'
  uvfits_read,hdr,params,layout,vis_arr,flag_arr,file_path_vis=fhd_vis_file,n_pol=2;, vis_time_average=8 ;slight modifier, n_freq was undefined at one point
  
  vis_XX = GETVAR_SAVEFILE('/fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/vis_data/1061316296_vis_XX.sav', 'vis_ptr') ;restore array of calibrated visibilities
  vis_YY = GETVAR_SAVEFILE('/fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/vis_data/1061316296_vis_YY.sav', 'vis_ptr')
  vis_arr=[vis_XX,vis_YY]
  vis_average,vis_arr,flag_arr,params,hdr,vis_time_average=4,n_freq=n_freq
  
quick_compare=1
if keyword_set(quick_compare) then begin
  flags_0 = where(*flag_arr[0] EQ 0,n_count)
  flags_1 = where(*flag_arr[1] EQ 0,n_count2)
  if n_count GT 0 then (*vis_arr[0])[flags_0] = 0 
  if n_count2 GT 0 then (*vis_arr[1])[flags_1] = 0 
  
  flags_0 = where(*flag_rts[0] EQ 0,n_count)
  flags_1 = where(*flag_rts[1] EQ 0,n_count2)
  if n_count GT 0 then (*vis_rts[0])[flags_0] = 0 
  if n_count2 GT 0 then (*vis_rts[1])[flags_1] = 0 

  obs = fhd_struct_init_obs(rts_vis_file,hdr_rts,params_rts,dimension=2048,elements=2048,n_pol=2) 
  vis_noise_calc, obs, vis_arr, flag_arr, noise_arr = noise_fhd
  fhd_noise_prt = Pointer_copy(obs.vis_noise)
  vis_noise_calc, obs, vis_rts, flag_rts, noise_arr = noise_rts 
  rts_noise_prt = Pointer_copy(obs.vis_noise)

endif
stop
  vis_diff_xx=(abs((*vis_arr[0])[368:383,*])-abs(*vis_rts[0]))/abs(*vis_rts[0])
  vis_diff_yy=(abs((*vis_arr[1])[368:383,*])-abs(*vis_rts[1]))/abs(*vis_rts[1])
  
  undefine, vis_diff_xx_long,vis_diff_yy_long
  
  for freq_i=0, (size(vis_diff_xx))[1]-1 do begin
    if keyword_set(vis_diff_xx_long) then vis_diff_xx_long=[vis_diff_xx_long,vis_diff_xx[freq_i,*]] else vis_diff_xx_long=vis_diff_xx[freq_i,*]
    if keyword_set(vis_diff_yy_long) then vis_diff_yy_long=[vis_diff_yy_long,vis_diff_yy[freq_i,*]] else vis_diff_yy_long=vis_diff_yy[freq_i,*]
  endfor
  
  ;for freq_i=0, (size(vis_diff_yy))[1]-1 do if keyword_set(vis_diff_yy_long) then vis_diff_yy_long=[vis_diff_yy_long,vis_diff_yy[freq_i,*]] else vis_diff_yy_long=vis_diff_yy[freq_i,*]
  ;for freq_i=0, (size(vis_diff_xx))[1]-1 do if keyword_set(vis_diff_xx_long) then vis_diff_xx_long=[vis_diff_xx_long,vis_diff_xx[freq_i,*]] else vis_diff_xx_long=vis_diff_xx[freq_i,*]
  
  ;The binsize is arbitrary, but it does give a pretty plot. Could be increased. Must be even.
  binsize=.05
  
  ;Histogram the uv baselines to the binsize provided above, and record that in the reverse indices for applying
  ;to the visibility difference later. Convert the uv from light speed time to wavelengths.
  ind_xx=where(vis_diff_xx_long EQ -1, xx_count)
  ind_yy=where(vis_diff_yy_long EQ -1, yy_count)
  If xx_count GT 0 then vis_diff_xx_long[ind_xx]=-2.
  If yy_count GT 0 then vis_diff_yy_long[ind_yy]=-2.
  
  
  result=histogram([vis_diff_xx_long,vis_diff_yy_long],binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
  
  ;Make the x input and y input pretty for the plot
  y_arr=[0, result , 0]
  x_arr=[locations[0],locations+binsize/2,omax]
  
  hostmean=Total(y_arr*x_arr)/(N_elements(vis_diff_xx_long)+N_Elements(vis_diff_yy_long))
  print, hostmean
  
  stop
  
  RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316296_obs.sav' ;restore obs structure
  ;bin_i=split_vis_flags_vissigma_timestep(obs,flag_arr,bi_use=bi_use,/preserve_flags,time=time_i,nt=7)
  data_diff_fhd =Imaginary( (*vis_arr[pol_i])[368:383,*bi_use[0]])-Imaginary((*vis_arr[pol_i])[368:383,*bi_use[1]])
  sigma_fhd=stddev(data_diff_fhd[*])/sqrt(2.)
  
  ;bin_i=split_vis_flags_vissigma_timestep(obs,flag_arr,bi_use=bi_use,/preserve_flags,time=time_i,nt=7)
  data_diff_rts =Imaginary( (*vis_arr_rts[pol_i])[*,*bi_use[0]])-Imaginary((*vis_arr_rts[pol_i])[*,*bi_use[1]])
  sigma_rts=stddev(data_diff_rts[*])/sqrt(2.)
  ;vis_noise_calc_sigmahistograms,vis_arr=vis_arr,flag_arr=flag_arr,nt=7
  
  stop
  
  filename='/nfs/eor-00/h1/nbarry/rts_fhd_vis_diff.png'
  cgPS_Open,filename,scale_factor=2,/quiet,/nomatch
  
  cgplot, x_arr, y_arr, psym=10,xrange=[-1,1],yrange=[0,8.*10^4.],title='Visibility amplitude fractional difference (V!Ifhd!N-V!Irts!N)/V!Irts!N', charsize=1,ytitle='density',xtitle='fractional amplitude difference'
  
  cgPS_Close,/png,Density=75,Resize=100.,/allow_transparent,/nomessage
  stop
  
  n_pol=2
  n_freq=384
  freq_start=1.9655500*10^8. ;index 368
  freq_end=1.9775500*10^8. ;index 383
  obs=fhd_struct_init_obs(vis_file_list,hdr,params,n_pol=n_pol,dft_threshold=dft_threshold,_Extra=extra)
  flag_arr=vis_flag_basic(flag_arr,obs,params,n_pol=n_pol,n_freq=n_freq,freq_start=freq_start,freq_end=freq_end,tile_flag_list=tile_flag_list,vis_ptr=vis_arr,_Extra=extra)
  stop
;vis_flag_update,flag_arr,obs,psf,params,_Extra=extra
  
end
