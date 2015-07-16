FUNCTION input_gain_plotter,cal,obs,output_file_path=output_file_path,tile_index=tile_index,tile_name=tile_name,obs_use_binary=obs_use_binary,$
    auto_input=auto_input, autoAug27=autoAug27,auto_params_Aug23=auto_params_Aug23,auto_params_Aug27=auto_params_Aug27,$
    gainAug27_forratio=gainAug27_forratio,_Extra=extra
  ;This function is the subprocedure of input_gain_plot_runner.pro. It makes panelled quick_image plots of the input gains, reconstructed in input_gain_plot_runner,
  ;both of the cross correlated and the auto correlated calculated gains.
    
  auto=auto_input
  
  ;Extract needed elements from the input structures
  gain_arr_ptr=cal.gain
  n_pol=cal.n_pol
  n_freq=cal.n_freq
  n_tile=62
  n_obs=62
  
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_std_May2015/metadata/1061316296_obs.sav'
  
  ;Find unflagged frequencies, but in an index array
  IF N_Elements(obs) GT 0 THEN freq_use=where((*obs.baseline_info).freq_use) ELSE freq_use=lindgen(n_freq)
  
  ;Build an index array from a binary array
  obs_use=where(obs_use_binary)
  freq_arr=cal.freq
  
  ;Find element numbers for loops later
  nf_use=N_Elements(freq_use)
  n_pol=N_Elements(gain_arr_ptr)
  
  ;Initialize input cross correlated gain arrays
  gain_arr_ptr2=Ptrarr(n_pol,/allocate)
  
  ;Get the cable lengths by accessing the cable reflection txt file.  This path could change, so be aware.
  mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  ;Find the cable associated with the current tile
  cable_tile=UINT(cable_len[tile_index])
  cable_tile_name=strtrim(string(cable_tile),2)
  
  
  ;Define pointing lines, stolen from bp_comparisons
  pointing_line_y=[0,336]
  ;obs_num_per_pointing=[16,15,15,15,15,18] ;Aug23 golden set
  obs_num_per_pointing=[16,15,15,15,15,18] ;Aug27 golden set
  ;obs_num_per_pointing=[13,15,14,13,6,1] ;Aug27 longrun set
  minustwo_line_x=[obs_num_per_pointing[0],obs_num_per_pointing[0]]
  minusone_line_x=[total(obs_num_per_pointing[0:1]),total(obs_num_per_pointing[0:1])]
  zenith_line_x=[total(obs_num_per_pointing[0:2]),total(obs_num_per_pointing[0:2])]
  plusone_line_x=[total(obs_num_per_pointing[0:3]),total(obs_num_per_pointing[0:3])]
  plustwo_line_x=[total(obs_num_per_pointing[0:4]),total(obs_num_per_pointing[0:4])]
  
  ;Begin pol loop for quick_image pane graph
  for pol_i=0, n_pol-1 do begin
  
    ;treat the flagged obs, make them deep blue (negative)
    ;flagged_obs=where(obs_use_binary EQ 0,n_flagged_obs)
    ;IF n_flagged_obs NE 0 then begin
    ;  input_nightave[*,flagged_obs]=-1.
    ;  input_pointing[*,flagged_obs]=-1.
    ;endif
    ;n_unflagged_obs=n_obs-n_flagged_obs
    if pol_i EQ 0 then pol_name='xx'
    if pol_i EQ 1 then pol_name='yy'
    
    
    ;make_panelled_plots=1
    If keyword_set(make_panelled_plots) then begin
    
      ;Set up multi pane parameters
      If keyword_set(autoAug27) then ncol=3 else ncol=2
      nrow=2
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
      window_num=1
      
      ;if pol_i EQ 0 then position_range=[0,1,2] else position_range=[3,4,5]
      if pol_i EQ 0 AND ~keyword_set(autoAug27) then position_range=[0,1]
      if pol_i EQ 1 AND ~keyword_set(autoAug27) then position_range=[2,3]

      if pol_i EQ 0 AND keyword_set(autoAug27) then position_range=[0,1,2]
      if pol_i EQ 1 AND keyword_set(autoAug27) then position_range=[3,4,5]
      if pol_i EQ 1 then undefine, start_multi_params
      ;save_file_name='/nfs/eor-00/h1/nbarry/Aug23_pertile/'+tile_name+'_inputgains_comparison_variscale'
      save_file_name=output_file_path+tile_name+'_fullautocompare'

      ;temp
      n_unflagged_obs=94
      
      IF pol_i EQ 0 then begin
        quick_image,transpose(abs((*gain_arr_ptr[pol_i])[freq_use,*])),multi_pos=positions,start_multi_params=start_multi_params,window_num=window_num,$
          savefile=save_file_name,png=1, data_range=data_range,$
          title=tile_name+' '+pol_name+' Aug27 g!Io!N',xtitle='Observation index',ytitle='Unflagged freqency index',$
          note=cable_tile_name+'m cable, ('+strtrim(string(n_unflagged_obs),2)+')'
        cgoplot, minustwo_line_x, pointing_line_y,Linestyle=1
        cgoplot, minusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, zenith_line_x, pointing_line_y, Linestyle=1
        cgoplot, plusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, plustwo_line_x, pointing_line_y, Linestyle=1
        
        ;scale_factor=((data_range[1]-data_range[0]))/auto_params_Aug23[0,1]
        ;scale_factor=((data_range[1]-data_range[0]))
        ;data_range_modified=[min(abs((*auto[pol_i])[freq_use,*])),min(abs((*auto[pol_i])[freq_use,*]))+scale_factor]
        data_range_modified=[min(abs((*auto[pol_i])[freq_use,*])),max(abs((*auto[pol_i])[freq_use,0:50]))]
        
        quick_image,transpose(abs((*auto[pol_i])[freq_use,*])),multi_pos=positions[*,1],window_num=window_num,$
          savefile=save_file_name,png=1,data_range=data_range_modified,$
          title=tile_name+' '+pol_name+' Aug27 unscaled g!Iauto!N',xtitle='Observation index',ytitle='Unflagged freqency index',$
          note=cable_tile_name+'m cable, ('+strtrim(string(n_unflagged_obs),2)+')'
        cgoplot, minustwo_line_x, pointing_line_y,Linestyle=1
        cgoplot, minusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, zenith_line_x, pointing_line_y, Linestyle=1
        cgoplot, plusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, plustwo_line_x, pointing_line_y, Linestyle=1
        
        
        If keyword_set(autoAug27) then begin
          quick_image,transpose(abs((*autoAug27[pol_i])[freq_use,*])),multi_pos=positions[*,2],window_num=window_num,$
            savefile=save_file_name,png=1,data_range=data_range_modified,$
            title=tile_name+' '+pol_name+' Aug27 unscaled g!Iauto!N',xtitle='Observation index',ytitle='Unflagged freqency index',$
            note=cable_tile_name+'m cable, ('+strtrim(string(n_unflagged_obs),2)+')'
          cgoplot, minustwo_line_x, pointing_line_y,Linestyle=1
          cgoplot, minusone_line_x, pointing_line_y, Linestyle=1
          cgoplot, zenith_line_x, pointing_line_y, Linestyle=1
          cgoplot, plusone_line_x, pointing_line_y, Linestyle=1
          cgoplot, plustwo_line_x, pointing_line_y, Linestyle=1
        endif
        undefine, data_range
      endif else begin
        quick_image,transpose(abs((*gain_arr_ptr[pol_i])[freq_use,*])),multi_pos=positions[*,position_range[0]],window_num=window_num,$
          savefile=save_file_name,png=1,data_range=data_range,$
          title=tile_name+' '+pol_name+' Aug27 g!Io!N',xtitle='Observation index',ytitle='Unflagged freqency index'
        cgoplot, minustwo_line_x, pointing_line_y,Linestyle=1
        cgoplot, minusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, zenith_line_x, pointing_line_y, Linestyle=1
        cgoplot, plusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, plustwo_line_x, pointing_line_y, Linestyle=1
        
        data_range_modified=[min(abs((*auto[pol_i])[freq_use,*])),max(abs((*auto[pol_i])[freq_use,0:50]))]
        
        quick_image,transpose(abs((*auto[pol_i])[freq_use,*])),multi_pos=positions[*,position_range[1]],window_num=window_num,$
          savefile=save_file_name,png=1,data_range=data_range_modified,$
          title=tile_name+' '+pol_name+' Aug27 unscaled g!Iauto!N',xtitle='Observation index',ytitle='Unflagged freqency index'
        cgoplot, minustwo_line_x, pointing_line_y,Linestyle=1
        cgoplot, minusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, zenith_line_x, pointing_line_y, Linestyle=1
        cgoplot, plusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, plustwo_line_x, pointing_line_y, Linestyle=1
        
        If keyword_set(autoAug27) then begin
          quick_image,transpose(abs((*autoAug27[pol_i])[freq_use,*])),multi_pos=positions[*,position_range[2]],window_num=window_num,$
            savefile=save_file_name,png=1,data_range=data_range_modified,$
            title=tile_name+' '+pol_name+' Aug27 unscaled g!Iauto!N',xtitle='Observation index',ytitle='Unflagged freqency index',$
            note=cable_tile_name+'m cable, ('+strtrim(string(n_unflagged_obs),2)+')'
          cgoplot, minustwo_line_x, pointing_line_y,Linestyle=1
          cgoplot, minusone_line_x, pointing_line_y, Linestyle=1
          cgoplot, zenith_line_x, pointing_line_y, Linestyle=1
          cgoplot, plusone_line_x, pointing_line_y, Linestyle=1
          cgoplot, plustwo_line_x, pointing_line_y, Linestyle=1
        endif
        
      endelse
      
      If pol_i EQ 1 then cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
    ENDIF
    
    ;make_gif_plots=1
    IF keyword_set(make_gif_plots) then begin
      save_file_name=output_file_path+tile_name+'_obs'
      
      for obs_i=0,93 do begin
        ;cgPS_Open,save_file_path+strtrim(string(tile_i, FORMAT='(I03)'),1)+'.png',/quiet,/nomatch
        cgplot,((*auto[0])[*,obs_i]), yrange=[min(abs((*auto[pol_i])[freq_use,*])),max(abs((*auto[pol_i])[freq_use,0:50]))], $
          title='Aug23 unscaled auto gains, xx',ytitle='Unscaled gain value',xtitle='Unflagged freqency index'
          stop
        ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      endfor
      
    endif
    
    
    
    make_ratio_plots=1
    If keyword_set(make_ratio_plots) then begin
    
      ;Set up multi pane parameters
      If keyword_set(gainAug27_forratio) then ncol=2 else ncol=1
      nrow=2
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
      window_num=1
      
      ;if pol_i EQ 0 then position_range=[0,1,2] else position_range=[3,4,5]
      if pol_i EQ 0 AND ~keyword_set(gainAug27_forratio) then position_range=[0]
      if pol_i EQ 1 AND ~keyword_set(gainAug27_forratio) then position_range=[1]

      if pol_i EQ 0 AND keyword_set(gainAug27_forratio) then position_range=[0,1]
      if pol_i EQ 1 AND keyword_set(gainAug27_forratio) then position_range=[2,3]
      if pol_i EQ 1 then undefine, start_multi_params
      ;save_file_name='/nfs/eor-00/h1/nbarry/Aug23_pertile/'+tile_name+'_inputgains_comparison_variscale'
      save_file_name=output_file_path+tile_name+'_ratio'

      ;temp
      n_unflagged_obs=94
      
      IF pol_i EQ 0 then begin
        quick_image,transpose(abs((*gain_arr_ptr[pol_i])[freq_use,*]))/transpose(abs((*auto[pol_i])[freq_use,*])),multi_pos=positions,start_multi_params=start_multi_params,window_num=window_num,$
          savefile=save_file_name,png=1,$
          title=tile_name+' '+pol_name+' Aug23 g!Io!N/g!Iauto!N unscaled',xtitle='Observation index',ytitle='Unflagged freqency index',$
          note=cable_tile_name+'m cable, ('+strtrim(string(n_unflagged_obs),2)+')'
        cgoplot, minustwo_line_x, pointing_line_y,Linestyle=1
        cgoplot, minusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, zenith_line_x, pointing_line_y, Linestyle=1
        cgoplot, plusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, plustwo_line_x, pointing_line_y, Linestyle=1       
        
        If keyword_set(gainAug27_forratio) then begin
          quick_image,transpose(abs((*gainAug27_forratio[pol_i])[freq_use,*]))/transpose(abs((*autoAug27[pol_i])[freq_use,*])),multi_pos=positions[*,1],window_num=window_num,$
            savefile=save_file_name,png=1,$
            title=tile_name+' '+pol_name+' Aug27 g!Io!N/g!Iauto!N unscaled',xtitle='Observation index',ytitle='Unflagged freqency index',$
            note=cable_tile_name+'m cable, ('+strtrim(string(n_unflagged_obs),2)+')'
          cgoplot, minustwo_line_x, pointing_line_y,Linestyle=1
          cgoplot, minusone_line_x, pointing_line_y, Linestyle=1
          cgoplot, zenith_line_x, pointing_line_y, Linestyle=1
          cgoplot, plusone_line_x, pointing_line_y, Linestyle=1
          cgoplot, plustwo_line_x, pointing_line_y, Linestyle=1
        endif

      endif else begin
        quick_image,transpose(abs((*gain_arr_ptr[pol_i])[freq_use,*]))/transpose(abs((*auto[pol_i])[freq_use,*])),multi_pos=positions[*,position_range[0]],window_num=window_num,$
          savefile=save_file_name,png=1,$
          title=tile_name+' '+pol_name+' Aug23 g!Io!N/g!Iauto!N unscaled',xtitle='Observation index',ytitle='Unflagged freqency index'
        cgoplot, minustwo_line_x, pointing_line_y,Linestyle=1
        cgoplot, minusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, zenith_line_x, pointing_line_y, Linestyle=1
        cgoplot, plusone_line_x, pointing_line_y, Linestyle=1
        cgoplot, plustwo_line_x, pointing_line_y, Linestyle=1
        
        If keyword_set(autoAug27) then begin
          quick_image,transpose(abs((*gainAug27_forratio[pol_i])[freq_use,*]))/transpose(abs((*autoAug27[pol_i])[freq_use,*])),multi_pos=positions[*,position_range[1]],window_num=window_num,$
            savefile=save_file_name,png=1,$
            title=tile_name+' '+pol_name+' Aug27 g!Io!N/g!Iauto!N unscaled',xtitle='Observation index',ytitle='Unflagged freqency index',$
            note=cable_tile_name+'m cable, ('+strtrim(string(n_unflagged_obs),2)+')'
          cgoplot, minustwo_line_x, pointing_line_y,Linestyle=1
          cgoplot, minusone_line_x, pointing_line_y, Linestyle=1
          cgoplot, zenith_line_x, pointing_line_y, Linestyle=1
          cgoplot, plusone_line_x, pointing_line_y, Linestyle=1
          cgoplot, plustwo_line_x, pointing_line_y, Linestyle=1
        endif
        
      endelse
      
      If pol_i EQ 1 then cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
    ENDIF
    
    
    
    
  endfor
  
END
