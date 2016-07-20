PRO vis_res_seti
  ; Script to plot visibility sigmas given various time interleaving or 2 second time visibilities differenced at various time steps

  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  
  vis_files=findfile(dir+'vis_data/106*_vis_XX.sav') ; use this as a proxy to get the obsids
  obsids=file_basename(vis_files,'_vis_XX.sav')
  
  ;****Restore the necessary information from the standard run to run this script outside of FHD.
  
  for obs_i=0, N_elements(obsids)-1 do begin
  ;for obs_i=0, 3 do begin
  
    vis_XX = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_XX.sav', 'vis_ptr') ;restore array of calibrated visibilities
    ;vis_YY= GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_YY.sav', 'vis_ptr')
    
    vis_XX_model = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_model_XX.sav', 'vis_model_ptr') ;restore array of model visibilities
    ;vis_YY_model= GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_model_YY.sav', 'vis_model_ptr')
    
    ;RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/vis_data/1061316296_flags.sav'
    ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/calibration/1061316296_cal.sav'
    ;****End of restore the necessary information from the standard run to run this script outside of FHD.
    
    
    vis_XX_res=*vis_XX-*vis_XX_model
    ;vis_YY_res=*vis_YY-*vis_YY_model
    
    ;Take the difference between neighbooring frequencies
    if keyword_set(freq_diff) then begin
      for freq_i=0,383,2 do begin
        If freq_i EQ 0 then vis_XX_freq_diff=abs(vis_XX_res[freq_i,*]-vis_XX_res[freq_i+1,*])
        vis_XX_freq_diff=[vis_XX_freq_diff,abs(vis_XX_res[freq_i,*]-vis_XX_res[freq_i+1,*])]
        
        If freq_i EQ 0 then vis_YY_freq_diff=abs(vis_YY_res[freq_i,*]-vis_YY_res[freq_i+1,*])
        vis_YY_freq_diff=[vis_YY_freq_diff,abs(vis_YY_res[freq_i,*]-vis_YY_res[freq_i+1,*])]
      endfor
      
      binsize=1.
      result=histogram(abs(vis_XX_freq_diff),binsize=binsize,locations=locations,omax=omax, /NAN)
      x_arr=[locations[0],locations+binsize/2,omax]
      y_arr=[result[0], result, result[N_elements(result)-1]]
    endif
    
    binsize=1.
    result=histogram(abs(vis_XX_res),binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
    ;x_arr=[locations[0],locations+binsize/2,omax]
    ;y_arr=[result[0], result, result[N_elements(result)-1]]
    
    IF obs_i EQ 0 then begin
      binned_diff=result
    endif else begin
      loc_size=(size(locations))[1]
      binned_size=(size(binned_diff))[1]
      If loc_size LE binned_size then begin
        binned_diff[0:loc_size-1] = binned_diff[0:loc_size-1]+result
      endif else begin
        binned_diff[0:binned_size-1]=binned_diff[0:binned_size-1]+result[0:binned_size-1]
        binned_diff=[binned_diff,result[binned_size:loc_size-1]]
      endelse
    endelse
    
    ;largest_bin=locations[N_elements(locations)-1]+.5
    outliers=ri[ri[650]:ri[N_elements(result)]-1]
    s = SIZE(vis_XX_res)
    ncol = s[1]
    col = outliers MOD ncol
    row = outliers / ncol
    
    IF obs_i EQ 0 then all_col = col else all_col=[col,all_col]
    IF obs_i EQ 0 then all_row = row else all_row=[row,all_row]
    obs_to_fill = strarr((N_elements(all_row)))
    obs_to_fill[*]=obsids[obs_i]
    If obs_i EQ 0 then all_obs = obs_to_fill else all_obs=[obs_to_fill,all_obs]
    
  endfor
  
  
  x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
  y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
  save, binned_diff, filename='/nfs/eor-00/h1/nbarry/seti_binned_diff.sav'
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/seti_res_vis.png',/quiet,/nomatch
  cgplot, x_arr_full,y_arr_full, xrange=[1,500], xtitle='Residual Visibility Amplitude', ytitle='Binned Result', title='Fall 2013 Semester Residual Visibilties'
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  
  obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  freq_arr=(*obs.baseline_info).freq/10.^6.
  binsize=1
  outlier_result=histogram(all_col,binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
  x_arr=[freq_arr[locations[0]],freq_arr[locations[0]],freq_arr[locations]+.080*binsize/2,freq_arr[omax]+.080*binsize,freq_arr[omax]+.080*binsize]
  y_arr=[0,outlier_result[0], outlier_result, outlier_result[N_elements(outlier_result)-1],0]
  
  save, all_col, filename='/nfs/eor-00/h1/nbarry/seti_all_col.sav'
  save, all_row, filename='/nfs/eor-00/h1/nbarry/seti_all_row.sav'
  save, all_obs, filename='/nfs/eor-00/h1/nbarry/seti_all_obs.sav'
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/seti_outlier_freq.png',/quiet,/nomatch
  cgplot, x_arr, y_arr,  xtitle='Frequency (MHz)', ytitle='Binned Result', title='Fall 2013 Semester 20$\sigma$ Visibility Outliers',psym=10
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
end