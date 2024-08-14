PRO vis_res_seti_2D,save_point=save_point
  ; Script to plot visibility sigmas given various time interleaving or 2 second time visibilities differenced at various time steps

save_point = '650'

  ;Would like pointing information, LST information, row info

  ;Long run directory
  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  
  ;Use this as a proxy to get the obsids that were run
  vis_files=findfile(dir+'vis_data/106*_vis_XX.sav')
  obsids=file_basename(vis_files,'_vis_XX.sav')
  
  ;Cut obsids from long run, 32 hours
  ;filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/beardsley_thesis_list.txt'
  ;readcol, filename, obsids, format='A', /silent
  
  ;Need the frequency array, which is unchanging
  obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  freq_arr=(*obs.baseline_info).freq/10.^6.
  
  ;Set up histogram arrays
  all_col=LONG(INTARR(384))
  all_row=LONG(INTARR(500000))
  all_2D_sec=LONG(INTARR(384,56))
  all_obs=LONG(INTARR(3000))
  
  ;Autos already flagged, but make sure for safety
  tiles = (*obs.baseline_info).tile_a - (*obs.baseline_info).tile_b
  autos = where(tiles EQ 0)
  
  undefine, tiles
  
  ;Read in histograms from a save point
  If keyword_set(save_point) then begin
  
    all_col=getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/seti_all_col_40sig_'+save_point+'.sav','all_col')
    all_row=getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/seti_all_row_40sig_'+save_point+'.sav','all_row')
    all_2D_sec=getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/seti_all_2D_sec_40sig_'+save_point+'.sav','all_2D_sec')
    all_obs=getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/seti_all_obs_40sig_'+save_point+'.sav','all_obs')
    binned_diff=getvar_savefile('/nfs/eor-00/h1/nbarry/vis_res/seti_binned_diff_40sig_'+save_point+'.sav','binned_diff')
    obs_i_start=UINT(save_point)+1
    
    IF ~keyword_set(all_col) OR ~keyword_set(all_row) OR ~keyword_set(binned_diff) then message, 'Uh oh'
    
  endif else obs_i_start=0
  
  ;****Restore the necessary information from the standard run to run this script outside of FHD.
  for obs_i=obs_i_start, N_elements(obsids)-1 do begin
    If obs_i EQ 1393 then continue
    print, obs_i
    ;for obs_i=0, 3 do begin
    
    vis_XX = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_XX.sav', 'vis_ptr') ;restore array of calibrated visibilities
    
    vis_XX_model = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_model_XX.sav', 'vis_model_ptr') ;restore array of model visibilities
    
    vis_XX_res=*vis_XX-*vis_XX_model
    vis_XX_res[autos]=0
    
    undefine, vis_XX, vis_XX_model
    
    ;Split up histogram by seconds
    ;(*obs.baseline_info).tile_a
    
    ;Number of seconds
    ;num_baselines = 8255
    ;num_sec=(size(vis_XX_res))[2]/8255

    ;For time_i=1, (size(vis_XX_res))[1]/((*obs.baseline_info).tile_a)[1] - 1
    
    binsize=1.
    result=histogram(abs(vis_XX_res),binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
    ;0bs_i = 301, 681, 1393 not recorded for 40sig test
    ;302, 681, 1393 not recorded for 2D test
    
    IF obs_i EQ 0 then begin
      binned_diff=uLong64(result)
    endif else begin
      loc_size=(size(locations))[1]
      binned_size=(size(binned_diff))[1]
      If loc_size LE binned_size then begin
        binned_diff[0:loc_size-1] = binned_diff[0:loc_size-1]+uLong64(result)
      endif else begin
        binned_diff[0:binned_size-1]=binned_diff[0:binned_size-1]+uLong64(result[0:binned_size-1])
        binned_diff=[binned_diff,uLong64(result[binned_size:loc_size-1])]
      endelse
    endelse
    
    ;largest_bin=locations[N_elements(locations)-1]+.5
    outlier_min = 650 ; ~20sigma
    ;outlier_min = 1230 ; ~40sigma
    ;outlier_min=2
    
    IF N_elements(result)-1 GE outlier_min then begin
      outliers=ri[ri[outlier_min]:ri[N_elements(result)]-1]
      s = SIZE(vis_XX_res)
      ncol = s[1]
      col = outliers MOD ncol
      row = outliers / ncol

      ;binsize=1
      
      ;outlier_result=histogram(col,binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
      ;outlier_row=histogram(row,binsize=binsize,locations=locations_row,omax=omax, /NAN,reverse_indices=ri)
      
      For loc_i = 0, N_elements(col)-1 do begin
        all_col[col[loc_i]]=all_col[col[loc_i]]+uLong64(1)
        all_row[row[loc_i]]=all_row[row[loc_i]]+uLong64(1)
        all_2D_sec[col[loc_i],row[loc_i] mod 56]=all_2D_sec[col[loc_i],row[loc_i] mod 56]+uLong64(1)
      endfor
      
      all_obs[obs_i]=N_elements(col)
      

    endif
    
    if (obs_i EQ 5) then stop
    
    if (obs_i EQ 650) OR (obs_i EQ 500) OR (obs_i EQ 750) OR (obs_i EQ 1000) OR (obs_i EQ 1250) OR (obs_i EQ 2000) then begin
      ;x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
      ;y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
      
      save, binned_diff, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_binned_diff_40sig_'+strtrim(STRING(obs_i),2)+'.sav'
      save, all_col, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_all_col_40sig_'+strtrim(STRING(obs_i),2)+'.sav'
      save, all_row, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_all_row_40sig_'+strtrim(STRING(obs_i),2)+'.sav'
      save, all_2D_sec, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_all_2D_sec_40sig_'+strtrim(STRING(obs_i),2)+'.sav'
      save, all_obs, filename='/nfs/eor-00/h1/nbarry/vis_res/seti_all_obs_40sig_'+strtrim(STRING(obs_i),2)+'.sav'
      
      ;cgPS_Open,'/nfs/eor-00/h1/nbarry/seti_res_vis_40sig_'+strtrim(STRING(obs_i),2)+'.png',/quiet,/nomatch
      ;cgplot, x_arr_full,y_arr_full, xrange=[1,500], xtitle='Residual Visibility Amplitude (Jy)', ytitle='Binned Result', title='Fall 2013 (Obs '+strtrim(STRING(obs_i),2)+') Semester Residual Visibilties', charsize=1
      ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
      ;binsize=1
      ;outlier_result=histogram(all_col,binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
      
      ;min_index = min(where(all_col NE 0))
      ;max_index = max(where(all_col NE 0))
      ;x_arr=[freq_arr[min_index],freq_arr[min_index],freq_arr[min_index:max_index]+.080*binsize/2,freq_arr[max_index]+.080*binsize,freq_arr[max_index]+.080*binsize]
      ;y_arr=[0,all_col[min_index], all_col[min_index:max_index], all_col[max_index],0]
      
      
      
      ;save, all_obs, filename='/nfs/eor-00/h1/nbarry/seti_all_obs.sav'
      
      ;cgPS_Open,'/nfs/eor-00/h1/nbarry/seti_outlier_freq_40sig_'+strtrim(STRING(obs_i),2)+'.png',/quiet,/nomatch
      ;cgplot, x_arr, y_arr,  xtitle='Frequency (MHz)', ytitle='Binned Result', title='Fall 2013 (Obs '+strtrim(STRING(obs_i),2)+') Semester 40$\sigma$ Visibility Outliers',psym=10,charsize=1
      ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    endif
    
  endfor
  
  
  x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
  y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
  save, binned_diff, filename='/nfs/eor-00/h1/nbarry/seti_binned_diff_40sig.sav'
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/seti_res_vis_40sig.png',/quiet,/nomatch
  cgplot, x_arr_full,y_arr_full, xrange=[1,500], xtitle='Residual Visibility Amplitude (Jy)', ytitle='Binned Result', title='Fall 2013 Semester Residual Visibilties', charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  ;binsize=1
  ;outlier_result=histogram(all_col,binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
  
  min_index = min(where(all_col NE 0))
  max_index = max(where(all_col NE 0))
  x_arr=[freq_arr[min_index],freq_arr[min_index],freq_arr[min_index:max_index]+.080*binsize/2,freq_arr[max_index]+.080*binsize,freq_arr[max_index]+.080*binsize]
  y_arr=[0,all_col[min_index], all_col[min_index:max_index], all_col[max_index],0]
  
  
  save, all_col, filename='/nfs/eor-00/h1/nbarry/seti_all_col_40sig.sav'
  save, all_row, filename='/nfs/eor-00/h1/nbarry/seti_all_row_40sig.sav'
  ;save, all_obs, filename='/nfs/eor-00/h1/nbarry/seti_all_obs.sav'
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/seti_outlier_freq_40sig.png',/quiet,/nomatch
  cgplot, x_arr, y_arr,  xtitle='Frequency (MHz)', ytitle='Binned Result', title='Fall 2013 Semester 40$\sigma$ Visibility Outliers',psym=10,charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
end