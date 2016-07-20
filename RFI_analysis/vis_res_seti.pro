PRO vis_res_seti,save_point=save_point, quick_run=quick_run
  ; Script to plot visibility sigmas given various time interleaving or 2 second time visibilities differenced at various time steps
  ;Modified to support another outlier cut

  ;parallelizing
  paralleled = 1
  pointing=1
  quick_run=1
  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/beardsley_thesis_list.txt'
  readcol, filename, obsids, format='A', /silent
  
  if keyword_set(paralleled) then begin
    compile_opt strictarr
    args = Command_Line_Args(count=nargs)
    if keyword_set(chunk) then begin
      obs_id_chunk = args[0]
      if (100*obs_id_chunk) LT N_elements(obsids) then obsids =obsids[(100*obs_id_chunk)-100:(100*obs_id_chunk)-1] else obsids =obsids[(100*obs_id_chunk)-100:N_elements(obsids)-1]
    endif
    if keyword_set(pointing) then begin
      obs_id_pointing=args[0]
      pointing_names=['minustwo','minusone','zenith','plusone','plustwo','plusthree']
      longrun_names_match,obs_names=obs_names
      where_pointing = where(obs_names[*,2] EQ pointing_names[obs_id_pointing],n_count)
      if n_count LE 0 then message, 'n_count is zero'
      obsids = obsids[where_pointing]
    endif
  endif
  
  
  dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  
  ;if (100*obs_id_chunk) LT N_elements(obsids) then obsids =obsids[(100*obs_id_chunk)-100:(100*obs_id_chunk)-1] else obsids =obsids[(100*obs_id_chunk)-100:N_elements(obsids)-1]
  
  ;vis_files=findfile(dir+'vis_data/106*_vis_XX.sav') ; use this as a proxy to get the obsids
  ;obsids=file_basename(vis_files,'_vis_XX.sav')
  
  obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  freq_arr=(*obs.baseline_info).freq/10.^6.
  bin_offset = (*obs.baseline_info).bin_offset
  
  all_3D_20=LONG(INTARR(384,56,3000))
  all_3D_40=LONG(INTARR(384,56,3000))
  tiles_20 = LONG(INTARR(129,129,3000)) ;tile name indexed
  tiles_40 = LONG(INTARR(129,129,3000))
  
  
  tiles = (*obs.baseline_info).tile_a - (*obs.baseline_info).tile_b
  autos = where(tiles EQ 0)
  
  undefine, vis_files, tiles
  
  
  If keyword_set(save_point) then begin
  
    all_col=getvar_savefile('/nfs/eor-00/h1/nbarry/seti_all_col_40sig_'+save_point+'.sav','all_col')
    all_row=getvar_savefile('/nfs/eor-00/h1/nbarry/seti_all_row_40sig_'+save_point+'.sav','all_row')
    binned_diff=getvar_savefile('/nfs/eor-00/h1/nbarry/seti_binned_diff_40sig_'+save_point+'.sav','binned_diff')
    obs_i_start=UINT(save_point)+2
    
    IF ~keyword_set(all_col) OR ~keyword_set(all_row) OR ~keyword_set(binned_diff) then message, 'Uh oh'
    
  endif else obs_i_start=0
  
  ;****Restore the necessary information from the standard run to run this script outside of FHD.
  for obs_i=obs_i_start, N_elements(obsids)-1 do begin
    ;for obs_i=0, 3 do begin
  
    print, obs_i
    
    vis_XX = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_XX.sav', 'vis_ptr') ;restore array of calibrated visibilities
    
    vis_XX_model = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_model_XX.sav', 'vis_model_ptr') ;restore array of model visibilities
    
    vis_XX_res=*vis_XX-*vis_XX_model
    vis_XX_res[autos]=0
    
    undefine, vis_XX, vis_XX_model
    
    even_odd=1
    If keyword_set(even_odd) then begin
    
      obs = GETVAR_SAVEFILE(dir+'metadata/'+obsids[obs_i]+'_obs.sav', 'obs')
      flag_arr = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_flags.sav', 'flag_arr')
      bin_i=split_vis_flags(obs,flag_arr,bi_use=bi_use,/preserve_flags)
      
      ;bin_i = N_elements(vis_XX_res[0,*])
      ;even_bins = where(bin_i mod (2.*bin_offset[1]) EQ 0, n_even)
      ;odd_bins = where(bin_i mod (2.*bin_offset) EQ bin_offset, n_odd)
      
      ;If (n_even GT 0) and (n_odd GT 0) then begin
      vis_XX_res = vis_XX_res[*,*bi_use[0]] - vis_XX_res[*,*bi_use[1]]
      ;endif
      
      tile_a = ((*obs.baseline_info).tile_a)[*bi_use[0]] ; tile_a iformation with the even-odd formulism
      tile_b = ((*obs.baseline_info).tile_b)[*bi_use[0]]
      
    endif
    
    binsize=1.
    result=histogram(abs(vis_XX_res),binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
    ;0bs_i = 301, 681, 1393 not recorded for 40sig test
    
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
    
    If ~keyword_set(quick_run) then begin
    
      ;largest_bin=locations[N_elements(locations)-1]+.5
      ;outlier_min = 650 ; ~20sigma
      ;outlier_min = 1230 ; ~40sigma
      outlier_min_arr = [300,650] ;beginning of tail and extended part of tail
      
      for out_i=0, N_elements(outlier_min_arr) - 1 do begin
        outlier_min = outlier_min_arr[out_i]
        
        IF N_elements(result)-1 GE outlier_min then begin
          outliers=ri[ri[outlier_min]:ri[N_elements(result)]-1]
          s = SIZE(vis_XX_res)
          ncol = s[1]
          col = outliers MOD ncol
          row = outliers / ncol
          
          binsize=1
          ;outlier_result=histogram(col,binsize=binsize,locations=locations,omax=omax, /NAN)
          ;outlier_row=histogram(row,binsize=binsize,locations=locations_row,omax=omax, /NAN)
          
          
          
          ;tile_a_outliers=tile_a[row]
          ;tile_b_outliers=tile_b[row]
          his_tile_a_outliers=histogram(tile_a[row],binsize=binsize,locations=locations_tile_a,omax=omax,reverse_indices=ri_a, /NAN) ;index it from 0
          his_tile_b_outliers=histogram(tile_b[row],binsize=binsize,locations=locations_tile_b,omax=omax, /NAN)
          his_tile_a_outliers_full = LONG(INTARR(129))
          his_tile_b_outliers_full = LONG(INTARR(129))
          his_tile_a_outliers_full[locations_tile_a] = his_tile_a_outliers
          his_tile_b_outliers_full[locations_tile_b] = his_tile_b_outliers
          
          
          
          
          if keyword_set(even_odd) then mod_num=28 else mod_num=56 ;number of time samples
          
          If out_i EQ 0 then begin
            For tile_i=0, N_elements(his_tile_a_outliers) - 1 do begin
              if ri_a[tile_i] ne ri_a[tile_i+1] then tile_subset=ri_a[ri_a[tile_i]:ri_a[tile_i+1]-1] else continue
              for tile_j=0, N_elements(tile_subset)-1 do begin
                tiles_20[locations_tile_a[tile_i],tile_b[row[tile_subset[tile_j]]],obs_i] += 1
              endfor
            endfor
            For loc_i = 0, N_elements(col)-1 do begin
              all_3D_20[col[loc_i],row[loc_i] mod mod_num,obs_i]=all_3D_20[col[loc_i],row[loc_i] mod mod_num,obs_i]+uLong64(1)
            endfor
          endif else begin
            tiles_40[*,*,obs_i] = [his_tile_a_outliers_full,his_tile_b_outliers_full]
            For loc_i = 0, N_elements(col)-1 do begin
              all_3D_40[col[loc_i],row[loc_i] mod mod_num,obs_i]=all_3D_40[col[loc_i],row[loc_i] mod mod_num,obs_i]+uLong64(1)
            endfor
          endelse
          stop
        ;obs_to_fill = strarr((N_elements(all_row)))
        ;obs_to_fill[*]=obsids[obs_i]
        ;If obs_i EQ 0 then all_obs = obs_to_fill else all_obs=[obs_to_fill,all_obs]
        endif
        
      endfor
      
    Endif
    
    
    if (obs_i EQ 30) OR (obs_i EQ 100) OR (obs_i EQ 300) OR (obs_i EQ 750) OR (obs_i EQ 1000) OR (obs_i EQ 500) then begin
      x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
      y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
      
      save, binned_diff, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_binned_diff_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.sav'
      
      If ~keyword_set(quick_run) then begin
        save, all_3D_20, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_20_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.sav'
        save, all_3D_40, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_40_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.sav'
        save, tiles_20, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_tiles_20_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.sav'
        save, tiles_40, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_tiles_40_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.sav'
        
      ;cgPS_Open,'/nfs/eor-00/h1/nbarry/seti_res_vis_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.png',/quiet,/nomatch
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
    endif
    
    
  endfor
  
  
  if ~keyword_set(paralleled) then begin
    ;Save up front due to weirdness with plot making sometimes causing a stop
    save, binned_diff, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_binned_diff_thesis_evenodd2.sav'
    save, all_3D_20, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_20_thesis_evenodd2.sav'
    save, all_3D_40, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_40_thesis_evenodd2.sav'
    save, tiles_20, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_tiles_20_thesis_evenodd2.sav'
    save, tiles_40, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_tiles_40_thesis_evenodd2.sav'
    
    
    x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
    y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
    
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_res_vis_thesis_evenodd2.png',/quiet,/nomatch
    cgplot, x_arr_full,y_arr_full, xrange=[1,500], xtitle='Residual Visibility Amplitude (Jy)', ytitle='Binned Result', title='Fall 2013 Semester Residual Visibilties', charsize=1
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    
  endif else begin
  
    if keyword_set(obs_id_chunk) then begin
      save, binned_diff, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_binned_diff_thesis_evenodd2_'+strtrim(STRING(obs_id_chunk),2)+'.sav'
      save, all_3D_20, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_20_thesis_evenodd2_'+strtrim(STRING(obs_id_chunk),2)+'.sav'
      save, all_3D_40, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_40_thesis_evenodd2_'+strtrim(STRING(obs_id_chunk),2)+'.sav'
    endif
    if keyword_set(obs_id_pointing) then begin
      save, binned_diff, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_binned_diff_thesis_evenodd2_'+pointing_names[obs_id_pointing]+'.sav'
      if ~keyword_set(quick_run) then begin
        save, all_3D_20, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_20_thesis_evenodd2_'+pointing_names[obs_id_pointing]+'.sav'
        save, all_3D_40, filename='/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_40_thesis_evenodd2_'+pointing_names[obs_id_pointing]+'.sav'
      endif
    endif
  endelse
  
  stop
  
  ;binsize=1
  ;outlier_result=histogram(all_col,binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
  
  min_index = min(where(all_col NE 0))
  max_index = max(where(all_col NE 0))
  x_arr=[freq_arr[min_index],freq_arr[min_index],freq_arr[min_index:max_index]+.080*binsize/2,freq_arr[max_index]+.080*binsize,freq_arr[max_index]+.080*binsize]
  y_arr=[0,all_col[min_index], all_col[min_index:max_index], all_col[max_index],0]
  
  
  
  ;save, all_obs, filename='/nfs/eor-00/h1/nbarry/seti_all_obs.sav'
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/seti_outlier_freq_40sig.png',/quiet,/nomatch
  cgplot, x_arr, y_arr,  xtitle='Frequency (MHz)', ytitle='Binned Result', title='Fall 2013 Semester 40$\sigma$ Visibility Outliers',psym=10,charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  
  
  ;Little bit of code to unravel the parallelized save files. To be run on command line
  arr=0
  for i=1,11 do begin
    restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_binned_diff_thesis_evenodd2_'+strtrim(STRING(i),2)+'.sav'
    arr = [arr, N_elements(binned_diff)]
  endfor
  arr_length = max(arr)
  binned_diff_total = ULONG64(INTARR(arr_length))
  for i=1,11 do begin
    restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_binned_diff_thesis_evenodd2_'+strtrim(STRING(i),2)+'.sav'
    binned_diff_total[0:N_elements(binned_diff)-1] = binned_diff_total[0:N_elements(binned_diff)-1]+binned_diff
  endfor
  save, binned_diff_total, filename = '/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_binned_diff_thesis_evenodd2_total.sav'
  
  all_3D_20_total = LONG(INTARR(384,56,3000))
  all_3D_40_total = LONG(INTARR(384,56,3000))
  
  for i=1,11 do begin
    restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_20_thesis_evenodd2_'+strtrim(STRING(i),2)+'.sav'
    all_3D_20_total = all_3D_20_total + all_3D_20
    restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_40_thesis_evenodd2_'+strtrim(STRING(i),2)+'.sav'
    all_3D_40_total = all_3D_40_total + all_3D_40
  endfor
  save, all_3D_20_total, filename = '/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_20_thesis_evenodd2_total.sav'
  save, all_3D_40_total, filename = '/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_all_3D_40_thesis_evenodd2_total.sav'
  
  
end