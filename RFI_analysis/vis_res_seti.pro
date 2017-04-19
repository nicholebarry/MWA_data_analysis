PRO vis_res_seti,save_point=save_point, quick_run=quick_run, even_odd=even_odd
  ; Script to plot visibility sigmas given various time interleaving or 2 second time visibilities differenced at various time steps
  ;Modified to support another outlier cut
  ;Script to find statistics regarding visibility residuals.
  ;Will return a variety of outputs. Will return the binned visibility residuals over season 1 from the thesis cut
  ;Can make various/multiple sigma cuts on the data

outdir = '/nfs/eor-00/h1/nbarry/vis_res/Oct23EoR1/'

  ;Things that can be set in the code below to make grid engine compatibility:
  ;parallelizing for Grid Engine. Unset if testing code on head node
;  paralleled=1
  ;pointing
  ;chunk=1
  quick_run=0
  ;save_point
  even_odd=1
  ;hot_baselines=1
  ;outlier_min_arr = [300,650] ;beginning of tail and extended part of tail for even-odd
  ;outlier_min_arr = [650,1230] ; ~20 and ~40 sigma for residual visibilities
  outlier_min_arr = [300,600] ;
  
  ;Looks like I can either run the parallelization by chunking up the obsids sequentially or by pointing
  ;pointing=1
  ;chunk=1
  
  ;Modify the number of time samples if even-odd
  if keyword_set(even_odd) then mod_num=28 else mod_num=56 ;number of time samples
  
  print, "Running visibility statistics"
  
  ;Read in the obsids
  ;filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/beardsley_thesis_list.txt'
  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/EoR1/Oct23.txt
  readcol, filename, obsids, format='A', /silent
  
  if keyword_set(paralleled) then begin
    print, "Begin parallelization"
    compile_opt strictarr
    args = Command_Line_Args(count=nargs)
    if keyword_set(chunk) then begin
      obs_id_chunk = args[0]
      if (100*ULONG(obs_id_chunk)) LT N_elements(obsids) then obsids =obsids[(100*ULONG(obs_id_chunk))-100:(100*ULONG(obs_id_chunk))-1] $
      else obsids =obsids[(100*ULONG(obs_id_chunk))-100:N_elements(obsids)-1]
      print, "Will analyze obsids from " + strtrim(obsids[0],2) + " to " + strtrim(obsids[N_elements(obsids)-1],2) + $
        ". There are " + strtrim(N_elements(obsids),2) + " obsids in this chunk."
    endif
    if keyword_set(pointing) then begin
      obs_id_pointing=args[0]
      pointing_names=['minustwo','minusone','zenith','plusone','plustwo','plusthree']
      longrun_names_match,obs_names=obs_names
      where_pointing = where(obs_names[*,2] EQ pointing_names[obs_id_pointing],n_count)
      if n_count LE 0 then message, 'n_count is zero'
      obsids = obsids[where_pointing]
      print, "Will analyze obsids as a function of pointing."
    endif
  endif
  
  ;Setup dir, a sample obs structure, frequency array, bin offset
  ;dir='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/'
  dir='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_presidelobe_Oct23_EoR1/'
  obs = GETVAR_SAVEFILE(dir+'metadata/1066586152_settings.txt','obs')
  ;obs = GETVAR_SAVEFILE(dir+'metadata/1061316296_obs.sav', 'obs')
  freq_arr=(*obs.baseline_info).freq/10.^6.
  bin_offset = (*obs.baseline_info).bin_offset
  
  ;if (100*obs_id_chunk) LT N_elements(obsids) then obsids =obsids[(100*obs_id_chunk)-100:(100*obs_id_chunk)-1] else obsids =obsids[(100*obs_id_chunk)-100:N_elements(obsids)-1]
  
  ;vis_files=findfile(dir+'vis_data/106*_vis_XX.sav') ; use this as a proxy to get the obsids
  ;obsids=file_basename(vis_files,'_vis_XX.sav')
  
  all_3D_20=LONG(INTARR(384,56,3000))
  all_3D_40=LONG(INTARR(384,56,3000))
  tiles_20 = LONG(INTARR(129,129,3000)) ;tile name indexed
  tiles_40 = LONG(INTARR(129,129,3000))
  
  ;Find where the autos are in the visibility data
  tiles = (*obs.baseline_info).tile_a - (*obs.baseline_info).tile_b
  autos = where(tiles EQ 0)
  
  undefine, vis_files, tiles
  
  ;Read-in from a save point if set
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
    
    ;Restore array of calibrated and model visibilities to generate residual visibilities (and make sure the autos are 0!)
    vis_XX = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_XX.sav', 'vis_ptr')
    vis_XX_model = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_vis_model_XX.sav', 'vis_model_ptr')
    vis_XX_res=*vis_XX-*vis_XX_model
    vis_XX_res[autos]=0
    
    ;option to select only hot baselines seen from a previous run
    if keyword_set(hot_baselines) then begin
      params = GETVAR_SAVEFILE(dir+'metadata/1061316296_params.sav', 'params')
      hot = getvar_savefile('/nfs/eor-00/h1/nbarry/hot_baselines.sav','hot')
      hot_col = hot mod 129
      hot_row = hot / 129
      vis_mask = FLTARR(N_elements(vis_XX))
      for hot_i=0, N_elements(hot)-1 do begin
        where_row = where(params.antenna1 EQ hot_row[hot_i])
        where_col = where(params.antenna2[where_row] EQ hot_col[hot_i])
        vis_mask[params.antenna2[where_row[where_col]]] = 1
      endfor
      vis_XX_res = vis_XX_res * vis_mask
      undefine, vis_mask
    endif
    
    undefine, vis_XX, vis_XX_model
    
    ;Take the difference between even and odd time samples in the visibilities to remove sky signal
    If N_elements(even_odd) EQ 0 then even_odd=1
    If keyword_set(even_odd) then begin
    
      ;Replace sample obs and flag_arr with the real versions for even-odd run
      obs = GETVAR_SAVEFILE(dir+'metadata/'+obsids[obs_i]+'_obs.sav', 'obs')
      flag_arr = GETVAR_SAVEFILE(dir+'vis_data/'+obsids[obs_i]+'_flags.sav', 'vis_weights')
      
      ;Run FHD util to find the right bin indicies for even and odd samples, then subtract
      bin_i=split_vis_weights(obs,flag_arr,bi_use=bi_use,/preserve_flags)
      vis_XX_res = vis_XX_res[*,*bi_use[0]] - vis_XX_res[*,*bi_use[1]]
      
      ;Bins refer to time samples, so an arbitrary one will work to pick out the right tiles
      tile_a = ((*obs.baseline_info).tile_a)[*bi_use[0]] ; tile_a information with the even-odd formulism
      tile_b = ((*obs.baseline_info).tile_b)[*bi_use[0]]
      
    endif
    
    ;Set the binsize and histogram the visibilities (at this point, either residuals or even-odd residuals)
    binsize=1.
    result=histogram(abs(vis_XX_res),binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
    
    ;binned_diff is the binned visibility residuals over all obsids, and result is the binned visibility residuals for this particular obsid
    ;BUILT-IN ASSUMPTION: There are at least some in the 0th bin. It's a fantastic assumption; most are 0.
    ;If this is the first obsid to be binned, set that equal to the all-obs binning.
    IF obs_i EQ obs_i_start then begin
      binned_diff=uLong64(result)
    endif else begin
    
      ;Get the amount of bins of the new obsid and compare that to the amount of bins of the total
      loc_size=(size(locations))[1]
      binned_size=(size(binned_diff))[1]
      
      ;If the amount of new bins is less than the total, than it is fully encompasses by the total bin amount. Straight addition can be done
      If loc_size LE binned_size then begin
        binned_diff[0:loc_size-1] = binned_diff[0:loc_size-1]+uLong64(result)
        
      ;Otherwise, do the straight addition and then tack on the leftover amount of bins to the end of the total array
      endif else begin
        binned_diff[0:binned_size-1]=binned_diff[0:binned_size-1]+uLong64(result[0:binned_size-1])
        binned_diff=[binned_diff,uLong64(result[binned_size:loc_size-1])]
      endelse
      
    endelse ;end else for obs=0 or not
    
    ;A quick run will skip any processing other than binning the residual visibilities. Otherwise, enter loop to begin
    ;extra processing
    If ~keyword_set(quick_run) then begin
    
      ;Sort processing by various "sigma" cuts
      for out_i=0, N_elements(outlier_min_arr) - 1 do begin
        outlier_min = outlier_min_arr[out_i]
        undefine, locations_tile_a, locations_tile_b, ri_a
        
        ;If their are more bins for this obsid than that of the specified cut, there are visibilities to process
        ;BUILT-IN ASSUMPTION: There are at least some in the 0th bin. It's a fantastic assumption; most are 0.
        IF N_elements(result)-1 GE outlier_min then begin
          ;Grab were the outliers occured in the visibility data using reverse indicies, and get their corresponding row and column numbers
          outliers=ri[ri[outlier_min]:ri[N_elements(result)]-1]
          s = SIZE(vis_XX_res)
          ncol = s[1]
          col = outliers MOD ncol
          row = outliers / ncol
          
          binsize=1
          
          ;dependent on even-odd
          
          If keyword_set(even_odd) then begin
            ;Histograms the outliers with tile information. The tile names corresponding to the outlier exist in locations
            his_tile_a_outliers=histogram(tile_a[row],binsize=binsize,locations=locations_tile_a,omax=omax,reverse_indices=ri_a, /NAN) ;index it from 0
            his_tile_b_outliers=histogram(tile_b[row],binsize=binsize,locations=locations_tile_b,omax=omax, /NAN)
            
            ;At each of the tile names, insert how many times it contributed to the outliers
            his_tile_a_outliers_full = LONG(INTARR(129))
            his_tile_b_outliers_full = LONG(INTARR(129))
            his_tile_a_outliers_full[locations_tile_a] = his_tile_a_outliers
            his_tile_b_outliers_full[locations_tile_b] = his_tile_b_outliers
            
            ;For each tile in tile_a's outlier histogram
            For tile_i=0, N_elements(his_tile_a_outliers) - 1 do begin
            
              ;If the i-vector does not change between the bins, then no outliers were present so continue the for loop
              ;otherwise, pick out the row values for that tile_a
              if ri_a[tile_i] ne ri_a[tile_i+1] then tile_subset=ri_a[ri_a[tile_i]:ri_a[tile_i+1]-1] else continue
              
              ;For every row value for that tile_a (which is every tile_b), pick out the name that corresponds to tile a
              ;and the name that corresponds to tile b, and add 1 to that bin
              ;OOPS, locations_tile_a is indexed from 0, tile_b is indexed from 0. Saved data in July 2016 is off by 1 bin
              If out_i EQ 0 then begin
                for tile_j=0, N_elements(tile_subset)-1 do tiles_20[locations_tile_a[tile_i],tile_b[row[tile_subset[tile_j]]],obs_i] += 1
              endif else for tile_j=0, N_elements(tile_subset)-1 do tiles_40[locations_tile_a[tile_i],tile_b[row[tile_subset[tile_j]]],obs_i] += 1
              
            endfor
            
          endif
          
          ;Find the times of ....
          If out_i EQ 0 then begin
            For loc_i = 0, N_elements(col)-1 do all_3D_20[col[loc_i],row[loc_i] mod mod_num,obs_i]=all_3D_20[col[loc_i],row[loc_i] mod mod_num,obs_i]+uLong64(1)
          endif else For loc_i = 0, N_elements(col)-1 do all_3D_40[col[loc_i],row[loc_i] mod mod_num,obs_i]=all_3D_40[col[loc_i],row[loc_i] mod mod_num,obs_i]+uLong64(1)
          
        endif
        
      endfor
      
    Endif
    
    If ~keyword_set(obs_id_chunk) then begin
      if (obs_i EQ 30) OR (obs_i EQ 100) OR (obs_i EQ 300) OR (obs_i EQ 750) OR (obs_i EQ 1000) OR (obs_i EQ 500) then begin
        x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
        y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
        
        save, binned_diff, filename=outdir + 'seti_binned_diff_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.sav'
        
        If ~keyword_set(quick_run) then begin
          save, all_3D_20, filename=outdir + 'seti_all_3D_20_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.sav'
          save, all_3D_40, filename=outdir + 'seti_all_3D_40_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.sav'
          save, tiles_20, filename=outdir + 'seti_tiles_20_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.sav'
          save, tiles_40, filename=outdir + 'seti_tiles_40_thesis_evenodd2_'+strtrim(STRING(obs_i),2)+'.sav'
          
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
    endif
    
  endfor
  
  
  if ~keyword_set(paralleled) then begin
    ;Save up front due to weirdness with plot making sometimes causing a stop
    save, binned_diff, filename=outdir + 'seti_binned_diff_thesis_evenodd2.sav'
    save, all_3D_20, filename=outdir + 'seti_all_3D_20_thesis_evenodd2.sav'
    save, all_3D_40, filename=outdir + 'seti_all_3D_40_thesis_evenodd2.sav'
    save, tiles_20, filename=outdir + 'seti_tiles_20_thesis_evenodd2.sav'
    save, tiles_40, filename=outdir + 'seti_tiles_40_thesis_evenodd2.sav'
    
    
    x_arr_full=[0,FINDGEN(N_elements(binned_diff))+.5,N_elements(binned_diff)-.5]
    y_arr_full=[binned_diff[0],binned_diff,binned_diff[N_elements(binned_diff)-1]]
    
    
    cgPS_Open,outdir + 'seti_res_vis_thesis_evenodd2.png',/quiet,/nomatch
    cgplot, x_arr_full,y_arr_full, xrange=[1,500], xtitle='Residual Visibility Amplitude (Jy)', ytitle='Binned Result', title='Fall 2013 Semester Residual Visibilties', charsize=1
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    
  endif else begin
  
    if keyword_set(obs_id_chunk) then begin
      save, binned_diff, filename=outdir + 'seti_binned_diff_thesis_residual_'+strtrim(STRING(obs_id_chunk),2)+'.sav'
      save, all_3D_20, filename=outdir + 'seti_all_3D_20_thesis_residual_'+strtrim(STRING(obs_id_chunk),2)+'.sav'
      save, all_3D_40, filename=outdir + 'seti_all_3D_40_thesis_residual_'+strtrim(STRING(obs_id_chunk),2)+'.sav'
      save, tiles_20, filename=outdir + 'seti_tiles20_thesis_residual_'+strtrim(STRING(obs_id_chunk),2)+'.sav'
      save, tiles_40, filename=outdir + 'seti_tiles40_thesis_residual_'+strtrim(STRING(obs_id_chunk),2)+'.sav'
    endif
    if keyword_set(obs_id_pointing) then begin
      save, binned_diff, filename=outdir + 'seti_binned_diff_thesis_evenodd2_'+pointing_names[obs_id_pointing]+'.sav'
      if ~keyword_set(quick_run) then begin
        save, all_3D_20, filename=outdir + 'seti_all_3D_20_thesis_evenodd2_'+pointing_names[obs_id_pointing]+'.sav'
        save, all_3D_40, filename=outdir + 'seti_all_3D_40_thesis_evenodd2_'+pointing_names[obs_id_pointing]+'.sav'
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
  
  ;*********For even-odd, paralleled, time analysis
  all_3D_20_total = LONG(INTARR(384,56,3000))
  all_3D_40_total = LONG(INTARR(384,56,3000))
  
  for i=1,11 do begin
    restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/hot_baselines/seti_all_3D_20_thesis_evenodd2_'+strtrim(STRING(i),2)+'.sav'
    all_3D_20_total[*,*,(i-1)*100:(i*100-1)] = all_3D_20[*,*,0:99]
    restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/hot_baselines/seti_all_3D_40_thesis_evenodd2_'+strtrim(STRING(i),2)+'.sav'
    all_3D_40_total[*,*,(i-1)*100:(i*100-1)] = all_3D_40[*,*,0:99]
  endfor
  save, all_3D_20_total, filename = '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/hot_baselines/seti_all_3D_20_thesis_evenodd2_total.sav'
  save, all_3D_40_total, filename = '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/hot_baselines/seti_all_3D_40_thesis_evenodd2_total.sav'
  
  ;***********For even-odd, paralleled, antenna analysis
  tiles20_total = LONG(INTARR(129,129,3000))
  tiles40_total = LONG(INTARR(129,129,3000))
  
  for i=1,11 do begin
    restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/seti_tiles20_thesis_evenodd2_'+strtrim(STRING(i),2)+'.sav'
    tiles20_total[*,*,(i-1)*100:(i*100-1)] = tiles_20[*,*,0:99]
    restore, '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/seti_tiles40_thesis_evenodd2_'+strtrim(STRING(i),2)+'.sav'
    tiles40_total[*,*,(i-1)*100:(i*100-1)] = tiles_40[*,*,0:99]
  endfor
  save, tiles20_total, filename = '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/seti_tiles20_thesis_evenodd2_total.sav'
  save, tiles40_total, filename = '/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/seti_tiles40_thesis_evenodd2_total.sav'
  
  
end