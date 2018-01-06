pro pfb_edge_compare

  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/beardsley_thesis_list.txt'
  readcol, filename, obsids, format='A', /silent
  
  longrun_names_match, obs_names=obs_names
  day_name_array=['Aug23','Aug27','Sep02','Sep04','Sep06','Sep09','Sep11','Sep13','Sep17','Sep19','Sep30','Oct02','Oct04','Oct08','Oct10','Oct15','Oct23','Oct25','Oct28','Oct29','Oct31','Nov17','Nov18','Nov29','Nov30']
  pointing_name_array=['minustwo','minusone','zenith','plusone','plustwo']
  pfb_set = FLTARR(N_elements(day_name_array),2,336,6)
  
  amp_degree=2
  
  obs = getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/metadata/1061316296_obs.sav','obs')
  freq_use=where((*obs.baseline_info).freq_use,nf_use)
  anti_freq_use=where((*obs.baseline_info).freq_use EQ 0,nf_use)
  tile_names=ULONG((*obs.baseline_info).tile_names)
  
  ;Using preexisting file to extract information about which tiles have which cable length
  mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  ;Taking tile information and cross-matching it with the nonflagged tiles array, resulting in nonflagged tile arrays
  ;grouped by cable length
  cable_length_ref=cable_len[Uniq(cable_len,Sort(cable_len))]
  n_cable=N_Elements(cable_length_ref)
  tile_use_arr=Ptrarr(n_cable)
  tile_arr=FINDGEN(128)
  FOR cable_i=0,n_cable-1 DO tile_use_arr[cable_i]=Ptr_new(where(cable_len EQ cable_length_ref[cable_i]))
  tile_cable_use = [(*tile_use_arr[0])[0],(*tile_use_arr[1])[0],(*tile_use_arr[2])[0],(*tile_use_arr[3])[0],(*tile_use_arr[4])[0],(*tile_use_arr[5])[0]]
  
  subband_low =  where((findgen(384) mod 16) EQ 0) + 1
  subband_high =  where((findgen(385) mod 16) EQ 0) -2
  
  ;****************************Read in and parse temperature data
  
  ;Find the smaller temperature data file with just the right obs to get temperature data.
  ;To get this file, run this code with make_obs_to_query set and then run query.sh
  filename='/nfs/eor-00/h1/nbarry/longrun_poly_bftemp/obs_queried_v2.txt'
  
  ;Read out temperature data in the form of the string due to the funky format
  textfast,data_array,/read,file_path=filename,string=1
  
  temperature_array=FLTARR(N_elements(day_name_array),128) ;Want day x tile
  
  tile_temp=STRARR(8)
  
  for day_i=0,N_elements(day_name_array)-1 do begin
    day_inds = where((obs_names[*,1] EQ day_name_array[day_i]) AND (obs_names[*,2] EQ 'zenith'),n_count)
    if n_count EQ 0 then continue
    temp_index=-1
    for inday_i=0, n_count-1 do begin
      day_index=where(strmatch(data_array,'*'+obs_names[day_inds[inday_i],0]+'*') EQ 1,n_count2)
      if n_count2 GT 0 then temp_index = day_index
    endfor
    if temp_index[0] LT 0 then continue
    split_temp_array=STRARR(N_elements(temp_index),3)
    for i=0,N_elements(temp_index)-1 do begin
      split_temp_array[i,*]=strsplit(data_array[temp_index[i]], '|',/EXTRACT)  ;Split line in data file into receiver, obsid full, and temps per tile
      split_temp_array[i,1]=(strsplit(split_temp_array[i,1], '|',/EXTRACT))[0]  ;Remove extra ms from end of obsid
      split_temp_array[i,0]=STRING(ULONG(split_temp_array[i,0])*10)  ;Add zero to end of receiver for getting tile names easier
      split_temp_array[i,2]=strmid(split_temp_array[i,2],1) ;Remove '{' from beginning of temp data
      str_length=strlen(split_temp_array[i,2]) ;Find length of string for next command
      split_temp_array[i,2]=strmid(split_temp_array[i,2],0,str_length-1) ;Remove '}' from end of temp data
      tile_temp[*]=strsplit(split_temp_array[i,2], ',',/EXTRACT) ;Split up temp data into 8 discrete temperatures per receiver
      
      
      
      for tile_i=0,7 do begin
        tile_index=where(tile_names EQ (ULONG(split_temp_array[i,0])+tile_i+1))
        tile_temp_float=Convert_To_Type(tile_temp[tile_i],4)
        temperature_array[day_i,tile_index]=tile_temp_float ;Add per tile data into array as a float
      endfor
    endfor
    
  endfor
  
  ;****************************End of read in and parse temperature data
  
  dig_jump_edges = FLTARR(2,N_elements(day_name_array),n_cable)
  digjump_amount_byday = FLTARR(2,N_elements(day_name_array),n_cable)
  edges = FLTARR(2,N_elements(day_name_array),24,n_cable)
  
  for day_i=0, N_elements(day_name_array)-1 do begin
    print, day_i
    obsid_inds = where((obs_names[*,1] EQ day_name_array[day_i]) AND (obs_names[*,2] EQ 'zenith'),n_count)
    if n_count GT 0 then obsids = reform(obs_names[obsid_inds,0]) else continue
    
    ;filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/Oct10_zenith.txt'
    ;readcol, filename, obsids, format='A', /silent
    n_obs = N_elements(obsids)
    
    
    all_subbands = FLTARR(2,n_obs,24,14,6) ;pol,obs,subband num, num of freq channels in subband, cable
    
    gain = (FLTARR(2,n_obs,384,128))
    digjump_amount = (FLTARR(2,n_obs,128))
    digjump_amount_cableavg = (FLTARR(2,n_obs,6))
    ave_gain = (FLTARR(2,n_obs,384,n_cable))
    for pol_i=0,0 do begin
      for obs_i=0, n_obs-1 do begin
        cal = getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/calibration/'+strtrim(obsids[obs_i],2)+'_cal.sav','cal')
        (*cal.gain[pol_i]) = (*cal.gain[pol_i]) + (*cal.gain_residual[pol_i])
        
        ;        for tile_j=0, N_elements(*tile_use_arr[1])-1 do begin
        ;          tile_i = (*tile_use_arr[1])[tile_j]
        ;          if (cal.mode_params[pol_i,tile_i] NE !NULL) then begin
        ;            mode_i = (*cal.mode_params[pol_i,tile_i])[0]
        ;            amp_use = (*cal.mode_params[pol_i,tile_i])[1]
        ;            phase_use = (*cal.mode_params[pol_i,tile_i])[2]
        ;            gain_mode_fit=amp_use*exp(-Complex(0,1)*2.*!Pi*(mode_i*findgen(384)/384)+Complex(0,1)*phase_use)
        ;            (*cal.gain[pol_i])[*,tile_i] -= gain_mode_fit
        ;          endif
        ;        endfor
        
        gain_fit = FLTARR(384,128)
        ;        gain_ref = abs((*cal.gain[pol_i])[freq_use,*])
        ;
        ;        for tile_i=0,127 do begin
        ;
        ;          fit_params1=poly_fit(freq_use[0:223],gain_ref[0:223,tile_i],1)
        ;          fit_params2=poly_fit(freq_use[224:335],gain_ref[224:335,tile_i],1)
        ;          gain_fit[freq_use[0]:freq_use[223],tile_i] = fit_params1[0]*findgen(freq_use[223])^0 +fit_params1[1]*findgen(freq_use[223])^1
        ;          gain_fit[freq_use[224]:freq_use[335],tile_i] = fit_params2[0]*(findgen(freq_use[335] - freq_use[224]+1) + freq_use[224])^0 +fit_params2[1]*(findgen(freq_use[335] - freq_use[224]+1) + freq_use[224])^1
        ;        endfor
        
        gain[pol_i,obs_i,*,*] = abs(*cal.gain[pol_i]); / gain_fit
        
        subband_mean_freq = FLTARR(14,6)
        for subband_i=0, 23 do begin
          for tile_i=0, 127 do begin
            fit_params_sub=poly_fit(INDGEN(14),reform(abs(gain[pol_i,obs_i,subband_low[subband_i]:subband_high[subband_i+1],tile_i])),1)
            gain_fit[subband_low[subband_i]:subband_high[subband_i+1],tile_i] = fit_params_sub[0]*findgen(14)^0 + fit_params_sub[1]*findgen(14)^1
          endfor
          for cable_i=0,5 do begin
            ;subband_mean = mean(reform(abs(gain[pol_i,obs_i,subband_low[subband_i]:subband_high[subband_i+1],*tile_use_arr[cable_i]]) / gain_fit[subband_low[subband_i]:subband_high[subband_i+1],*tile_use_arr[cable_i]] ),dim = 1)
            ;for freq_i=0, 13 do subband_mean_freq[freq_i,*] = subband_mean
            all_subbands[pol_i,obs_i,subband_i,*,cable_i]=reform(mean(abs(gain[pol_i,obs_i,subband_low[subband_i]:subband_high[subband_i+1],*tile_use_arr[cable_i]])/ gain_fit[subband_low[subband_i]:subband_high[subband_i+1],*tile_use_arr[cable_i]],dim=4,/NAN)); / subband_mean_freq
          endfor
        endfor
        
        digjump_amount[pol_i,obs_i,*]=reform(gain[pol_i,obs_i,freq_use[223],*] - gain[pol_i,obs_i,freq_use[224],*])
        
        for cable_i=0,5 do begin
          digjump_amount_cableavg[pol_i,obs_i,cable_i] = mean(reform(digjump_amount[pol_i,obs_i,*tile_use_arr[cable_i]]),/NAN)
        endfor
        
      endfor
    endfor
    
    digjump_amount_byday[*,day_i,*] = mean(reform(digjump_amount_cableavg),dim=2,/NAN)
    all_sub_obsave = mean(all_subbands,dim=2,/NAN)
    for pol_i=0,1 do begin
      for cable_i=0, n_cable - 1 do begin
        bart_plot = FLTARR(384)
        for subband_i=0, 23 do begin
          if subband_i EQ 0 then bart_plot = reform(all_sub_obsave[pol_i,subband_i,*,cable_i]) else bart_plot = [bart_plot,reform(all_sub_obsave[pol_i,subband_i,*,cable_i])]
        endfor
        pfb_set[day_i,pol_i,*,cable_i] = bart_plot
      endfor
      dig_jump_coarse = mean(reform(all_sub_obsave[pol_i,[16,17,18,19,20,21,22,23],*,*]),dim=1) ;after dig jump = [16,17,18,19,20,21,22,23] ;before dig jump = [10,11,12,13,14,15]
      dig_jump_edges[pol_i,day_i,*] = mean(dig_jump_coarse[[0,14],*], dim=1) ;choose first and last frequency channels for each cable type
    ;edges[pol_i,day_i,*,*] =  mean(reform(all_sub_obsave[pol_i,*,[0,14],*]),dim=2)
    ;digjump_amount_obsavg
    endfor
    make_plots=0
    if keyword_set(make_plots) then begin
      freq_arr = (*obs.baseline_info).freq
      cgps_open, '/nfs/eor-00/h1/nbarry/pfb/gain_residuals/'+string(strtrim(day_i,2),FORMAT='(I02)')+'pfb_'+day_name_array[day_i]+'zenith.png',/quiet,/nomatch
      cgplot, freq_arr[freq_use], pfb_set[day_i,0,*,3], yrange=[.97,1.03], xrange = [freq_arr[freq_use[0]], freq_arr[freq_use[335]]], $
        Position=[0.10, 0.85, 0.9, 0.95], XTICKFORMAT="(A1)", charsize=0.7, title = 'LMR400_320'
      cgplot, freq_arr[freq_use], pfb_set[day_i, 0,*,5], yrange=[.97,1.03], xrange = [freq_arr[freq_use[0]], freq_arr[freq_use[335]]], $
        Position=[0.10, 0.7, 0.9, 0.80],/NoErase, XTICKFORMAT="(A1)", charsize=0.7, title='LMR400_524'
      cgplot, freq_arr[freq_use], pfb_set[day_i,0,*,2], yrange=[.97,1.03], xrange = [freq_arr[freq_use[0]], freq_arr[freq_use[335]]], $
        Position=[0.10, 0.55, 0.9, 0.65],/NoErase, XTICKFORMAT="(A1)", charsize=0.7, title='RG6_230'
      cgplot, freq_arr[freq_use], pfb_set[day_i,0,*,0], yrange=[.97,1.03], xrange = [freq_arr[freq_use[0]], freq_arr[freq_use[335]]], $
        Position=[0.10, 0.4, 0.9, 0.50],/NoErase, XTICKFORMAT="(A1)", charsize=0.7, title='RG6_90'
      cgplot, freq_arr[freq_use], pfb_set[day_i,0,*,4], yrange=[.97,1.03], xrange = [freq_arr[freq_use[0]], freq_arr[freq_use[335]]], $
        Position=[0.10, 0.25, 0.9, 0.35],/NoErase, XTICKFORMAT="(A1)", charsize=0.7, title='LMR400_400'
      cgplot, freq_arr[freq_use], pfb_set[day_i,0,*,2], yrange=[.97,1.03], xrange = [freq_arr[freq_use[0]], freq_arr[freq_use[335]]], $
        Position=[0.10, 0.1, 0.9, 0.20],/NoErase, charsize=0.7, title='RG6_150', xtitle='Frequency'
      cgtext,.1,.96,day_name_array[day_i] + ' zenith pointing',/normal, charsize=1
      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    endif
  endfor
  
  ;rgbcolors = [[59,245,7],[49,219,26],[29,181,69],[11,158,107],[10,245,174],[10,196,178],[5,215,247],[6,145,191],[9,104,237],[8,89,201],$
  ;  [55,87,230],[19,38,242],[74,44,245],[108,69,247],[59,12,247],[160,78,237],[122,5,247]]
  
  rgbcolors= $
    [[177,22,51],$
    [255,33,21],$
    [201,104,0],$
    [255,223,48],$
    [212,223,0],$
    [200,248,78],$
    [211,244,126],$
    [85,161,0],$
    [53,126,0],$
    [29,132,52],$
    [109,255,137],$
    [151,253,166],$
    [0,190,101],$
    [0,100,211],$
    [0,84,238],$
    [0,56,185],$
    [132,115,255]]
  rgbcolors=reverse(rgbcolors)
  
  
  ;rgbcolors=[[69,190,207],[144,23,3],[240,87,249],[45,165,30],[42,25,77],[1,53,8],[167,138,28],[251,193,236],[52,86,193],[246,229,194],[148,33,145],[253,88,89],[239,233,110],$
  ;  [14,128,63],[96,50,14],[176,247,108],[74,126,157],[47,14,36],[212,120,250],[202,29,168]]
  color_num = [10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48]
  
  nonzero = where((dig_jump_edges[0,*,0] NE 0) AND (temperature_array[*,0] NE 0))
  for n_i=0, N_elements(nonzero)-1 do $
    TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]
    
  color_byte = [10B,12B,14B,16B,18B,20B,22B,24B,26B,28B,30B,32B,34B,36B,38B,40B,42B];,44B,46B,48B]
  temperature_array2 = temperature_array[nonzero,*]
  stop
  dig_jump_plot=1
  if keyword_set(dig_jump_plot) then begin
    cgps_open,'/nfs/eor-00/h1/nbarry/pfb/gain_residuals/temp_season_digjump.pdf',/quiet,/nomatch
    x_range= [.01,.1]
    dig_jump_edges = digjump_amount_byday
  endif else begin
    cgps_open,'/nfs/eor-00/h1/nbarry/pfb/gain_residuals/temp_season_postdig.png',/quiet,/nomatch
    
  ;x_range= [.984,.988] ;predig
    x_range= [.985,1.005] ;postdig
  endelse
  
  cgplot, dig_jump_edges[0,nonzero,0], mean(temperature_array2[*,(*tile_use_arr[0])],dim=2), xrange=x_range, yrange=[10,40], psym=16, charsize=.9, Position=[0.10, 0.70, 0.45, 0.92], /NoErase,/NoData, title='RG6_90'
  for n_i=0, N_elements(nonzero)-1 do cgplot, dig_jump_edges[0,nonzero[n_i],0], mean(temperature_array2[n_i,(*tile_use_arr[0])],dim=2), color=color_byte[n_i],$
    xrange=x_range, yrange=[10,40], psym=16, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", Position=[0.10, 0.70, 0.45, 0.920],/NoErase,xstyle=4,ystyle=4
    
  cgplot, dig_jump_edges[0,nonzero,1], mean(temperature_array2[*,(*tile_use_arr[1])],dim=2), xrange=x_range, yrange=[10,40], psym=16, charsize=.9, Position=[0.5, 0.70, 0.85, 0.92],/NoErase,/NoData,title='RG6_150'
  for n_i=0, N_elements(nonzero)-1 do cgplot, dig_jump_edges[0,nonzero[n_i],1], mean(temperature_array2[n_i,(*tile_use_arr[1])],dim=2), color=color_byte[n_i],$
    xrange=x_range, yrange=[10,40], psym=16, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", Position=[0.5, 0.70, 0.85, 0.920],/NoErase,xstyle=4,ystyle=4
    
  cgplot, dig_jump_edges[0,nonzero,2], mean(temperature_array2[*,(*tile_use_arr[2])],dim=2), xrange=x_range, yrange=[10,40], psym=16, title='RG6_230' ,charsize=.9, Position=[0.10, 0.4, 0.45, 0.620],/NoErase,/NoData
  for n_i=0, N_elements(nonzero)-1 do cgplot, dig_jump_edges[0,nonzero[n_i],2], mean(temperature_array2[n_i,(*tile_use_arr[2])],dim=2), color=color_byte[n_i],$
    xrange=x_range, yrange=[10,40], psym=16, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", Position=[0.10, 0.4, 0.45, 0.620],/NoErase,xstyle=4,ystyle=4
    
  cgplot, dig_jump_edges[0,nonzero,3], mean(temperature_array2[*,(*tile_use_arr[3])],dim=2), xrange=x_range, yrange=[10,40], psym=16, charsize=.9, Position=[0.5, 0.4, 0.85, 0.620],/NoErase,/NoData, title='LMR400_320'
  for n_i=0, N_elements(nonzero)-1 do cgplot, dig_jump_edges[0,nonzero[n_i],3], mean(temperature_array2[n_i,(*tile_use_arr[3])],dim=2), color=color_byte[n_i],$
    xrange=x_range, yrange=[10,40], psym=16, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", Position=[0.5, 0.4, 0.85, 0.620],/NoErase,xstyle=4,ystyle=4
    
  cgplot, dig_jump_edges[0,nonzero,4], mean(temperature_array2[*,(*tile_use_arr[4])],dim=2), xrange=x_range, yrange=[10,40], psym=16, charsize=.9, Position=[0.10, 0.1, 0.45, 0.320],/NoErase,/NoData, title='LMR400_400'
  for n_i=0, N_elements(nonzero)-1 do cgplot, dig_jump_edges[0,nonzero[n_i],4], mean(temperature_array2[n_i,(*tile_use_arr[4])],dim=2), color=color_byte[n_i],$
    xrange=x_range, yrange=[10,40], psym=16, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", Position=[0.10, 0.1, 0.45, 0.320],/NoErase,xstyle=4,ystyle=4
    
  cgplot, dig_jump_edges[0,nonzero,5], mean(temperature_array2[*,(*tile_use_arr[5])],dim=2), xrange=x_range, yrange=[10,40], psym=16, charsize=.9, Position=[0.5, 0.1, 0.85, 0.320],/NoErase,/NoData, title='LMR400_524'
  for n_i=0, N_elements(nonzero)-1 do cgplot, dig_jump_edges[0,nonzero[n_i],5], mean(temperature_array2[n_i,(*tile_use_arr[5])],dim=2), color=color_byte[n_i],$
    xrange=x_range, yrange=[10,40], psym=16, XTICKFORMAT="(A1)", YTICKFORMAT="(A1)", Position=[0.5, 0.1, 0.85, 0.320],/NoErase,xstyle=4,ystyle=4
    
  if keyword_set(dig_jump_plot) then cgtext, .335,.025, 'Average Digital Gain Jump',/Normal,charsize=1.2 else $
    cgtext, .17,.025, 'Average Edge Channel Gain (after digital gain jump)',/Normal,charsize=1.2
  cgtext, .05,.32, 'Average Temperature ('+cgSymbol('deg')+'C)',/Normal,charsize=1.2, orientation=90.
  
  cglegend, Title=day_name_array[nonzero], color=color_byte, location = [.885,.75], psym=16, charsize=1,Length=0.0
  
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  stop
  
end