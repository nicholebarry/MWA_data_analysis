pro simple_thesis_plots
  cable_bp=1
  ;saved_bp=1
  if keyword_set(cable_bp) then begin
    filename= '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/calibration/1061316296_cal.sav'
    restore, filename
    filename= '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/metadata/1061316296_obs.sav'
    restore, filename
    pol_i=0
    *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
    cal_bandpass=vis_cal_bandpass(cal,obs,cal_remainder=cal_remainder,cable_bandpass_fit=1)
    if keyword_set(saved_bp) then cal_bandpass=vis_cal_bandpass(cal,obs,cal_remainder=cal_remainder,saved_run_bp=1,cable_bandpass_fit=1)
    freq_arr = (*obs.baseline_info).freq
    
    cgps_open,'/nfs/eor-00/h1/nbarry/bp_cable.pdf',/quiet,/nomatch
    cgplot, freq_arr / 1E6, (*cal_bandpass.gain[0])[*,0], psym=16, yrange=[.8,1.2], xrange = [(freq_arr[0]) / 1E6 - 1., (freq_arr[383]) / 1E6 + 1.], $
      color='royal blue', ytitle='Normalized calibration amplitude', $
      title = 'Cable bandpass for zenith observation 8/23/2013', aspect=.5, charsize=1, symsize=.5, position=[.1,.4,.9,.95]
    cgoplot, freq_arr / 1E6, (*cal_bandpass.gain[0])[*,2], psym=16,color='firebrick', symsize=.5
    cgoplot, freq_arr / 1E6, (*cal_bandpass.gain[0])[*,37], psym=16,color='forest green', symsize=.5
    cgoplot, freq_arr / 1E6, (*cal_bandpass.gain[0])[*,35], psym=16,color='dark orchid', symsize=.5
    cgoplot, freq_arr / 1E6, (*cal_bandpass.gain[0])[*,33], psym=16,color='gold', symsize=.5
    cgoplot, freq_arr / 1E6, (*cal_bandpass.gain[0])[*,36], psym=16,color='turquoise', symsize=.5
    cglegend, title = ['90m','150m','230m','320m','400m','524m'], color = ['royal blue','firebrick','forest green','dark orchid','gold','turquoise'], location = [.75,.85], charsize=1, length=0, psym=16, symsize=.5
    ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    
    gain1 = (*cal.gain[0])
    
    tile_use = (*obs.baseline_info).tile_use
    freq_use = (*obs.baseline_info).freq_use
    freq_arr = (*obs.baseline_info).freq
    n_pol=2
    n_tile = 128
    n_freq=384
    temp = where(freq_use,n_count)
    nf_use = n_count
    temp = where(tile_use,n_count)
    nt_use = n_count
    
    bandpass_arr=Fltarr(n_pol+1,n_freq)
    bandpass_arr[0,*]=freq_arr
    temp2 = FLTARR(n_pol,n_freq,n_tile)
    pol_i=0
    gain_use=extract_subarray(gain1,where(freq_use),where(tile_use))
    amp=Abs(gain_use)
    phase=Atan(gain_use,/phase)
    amp2=fltarr(nf_use,nt_use)
    
    FOR tile_i=0,nt_use-1 DO BEGIN
      resistant_mean,amp[*,tile_i],2,res_mean
      IF res_mean NE 0 THEN amp2[*,tile_i]=amp[*,tile_i]/res_mean ELSE amp2[*,tile_i]=0.
    ENDFOR
    temp2[pol_i,where(freq_use),where(tile_use)] = amp2
    bandpass_single=Fltarr(nf_use)
    FOR f_i=0L,nf_use-1 DO BEGIN
      resistant_mean,amp2[f_i,*],2,res_mean
      bandpass_single[f_i]=res_mean
    ENDFOR
    bandpass_arr[pol_i+1,where(freq_use)] = bandpass_single
    
    ;cal_bandpass_cable=vis_cal_bandpass(cal,obs,cal_remainder=cal_remainder,cable_bandpass_fit=1)
    ;bandpass_arr = FLTARR(6, 384)
    ;bandpass_arr[0,*] = (*cal_bandpass_cable.gain[0])[*,0]
    ;bandpass_arr[1,*] = (*cal_bandpass_cable.gain[0])[*,2]
    ;bandpass_arr[2,*] = (*cal_bandpass_cable.gain[0])[*,37]
    ;bandpass_arr[3,*] = (*cal_bandpass_cable.gain[0])[*,35]
    ;bandpass_arr[4,*] = (*cal_bandpass_cable.gain[0])[*,33]
    ;bandpass_arr[5,*] = (*cal_bandpass_cable.gain[0])[*,36]
    
    ;cgplot, freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,0])/bandpass_arr[1,*], /noerase, position=[.1,.05,.9,.3],xrange = [(freq_arr[0]) / 1E6 - 1., (freq_arr[383]) / 1E6 + 1.],$
    ;  psym=16,color='navy', symsize=.5, charsize=1, ytitle='% difference to global',xtitle='Frequency (MHz)'
    cgplot, freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,0])/bandpass_arr[1,*], /noerase, position=[.1,.1,.9,.35],xrange = [(freq_arr[0]) / 1E6 - 1., (freq_arr[383]) / 1E6 + 1.],$
      color='royal blue', charsize=1, ytitle='% difference to global',xtitle='Frequency (MHz)', yrange=[-8,8]
    ;cgoplot, freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,0])/bandpass_arr[1,*],color='navy'
    ;cgoplot, freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,2])/bandpass_arr[1,*],psym=16,color='firebrick', symsize=.5
    cgoplot, freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,2])/bandpass_arr[1,*],color='firebrick'
    ;cgoplot,freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,37])/bandpass_arr[1,*], psym=16,color='forest green', symsize=.5
    cgoplot, freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,37])/bandpass_arr[1,*],color='forest green'
    ;cgoplot,freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,35])/bandpass_arr[1,*], psym=16,color='dark orchid', symsize=.5
    cgoplot, freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,35])/bandpass_arr[1,*],color='dark orchid'
    ;cgoplot,freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,33])/bandpass_arr[1,*], psym=16,color='gold', symsize=.5
    cgoplot, freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,33])/bandpass_arr[1,*],color='gold'
    ;cgoplot,freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,36])/bandpass_arr[1,*],psym=16,color='teal', symsize=.5
    cgoplot, freq_arr / 1E6, 100.*(bandpass_arr[1,*] - (*cal_bandpass.gain[0])[*,36])/bandpass_arr[1,*],color='turquoise'
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
    
  endif
  
  ;global_bp=1
  ;plot_grey_background=1
  ;Bandpass figure in chapter 2 and bandpass schematic in chapter 3
  ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if keyword_set(global_bp) then begin
    obs = getvar_savefile('/nfs/mwa-05/r1/EoRuvfits/EoR2013/fhd_nb_2013longrun_savedbp/metadata/1061316296_obs.sav','obs')
    cal = getvar_savefile('/nfs/mwa-05/r1/EoRuvfits/EoR2013/fhd_nb_2013longrun_savedbp/calibration/1061316296_cal.sav','cal')
    
    gain1 = (*cal.gain[0]) + (*cal.gain_residual[0])
    gain2 = (*cal.gain[1]) + (*cal.gain_residual[1])
    
    tile_use = (*obs.baseline_info).tile_use
    freq_use = (*obs.baseline_info).freq_use
    freq_arr = (*obs.baseline_info).freq
    n_pol=2
    n_tile = 128
    n_freq=384
    temp = where(freq_use,n_count)
    nf_use = n_count
    temp = where(tile_use,n_count)
    nt_use = n_count
    
    bandpass_arr=Fltarr(n_pol+1,n_freq)
    bandpass_arr[0,*]=freq_arr
    temp2 = FLTARR(n_pol,n_freq,n_tile)
    FOR pol_i=0,n_pol-1 DO BEGIN
      gain=(*cal.gain[pol_i]) + (*cal.gain_residual[pol_i])
      gain_use=extract_subarray(gain,where(freq_use),where(tile_use))
      amp=Abs(gain_use)
      phase=Atan(gain_use,/phase)
      amp2=fltarr(nf_use,nt_use)
      
      FOR tile_i=0,nt_use-1 DO BEGIN
        resistant_mean,amp[*,tile_i],2,res_mean
        IF res_mean NE 0 THEN amp2[*,tile_i]=amp[*,tile_i]/res_mean ELSE amp2[*,tile_i]=0.
      ENDFOR
      temp2[pol_i,where(freq_use),where(tile_use)] = amp2
      bandpass_single=Fltarr(nf_use)
      FOR f_i=0L,nf_use-1 DO BEGIN
        resistant_mean,amp2[f_i,*],2,res_mean
        bandpass_single[f_i]=res_mean
      ENDFOR
      bandpass_arr[pol_i+1,where(freq_use)] = bandpass_single
      
    ENDFOR
    
    cgps_open,'/nfs/eor-00/h1/nbarry/bp_global_example.pdf',/quiet,/nomatch
    cgplot, bandpass_arr[0,*] / 1E6,temp2[0,*,0],color = 'grey' ,xrange = [bandpass_arr[0,0] / 1E6 - 1.,bandpass_arr[0,383] / 1E6 + 1.], $
      yrange=[.8,1.2], ytitle='Normalized calibration amplitude', xtitle = 'Frequency (MHz)', title = 'Global bandpass for zenith observation 8/23/2013', $
      charsize=1, aspect = .5,/nodata
    plot_grey_background=0
    if keyword_set(plot_grey_background) then begin
      for tile_i=0,nt_use-1 do cgoplot, bandpass_arr[0,*] / 1E6,temp2[0,*,tile_i],color = 'grey'
      for tile_i=0,nt_use-1 do cgoplot, bandpass_arr[0,*] / 1E6,temp2[1,*,tile_i],color = 'grey'
    endif
    cgoplot, bandpass_arr[0,*] / 1E6,bandpass_arr[1,*], color='black', xrange = [bandpass_arr[0,0] / 1E6 - 1.,bandpass_arr[0,383] / 1E6 + 1.], $
      yrange=[.8,1.2], ytitle='Normalized calibration amplitude', xtitle = 'Frequency (MHz)', title = 'Global bandpass for zenith observation 8/23/2013', $
      charsize=1, aspect = .5
    ;cgoplot,bandpass_arr[0,*] / 1E6,bandpass_arr[2,*], color='firebrick'
    ;cglegend, title = ['XX','YY'], color=['blue','firebrick'], location=[.73,.73], charsize=1
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
    
  endif
  ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  auto_corr=1
  ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if keyword_set(auto_corr) then begin
    obs = getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_autos/metadata/1061316296_obs.sav','obs')
    cal = getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_autos/calibration/1061316296_cal.sav','cal')
    
    tile_use = (*obs.baseline_info).tile_use
    freq_use = (*obs.baseline_info).freq_use
    freq_arr = (*obs.baseline_info).freq
    n_pol=2
    n_tile = 128
    n_freq=384
    temp = where(freq_use,n_count)
    nf_use = n_count
    temp = where(tile_use,n_count)
    nt_use = n_count
    
    gain=abs(*cal.gain[0])
    gain_full = abs(*cal.gain[0] + *cal.gain_residual[0])
    mean_gain = mean(gain[freq_use,*], dimension=1)
    gain[where(freq_use EQ 0),*]=0
    gain_full[where(freq_use EQ 0),*]=0
    mean_gain_full = mean(gain_full[freq_use,*], dimension=1)
    gain_use=extract_subarray(gain,where(freq_use),where(tile_use))
    amp=Abs(gain_use)
    phase=Atan(gain_use,/phase)
    amp2=fltarr(nf_use,nt_use)
    
    cgps_open,'/nfs/eor-00/h1/nbarry/auto_corr.pdf',/quiet,/nomatch
    cgplot, freq_arr / 1E6,gain_full[*,3] / mean_gain_full[3],color = 'grey' ,xrange = [freq_arr[0] / 1E6 - 1.,freq_arr[383] / 1E6 + 1.], $
      yrange=[.85,1.25], ytitle='Normalized calibration amplitude', xtitle = 'Frequency (MHz)', title = 'Bandpass (Tile 3) for zenith observation 8/23/2013', $
      charsize=1, aspect = .5
    cgoplot, freq_arr / 1E6,gain[*,3] / mean_gain[3], color='black'
    ;cgoplot,bandpass_arr[0,*] / 1E6,bandpass_arr[2,*], color='firebrick'
    cglegend, title = ['Cross','Auto'], color=['grey','black'], location=[.71,.7], charsize=1
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
    
  endif
  ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  amplitude_poly=1
  if keyword_set(amplitude_poly) then begin
    day=['Aug23','Aug27','Sep02','Sep04','Sep06','Sep09','Sep11','Sep13','Sep17','Sep19','Sep30','Oct02','Oct04','Oct08','Oct10','Oct15','Oct23','Oct25','Oct29']
    day_num=N_elements(day)
    
    pointing_num=[-2,-1,0,1,2]
    pointing_name=['-2','-1','0','1','2']
    
    parsednames=STRARR(day_num,N_elements(pointing_num))
    for day_i=0,day_num-1 do parsednames[day_i,*]=[day[day_i]+'_minustwo',day[day_i]+'_minusone',day[day_i]+'_zenith',day[day_i]+'_plusone',day[day_i]+'_plustwo']
    
    restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/1061316296_obs.sav'
    tile_names=ULONG((*obs.baseline_info).tile_names)
    
    ;***************Loop to get the obsids of each day/pointing of the longrun
    obsid_day_pointing=PTRARR((size(pointing_num))[1],day_num,/allocate)
    obsid_count=INTARR((size(pointing_num))[1])
    beginning_obsid_count=INTARR((size(pointing_num))[1])
    bftemp_pointing=FLTARR((size(pointing_num))[1],day_num)
    
    FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
      undefine, obsid
      undefine, beginning_obsid
      
      FOR day_i=0,day_num-1 DO BEGIN
      
        filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/' + parsednames[day_i,j] + '.txt'
        obs_array='empty'
        readcol, filename, obs_array, format='A', /silent
        
        If obs_array[0] NE 'empty' Then begin
          textfast,obs_array,/read,file_path=filename,string=1
          
          If keyword_set(obsid) then obsid=[obsid,ULONG(obs_array[*])] else obsid=ULONG(obs_array[*])
          If keyword_set(beginning_obsid) then beginning_obsid=[beginning_obsid,ULONG(obs_array[0])] else beginning_obsid=ULONG(obs_array[0])
          *obsid_day_pointing[j, day_i]=(obs_array[*])
          if keyword_set(make_obs_to_query) then PRINTF, lun,(*obsid_day_pointing[j,day_i])[0]
          
        endif
        
      ENDFOR ; end day for
      
      obsid_count[j]=N_elements(obsid)
      beginning_obsid_count[j]=N_elements(beginning_obsid)
      
    endfor ;end pointing for
    ;***************End of loop to get the obsids of each day/pointing of the longrun
    
    ;****************************Read in the polyfit calibration solutions, and make the raw gains minus bp sol if specified
    ;gain=complex(DBLARR(2,day_num,384,128,(size(pointing_num))[1]))
    gain=complex(DBLARR(2,day_num,384,128,(size(pointing_num))[1],20))
    
    ;for j=0,(size(pointing_num))[1]-1 do begin
    FOR j=2,3 DO BEGIN
      print, 'Reading in pointing ' + pointing_name[j]
      for day_i=0,day_num-1 do begin
      
        If *obsid_day_pointing[j, day_i] NE !NULL then begin
          for obs_i=0, N_elements(*obsid_day_pointing[j, day_i])-1 do begin
          
            filename= '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/calibration/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_cal.sav'
            if ~file_test(filename) then continue
            restore,filename
            filename= '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/metadata/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_obs.sav'
            restore,filename
            for pol_i=0,1 do *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
            
            cal_bandpass=vis_cal_bandpass(cal,obs,cal_remainder=cal_remainder,saved_run_bp=1,cable_bandpass_fit=1)
            for pol_i=0,1 do *cal.gain[pol_i]=*cal_remainder.gain[pol_i]
            
            gain[0,day_i,*,*,j,obs_i]=*cal.gain[0]
            gain[1,day_i,*,*,j,obs_i]=*cal.gain[1]
          endfor
          
        endif else begin
          gain[0,day_i,*,*,j]=-1
          gain[1,day_i,*,*,j]=-1
        endelse
        
      endfor
    endfor ;end pointing for
    freq_use = where((*obs.baseline_info).freq_use[0:255])
    
    ;****************************Read in and parse temperature data
    
    ;Find the smaller temperature data file with just the right obs to get temperature data.
    ;To get this file, run this code with make_obs_to_query set and then run query.sh
    filename='/nfs/eor-00/h1/nbarry/longrun_poly_bftemp/obs_queried_v2.txt'
    
    ;Read out temperature data in the form of the string due to the funky format
    textfast,data_array,/read,file_path=filename,string=1
    
    temperature_array=FLTARR(day_num,128,8) ;Want day x tile x pointing
    
    tile_temp=STRARR(8)
    
    ;for j=0,(size(pointing_num))[1]-1 do begin
    FOR j=2,3 DO BEGIN
      for day_i=0,day_num-1 do begin
      
        If *obsid_day_pointing[j, day_i] NE !NULL then begin
        
          temp_index=where(strmatch(data_array,'*'+(*obsid_day_pointing[j, day_i])[0]+'*') EQ 1)
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
              temperature_array[day_i,tile_index,j]=tile_temp_float ;Add per tile data into array as a float
            endfor
          endfor
          
        endif
        
      endfor
      
    endfor
    
    ;****************************End of read in and parse temperature data
    
    pol_i=0
    tile_i=2
    j=2
    day_i=0
    ;Create skeleton of plot
    cgplot, mean(abs(gain[pol_i,0,freq_use,tile_i,j])),temperature_array[day_i,tile_i,j],xtitle='Average Gain (after bandpass removal)',ytitle='Beamformer Temperature (C)', $
      charsize=1.3, psym='Filled Circle', symsize=0.5,yrange=[10,40],xrange=[1.3,1.8],/NODATA, title='Temperature dependence of gain at zenith, tile 3'
      
    non_zero_days = where(temperature_array[*,tile_i,j] NE 0,n_count)
    
    undefine, color_array
    cgLoadCT, 25, clip=[0,190], Ncolors=day_num+1
    Device, Decomposed=0
    for day_i=0,N_elements(non_zero_days)-1 do cgoplot, mean(abs(gain[pol_i,non_zero_days[day_i],freq_use,tile_i,j])),temperature_array[non_zero_days[day_i],tile_i,j], color=day_i+1,psym='Filled Circle', symsize=1
    ;TVLCT, rgbcolors[0,day_i], rgbcolors[1,day_i], rgbcolors[2,day_i],10+day_i
    
    for day_i=0,N_elements(non_zero_days)-1 do if ~keyword_set(color_array) then color_array=day_i+1 else color_array=[color_array,day_i+1]
    
    cgLegend, Title=day[non_zero_days], $
      Color=color_array, Location=[0.8,0.87],charsize=1.2,VSpace=1.4,thick=2, $
      /center_sym,psyms=16, length=0
      
    Device, Decomposed=1
    
    amp_check=1
    if keyword_set(amp_check) then begin
      pol_i=0
      tile_i=2
      j=2
      day_i=0
      
      cgplot, abs(gain[pol_i,0,*,tile_i,j]),xtitle='Frequency index',ytitle='Polynomial amplitude value', $
        charsize=1.3, symsize=0.5,yrange=[1,4],/NODATA, title='tile 2'
        
      undefine, color_array
      cgLoadCT, 25, clip=[0,190], Ncolors=day_num+1
      Device, Decomposed=0
      for day_i=0,day_num-1 do cgoplot, reform(abs(gain[pol_i,day_i,*,tile_i,j])),color=day_i+1
      
      Device, Decomposed=1
    endif
    stop
    
  endif
  
end
