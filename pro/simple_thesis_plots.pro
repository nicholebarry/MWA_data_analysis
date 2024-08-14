pro simple_thesis_plots

  ;cable_bp=1
  ;saved_bp=1
  ;tile_bp=1
  ;auto_bp=1
  ;zeroed_bp=1
  if keyword_set(cable_bp) then begin
    filename= '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/calibration/1061316296_cal.sav'
    restore, filename
    filename= '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/metadata/1061316296_obs.sav'
    restore, filename
    pol_i=0
    *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
    cal_bandpass=vis_cal_bandpass(cal,obs,cal_remainder=cal_remainder,cable_bandpass_fit=1)
    if keyword_set(saved_bp) then cal_bandpass=vis_cal_bandpass(cal,obs,cal_remainder=cal_remainder,saved_run_bp=1,cable_bandpass_fit=1)
    if keyword_set(tile_bp) then begin
      restore, '/nfs/mwa-10/r1/EoRuvfits/analysis/ave_cals/full_ave_amp.sav' ;'amp_mean_obs3' 384, 5, 2, 128
      cal_bandpass = Pointer_copy(cal)
      ;cal_bandpass.gain[0] = reform(amp_mean_obs3[*,2,0,*]) ;choose zenith pointing, XX
      ;cal_bandpass.gain[1] = reform(amp_mean_obs3[*,2,1,*]) ;choose zenith pointing, XX
      
      cal.amp_degree = 4
      cal.mode_fit=0
      cal_polyfit = vis_cal_polyfit(cal,obs,digital_gain_jump_polyfit=1)
      amp_mean = FLTARR(384,128)
      freq_use=where((*obs.baseline_info).freq_use)
      
      for tile_i=0,127 do begin
        resistant_mean,abs((*cal.gain[0])[freq_use,tile_i]),2,res_mean
        IF res_mean NE 0 THEN amp_mean[freq_use,tile_i]=res_mean ELSE amp_mean[freq_use,tile_i]=0.
      endfor
      
      *cal_bandpass.gain[0] = abs(*cal_polyfit.gain[0])*abs(reform(amp_mean_obs3[*,2,0,*])) / amp_mean
      zeros = (where(amp_mean EQ 0))
      (*cal_bandpass.gain[0])[zeros]=0.
    ;*cal_bandpass.gain[1] = abs(*cal_polyfit.gain[1])*abs(reform(amp_mean_obs3[*,2,1,*]))*exp(Complex(0,1)* atan((*cal.gain[1]),/phase) )
    endif
    if keyword_set(zeroed_bp) then begin
      restore,'/nfs/mwa-10/r1/EoRuvfits/analysis/ave_cals/full_ave_amp_zeroed.sav' ;'amp_mean_obs3' 384, 5, 2, 128
      cal_bandpass = Pointer_copy(cal)
      ;cal_bandpass.gain[0] = reform(amp_mean_obs3[*,2,0,*]) ;choose zenith pointing, XX
      ;cal_bandpass.gain[1] = reform(amp_mean_obs3[*,2,1,*]) ;choose zenith pointing, XX
      
      cal.amp_degree = 4
      cal.mode_fit=0
      cal_polyfit = vis_cal_polyfit(cal,obs,digital_gain_jump_polyfit=1)
      amp_mean = FLTARR(384,128)
      freq_use=where((*obs.baseline_info).freq_use)
      
      for tile_i=0,127 do begin
        resistant_mean,abs((*cal.gain[0])[freq_use,tile_i]),2,res_mean
        IF res_mean NE 0 THEN amp_mean[freq_use,tile_i]=res_mean ELSE amp_mean[freq_use,tile_i]=0.
      endfor
      
      *cal_bandpass.gain[0] = abs(*cal_polyfit.gain[0])*abs(reform(amp_mean_obs3[*,2,0,*])) / amp_mean
      zeros = (where(amp_mean EQ 0))
      (*cal_bandpass.gain[0])[zeros]=0.
      
    endif
    if keyword_set(auto_bp) then begin
      restore, '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_Aug2017_autocal_wo_cable_w_digjump/calibration/1061316296_cal.sav'
      cal_bandpass = Pointer_copy(cal)
      for tile_i=0,127 do begin
        ;      amp_params = *cal.amp_params[0,tile_i]
        ;      fit = FLTARR(384)
        ;      fit[0:255] = amp_params[0,0] + amp_params[0,1]*FINDGEN(256)
        ;      fit[256:383] = amp_params[1,0] + amp_params[1,1]*(FINDGEN(128)+256.)\
        (*cal_bandpass.gain[0])[*,tile_i] = abs((*cal_bandpass.gain[0])[*,tile_i]) / mean(abs((*cal.gain[0])[*,tile_i]))
        zeros = where(reform((*obs.baseline_info).freq_use EQ 0))
        (*cal_bandpass.gain[0])[zeros,tile_i]=0.
      endfor
      
    endif
    
    freq_arr = (*obs.baseline_info).freq
    
    ;cgps_open,'/nfs/eor-00/h1/nbarry/bp_auto.pdf',/quiet,/nomatch
    cgplot, freq_arr / 1E6, (*cal_bandpass.gain[0])[*,0], psym=16, yrange=[.8,1.2], xrange = [(freq_arr[0]) / 1E6 - 1., (freq_arr[383]) / 1E6 + 1.], $
      color='royal blue', ytitle='Normalized calibration amplitude', $
      title = 'Auto bandpass for zenith observation 8/23/2013', aspect=.5, charsize=1, symsize=.5, position=[.1,.4,.9,.95]
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
    ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
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
  ;auto_corr=1
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
  
  
  ;amplitude_poly=1
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
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;oned_plot=1
  if keyword_set(oned_plot) then begin
    filename_24 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_bubble_only/ps/'
    filename_1 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_bubble_overfit/ps/'
    filename_2 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_hash_overfit/ps/'
    
    obsname='obs_id_6296'
    obsname_1='obs_id_6296'
    obsname_2='obs_id_6296'
    
    
    power_24 = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda1-100_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda1-100_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_24 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    power_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda1-100_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda1-100_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_1 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    power_2 = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda1-100_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_res_xx_averemove_bh_dencorr_kperplambda1-100_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_2 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    power_eor = getvar_savefile('/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave','power')
    k_centers_eor = getvar_savefile('/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave','k_centers')
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/thesis_images/1d_sim_ps.pdf',/quiet,/nomatch
    cgplot, k_24, (power_24),/ylog,/xlog,psym=10,xrange=[.01,1.2],charsize=1.3, xtitle = 'k (!8h!x Mpc$\exp-1$)',yrange=[10^1., 10^8.],$
      ytitle='P (mK$\exp2$ Mpc$\exp3$ !8h!x$\exp-3$)',title = '1D power, simulations', thick=4
    cgoplot, k_1, (power_1),/ylog,/xlog,psym=10,xrange=[.001,1.3],charsize=1, yrange=[10^2., 10^7.],color='royal blue',thick=3
    cgoplot, k_2, power_2,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='dark green',thick=3
    cgoplot, k_centers_eor, power_eor,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='purple',thick=4
    cglegend, title=['input Gaussian EoR','recovered Gaussian EoR','image-based EoR, no foregrounds','recovered image-based EoR'], $
      color=['purple','dark green','black','royal blue'], location=[.42,.85], charsize=1.1;, thick=[2,1,2,1]
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    
  endif
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;image_filter_1d=1
  if keyword_set(image_filter_1d) then begin
  
    filename_24 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_hash_overfit/ps/'
    filename_1 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_large_unflag_eorhash/ps/'
    filename_2 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_large_unflag_eorhash/ps/'
    filename_3 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_large_unflag_eorhash/ps/'
    filename_4 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_large_unflag_eorhash_pskspan100/ps/'
    obsname='obs_id_6296'
    obsname_1='obs_id_6296'
    obsname_2='obs_id_6296'
    obsname_3='obs_id_6296'
    obsname_4='obs_id_6296'
    
    cube_type='res'
    
    power_24 = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_'+cube_type+'_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_'+cube_type+'_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_24 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    power_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_large_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_large_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_1 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    power_2 = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_tk_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_tk_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_2 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    ;
    power_3 = getvar_savefile(filename_3 + 'Combined_obs_'+obsname_3+'_cubeXX__even_odd_joint_bh_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_3 + 'Combined_obs_'+obsname_3+'_cubeXX__even_odd_joint_bh_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_3 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    power_4 = getvar_savefile(filename_4 + 'Combined_obs_'+obsname_4+'_cubeXX__even_odd_joint_bh_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_4 + 'Combined_obs_'+obsname_4+'_cubeXX__even_odd_joint_bh_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_4 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    power_eor = getvar_savefile('/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave','power')
    k_centers_eor = getvar_savefile('/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave','k_centers')
    
    
    if ~keyword_set(delta) then ytitle = 'P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)' else ytitle = '$\Delta$$\up2$ (mK$\up2$)'
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/thesis_images/1d_sims_filters_3.pdf',/quiet,/nomatch
    cgplot, k_24, (power_24),/ylog,/xlog,psym=10,xrange=[.05,1.2],charsize=1.3, xtitle = 'k (!8h!x / Mpc)',yrange=[10^3., 10^9.],$
      ytitle=ytitle,title = '1D power, simulation residuals', color='dark green',thick=3
    cgoplot, k_1, (power_1),/ylog,/xlog,psym=10,color='navy',thick=3
    ;cgoplot, k_2, power_2,/ylog,/xlog,psym=10,color='firebrick', thick=3
    cgoplot, k_3, power_3,/ylog,/xlog,psym=10, color='orange',thick=3
    ;cgoplot, k_4, power_4,/ylog,/xlog,psym=10,color='turquoise', thick=4
    cgoplot, k_centers_eor, power_eor,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='purple',thick=3
    cglegend, title=['input Gaussian EoR','recovered Gaussian EoR','no filter','Blackman-Harris filter'],color=['purple','dark green','navy','orange'],location=[.42,.85], charsize=1.1
    ;cglegend, title=['input Gaussian EoR','recovered Gaussian EoR','no filter','Tukey filter','Blackman-Harris filter','Blackman-Harris filter, 100$\lambda$ extent'],$
    ;  color=['purple','dark green','navy','firebrick','orange','turquoise'],$
    ;  location=[.42,.85], charsize=1.1
    
    
    cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
    
    
  endif
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;uv_filter_1d=1
  if keyword_set(uv_filter_1d) then begin
  
    filename_24 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_hash_overfit/ps/'
    filename_1 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_large_unflag_eorhash/ps/'
    filename_2 = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_sim_large_unflag_eorhash_uvfilter/ps/'
    obsname='obs_id_6296'
    obsname_1='obs_id_6296'
    obsname_2='obs_id_6296'
    
    cube_type='res'
    
    power_24 = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_'+cube_type+'_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_24 + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_'+cube_type+'_xx_averemove_bh_dencorr_kperplambda10-50_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_24 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    power_1 = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_bh_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_1 + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_bh_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_1 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    power_2 = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_bh_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'power')
    k_edges = getvar_savefile(filename_2 + 'Combined_obs_'+obsname_2+'_cubeXX__even_odd_joint_bh_ch10-127_'+cube_type+'_xx_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_1dkpower.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_2 = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    power_eor = getvar_savefile('/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave','power')
    k_centers_eor = getvar_savefile('/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/eor_power_1d.idlsave','k_centers')
    
    
    if ~keyword_set(delta) then ytitle = 'P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)' else ytitle = '$\Delta$$\up2$ (mK$\up2$)'
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/thesis_images/1d_sims_uvfilters.pdf',/quiet,/nomatch
    cgplot, k_24, (power_24),/ylog,/xlog,psym=10,xrange=[.05,1.2],charsize=1.3, xtitle = 'k (!8h!x / Mpc)',yrange=[10^3., 10^9.],$
      ytitle=ytitle,title = '1D power, simulation residuals', color='dark green',thick=3
    cgoplot, k_1, (power_1),/ylog,/xlog,psym=10,color='light coral',thick=3
    cgoplot, k_2, power_2,/ylog,/xlog,psym=10,color='turquoise', thick=3
    cgoplot, k_centers_eor, power_eor,/ylog,/xlog,psym=10, yrange=[10^3., 10^7.],xrange=[.005,2.],charsize=1, color='purple',thick=4
    cglegend, title=['input Gaussian EoR','recovered Gaussian EoR','Blackman-Harris image filter','Blackman-Harris image/UV filter'], color=['purple','dark green','light coral','turquoise'],$
      location=[.42,.85], charsize=1.1
      
      
    cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
    
    
  endif
  
  
  
  ;phase_ave_plot=1
  if keyword_set(phase_ave_plot) then begin
    debug_phase_longrun='/nfs/mwa-10/r1/EoRuvfits/analysis/ave_cals/full_ave_phase_dayref.sav'
    longrun_phase = getvar_savefile(debug_phase_longrun,'phase_mean_pointing') ;384, 5, 2, 128
    filename= '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/metadata/1061316296_obs.sav'
    restore, filename
    
    pointing_num=mwa_get_pointing_number(obs,/string)
    poi_name = ['-2','-1','0','1','2']
    poi = where(pointing_num EQ poi_name)
    
    filename= '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/calibration/1061316296_cal.sav'
    cal = getvar_savefile(filename, 'cal')
    cal_base = getvar_savefile(filename, 'cal')
    *cal_base.gain[0] = *cal_base.gain[0] + *cal_base.gain_residual[0]
    zeros = where(atan(*cal_base.gain[0],/phase) EQ 0,n_count)
    if n_count GT 0 then (*cal_base.gain[0])[zeros] = !VALUES.F_NAN
    
    ;    longrun_names_match, obs_names=obs_names
    ;    wh_zen = where(obs_names[*,2] EQ 'zenith',n_count)
    ;    gain_base = complex(FLTARR(384,128,N_elements(wh_zen)))
    ;    pol_i=0
    ;    for obs_i=0, N_elements(wh_zen)-1 do begin
    ;      filename = '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/calibration/'+obs_names[wh_zen[obs_i],0]+'_cal.sav'
    ;      cal_base2 = getvar_savefile(filename, 'cal')
    ;
    ;      gain_base[*,*,obs_i] = *cal_base2.gain[pol_i] + *cal_base2.gain_residual[pol_i]
    ;    endfor
    ;    zeros = where(atan(gain_base,/phase) EQ 0,n_count)
    ;    if n_count GT 0 then gain_base[zeros] = !VALUES.F_NAN
    ;
    
    cal_base.mode_fit=0
    cal_polyfit1 = vis_cal_polyfit(cal_base,obs)
    
    *cal.gain[0] = abs(*cal.gain[0])*exp(Complex(0,1)*( atan((*cal_polyfit1.gain[0]),/phase)+reform(longrun_phase[*,poi,0,*])))
    ;*cal.gain[1] = abs(*cal.gain[1])*exp(Complex(0,1)*( atan((*cal_polyfit1.gain[1]),/phase)+reform(longrun_phase[*,poi,1,*])))
    
    freq_arr = (*obs.baseline_info).freq / 1e6
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/thesis_images/phase_example.pdf',/quiet,/nomatch
    
    cgplot, freq_arr,atan((*cal_base.gain[0])[*,0],/phase),color = 'grey' ,xrange = [freq_arr[0] - 1.,freq_arr[383] + 1.], yrange=[-3.14,3.14],thick=3,$
      xtitle='Frequency (MHz)', ytitle='Phase (radians)',title='Example raw and averaged phases',aspect=.5,charsize=1
    cgoplot, freq_arr,atan((*cal.gain[0])[*,0],/phase),color = 'blue' ,xrange = [freq_arr[0] - 1.,freq_arr[383] + 1.], yrange=[-3.14,3.14],thick=3
    
    cgoplot, freq_arr,atan((*cal_base.gain[0])[*,3],/phase),color = 'grey' ,xrange = [freq_arr[0] - 1.,freq_arr[383] + 1.], yrange=[-3.14,3.14],thick=3
    cgoplot, freq_arr,atan((*cal.gain[0])[*,3],/phase),color = 'purple' ,xrange = [freq_arr[0] - 1.,freq_arr[383] + 1.], yrange=[-3.14,3.14],thick=3
    
    cgoplot, freq_arr,atan((*cal_base.gain[0])[*,50],/phase),color = 'grey' ,xrange = [freq_arr[0] - 1.,freq_arr[383] + 1.], yrange=[-3.14,3.14],thick=3
    cgoplot, freq_arr,atan((*cal.gain[0])[*,50],/phase),color = 'dark green' ,xrange = [freq_arr[0] - 1.,freq_arr[383] + 1.], yrange=[-3.14,3.14],thick=3
    
    cgoplot, freq_arr,atan((*cal_base.gain[0])[*,60],/phase),color = 'grey' ,xrange = [freq_arr[0] - 1.,freq_arr[383] + 1.], yrange=[-3.14,3.14],thick=3
    cgoplot, freq_arr,atan((*cal.gain[0])[*,60],/phase),color = 'maroon' ,xrange = [freq_arr[0] - 1.,freq_arr[383] + 1.], yrange=[-3.14,3.14],thick=3
    
    
    cglegend, title=['tile 0'], color=['blue'], location=[.17,.3], charsize=1.1
    cglegend, title=['tile 3'], color=['purple'], location=[.342,.3], charsize=1.1
    cglegend, title=['tile 50'], color=['dark green'], location=[.522,.3], charsize=1.1
    cglegend, title=['tile 60'], color=['maroon'], location=[.7,.3], charsize=1.1
    
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif
  
  ;vis_stats_tv=1
  if keyword_Set(vis_stats_tv) then begin
    obsid='1061313128'
    dir = '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_vis_stats/'
    
    all_3D_20=LONG(INTARR(384,56))
    
    tiles_20 = LONG(INTARR(129,129)) ;tile name indexed
    
    
    ;Find where the autos are in the visibility data
    obs = GETVAR_SAVEFILE(dir+'metadata/'+obsid+'_obs.sav','obs')
    freq_arr = (*obs.baseline_info).freq / 1e6
    tiles = (*obs.baseline_info).tile_a - (*obs.baseline_info).tile_b
    autos = where(tiles EQ 0)
    
    
    ;Restore array of calibrated and model visibilities to generate residual visibilities (and make sure the autos are 0!)
    vis_XX = GETVAR_SAVEFILE(dir+'vis_data/'+obsid+'_vis_XX.sav', 'vis_ptr')
    vis_XX_model = GETVAR_SAVEFILE(dir+'vis_data/'+obsid+'_vis_model_XX.sav', 'vis_model_ptr')
    vis_XX_res=*vis_XX-*vis_XX_model
    vis_XX_res[autos]=0
    
    flag_arr = GETVAR_SAVEFILE(dir+'vis_data/'+obsid+'_flags.sav', 'vis_weights')
    ;Run FHD util to find the right bin indicies for even and odd samples, then subtract
    bin_i=split_vis_weights(obs,flag_arr,bi_use=bi_use,/preserve_flags)
    vis_XX_res = vis_XX_res[*,*bi_use[0]] - vis_XX_res[*,*bi_use[1]]
    
    ;Bins refer to time samples, so an arbitrary one will work to pick out the right tiles
    tile_a = ((*obs.baseline_info).tile_a)[*bi_use[0]] ; tile_a information with the even-odd formulism
    tile_b = ((*obs.baseline_info).tile_b)[*bi_use[0]]
    
    ;Set the binsize and histogram the visibilities (at this point, either residuals or even-odd residuals)
    binsize=1.
    result=histogram(abs(vis_XX_res),binsize=binsize,locations=locations,omax=omax, /NAN,reverse_indices=ri)
    
    ;binned_diff is the binned visibility residuals over all obsids, and result is the binned visibility residuals for this particular obsid
    ;BUILT-IN ASSUMPTION: There are at least some in the 0th bin. It's a fantastic assumption; most are 0.
    ;If this is the first obsid to be binned, set that equal to the all-obs binning.
    binned_diff=uLong64(result)
    
    outlier_min = 300
    mod_num=28
    
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
      even_odd=1
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
          
          for tile_j=0, N_elements(tile_subset)-1 do tiles_20[locations_tile_a[tile_i],tile_b[row[tile_subset[tile_j]]]] += 1
          
          
        endfor
        
      endif
      
      
      ;Find the times of ....
      
      For loc_i = 0, N_elements(col)-1 do all_3D_20[col[loc_i],row[loc_i] mod mod_num]=all_3D_20[col[loc_i],row[loc_i] mod mod_num]+uLong64(1)
      
      
      
      ;*******1D freq plots
      ;longrun_names_match, obs_names=obs_names
      ;minustwo_obs = where(obs_names[*,2] EQ 'minustwo')
      all_3D_20_final = total(all_3d_20,2)
      ;for obs_i=0,1028 do if (where(obs_i EQ minustwo_obs) NE -1) and (obs_names[obs_i,1] NE 'Oct31') then for freq_i=0,383 do all_3D_20_final[freq_i] = all_3D_20_final[freq_i] $
      ;  + total(all_3D_20_total[freq_i,*,obs_i])
      ;cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/antenna/remove_center_freq_pointing/plots/3D_20_evenodd2_freq_minustwo.png',/quiet,/nomatch
      cgplot, freq_arr,all_3D_20_final, ytitle='Counts', xtitle='Frequency (MHz)',xrange=[freq_arr[0],freq_arr[383]], title='-2 Fall 2013 Residual Even-Odd Visibilities, Center flagged',charsize=1.5
      ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      stop
    endif
    
  endif
  
  if keyword_set(delay_plot) then begin
    ;Call mkbp_updated
    quick_image, result[*,192:383], x_axis, y_axis[192:383], data_range=[0,25], data_aspect=1, window=2,xtitle='Calibration delay (real)', ytitle='Time (ns)', charsize=1, cb_title='counts', multi_pos = [0,0,.5,1],no_ps_close=1
    cgoplot, 384.*real_part(fft_amp_mean[192:383,2,0,tile_i]),y_axis[192:383]
  endif
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;slice_plots=1
  if keyword_set(slice_plots) then begin
  
    ;  ps_wrapper, '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/', 'beardsley_thesis_list',kpar_range_kperppower=[1.05,1.15],$
    ;  /exact_obsname, ps_foldername='ps_bh',/refresh_beam,plot_stdset=0, kperp_lambda_plot_range=[.5,100],kperp_range_lambda_1dave=[.5,200],$
    ;  /plot_kperp,image_window_name='Blackman-Harris',freq_ch_range=[0,127],pol_inc='yy'
  
    filename_kpar = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_bh/'
    filename_kperp = '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_bh/'
    
    obsname='beardsley_thesis_list'
    obsname_1='beardsley_thesis_list'
    
    averemove=''
    
    cube_type='res'
    kpar_range = '19-20'
    kperp_range='0.29-0.3'
    ;kperp_range='1.05-1.15'
    
    power_kpar = getvar_savefile(filename_kpar + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_maxuv100_bh_ch0-127_'+cube_type+'_yy_'+averemove+'swbh_dencorr_kperplambda'+kpar_range+'_kpar_power.idlsave', $
      'power')
    weights_kpar = getvar_savefile(filename_kpar + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_maxuv100_bh_ch0-127_'+cube_type+'_yy_'+averemove+'swbh_dencorr_kperplambda'+kpar_range+'_kpar_power.idlsave', $
      'weights')
    k_edges = getvar_savefile(filename_kpar + 'Combined_obs_'+obsname+'_cubeXX__even_odd_joint_maxuv100_bh_ch0-127_'+cube_type+'_yy_'+averemove+'swbh_dencorr_kperplambda'+kpar_range+'_kpar_power.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_kpar = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    hubble_param=.71
    k_kpar=k_kpar/hubble_param
    
    dsigma_kpar=weight_invert(sqrt(weights_kpar))
    
    power_kperp = getvar_savefile(filename_kperp + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_maxuv100_bh_ch0-127_'+cube_type+'_yy_'+averemove+'swbh_dencorr_kpar'+kperp_range+'_kperp_power.idlsave', $
      'power')
    weights_kperp = getvar_savefile(filename_kperp + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_maxuv100_bh_ch0-127_'+cube_type+'_yy_'+averemove+'swbh_dencorr_kpar'+kperp_range+'_kperp_power.idlsave', $
      'weights')
    k_edges = getvar_savefile(filename_kperp + 'Combined_obs_'+obsname_1+'_cubeXX__even_odd_joint_maxuv100_bh_ch0-127_'+cube_type+'_yy_'+averemove+'swbh_dencorr_kpar'+kperp_range+'_kperp_power.idlsave', $
      'k_edges')
    n_k=n_elements(k_edges)
    k_kperp = (k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
    
    k_kperp=k_kperp/hubble_param
    
    dsigma_kperp=weight_invert(sqrt(weights_kperp))
    
    ytitle = 'P (mK$\exp2$ Mpc$\exp3$ /  !8h!x$\exp3$)'
    perp = '!9' + String("136B) + '!X'
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/thesis_images/1d_kparslice.pdf',/quiet,/nomatch
    cgplot, k_kpar, (power_kpar),/ylog,/xlog,psym=10,xrange=[.0045,1.7],charsize=1.3, xtitle = 'k$\sub||$ (!8h!x / Mpc)',yrange=[10^3., 10^15.],$
      ytitle=ytitle, color='blue',thick=4,aspect=.4
    cgoplot, k_kpar, dsigma_kpar,/ylog,/xlog,psym=10, thick=3, color='blue',linestyle=2
    cglegend, title=['power','noise'], color=['blue','blue'],linestyle=[0,2],$
      location=[.45,.65], charsize=1.1,thick=4
    cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
    ;cgPS_Open,'/nfs/eor-00/h1/nbarry/thesis_images/1d_kperpslice.pdf',/quiet,/nomatch
    cgplot, k_kperp, abs(power_kperp),/ylog,/xlog,psym=10,xrange=[.0055,0.1],charsize=1.3, xtitle = 'k$\sub'+perp+'$ (!8h!x / Mpc)',yrange=[10^3., 10^15.],$
      ytitle=ytitle, thick=4,aspect=.4,color='red'
    cgoplot, k_kperp, dsigma_kperp,/ylog,/xlog,psym=10, thick=3, color='red',linestyle=2
    
    
    cglegend, title=['power','noise'], color=['red','red'],linestyle=[0,2],$
      location=[.45,.65], charsize=1.1,thick=4
    ;    cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
    
    
    ;cgPS_Close,/pdf,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
    
    
  endif
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  if keyword_set(contour_plot) then begin
    ps_wrapper, '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/', 'beardsley_thesis_list',$
      /exact_obsname, ps_foldername='ps_bh',image_window_name='Blackman-Harris', pol_inc='yy',plot_stdset=0,$
      freq_ch_range=*ch_bands[ch_i],plot_2d_masked=1
  endif
  
  
  ;limit_plot=1
  if keyword_set(limit_plot) then begin
  
    cut = 'beardsley_thesis_list'
    chans='ch10-127_'
    win='bh_'
    spec_window='bh_'
    basefile='/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_bh/Combined_obs_'+cut+'_cubeXX__even_odd_joint_maxuv100_'+win+chans
    cubes=['res']
    pols=['xx']
    averemove=''
    endfile='_'+averemove+'sw'+spec_window+'dencorr_no_110deg_wedge_kperplambda20-80_1dkpower.idlsave'
    
    ;limit_percent=0.9772 ;2 sigma
    limit_percent = 0.97725
    hubble_param=0.71
    
    ;color_array=['black','blue','green','purple']
    n_cubes=1
    color_num = [10,12,14,16]
    rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    ;rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101]];,[198,118,63]]
    
    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]
      
    color_array = [10B,12B,14B,16B]
    
    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 5 else thickness = 3
    
    for j=0,N_elements(pols)-1 do begin
    
      file=basefile+cubes+'_'+pols[j]+endfile
      restore,file
      ; find the best limit in this file
      n_k=n_elements(k_edges)
      k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
      delta=power*(k^3.)/(2.*!pi^2.)
      dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
      dsigma[0] = !Values.F_INFINITY
      ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
      limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf((delta)/dsigma/sqrt(2)))*sqrt(2))+delta
      limits_abs=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(abs(delta)/dsigma/sqrt(2)))*sqrt(2))+abs(delta)
      lim=min(limits,ind)
      header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut+' '+chans+' '+cubes+' '+win+' '+pols[j]+')'
      print,header
      
      k=k/hubble_param
      
      low_range=10^2.
      low_x_range=.15
      ;cgPS_Open,'/nfs/mwa-00/h1/nbarry/'+pols[j]+'_limit_noise.pdf',/quiet,/nomatch
      cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[low_x_range,1.7], yrange=[low_range,10^7.],ytitle='$\Delta$$\up2$ (mK$\up2$)', $
        xtitle='k (!8h!x Mpc$\up-1$)', charsize =1.25, color=color_array[0],thick=thickness,/nodata,aspect=.6
        
        
      error_low = (dsigma*2.); < (delta);*.99999
      lows = where(delta LT 0, n_count)
      if n_count GT 0 then error_low[lows]=abs(delta[lows])
      lows = where((abs(delta) - error_low) LT low_range,n_count)
      if n_count GT 0 then error_low[lows]=abs(delta[lows])-low_range
      error_high = dsigma*2.
      lows = where(k LT low_x_range,n_count)
      if n_count GT 0 then error_high[lows]=!Values.F_INFINITY
      if n_count GT 0 then limits[lows]=!Values.F_INFINITY
      delta_k=k[1]/2.
      for k_i=0, n_k-2 do $
        cgColorFill, [k[k_i]-delta_k, k[k_i]+delta_k, k[k_i]+delta_k,k[k_i]-delta_k], $
        [limits[k_i], limits[k_i], abs(delta[k_i])-error_low[k_i],abs(delta[k_i])-error_low[k_i]], $
        Color='grey',/checkforfinite
        
      cgoplot, k,(limits),/xlog,/ylog,psym=10,color=color_array[0],thick=thickness
      
      ;cgoplot, k,abs(delta),/xlog,/ylog,psym=10,color=color_array[1],thick=thickness,linestyle=2
      cgoplot, k,(delta),/xlog,/ylog,psym=10,color=color_array[1],thick=thickness
      
      cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2, color=color_array[1],thick=thickness
    ;cglegend, title=['upper limit','measured power','thermal noise'],color=[color_array[0],color_array[1],color_array[1]], location=[.15,.76], charsize=1, thick=5,linestyle=[0,0,2],length=.05
      
      
    endfor
    ;cglegend, title=['z=7.1','z=6.8','z=6.5','z=7'],color=color_array, location=[.2,.85], charsize=1, thick=3
    ;cglegend, title=['z=7 upper limit','1$\sigma$ noise','Beardsley 2016'],color=[color_array[3],color_array[3],color_array[2]], location=[.18,.85], charsize=1.25, thick=5, linestyle=[0,2,0]
    ;cglegend, title=['Beardsley z=7.1','z=7 upper limit, image BH, 100$\lambda$'],color=[color_array[2],color_array[1]], location=[.18,.85], charsize=1.25, thick=thickness, linestyle=[0,0,0]
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif
  
  limit_compare=1
  if keyword_set(limit_compare) then begin
  
    cut = 'beardsley_thesis_list'
    chans='ch10-127_'
    win='bh_'
    spec_window='bh_'
    basefile1='/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_bh/Combined_obs_'+cut+'_cubeXX__even_odd_joint_maxuv100_'+win+chans
    cubes=['res']
    pols=['yy']
    averemove='averemove_'
    endfile1='_'+averemove+'sw'+spec_window+'dencorr_no_110deg_wedge_kperplambda20-80_1dkpower.idlsave'
    
    cut='wedge_cut_plus_res_cut'
    chans='ch0-95_'
    win = ''
    basefile2='/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps_jan2016/Combined_obs_'+cut+'_cubeXX__even_odd_joint_'+win+chans
    cubes=['res']
    pols=['yy']
    endfile2='_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-60_kpar0.15-200_1dkpower.idlsave'
    
    ;limit_percent=0.9772 ;2 sigma
    limit_percent = 0.97725
    hubble_param=0.71
    
    ;color_array=['black','blue','green','purple']
    n_cubes=4
    color_num = [10,12,14,16]
    rgbcolors = [[137,117,202],[113,166,89],[203,86,131],[197,120,62]]
    ;rgbcolors=[[131,119,203],[114,166,89],[202,86,139],[91,169,101]];,[198,118,63]]
    
    for n_i=0, n_cubes-1 do $
      TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]
      
    color_array = [10B,12B,14B,16B]
    
    pdf_plot=1
    if keyword_set(pdf_plot) then thickness = 5 else thickness = 3
    
    for j=0,N_elements(pols)-1 do begin
    
      file=basefile1+cubes+'_'+pols[j]+endfile1
      restore,file
      ; find the best limit in this file
      n_k=n_elements(k_edges)
      k1=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
      delta=power*(k1^3.)/(2.*!pi^2.)
      dsigma1=(k1^3.)/(2.*!pi^2.)/sqrt(weights)
      dsigma1[0] = !Values.F_INFINITY
      ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
      limits1=dsigma1*(inverf(limit_percent-(1.-limit_percent)*erf((delta)/dsigma1/sqrt(2)))*sqrt(2))+delta
      limits_abs=dsigma1*(inverf(limit_percent-(1.-limit_percent)*erf(abs(delta)/dsigma1/sqrt(2)))*sqrt(2))+abs(delta)
      lim=min(limits1,ind)
      header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k1[ind]/hubble_param)+' h Mpc^-1 ('+cut+' '+chans+' '+cubes+' '+win+' '+pols[j]+')'
      print,header
      
      k1=k1/hubble_param
      
      file=basefile2+cubes+'_'+pols[j]+endfile2
      restore,file
      ; find the best limit in this file
      n_k=n_elements(k_edges)
      k2=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
      delta=power*(k2^3.)/(2.*!pi^2.)
      dsigma2=(k2^3.)/(2.*!pi^2.)/sqrt(weights)
      dsigma2[0] = !Values.F_INFINITY
      ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
      limits2=dsigma2*(inverf(limit_percent-(1.-limit_percent)*erf((delta)/dsigma2/sqrt(2)))*sqrt(2))+delta
      limits_abs=dsigma2*(inverf(limit_percent-(1.-limit_percent)*erf(abs(delta)/dsigma2/sqrt(2)))*sqrt(2))+abs(delta)
      lim=min(limits2,ind)
      header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k2[ind]/hubble_param)+' h Mpc^-1 ('+cut+' '+chans+' '+cubes+' '+win+' '+pols[j]+')'
      print,header
      
      k2=k2/hubble_param
      
      low_range=10^2.
      low_x_range=.15
      cgPS_Open,'/nfs/mwa-00/h1/nbarry/'+pols[j]+'_limit_compare_all.pdf',/quiet,/nomatch
      cgplot, k1,abs(limits1),/xlog,/ylog,psym=10, xrange=[low_x_range,1.7], yrange=[low_range,10^7.],ytitle='$\Delta$$\up2$ (mK$\up2$)', $
        xtitle='k (!8h!x Mpc$\up-1$)', charsize =1.25, color=color_array[0],thick=thickness,/nodata;,aspect=1
        
        
      ;      error_low = (dsigma*2.); < (delta);*.99999
      ;      lows = where(delta LT 0, n_count)
      ;      if n_count GT 0 then error_low[lows]=abs(delta[lows])
      ;      lows = where((abs(delta) - error_low) LT low_range,n_count)
      ;      if n_count GT 0 then error_low[lows]=abs(delta[lows])-low_range
      ;      error_high = dsigma*2.
      ;      lows = where(k1 LT low_x_range,n_count)
      ;      if n_count GT 0 then error_high[lows]=!Values.F_INFINITY
      ;      if n_count GT 0 then limits[lows]=!Values.F_INFINITY
      ;      delta_k=k1[1]/2.
      ;      for k_i=0, n_k-2 do $
      ;        cgColorFill, [k[k_i]-delta_k, k[k_i]+delta_k, k[k_i]+delta_k,k[k_i]-delta_k], $
      ;        [limits[k_i], limits[k_i], abs(delta[k_i])-error_low[k_i],abs(delta[k_i])-error_low[k_i]], $
      ;        Color='grey',/checkforfinite
        
      cgoplot, k1,(limits1),/xlog,/ylog,psym=10,color=color_array[0],thick=thickness
      cgoplot, k1,dsigma1,/xlog,/ylog,psym=10,linestyle=2, color=color_array[0],thick=thickness
      
      cgoplot, k2,(limits2),/xlog,/ylog,psym=10,color=color_array[1],thick=thickness
      cgoplot, k2,dsigma2,/xlog,/ylog,psym=10,linestyle=2, color=color_array[1],thick=thickness
      
      cglegend, title=['NB2018 z=7 upper limit','NB2018 z=7 thermal noise','AB2016 z=7.1 upper limit','AB2016 z=7.1 thermal noise'],color=[color_array[0],color_array[0],color_array[1],color_array[1]],$
        location=[.55,.275], charsize=1, thick=5,linestyle=[0,2,0,2],length=.05
        
        
    endfor
    ;cglegend, title=['z=7.1','z=6.8','z=6.5','z=7'],color=color_array, location=[.2,.85], charsize=1, thick=3
    ;cglegend, title=['z=7 upper limit','1$\sigma$ noise','Beardsley 2016'],color=[color_array[3],color_array[3],color_array[2]], location=[.18,.85], charsize=1.25, thick=5, linestyle=[0,2,0]
    ;cglegend, title=['Beardsley z=7.1','z=7 upper limit, image BH, 100$\lambda$'],color=[color_array[2],color_array[1]], location=[.18,.85], charsize=1.25, thick=thickness, linestyle=[0,0,0]
    
    add_others=1
    if keyword_set(add_others) then begin
      n_papers=7
      papers_array = PTRARR(n_papers,/allocate)
      name_array = STRARR(n_papers)
      
      ;Organized by paper: mK value, center of k bin/range, and redshift
      ali_2015 = transpose([22.4^2.,.325,8.4])
      *papers_array[0] = ali_2015
      name_array[0] = 'Ali, 2015'
      
      paciga_2013 = transpose([6.15E4,.5,8.6])
      *papers_array[1] = paciga_2013
      name_array[1] = 'Paciga, 2013'
      
      parsons_2014 = transpose([41^2.,.27,7.7])
      *papers_array[2] = parsons_2014
      name_array[2] = 'Parsons, 2014'
      
      patil_2017 = FLTARR(15,3)
      patil_2017[*,0] = [131.5^2.,242.1^2.,220.9^2.,337.4^2.,407.7^2.,86.4^2.,144.2^2.,184.7^2.,$
        296.1^2.,342.0^2.,79.6^2.,108.8^2.,148.6^2.,224.0^2.,366.1^2.]
      patil_2017[*,1] = [.053,.067,.083,.103,.128,.053,.067,.083,.103,.128,.053,.067,.083,.103,.128]
      patil_2017[*,2] = [8.3,8.3,8.3,8.3,8.3,9.15,9.15,9.15,9.15,9.15,10.1,10.1,10.1,10.1,10.1]
      patil_2017_short = FLTARR(6,3)
      patil_2017_short[*,0] = [131.5^2.,220.9^2.,86.4^2.,184.7^2.,79.6^2.,148.6^2.]
      patil_2017_short[*,1] = [.053,.083,.053,.083,.053,.083]
      patil_2017_short[*,2] = [8.3,8.3,9.15,9.15,10.1,10.1]
      *papers_array[3] = patil_2017_short
      name_array[3] = 'Patil, 2017'
      
      jacobs_2015 = FLTARR(4,3)
      jacobs_2015[*,0] = [1.64E4,6.9E3,6.9E3,2.4E3]
      jacobs_2015[*,1] = [.2,.1,.2,.2]
      jacobs_2015[*,2] = [10.29,8.54,7.94,7.55]
      *papers_array[4] = jacobs_2015
      name_array[4] = 'Jacobs, 2015'
      
      dillon_2014 = FLTARR(11,3)
      dillon_2014[*,0] = [2.6E7,1.16E6,8.64E5,6.7E5,1.3E6,1.28E7,5.26E7,5.67E8,4.58E6,2.93E8,6.92E8]
      dillon_2014[*,1] = [.058,.06,.063,.065,.0678,.0712,.0715,.149,.15,.15,.089]
      dillon_2014[*,2] = [11.68,10.868,10.153,9.518,8.444,7.985,7.57,7.1896,6.84,6.52,6.23]
      *papers_array[5] = dillon_2014
      name_array[5] = 'Dillon, 2014'
      
      dillon_2015 = FLTARR(3,3)
      dillon_2015[*,0] = [3.8E4,3.69E4,4.67E4]
      dillon_2015[*,1] = [.18,.18,.16]
      dillon_2015[*,2] = [6.4,6.8,7.25]
      *papers_array[6] = dillon_2015
      name_array[6] = 'Dillon, 2015'
      
;      beardsley_2016 = FLTARR(9,3)
;      beardsley_2016[*,0] = [3.67E4,2.70E4,3.56E4,3.02E4,4.70E4,3.22E4,3.2E4,2.6E4,2.5E4]
;      beardsley_2016[*,1] = [.231,.27,.24,.24,.20,.24,.16,.14,.14]
;      beardsley_2016[*,2] = [7.1,7.1,6.8,6.8,6.5,6.5,7.1,6.8,6.5]
;      beardsley_2016_short = FLTARR(6,3)
;      beardsley_2016_short[*,0] = [2.70E4,3.02E4,3.22E4,3.2E4,2.6E4,2.5E4]
;      beardsley_2016_short[*,1] = [.27,.24,.24,.16,.14,.14]
;      beardsley_2016_short[*,2] = [7.1,6.8,6.5,7.1,6.8,6.5]
;      *papers_array[7] = beardsley_2016_short
;      name_array[7] = 'Beardsley, 2016'
      
      sym_array = ['Filled Diamond','Filled Bowtie','Filled Up Triangle','Filled Lower Half Circle','Filled Standing Bar','Filled Circle','Filled Laying Bar','Filled Star']
      sym_num = [14,24,17,40,26,16,28,46]
      
      color_num = [10,12,14,16,18,20]
      ;rgbcolors = [[9,14,181],[70,104,206],[52,115,145],[60,158,149],[60,158,102],[52,175,63]]
      ;rgbcolors = [[60,158,190],[60,158,160],[60,158,145],[60,158,120],[60,158,102],[60,158,63]]
      rgbcolors = [[45,31,198],[52,116,173],[40,183,183],[39,158,96],[47,158,39]]
      
      ;for n_i=0, 4 do $
      ;  TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]
      
      ;color_array = [10B,12B,14B,16B,18B,20B]
      color_array = ['purple', 'teal', 'blue', 'firebrick','purple', 'teal', 'blue', 'firebrick']
      color_array = ['purple', 'teal', 'blue', 'firebrick','orange red', 'teal', 'firebrick']
      
      for papers_i=0, n_papers-1 do begin
        k_arr = *papers_array[papers_i]
        
        ;inds1 = where(k_arr[*,2] GT 6,n_count1)
        inds = where((k_arr[*,2] LE 8.0) AND (k_arr[*,2] GE 6.0),n_count)
        ;inds3 = where((k_arr[*,2] LT 8.0) AND (k_arr[*,2] GE 7.0),n_count3)
        ;inds4 = where((k_arr[*,2] LT 9.0) AND (k_arr[*,2] GE 8.0),n_count4)
        ;inds5 = where((k_arr[*,2] LT 10.0) AND (k_arr[*,2] GE 9.0),n_count5)
        ;inds5 = where(k_arr[*,2] GE 9.0,n_count5)
        
        if n_count GT 0 then cgoplot, k_arr[inds,1],k_arr[inds,0],psym=sym_array[papers_i],symsize = 1.25, color=color_array[papers_i]
      ;if n_count1 GT 0 then cgoplot, k_arr[inds1,1],k_arr[inds1,0],psym=sym_array[papers_i],thick=3, color=color_array[0]
      ;if n_count2 GT 0 then cgoplot, k_arr[inds2,1],k_arr[inds2,0],psym=sym_array[papers_i],thick=3, color=color_array[1]
      ;if n_count3 GT 0 then cgoplot, k_arr[inds3,1],k_arr[inds3,0],psym=sym_array[papers_i],thick=3, color=color_array[2]
      ;if n_count4 GT 0 then cgoplot, k_arr[inds4,1],k_arr[inds4,0],psym=sym_array[papers_i],thick=3, color=color_array[3]
      ;if n_count5 GT 0 then cgoplot, k_arr[inds5,1],k_arr[inds5,0],psym=sym_array[papers_i],thick=3, color=color_array[4]
      ;if n_count6 GT 0 then cgoplot, k_arr[inds6,1],k_arr[inds6,0],psym=sym_array[papers_i],thick=3, color=color_array[5]
        
        
      endfor
      cglegend, Title=name_array[[2,4,5,6]], color=color_array[[2,4,5,6]], location=[.73,.42], psym=sym_num[[2,4,5,6]], charsize=1,Length=0,symsize = 1.25
    endif
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif
  
end
