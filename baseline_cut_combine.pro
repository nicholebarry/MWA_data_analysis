pro baseline_cut_combine


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
  FOR cable_i=0,n_cable-1 DO tile_use_arr[cable_i]=Ptr_new(cable_len EQ cable_length_ref[cable_i])
  
  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23.txt'
  readcol, filename, obs_ids, format='A', /silent
  n_obs = N_Elements(obs_ids)
  
  cal_med_fft = complex(FLTARR(2,n_obs,384,128))
  cal_low_fft = complex(FLTARR(2,n_obs,384,128))
  cal_high_fft = complex(FLTARR(2,n_obs,384,128))
  cal_med_fft_cableave = complex(FLTARR(2,n_obs,384,n_cable))
  cal_low_fft_cableave = complex(FLTARR(2,n_obs,384,n_cable))
  cal_high_fft_cableave = complex(FLTARR(2,n_obs,384,n_cable))
  cal_med_fft_cableave_obsave = complex(FLTARR(2,6,384,n_cable))
  cal_low_fft_cableave_obsave = complex(FLTARR(2,6,384,n_cable))
  cal_high_fft_cableave_obsave = complex(FLTARR(2,6,384,n_cable))
  tile_use = PTRARR(n_obs,/allocate)
  
  for obs_i=0, n_obs-1 do begin
    ;cal_med = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/1061316296_cal_LT30GT10_GT150.sav','cal2')
    cal_med = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/selected_cal/GT7LT20_GT150/'+obs_ids[obs_i]+'_GT7LT20_GT150_cal.sav','cal2')
    cal_high = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/selected_cal/GT7_LT50/'+obs_ids[obs_i]+'_GT7_LT50_cal.sav','cal2')
    cal_low =getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/calibration/'+obs_ids[obs_i]+'_cal.sav','cal')
    obs =getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/'+obs_ids[obs_i]+'_obs.sav','obs')
    
    freq_use = (*obs.baseline_info).freq_use
    freq_arr = (*obs.baseline_info).freq
    freq_not_use = where(freq_use EQ 0)
    freq_use = where(freq_use)
    *tile_use[obs_i] = where((*obs.baseline_info).tile_use)
    
    for pol_i=0,1 do begin
      (*cal_low.gain[pol_i])= (*cal_low.gain[pol_i])+(*cal_low.gain_residual[pol_i])
      (*cal_med.gain[pol_i])[freq_not_use,*] = 0
      (*cal_low.gain[pol_i])[freq_not_use,*] = 0
      (*cal_high.gain[pol_i])[freq_not_use,*] = 0
      
      cal_med_fft[pol_i,obs_i,*,*] = fft(abs(*cal_med.gain[pol_i]),dim=1)
      cal_low_fft[pol_i,obs_i,*,*] = fft(abs(*cal_low.gain[pol_i]),dim=1)
      cal_high_fft[pol_i,obs_i,*,*] = fft(abs(*cal_high.gain[pol_i]),dim=1)
    endfor
  endfor
  
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  N=384.
  T= (80000.)
  X = FINDGEN((N - 1)/2) + 1
  x_axis = [0.0, X, N/2, -N/2 + X]/(N*T)*1E9
  
  selected_med_pos = where((x_axis LE 400.) AND (x_axis GE 150),nc)
  selected_med_neg = where((x_axis GE -400.) AND (x_axis LE -150),nc)
  selected_low_pos = where((x_axis LT 150.) AND (x_axis GE 0),nc)
  selected_low_neg = where((x_axis GT -150.) AND (x_axis LT 0),nc)
  selected_high_pos = where((x_axis GT 400.),nc)
  selected_high_neg = where((x_axis LT -400.),nc)
  
  cal_selected = complex(FLTARR(2,94,384,128))
  for pol_i=0,1 do $
    for obs_i=0, 93 do $
    cal_selected[pol_i,obs_i,*,*] = [reform(cal_low_fft[pol_i,obs_i,selected_low_pos,*]),reform(cal_med_fft[pol_i,obs_i,selected_med_pos,*]),reform(cal_high_fft[pol_i,obs_i,selected_high_pos,*]),$
    reform(cal_high_fft[pol_i,obs_i,selected_high_neg,*]),reform(cal_med_fft[pol_i,obs_i,selected_med_neg,*]),reform(cal_low_fft[pol_i,obs_i,selected_low_neg,*])]
    
  stop
  
  ;AVERAGING
  FOR cable_i=0,n_cable-1 DO BEGIN
    tile_arr_sub = tile_arr[where(*tile_use_arr[cable_i])]
    for obs_i=0, n_obs-1 do begin
      match, tile_arr_sub, *tile_use[obs_i], suba, subb, count=count
      if count GT 0 then begin
        cal_low_fft_cableave[*,obs_i,*,cable_i] = mean(cal_low_fft[*,obs_i,*,tile_arr_sub[suba]], dim=4)
        cal_med_fft_cableave[*,obs_i,*,cable_i] = mean(cal_med_fft[*,obs_i,*,tile_arr_sub[suba]], dim=4)
        cal_high_fft_cableave[*,obs_i,*,cable_i] = mean(cal_high_fft[*,obs_i,*,tile_arr_sub[suba]], dim=4)
      endif
    endfor
  endfor
  
  day='Aug23'
  parsednames=[day+'minustwo',day+'minusone',day+'zenith',day+'plusone',day+'plustwo',day+'plusthree']
  obs_ptr=PTRARR(10,/allocate_heap)
  pointing_num=[-2,-1,0,1,2,3]
  FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
    filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/' + parsednames[j] + '.txt'
    ;obs_temp='empty'
    readcol, filename, obs_pointing, format='A', /silent
    ;If j NE 0 then parsednumbers=[parsednumbers,N_elements(obs_temp)] else parsednumbers=N_elements(obs_temp)
    ;*obs_ptr[j]=obs_temp
    
    match, obs_ids, obs_pointing, suba, subb
    cal_high_fft_cableave_obsave[*,j,*,*] = mean(cal_high_fft_cableave[*,suba,*,*], dim=2)
    cal_med_fft_cableave_obsave[*,j,*,*] = mean(cal_med_fft_cableave[*,suba,*,*], dim=2)
    cal_low_fft_cableave_obsave[*,j,*,*] = mean(cal_low_fft_cableave[*,suba,*,*], dim=2)
    
    for cable_i=0, n_cable-1 do begin
      cgPS_Open,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/selected_cal/selected_delay_ave_plots/med_pointing_'+strtrim(pointing_num[j],2)+'_cable_'+strtrim(ULONG(cable_length_ref[cable_i]),2)+'.png',/quiet,/nomatch
      cgplot, x_axis, abs(cal_med_fft_cableave[0,subb[0],*,cable_i]), yrange=[0,.05], xrange=[30, 7000], color='grey', /xlog, charsize=1,$
        ytitle='Amp', xtitle='Delay (ns)', title='Pointing ' +strtrim(pointing_num[j],2)+ ', cable ' + strtrim(ULONG(cable_length_ref[cable_i]),2) + ', Med delays'
      for obs_i=1, N_elements(obs_pointing)-1 do begin
        cgoplot, x_axis, abs(cal_med_fft_cableave[0,subb[obs_i],*,cable_i]), yrange=[0,.01], xrange=[100, 7000], color='grey', /xlog
      endfor
      cgoplot, x_axis, abs(cal_med_fft_cableave_obsave[0,j,*,cable_i]), yrange=[0,.01], xrange=[100, 7000], color='blue', /xlog
      cgoplot, [400., 400.],[-1,1],  linestyle=2
      cgoplot, [150., 150.],[-1,1],  linestyle=2
      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
      cgPS_Open,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/selected_cal/selected_delay_ave_plots/low_pointing_'+strtrim(pointing_num[j],2)+'_cable_'+strtrim(ULONG(cable_length_ref[cable_i]),2)+'.png',/quiet,/nomatch
      cgplot, x_axis, abs(cal_low_fft_cableave[0,subb[0],*,cable_i]), yrange=[0,.05], xrange=[30, 7000], color='grey', /xlog, charsize=1,$
        ytitle='Amp', xtitle='Delay (ns)', title='Pointing ' +strtrim(pointing_num[j],2)+ ', cable ' + strtrim(ULONG(cable_length_ref[cable_i]),2) + ', Low delays'
      for obs_i=1, N_elements(obs_pointing)-1 do begin
        cgoplot, x_axis, abs(cal_low_fft_cableave[0,subb[obs_i],*,cable_i]), yrange=[0,.01], xrange=[100, 7000], color='grey', /xlog
      endfor
      cgoplot, x_axis, abs(cal_low_fft_cableave_obsave[0,j,*,cable_i]), yrange=[0,.01], xrange=[100, 7000], color='blue', /xlog
      ;cgoplot, [400., 400.],[-1,1],  linestyle=2
      cgoplot, [150., 150.],[-1,1],  linestyle=2
      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
      cgPS_Open,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/selected_cal/selected_delay_ave_plots/high_pointing_'+strtrim(pointing_num[j],2)+'_cable_'+strtrim(ULONG(cable_length_ref[cable_i]),2)+'.png',/quiet,/nomatch
      cgplot, x_axis, abs(cal_high_fft_cableave[0,subb[0],*,cable_i]), yrange=[0,.05], xrange=[30, 7000], color='grey', /xlog, charsize=1,$
        ytitle='Amp', xtitle='Delay (ns)', title='Pointing ' +strtrim(pointing_num[j],2)+ ', cable ' + strtrim(ULONG(cable_length_ref[cable_i]),2) + ', High delays'
      for obs_i=1, N_elements(obs_pointing)-1 do begin
        cgoplot, x_axis, abs(cal_high_fft_cableave[0,subb[obs_i],*,cable_i]), yrange=[0,.01], xrange=[100, 7000], color='grey', /xlog
      endfor
      cgoplot, x_axis, abs(cal_high_fft_cableave_obsave[0,j,*,cable_i]), yrange=[0,.01], xrange=[100, 7000], color='blue', /xlog
      cgoplot, [400., 400.],[-1,1],  linestyle=2
      ;cgoplot, [150., 150.],[-1,1],  linestyle=2
      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    endfor
    
  ENDFOR
  
  cal_selected = complex(FLTARR(2,6,384,6))
  for pol_i=0,1 do begin
    for cable_i=0, 5 do begin
      for poi_i=0,5 do begin
        cal_selected[pol_i,cable_i,*,poi_i] = [reform(cal_low_fft_cableave_obsave[pol_i,cable_i,selected_low_pos,poi_i]),reform(cal_med_fft_cableave_obsave[pol_i,cable_i,selected_med_pos,poi_i]),reform(cal_high_fft_cableave_obsave[pol_i,cable_i,selected_high_pos,poi_i]),$
          reform(cal_high_fft_cableave_obsave[pol_i,cable_i,selected_high_neg,poi_i]),reform(cal_med_fft_cableave_obsave[pol_i,cable_i,selected_med_neg,poi_i]),reform(cal_low_fft_cableave_obsave[pol_i,cable_i,selected_low_neg,poi_i])]
      endfor
    endfor
  endfor
  
  j=5
  colors=['blue','red','green','purple','yellow','black']
  cgPS_Open,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/selected_cal/selected_delay_ave_plots/allcables_pointing_'+strtrim(pointing_num[j],2)+'.png
  cgplot, x_axis, abs(cal_selected[0,j,*,0]), yrange=[0,.01], xrange=[100, 7000], color=colors[0], /xlog,  charsize=1,$
    ytitle='Amp', xtitle='Delay (ns)', title='Selected cal, pointing ' +strtrim(pointing_num[j],2)
  for cable_i=1,5 do cgoplot,  x_axis, abs(cal_selected[0,j,*,cable_i]), color=colors[cable_i]
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  ;updated_cal = fft(cal_selected,/inverse,dim=1)
  
  denom = (80000.*16.)^2. + (80000.)^2.
  
  for harmonic_i=1,7 do begin
    high_line = ( ((80000.*16.*float(harmonic_i))/denom)-(80000./denom) )*1E9
    low_line = ( ((80000.*16.*float(harmonic_i))/denom)+(80000./denom) )*1E9
    ;cgoplot, [high_line, high_line],[-1,1],  linestyle=2
    ;cgoplot, [low_line, low_line],[-1,1],  linestyle=2
    coarse_band_wh = where((abs(x_axis) LT low_line) AND (abs(x_axis) GT high_line),nc)
    if harmonic_i EQ 1 then coarse_bands_wh = coarse_band_wh else coarse_bands_wh = [coarse_bands_wh,coarse_band_wh]
    if harmonic_i EQ 1 then first_coarse = where(abs(x_axis) GT low_line)
  endfor
  
  cal_selected_zero=cal_selected
  cal_selected_zero[first_coarse,*]=0
  cal_selected_zero[coarse_bands_wh,*]=cal_selected[coarse_bands_wh,*]
  updated_cal = fft(cal_selected_zero,/inverse,dim=1)
  stop
  stop
  
  for tile_i=0,127 do begin
    cgPS_Open,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/selected_cal/amp'+string(strtrim(tile_i,2),FORMAT='(I03)')+'.png',/quiet,/nomatch
    cgplot, freq_arr[freq_use]/1E6, abs(updated_cal[freq_use,tile_i]) - abs((*cal_low.gain[0])[freq_use,tile_i]), $
      xtitle='Frequency (MHz)',ytitle='Amplitude difference',charsize=1, title='Selected cal minus per-freq cal, ' + strtrim(tile_i,2)
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  endfor
  
end