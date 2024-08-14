pro longrun_gain_ave

  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/beardsley_thesis_list.txt'
  readcol, filename, obsids, format='A', /silent
  
  longrun_names_match, obs_names=obs_names
  good_days = where((obs_names[*,1] NE 'Oct15') AND (obs_names[*,1] NE 'Oct31'),nc)
  obsids = obsids[good_days]
  
  n_obs = N_elements(obsids)
  
  gain = complex(FLTARR(2,n_obs,384,128))
  amp_degree = 2
  phase_degree = 1
  obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/1061316296_obs.sav','obs')
  freq_use=where((*obs.baseline_info).freq_use,nf_use)
  
  saved=1
  if ~keyword_set(saved) then begin
  
    for pol_i=0, 1 do begin
      for obs_i=0, n_obs-1 do begin
        print, obs_i
        cal = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_notileref/calibration/'+strtrim(obsids[obs_i],2)+'_cal.sav','cal')
        if ~keyword_set(cal) then begin
          print, 'Obs_i ' + strtrim(obs_i,2) + ' not defined.'
          continue
        endif
        gain_fit = fltarr(384,128)
        phase_fit = fltarr(384,128)
        
        
        gain_ref = reform(abs( (*cal.gain[pol_i])[freq_use,*] + (*cal.gain_residual[pol_i])[freq_use,*] ))
        for tile_i=0,127 do begin
        
          if (cal.amp_params[pol_i,tile_i] NE !NULL) then $
            FOR di=0L,amp_degree DO gain_fit[*,tile_i]+=(*cal.amp_params[pol_i,tile_i])[di]*findgen(384)^di
          ;fit_params1=poly_fit(freq_use[0:223],gain_ref[0:223,tile_i],1)
          ;fit_params2=poly_fit(freq_use[224:335],gain_ref[224:335,tile_i],1)
          ;gain_fit[freq_use[0]:freq_use[223],tile_i] = fit_params1[0]*findgen(freq_use[223])^0 +fit_params1[1]*findgen(freq_use[223])^1
          ;gain_fit[freq_use[224]:freq_use[335],tile_i] = fit_params2[0]*(findgen(freq_use[335] - freq_use[224]+1) + freq_use[224])^0 +fit_params2[1]*(findgen(freq_use[335] - freq_use[224]+1) + freq_use[224])^1
          if (cal.phase_params[pol_i,tile_i] NE !NULL) then $
            FOR di=0L,phase_degree DO phase_fit[*,tile_i]+=(*cal.phase_params[pol_i,tile_i])[di]*findgen(384)^di
        endfor
        
        unwrapped_phase = phunwrap(atan((*cal.gain[pol_i]) + (*cal.gain_residual[pol_i]),/phase))
        new_phase = atan((exp(Complex(0,1)*unwrapped_phase)) / (exp(Complex(0,1)*phase_fit)),/phase)
        
        gain[pol_i,obs_i,*,*] = abs((*cal.gain[pol_i]) + (*cal.gain_residual[pol_i])) / gain_fit * $
          exp(Complex(0,1)* new_phase)
          
        undefine_fhd, cal
      endfor
    endfor
    
  endif else restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/longrun_gain_dig_poi_refave_extras.sav'
  
  longrun_gain = complex(FLTARR(2,5,384,128))
  poi_name=['minustwo','minusone','zenith','plusone','plustwo']
  day_name_array=['Aug23','Aug27','Sep02','Sep04','Sep06','Sep09','Sep11','Sep13','Sep17','Sep19','Sep30','Oct02','Oct04','Oct08','Oct10','Oct23','Oct25','Oct28','Oct29','Nov18','Nov29']
  Result = FILE_SEARCH('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_notileref/calibration/', '*_cal.sav',count=cal_count)
  obs_use = file_basename(result, '_cal.sav')
  match, obs_use, reform(obs_names[good_days,0]), suba, subb
  
  indices = PTRARR(N_elements(day_name_array),/allocate)
  for day_i=0, N_elements(day_name_array)-1 do begin
    temp = where(obs_names[good_days[subb],1] EQ day_name_array[day_i],n_count)
    *indices[day_i] = temp
    if day_i EQ 0 then n_count_days = n_count else n_count_days = [n_count_days,n_count]
  endfor
  
  if keyword_set(ref_by_obs) then begin
    ref_avg = FLTARR(2,n_obs,384)
    for pol_i=0,1 do begin
      for obs_i=0,n_obs-1 do begin
        for freq_i=0, 383 do begin
          resistant_mean, reform(atan(gain[pol_i,obs_i,freq_i,*],/phase)), 2, phase_mean, /silent
          ref_avg[pol_i,obs_i,freq_i] = phase_mean
        endfor
      endfor
    endfor
    
    reref_phase = atan(gain,/phase)
    for tile_i=0,127 do reref_phase[*,*,*,tile_i] = atan(gain[*,*,*,tile_i],/phase) - ref_avg
  endif
  ref_by_day=1
  if keyword_set(ref_by_day) then begin
    reref_phase = atan(gain,/phase)
    ref_avg = FLTARR(2,N_elements(day_name_array),384)
    for pol_i=0,1 do begin
      for day_i=0,N_elements(day_name_array)-1 do begin
        day_inds = where(obs_names[good_days,1] EQ day_name_array[day_i])
        for freq_i=0, 383 do begin
          resistant_mean, reform(atan(gain[pol_i,day_inds,freq_i,*],/phase)), 2, phase_mean, /silent
          ref_avg[pol_i,day_i,freq_i] = phase_mean
        endfor
        for tile_i=0,127 do reref_phase[pol_i,day_inds,*,tile_i] = atan(gain[pol_i,day_inds,*,tile_i],/phase) - rebin(ref_avg[pol_i,day_i,*],N_elements(day_inds),384)
      endfor
    endfor
  endif
  stop
  
  
  for pol_i=0, 1 do begin
    for poi_i=0,4 do begin
      for freq_i=0, 383 do begin
        for tile_i=0,127 do begin
          poi_tiles = where(obs_names[*,2] EQ poi_name[poi_i])
          resistant_mean, reform(abs(gain[pol_i,poi_tiles,freq_i,tile_i])), 2, res_mean, /silent
          resistant_mean, reform(reref_phase[pol_i,poi_tiles,freq_i,tile_i]), 2, phase_mean, /silent
          longrun_gain[pol_i,poi_i,freq_i,tile_i] = res_mean * exp(Complex(0,1) * phase_mean)
        endfor
      endfor
    endfor
  endfor
  
  stop
  stop
  mean_phase=FLTARR(214,336) ; 64,60,23,67
  for obs_i=0, 213 do mean_phase[obs_i,*] = mean(reform(atan(gain[0,obs_i,freq_use,*],/phase)),dimension=2,/NAN)
  quick_image, mean_phase
  
  
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/norm_gain_plus_phase_poi.sav'
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/longrun_gain_plus_phase_poi.sav'
  poi_tiles_minustwo = where(obs_names[*,2] EQ poi_name[0])
  
  mean_val = mean(reform(abs(gain[0,poi_tiles_minustwo,*,*])),dimension = 1)
  err = FLTARR(384,128)
  for freq_i=0, 383 do for tile_i=0, 127 do $
    err[freq_i,tile_i] =  sqrt(total((mean_val[freq_i,tile_i] - reform(abs(gain[0,poi_tiles_minustwo,freq_i,tile_i])))^2.,1)/N_elements(poi_tiles_minustwo))
    
    
  tile_i=0
  cgplot, freq_arr, abs(gain[0,0,*,tile_i]), color='grey', xrange=[freq_arr[freq_use[0]],freq_arr[freq_use[335]]], yrange=[.8,1.2], xtitle='Freq (MHz)', ytitle='Normalized Amp',/NoData
  for obs_i=0,N_elements(poi_tiles_minustwo)-1 do cgoplot,freq_arr, abs(gain[0,poi_tiles_minustwo[obs_i],*,tile_i]), color='grey', xrange=[freq_arr[0],freq_arr[383]]
  cgoplot, freq_arr[*], abs(longrun_gain[0,0,*,tile_i]), color='blue'
  cgoplot, freq_arr, reform(abs(longrun_gain[0,0,*,tile_i])) - err[*,tile_i], color='purple'
  cgoplot, freq_arr, reform(abs(longrun_gain[0,0,*,tile_i])) + err[*,tile_i], color='purple'
  
  
  mean_val = mean(reform((atan(gain[0,poi_tiles_minustwo,*,*],/phase))) ,dimension = 1, /NAN)
  phase_err = FLTARR(384,128)
  for freq_i=0, 383 do for tile_i=0, 127 do $
    phase_err[freq_i,tile_i] =  sqrt(total((mean_val[freq_i,tile_i] - reform(phunwrap(atan(gain[0,poi_tiles_minustwo,freq_i,tile_i],/phase))))^2.,1)/N_elements(poi_tiles_minustwo))
    
  freq_not_use=where((*obs.baseline_info).freq_use EQ 0,nf_use)
  longrun_gain[*,*,freq_not_use,*]=!VALUES.F_NAN
  gain[*,*,freq_not_use,*]=!VALUES.F_NAN
  phase_err[freq_not_use,*]=!VALUES.F_NAN
  
  tile_i=30
  cgps_open,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/dig_phase_'+strtrim(tile_i)+'.png',/quiet,/nomatch
  cgplot, freq_arr, atan(gain[0,0,*,tile_i],/phase), color='grey', xrange=[freq_arr[freq_use[0]],freq_arr[freq_use[335]]], yrange=[-.5,.5], xtitle='Freq (MHz)', ytitle='Normalized Phase',$
    title='Tile ' + strtrim(tile_i,2),/NoData
  for obs_i=0,N_elements(poi_tiles_minustwo)-1 do cgoplot,freq_arr, atan(gain[0,poi_tiles_minustwo[obs_i],*,tile_i],/phase), color='grey', xrange=[freq_arr[0],freq_arr[383]]
  cgoplot, freq_arr[*], atan(longrun_gain[0,0,*,tile_i],/phase), color='blue'
  cgoplot, freq_arr, reform(atan(longrun_gain[0,0,*,tile_i],/phase)) - phase_err[*,tile_i], color='purple'
  cgoplot, freq_arr, reform(atan(longrun_gain[0,0,*,tile_i],/phase)) + phase_err[*,tile_i], color='purple'
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  
  obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/1061316296_obs.sav','obs')
  freq_arr=(*obs.baseline_info).freq/10.^6.
  freq_use = (*obs.baseline_info).freq_use
  freq_use = where(freq_use)
  
  
  for tile_i=0, 127 do begin
    cgPS_Open,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/' + string(strtrim(tile_i,2),FORMAT='(I03)') + 'png',/quiet,/nomatch
    cgplot, freq_arr[freq_use], abs(gain[0,0,freq_use,tile_i]), color='grey', xrange=[freq_arr[freq_use[0]],freq_arr[freq_use[335]]], yrange=[.5,1.5], xtitle='Freq (MHz)', ytitle='Normalized Amp', title='Longrun norm ave, xx, tile ' + strtrim(tile_i,2)
    for obs_i=1,1028 do cgoplot,freq_arr[freq_use], abs(gain[0,obs_i,freq_use,tile_i]), color='grey', xrange=[freq_arr[freq_use[0]],freq_arr[freq_use[335]]]
    cgoplot, freq_arr[freq_use], longrun_gain[0,freq_use,tile_i], color='blue'
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  endfor
  
end