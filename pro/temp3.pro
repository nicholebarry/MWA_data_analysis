pro temp3

  gain  = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/norm_all_gain_plus_phase.sav', 'gain')
  n_obs = (size(gain))[2]

  longrun_gain  = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/longrun_gain_plus_phase.sav', 'longrun_gain')
  
  obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/1061316296_obs.sav','obs')
  freq_arr=(*obs.baseline_info).freq/10.^6.
  freq_use = (*obs.baseline_info).freq_use
  freq_not_use = where(freq_use EQ 0)
  freq_use = where(freq_use)
  
  gain[*,*,freq_not_use,*]=0
  longrun_gain[*,freq_not_use,*]=0
  
  
  for tile_i=0, 127 do begin
    cgPS_Open,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/amp_' + string(strtrim(tile_i,2),FORMAT='(I03)') + '.png',/quiet,/nomatch
    cgplot, freq_arr, abs(gain[0,0,*,tile_i]), color='grey', xrange=[freq_arr[freq_use[0]],freq_arr[freq_use[335]]], yrange=[.5,1.5], xtitle='Freq (MHz)', ytitle='Normalized Amp', title='Longrun norm ave, xx, tile ' + strtrim(tile_i,2)
    for obs_i=1,n_obs-1 do cgoplot,freq_arr, abs(gain[0,obs_i,*,tile_i]), color='grey', xrange=[freq_arr[freq_use[0]],freq_arr[freq_use[335]]]
    cgoplot, freq_arr, abs(longrun_gain[0,*,tile_i]), color='blue'
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    
        cgPS_Open,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/phase_' + string(strtrim(tile_i,2),FORMAT='(I03)') + '.png',/quiet,/nomatch
    cgplot, freq_arr, atan(gain[0,0,*,tile_i],/phase), color='grey', xrange=[freq_arr[freq_use[0]],freq_arr[freq_use[335]]], yrange=[-.5,.5], xtitle='Freq (MHz)', ytitle='Normalized Amp', title='Longrun norm ave, xx, tile ' + strtrim(tile_i,2)
    for obs_i=1,n_obs-1 do cgoplot,freq_arr, atan(gain[0,obs_i,*,tile_i],/phase), color='grey', xrange=[freq_arr[freq_use[0]],freq_arr[freq_use[335]]]
    cgoplot, freq_arr, atan(longrun_gain[0,*,tile_i],/phase), color='blue'
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  endfor
  
  
end