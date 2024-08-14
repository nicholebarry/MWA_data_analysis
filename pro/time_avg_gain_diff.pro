pro time_avg_gain_diff
  cal2 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_notimeavg/calibration/1061316176_cal.sav', 'cal')
  gain_arr = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_notimeavg_delaymodel/calibration/1061316176_cal.sav', 'cal')
  
  ;for tile_i=0,127 do begin
  ;  cgPS_Open,'/nfs/mwa-00/h1/nbarry/time_avg_gain/delaymodel/diff_'+string(strtrim(tile_i,2),FORMAT='(I03)')+'.png',/quiet,/nomatch
  ;  cgplot, abs((*cal2.gain[0])[*,tile_i]+(*cal2.gain_residual[0])[*,tile_i]) - abs((*gain_arr.gain[0])[*,tile_i]+(*gain_arr.gain_residual[0])[*,tile_i]), yrange=[-.01,.01]
  ;  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  ;endfor
  ;stop
  gain_arr_tau = fft(((*gain_arr.gain[0])+(*gain_arr.gain_residual[0])),dim=1)
  cal2_tau = fft(((*cal2.gain[0])+(*cal2.gain_residual[0])),dim=1)
  
  gain_arr_tau_fold = (gain_arr_tau[0:384/2-1,*])+(reverse(gain_arr_tau[384/2:383,*]))
  cal2_tau_fold = (cal2_tau[0:384/2-1,*])+(reverse(cal2_tau[384/2:383,*]))
  
  diff_tau_fold = (gain_arr_tau_fold-cal2_tau_fold)
  
  x_axis = FINDGEN(384/2)*(1./80000.)*1E6
  
  for tile_i=0,127 do begin
  
    cgPS_Open,'/nfs/mwa-00/h1/nbarry/time_avg_gain/delaymodel/fft_diff_'+string(strtrim(tile_i,2),FORMAT='(I03)')+'.png',/quiet,/nomatch
    cgplot, x_axis, diff_tau_fold[*,tile_i]*384., yrange = [-0.3,0.3], ytitle = 'arbitrary amplitude',xtitle='ns',title='Diff of fft of amp in no delay filter and delay filter gains, tile ' + strtrim(tile_i,2), charsize=1
    
    denom = (80000.*16.)^2. + (80000.)^2.
    
    for harmonic_i=1,4 do begin
    
      high_line = ( ((80000.*16.*float(harmonic_i))/denom)-(80000./denom) )*1E9
      low_line = ( ((80000.*16.*float(harmonic_i))/denom)+(80000./denom) )*1E9
      cgoplot, [high_line, high_line],[-1,1],  linestyle=2
      cgoplot, [low_line, low_line],[-1,1],  linestyle=2
    endfor
    
    
    
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  endfor
end