pro temp2

  cal2 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/calibration/1061316176_cal.sav', 'cal')
  gain_arr = *cal2.gain[0]; + *cal2.gain_residual[0]
  cal_return = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/gain_arr_notimeavg2.sav', 'gain_arr')
  gain_arr_tau = cal_return
  
  obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/1061316176_obs.sav', 'obs')
  freq = (*obs.baseline_info).freq
  freq_use = (*obs.baseline_info).freq_use
  freq_use[*]=1
  freq_inds = where(freq_use GT 0)
  
  for tile_i=0,127 do begin
    cgPS_Open,'/nfs/mwa-00/h1/nbarry/tau_cal/sim_phase/'+string(strtrim(tile_i,2),FORMAT='(I03)')+'.png',/quiet,/nomatch
    cgplot, atan(gain_arr[freq_inds,tile_i],/phase), ytitle='Gain phase',xtitle='Freq channel', title='Tile ' + strtrim(tile_i,2) + ' Phase', yrange=[-.01,.01], xrange=[0,N_elements(freq_inds)], charsize=1,color='green'
    cgoplot, atan(gain_arr_tau[freq_inds,tile_i],/phase), color='black'
    cglegend, Title=['Freq gains','Freq gains via delay gains'], $
      Color=['green','black'],Length=.03,charsize=1,$;Psym=[2,2,2,2,2,2],Length=0.0 $
      Location=[0.56,0.87]
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  endfor
  
end