pro baseline_plotter

  restore, '/nfs/eor-00/h1/nbarry/vis_cal_diff.sav'
  params=getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_beamperchannel_novisflagbasic_modelnoflag_noeor_beamperchannelforreal_test/metadata/1061316176_params.sav','params')
  obs_cal=getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_beamperchannel_novisflagbasic_modelnoflag_noeor_beamperchannelforreal_test/metadata/1061316176_obs.sav','obs')
  
  ;baseline_chosen=FLTARR(384,100)
  
  for freq_i=0,383 do begin
    wh_cal=where(abs(vis_cal_diff[freq_i,*]) GT 60,n_count)
    
    if n_count GT 0 then begin
      uv_hypotenuse=sqrt(params.uu[wh_cal]^2+params.vv[wh_cal]^2)
      baseline_chosen_temp=(*obs_cal.baseline_info).freq[freq_i]*uv_hypotenuse
      ;baseline_chosen[freq_i, 0:(size(baseline_chosen_temp))[1]-1]=baseline_chosen_temp
      if freq_i EQ 0 then baseline_chosen=baseline_chosen_temp else baseline_chosen=[baseline_chosen,baseline_chosen_temp]
    ;if freq_i EQ 0 then uu_chosen=params.uu[wh_cal] else
    endif
  endfor
  
  
  binsize=5.
  result=histogram(baseline_chosen, binsize=binsize,/NAN, reverse_indices=ri, locations=locations, omax=omax)
  y_arr=[result[0], result, result[N_elements(result)-1]]
  x_arr=[locations[0],locations+binsize/2,omax]
  cgplot, x_arr,y_arr,psym=10
  stop
  
end