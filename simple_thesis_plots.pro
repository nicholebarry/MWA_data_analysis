pro simple_thesis_plots

  global_bp=1
  if keyword_set(global_bp) then begin
    obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/1061316296_obs.sav','obs')
    cal = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/calibration/1061316296_cal.sav','cal')
    
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
    
    cgps_open,'/nfs/eor-00/h1/nbarry/test_global_bp_plot2.pdf',/quiet,/nomatch
    cgplot, bandpass_arr[0,*] / 1E6,temp2[0,*,0],color = 'grey' ,xrange = [bandpass_arr[0,0] / 1E6 - 1.,bandpass_arr[0,383] / 1E6 + 1.], $
      yrange=[.8,1.2], ytitle='Normalized calibration amplitude', xtitle = 'Frequency (MHz)', title = 'Global bandpass for zenith observation 8/23/2013', $
      charsize=1, aspect = .5
    for tile_i=1,nt_use-1 do cgoplot, bandpass_arr[0,*] / 1E6,temp2[0,*,tile_i],color = 'grey'
    for tile_i=0,nt_use-1 do cgoplot, bandpass_arr[0,*] / 1E6,temp2[1,*,tile_i],color = 'grey'
    cgoplot, bandpass_arr[0,*] / 1E6,bandpass_arr[1,*], color='blue', xrange = [bandpass_arr[0,0] / 1E6 - 1.,bandpass_arr[0,383] / 1E6 + 1.], $
      yrange=[.8,1.2], ytitle='Normalized calibration amplitude', xtitle = 'Frequency (MHz)', title = 'Global bandpass for zenith observation 8/23/2013', $
      charsize=1, aspect = .5
    cgoplot,bandpass_arr[0,*] / 1E6,bandpass_arr[2,*], color='firebrick'
    cglegend, title = ['XX','YY'], color=['blue','firebrick'], location=[.73,.73], charsize=1
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
    
  endif
  
end
