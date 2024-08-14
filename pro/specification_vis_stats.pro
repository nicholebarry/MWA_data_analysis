pro specification_vis_stats
  ;Script to try to understand how deviations in the calibration propogate to the power spectrum.

  ;******************Data read-in

  ;Simulation visibility data: Traditional cal, Averaged Cal, Perfect cal (Just xx)
  vis_overfit = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/vis_data/1061316176_vis_XX.sav','vis_ptr')
  vis_saved = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_savedfit_cal_eor_ones_short_baselines_included_zenithpointing/vis_data/1061316176_vis_XX.sav','vis_ptr')
  vis_perfect = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included/vis_data/1061316176_vis_XX.sav','vis_ptr')
  
  ;Simulation calibration data: Traditional cal, averaged cal, perfect cal (just xx)
  cal_overfit = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included/calibration/1061316176_cal.sav','cal')
  cal_saved = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_savedfit_cal_eor_ones_short_baselines_included_zenithpointing/calibration/1061316176_cal.sav','cal')
  cal_perfect = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included/calibration/1061316176_cal.sav','cal')
  
  ;Uh oh found that flagged tile in data creeping in
  (*cal_overfit.gain[0])[*,76] = 0
  ;******************Data read-in done
  
  ;******************AMP
  ;
  ;******Visibilities
  ; - [ deltab* + deltaa + deltaa deltab* ]
  vis_perf_over_diff = (*vis_perfect[0]-*vis_overfit[0])/(*vis_perfect[0])
  
  ;Get histograms of the amplitude gain deviations
  binsize = .0001
  his_perf_over_diff = histogram(abs(vis_perf_over_diff), binsize=binsize, /NAN, locations=locations )
  
  ;Display amp histogram
  x_arr_full_amp=[0,FINDGEN(N_elements(his_perf_over_diff))*binsize+binsize/2.,N_elements(his_perf_over_diff)*binsize-binsize/2.]
  y_arr_full_amp=[his_perf_over_diff[0],his_perf_over_diff,his_perf_over_diff[N_elements(his_perf_over_diff)-1]]
  
  ;Make the plot
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/specs/trad_amp_all.png',/quiet,/nomatch
  cgplot, x_arr_full_amp,y_arr_full_amp, xrange=[0,.02], xtitle='Gain deviations', ytitle='Binned Result', $
    title='Deviation Amplitudes from Traditional Cal', charsize=1,/ylog,yrange=[1,10^7.]
    
  mean_string = strtrim(string(mean(abs(vis_perf_over_diff),/NAN),Format='(F10.4)'),2)
  stddev_string = strtrim(string(stddev(abs(vis_perf_over_diff),/NAN),Format='(F10.4)'),2)
  max_number = where(his_perf_over_diff EQ max(his_perf_over_diff))
  max_string = strtrim(string(locations[max_number],Format='(F10.4)'),2)
  
  cglegend, title=['Vis deviations ($\Delta$g$\downb$$\up*$ + $\Delta$g$\downa$ + $\Delta$g$\downa$$\Delta$g$\downb$$\up*$)','mean='+ mean_string, $
    '$\sigma$='+ stddev_string, 'max=' + max_string], tcolors=['black','black','black', 'black'],length=0,location=[.5,.85],charsize=1
  ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  ;******End of Visibilities
    
  ;******Cal
  cal_perf_over_gain = reform(abs(*cal_perfect.gain[0])-abs(*cal_overfit.gain[0]), 384.*128.)
  
  ;Histogram the cal data
  binsize_cal = .0001
  his_cal_perf_over_diff = histogram(abs(cal_perf_over_gain), binsize=binsize_cal, /NAN, locations=locations )
  
  x_arr_full_calamp=[0,FINDGEN(N_elements(his_cal_perf_over_diff))*binsize+binsize/2.,N_elements(his_cal_perf_over_diff)*binsize-binsize/2.]
  y_arr_full_calamp=[his_cal_perf_over_diff[0],his_cal_perf_over_diff,his_cal_perf_over_diff[N_elements(his_cal_perf_over_diff)-1]]
  
  ;Make the plot
  ;cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/specs/trad_amp_cal.png',/quiet,/nomatch
  ;cgplot, x_arr_full_calamp,y_arr_full_calamp, xrange=[0,1], xtitle='Gain deviations', ytitle='Binned Result', $
  ;  title='Amplitude of $\Delta$g for Traditional Cal', charsize=1
  cgoplot, x_arr_full_calamp,y_arr_full_calamp, xrange=[0,1], xtitle='Gain deviations', ytitle='Binned Result', $
    title='Amplitude of $\Delta$g for Traditional Cal', charsize=1, /ylog, color='blue'
    
  resistant_mean, abs(cal_perf_over_gain), 4, res_mean
  mean_string = strtrim(string(res_mean,Format='(F10.4)'),2)
  stddev_string = strtrim(string(stddev(abs(cal_perf_over_gain),/NAN),Format='(F10.4)'),2)
  max_number = where(his_cal_perf_over_diff EQ max(his_cal_perf_over_diff))
  max_string = strtrim(string(locations[max_number],Format='(F10.4)'),2)
  
  cglegend, title=['Cal deviations ($\Delta$g)','mean='+ mean_string, '$\sigma$='+ stddev_string, 'max=' + max_string], $
    tcolors=['blue','blue','blue','blue'],length=0,location=[.5,.7],charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  
  ;
  ;*******Average
  ;
  ;*******Visibilities
  ; - [ deltab* + deltaa + deltaa deltab* ]
  vis_perf_ave_diff = (*vis_perfect[0]-*vis_saved[0])/(*vis_perfect[0])
  
  ;Get histograms of the amplitude and phase of the gain deviations
  binsize = .000005
  his_perf_ave_diff = histogram(abs(vis_perf_ave_diff), binsize=binsize, /NAN, locations=locations )
  
  ;Display amp histogram
  x_arr_full_amp=[0,FINDGEN(N_elements(his_perf_ave_diff))*binsize+binsize/2.,N_elements(his_perf_ave_diff)*binsize-binsize/2.]
  y_arr_full_amp=[his_perf_ave_diff[0],his_perf_ave_diff,his_perf_ave_diff[N_elements(his_perf_ave_diff)-1]]
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/specs/ave_amp_all.png',/quiet,/nomatch
  cgplot, x_arr_full_amp,y_arr_full_amp, xrange=[0,.0003], xtitle='Gain deviations', ytitle='Binned Result', $
    title='Deviation Amplitude from Average Cal', charsize=1, ystyle=8;,/ylog, yrange = [1,10^7.]
    
  mean_string = strtrim(string(mean(abs(vis_perf_ave_diff),/NAN),Format='(F10.6)'),2)
  stddev_string = strtrim(string(stddev(abs(vis_perf_ave_diff),/NAN),Format='(F10.6)'),2)
  max_number = where(his_perf_ave_diff EQ max(his_perf_ave_diff))
  max_string = strtrim(string(locations[max_number],Format='(F10.6)'),2)
  
  cglegend, title=['Vis deviations ($\Delta$g$\downb$$\up*$ + $\Delta$g$\downa$ + $\Delta$g$\downa$$\Delta$g$\downb$$\up*$)','mean='+ mean_string, $
    '$\sigma$='+ stddev_string, 'max=' + max_string], tcolors=['black','black','black', 'black'],length=0,location=[.5,.85],charsize=1
  ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  ;*******
    
  ;*******Saved Cal
  cal_perf_ave_gain = reform(abs(*cal_perfect.gain[0])-abs(*cal_saved.gain[0]), 384.*128.)
  binsize = .000005
  
  ;Had to modify the equation for the phase, since perfect phase is 0, and dividing by 0 was not helpful
  his_cal_perf_ave_diff = histogram(abs(cal_perf_ave_gain), binsize=binsize, /NAN, locations=locations )
  
  x_arr_full_amp=[0,FINDGEN(N_elements(his_cal_perf_ave_diff))*binsize+binsize/2.,N_elements(his_cal_perf_ave_diff)*binsize-binsize/2.]
  y_arr_full_amp=[his_cal_perf_ave_diff[0],his_cal_perf_ave_diff,his_cal_perf_ave_diff[N_elements(his_cal_perf_ave_diff)-1]]
  
  ;cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/specs/ave_amp_cal.png',/quiet,/nomatch
  cgaxis, yaxis=1, yrange=[0,10000],color='blue',/save
  cgoplot, x_arr_full_amp,y_arr_full_amp, xtitle='Gain deviations', ytitle='Binned Result',xrange=[0,.0003], $
    title='Amplitude of $\Delta$g for Average Cal', charsize=1, color='blue'
    
  mean_string = strtrim(string(mean(abs(cal_perf_ave_gain),/NAN),Format='(F10.6)'),2)
  stddev_string = strtrim(string(stddev(abs(cal_perf_ave_gain),/NAN),Format='(F10.6)'),2)
  max_number = where(his_cal_perf_ave_diff EQ max(his_cal_perf_ave_diff))
  max_string = strtrim(string(locations[max_number],Format='(F10.6)'),2)
  
  cglegend, title=['Cal deviations ($\Delta$g)','mean='+ mean_string, '$\sigma$='+ stddev_string, 'max=' + max_string], $
    tcolors=['blue','blue','blue','blue'],length=0,location=[.5,.7],charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  ;********
  
  ;******************END OF AMP
  
  ;******************PHASE
  
    ;***Vis phase
  binsize = .0001
  ;Had to modify the equation for the phase, since perfect phase is 0, and dividing by 0 was not helpful
  his_perf_over_diff_phase = histogram(atan((*vis_perfect[0]-*vis_overfit[0]),/phase), binsize=binsize, /NAN, locations=locations )
  
  ;Display phase histogram
  ;x_arr_full_phase=[0,FINDGEN(N_elements(his_perf_over_diff_phase))*binsize+binsize/2.,N_elements(his_perf_over_diff_phase)*binsize-binsize/2.]
  ;y_arr_full_phase=[his_perf_over_diff_phase[0],his_perf_over_diff_phase,his_perf_over_diff_phase[N_elements(his_perf_over_diff_phase)-1]]
  ;cgPS_Open,'/nfs/eor-00/h1/nbarry/vis_res/thesis/seti_res_vis_thesis_log_evenodd2.png',/quiet,/nomatch
  ;cgplot, x_arr_full_phase,y_arr_full_phase, xrange=[0,1], xtitle='Gain deviations', ytitle='Binned Result', $
  ;  title='Phase of $\Delta$g$\downb$$\up*$ + $\Delta$g$\downa$ + $\Delta$g$\downa$$\Delta$g$\downb$$\up*$ for Traditional Cal', charsize=1
  ;mean_string = strtrim(string(mean(atan(vis_perf_over_diff,/phase),/NAN),Format='(F10.4)'),2)
  ;stddev_string = strtrim(string(stddev(atan(vis_perf_over_diff,/phase),/NAN),Format='(F10.4)'),2)
  ;cglegend, title=['mean='+ mean_string, '$\sigma$='+ stddev_string], color=['black','black'],length=0,location=[.7,.8],charsize=1
  ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  ;*******
  
end