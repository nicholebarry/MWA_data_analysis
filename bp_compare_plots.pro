pro bp_compare_plots

  ;************Levine memo
  ;Setup is for the Levine memo
  bp_filename = filepath('bandpass_Levine.txt',root=rootdir('FHD'),subdir='instrument_config')
  readcol, bp_filename, memo_bp_gains
  
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_patti_catalog/metadata/1061316296_obs.sav'
  n_freq=384
  n_pol=2
  n_tile=128
  
  ;The memo gains have 128 fine frequency channels, so break the fine frequencies into sets that fit the data resolution
  range_memo_low = where((findgen(128) mod ULONG(obs.freq_res / 10000.)) EQ 0)
  range_memo_high = where((findgen(128) mod ULONG(obs.freq_res / 10000.)) EQ ULONG(obs.freq_res / 10000.)-1)
  
  ;Average the memo gains in sets to create the correct freq resolution
  num_bp_channels = ULONG(1.28E6 / obs.freq_res)
  freq_use = (*obs.baseline_info).freq_use
  ave_bp_gains = FLTARR(num_bp_channels)
  for channel_i=0,num_bp_channels-1 do $
    ave_bp_gains[channel_i] = median(memo_bp_gains[range_memo_low[channel_i]:range_memo_high[channel_i]])
  ;Create bandpass effects over the whole band
  ave_bp_gains_fullband = ave_bp_gains
  for band_i=1,n_freq/num_bp_channels-1 do ave_bp_gains_fullband = [ave_bp_gains_fullband,ave_bp_gains]
  ;************end of Levine memo
  
  subband_low =  where((findgen(384) mod 16) EQ 0) + 1
  subband_high =  where((findgen(385) mod 16) EQ 0) -2
  
  filenames = ['/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_bandpass_division_test/calibration/1061316296_cal.sav',$
    '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_Aug24_2014/calibration/1092938520_cal.sav',$
    '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_Sep14_2014/calibration/1094747968_cal.sav',$
    '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_Sep10_2015/calibration/1125939344_cal.sav']
  bp_day = ['Aug23_2013','Aug24_2014', 'Sep14_2014','Sep10_2015']
  colors=['green','cyan']
  mean_subband_all=FLTARR(N_elements(filenames),14*24)
  
  
  
  for bp_comparison_i=0,N_elements(filenames)-1 do begin
    cal_bp = getvar_savefile(filenames[bp_comparison_i],'cal')
    for pol_i=0,1 do *cal_bp.gain[pol_i]=*cal_bp.gain[pol_i]+*cal_bp.gain_residual[pol_i]
    
    pol=0
    
    subband_mean = mean(abs((*cal_bp.gain[pol])[subband_low[0]+3:subband_high[1]-3,0]))
    ;if bp_comparison_i EQ 0 then cgplot, abs((*cal_bp.gain[pol])[subband_low[0]:subband_high[1],0])/subband_mean,yrange=[.6,1.2], xrange=[0,13], $
    ;  title='Bandpass shape comparison, XX',ytitle='Normalized gain', xtitle='Frequency channel index',/NODATA
    
    all_subbands = FLTARR(24,128,14)
    mean_subband=FLTARR(24,14)
    
    ;!Y.OMargin = [2, 8]
    ;!X.OMargin = [2, 6]
    ;!P.Multi = [0, 6, 4]
    
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/bp_compare_'+bp_day[bp_comparison_i]+'.png',/quiet,/nomatch
    cgDisplay
    pos=cgLayout([6,4],XGap=3, YGap=3,OXMargin=[5,5], OYMargin=[5,7])
    
    for subband_i=0, 23 do begin
    
      for tile_i=0, 127 do begin
        subband_mean = mean(abs((*cal_bp.gain[pol])[subband_low[subband_i]+2:subband_high[subband_i+1]-2,tile_i]))
        ;cgoplot, abs((*cal_bp.gain[pol])[subband_low[subband_i]:subband_high[subband_i+1],tile_i])/subband_mean
        ;if tile_i EQ 2 then stop
        all_subbands[subband_i,tile_i,*]=abs((*cal_bp.gain[pol])[subband_low[subband_i]:subband_high[subband_i+1],tile_i])/subband_mean
        if max(all_subbands[subband_i,tile_i,*]) EQ 0 then all_subbands[subband_i,tile_i,*] = !Values.F_NAN
      endfor
      
      
      for freq_i=0,13 do begin
        resistant_mean,all_subbands[subband_i,*,freq_i],2, res_mean
        mean_subband[subband_i,freq_i]=res_mean
      endfor
      
      residual_gains = reform(all_subbands[subband_i,0,*]-mean_subband[subband_i,*])
      for tile_i=1, 127 do residual_gains = [residual_gains,reform(all_subbands[subband_i,tile_i,*]-mean_subband[subband_i,*])]
      
      binsize=.005
      his_subband=histogram(residual_gains,binsize=binsize,locations=locations,omax=omax,reverse_indices=ri_a, /NAN)
      x_arr=[locations[0],locations,locations[N_elements(locations)-1]]
      y_arr=[his_subband[0],his_subband,his_subband[N_elements(his_subband)-1]]
      
      ;If subband_i EQ 0 then $
      cgplot, x_arr, y_arr,/ylog, yrange=[1,10^3.], xrange=[-.2,.2], position=pos[*,subband_i], charsize=0.6,NoErase=subband_i NE 0; $
      ;else $
      ;  cgoplot, x_arr, y_arr,/ylog, yrange=[1,10^3.], xrange=[-.2,.2], position=pos[*,subband_i], charsize=1
      
      yfit = GAUSSFIT(x_arr, y_arr, coeff, NTERMS=3)
      
      cgoplot, x_arr, coeff[0]*exp(-(x_arr-coeff[1])^2./2./coeff[2]^2.), color= 'blue'
      
    endfor
    
    cgText, 0.5, 0.95, Alignment=0.5, /Normal, 'Zenith, ' + bp_day[bp_comparison_i]
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    !P.Multi = 0
    
    
    for subband_i=0,23 do if subband_i EQ 0 then mean_subband_perobs = reform(mean_subband[subband_i,*]) else mean_subband_perobs = [mean_subband_perobs,reform(mean_subband[subband_i,*])]
    mean_subband_all[bp_comparison_i,*] = mean_subband_perobs
    
  endfor
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/all_bp_compare.png',/quiet,/nomatch
  cgplot, mean_subband_all[0,*], yrange=[.8,1.65], xrange=[0,336], ytitle = 'Normalized gain, staggered', xtitle='Frequency channel', title ='Normalized unfit gains for zenith observations', charsize=1
  cgoplot, mean_subband_all[1,*]+.2,color='blue'
  cgoplot, mean_subband_All[2,*]+.4, color='forest green'
  cgoplot, mean_subband_All[3,*]+.6, color='purple'
  
  cgText, 0.8, 0.15, Alignment=0.5, /Normal, bp_day[0]
  cgText, 0.8, 0.32, Alignment=0.5, /Normal, bp_day[1], color='blue'
  cgText, 0.8, 0.5, Alignment=0.5, /Normal, bp_day[2], color='forest green'
  cgText, 0.8, 0.685, Alignment=0.5, /Normal, bp_day[3], color='purple'
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  
  ;   cgPS_Open,'/nfs/eor-00/h1/nbarry/all_bp_compare_residuals.png',/quiet,/nomatch
  mean_subband_all[0,*]=mean_subband_all[0,*]-ave_bp_gains_fullband[where(freq_use)]
  mean_subband_all[1,*]=mean_subband_all[1,*]-ave_bp_gains_fullband[where(freq_use)]
  mean_subband_all[2,*]=mean_subband_all[2,*]-ave_bp_gains_fullband[where(freq_use)]
  mean_subband_all[3,*]=mean_subband_all[3,*]-ave_bp_gains_fullband[where(freq_use)]
  
  large_trans = findgen(384)
  large_trans[*]=-1
  large_trans[where(freq_use)] = findgen(336)
  
  tails_low=[subband_low + 1];,subband_high[1:24] - 1]
  tails_low=tails_low[sort(tails_low)]
  tails_high=[subband_high[1:24] - 1]
  tails_high=tails_high[sort(tails_high)]
  
  tails_one_more=[subband_low + 2,subband_high[1:24] - 2]
  tails_one_more=tails_one_more[sort(tails_one_more)]
  
  mean_subband_tails=FLTARR(2,4,24)
  mean_subband_tails_one_more=FLTARR(4,48)
  
  mean_subband_tails[0,0,*]= mean_subband_all[ 0,large_trans[tails_low] ]
  mean_subband_tails[0,1,*]= mean_subband_all[ 1,large_trans[tails_low] ]
  mean_subband_tails[0,2,*]= mean_subband_all[ 2,large_trans[tails_low] ]
  mean_subband_tails[0,3,*]= mean_subband_all[ 3,large_trans[tails_low] ]
  mean_subband_tails[1,0,*]= mean_subband_all[ 0,large_trans[tails_high] ]
  mean_subband_tails[1,1,*]= mean_subband_all[ 1,large_trans[tails_high] ]
  mean_subband_tails[1,2,*]= mean_subband_all[ 2,large_trans[tails_high] ]
  mean_subband_tails[1,3,*]= mean_subband_all[ 3,large_trans[tails_high] ]
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/all_bp_compare_residuals_start.png',/quiet,/nomatch
  cgplot, mean_subband_tails[0,0,*],ytitle = 'Normalized residual gain', xtitle='Frequency channel', $
    title ='Fractional difference, unfit gains for zenith observations, first channel', charsize=1, linestyle='none',psym=1
  cgoplot, mean_subband_tails[0,1,*],color='blue',linestyle='none',psym=1
  cgoplot, mean_subband_tails[0,2,*], color='forest green',linestyle='none',psym=1
  cgoplot, mean_subband_tails[0,3,*], color='purple',linestyle='none',psym=1
  cglegend, title=[bp_day[0],bp_day[1],bp_day[2],bp_day[3]], color=['black','blue','forest green','purple'],length=0,psyms=[1,1,1,1],location=[.15,.85]
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
    cgPS_Open,'/nfs/eor-00/h1/nbarry/all_bp_compare_residuals_end.png',/quiet,/nomatch
  cgplot, mean_subband_tails[1,0,*],ytitle = 'Normalized residual gain', xtitle='Frequency channel', $
    title ='Fractional difference, unfit gains for zenith observations, last channel', charsize=1, linestyle='none',psym=1
  cgoplot, mean_subband_tails[1,1,*],color='blue',linestyle='none',psym=1
  cgoplot, mean_subband_tails[1,2,*], color='forest green',linestyle='none',psym=1
  cgoplot, mean_subband_tails[1,3,*], color='purple',linestyle='none',psym=1
  cglegend, title=[bp_day[0],bp_day[1],bp_day[2],bp_day[3]], color=['black','blue','forest green','purple'],length=0,psyms=[1,1,1,1],location=[.15,.35]
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  
  cgoplot, mean_subband_all[0,*], color=colors[0]
  cgoplot, mean_subband_all[1,*], color=colors[1]
  
  cgoplot, ave_bp_gains[1:14], color='yellow'
  cglegend, title=['Mean of raw data, 2013','Mean of raw data, 2014','Levine memo'], color=[colors[0],colors[1],'yellow'],length=.05,location=[.35,.25]
  ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  
  stop
end