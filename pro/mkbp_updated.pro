pro mkbp_updated

  half_band=1
  if keyword_set(half_band) then n_freq=256 else n_freq=384
  
  day=['Aug23','Aug27','Sep02','Sep04','Sep06','Sep09','Sep11','Sep13','Sep17','Sep19','Sep30','Oct02','Oct04','Oct08','Oct10','Oct15','Oct23','Oct25','Oct29','Oct30','Nov17','Nov18','Nov29','Nov30']
  day_num=N_elements(day)
  
  ;Oct31 is bad.
  
  pointing_num=[-2,-1,0,1,2]
  pointing_name=['-2','-1','0','1','2']
  
  parsednames=STRARR(day_num,N_elements(pointing_num))
  for day_i=0,day_num-1 do parsednames[day_i,*]=[day[day_i]+'_minustwo',day[day_i]+'_minusone',day[day_i]+'_zenith',day[day_i]+'_plusone',day[day_i]+'_plustwo']
  
  restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/1061316296_obs.sav'
  tile_names=ULONG((*obs.baseline_info).tile_names)
  
  amp_mean_obs = FLTARR(n_freq, 17, N_elements(pointing_num), 2, 128)
  phase_mean_obs = FLTARR(n_freq, 17, N_elements(pointing_num), 2, 128)
  
  ;***************Loop to get the obsids of each day/pointing of the longrun
  obsid_day_pointing=PTRARR((size(pointing_num))[1],day_num,/allocate)
  obsid_day_pointing_match=PTRARR((size(pointing_num))[1],day_num,/allocate)
  obsid_count=INTARR((size(pointing_num))[1])
  
  FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
  
    FOR day_i=0,day_num-1 DO BEGIN
    
      filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/' + parsednames[day_i,j] + '.txt'
      obs_array='empty'
      readcol, filename, obs_array, format='A', /silent
      filename_uncut='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/uncut/' + parsednames[day_i,j] + '.txt'
      
      If obs_array[0] NE 'empty' Then begin
        textfast,obs_array,/read,file_path=filename,string=1
        *obsid_day_pointing[j, day_i]=(obs_array[*])
        
        ;Also read in an uncut set for the pointing and match them. This will LST match obsids later for diagnostic testing
        textfast,obs_array_uncut,/read,file_path=filename_uncut,string=1
        match, obs_array, obs_array_uncut, suba, subb
        *obsid_day_pointing_match[j,day_i] = subb
      endif
      
    ENDFOR ; end day for
    
  endfor ;end pointing for
  ;***************End of loop to get the obsids of each day/pointing of the longrun
  
  ;****************************Read in the polyfit calibration solutions, and make the raw gains minus bp sol if specified
  ;gain=complex(DBLARR(2,day_num,384,128,(size(pointing_num))[1]))
  gain=complex(DBLARR(2,day_num,n_freq,128,(size(pointing_num))[1],20))
  polygain=complex(DBLARR(2,day_num,n_freq,128,(size(pointing_num))[1],20))
  phase_mean_array=complex(DBLARR(2,day_num,n_freq,128,(size(pointing_num))[1],20))
  if keyword_set(half_band) then amp_fine_struct_match=FLTARR(2,day_num,224,128,(size(pointing_num))[1],20) else amp_fine_struct_match=FLTARR(2,day_num,336,128,(size(pointing_num))[1],20)
  
  ;restore_poly=1
  if keyword_set(restore_poly) then restore, '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/polygain.sav'
  
  for j=0,(size(pointing_num))[1]-1 do begin
    ;FOR j=2,2 DO BEGIN
    print, 'Reading in pointing ' + pointing_name[j]
    for day_i=0,day_num-1 do begin
    
      If *obsid_day_pointing[j, day_i] NE !NULL then begin
        for obs_i=0, N_elements(*obsid_day_pointing[j, day_i])-1 do begin
        
          filename= '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/calibration/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_cal.sav'
          if ~file_test(filename) then continue
          if keyword_set(restored_obs) then restored_obs = [restored_obs,(*obsid_day_pointing[j, day_i])[obs_i]] else restored_obs = (*obsid_day_pointing[j, day_i])[obs_i]
          restore,filename
          filename= '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/metadata/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_obs.sav'
          restore,filename
          
          if keyword_set(auto_cal) AND day_i EQ 0 then begin
            filename= '/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_Aug2017_autocal2/calibration/'+(*obsid_day_pointing[j, day_i])[obs_i]+'_cal.sav'
          endif
          
          for pol_i=0,1 do begin
            *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
            if keyword_set(half_band) then *cal.gain[pol_i] = (*cal.gain[pol_i])[0:255,*]
          endfor
          
          ;Only analyze the lower part of the band, avoiding all data after the digital gain jump
          if keyword_set(half_band) then begin
            obs.n_freq = 256
            cal.n_freq = 256
            cal.freq = cal.freq[0:255]
            (*obs.baseline_info).freq_use[256:383] = 0
          endif
          
          if keyword_set(restore_poly) then begin
            for pol_i=0,1 do *cal.gain[pol_i] = *cal.gain[pol_i] * weight_invert(polygain[pol_i,day_i,*,*,j,obs_i])
          endif
          
          ;cal_bandpass=vis_cal_bandpass(cal,obs,cal_remainder=cal_remainder,saved_run_bp=0,cable_bandpass_fit=1)
          if ~keyword_set(restore_poly) then begin
            cal.amp_degree = 4
            cal_polyfit = vis_cal_polyfit(cal,obs);, jump_longrun=1) ;(cal_remainder,obs,jump_longrun=1)
            ;cal=vis_cal_combine(cal_polyfit,cal_bandpass)
            polygain[0,day_i,*,*,j,obs_i]=*cal_polyfit.gain[0]
            polygain[1,day_i,*,*,j,obs_i]=*cal_polyfit.gain[1]
          endif else cal = cal_bandpass
          
          gain[0,day_i,*,*,j,obs_i]=*cal.gain[0]
          gain[1,day_i,*,*,j,obs_i]=*cal.gain[1]
        endfor
        
      endif
      
    endfor
    
    
    freq_use = where((*obs.baseline_info).freq_use)
    
    amp=1
    if keyword_set(amp) then begin
      ;Create just fine structure (traditional polyfitting) for averaging purposes, amplitude only
      amp_fine_struct = abs(gain * weight_invert(polygain))
      zeros_amp = where(amp_fine_struct EQ 0, n_count)
      if n_count GT 0 then amp_fine_struct[zeros_amp] = !VALUES.F_NAN ;resistant_mean doesn't use NaN, force 0's to be NaN
      
      ;cgPS_Open,'/nfs/eor-00/h1/nbarry/obscheck.png',/quiet,/nomatch
      
      ;Create template for the plot
      obs_num=INTARR(day_num)
      for day_i=0, day_num-1 do obs_num[day_i] = N_elements(*obsid_day_pointing[0,day_i])
      ;cgplot, amp_fine_struct[0,0,freq_use,2,0,0], yrange=[-2.2,1.2], /nodata, $
      ;  title = 'Tile 2, all day averages sep by obsid count', xrange=[0,N_elements(freq_use)], $
      ;  xtitle='Unflagged freq index, no dig jump', ytitle='Norm fine struct'
      
      ;Loop over all days and plot fine structure with matched LSTs in a group. This explores
      ;how consistent obsids are and how they differ by LST. There should be differences due to
      ;differing modulation contamination.
      for obs_i=0, max(obs_num) - 1 do begin
      
        for day_i=0, day_num-1 do begin
          match_logic = where(*obsid_day_pointing_match[j,day_i] EQ obs_i,n_count)
          
          if n_count GT 0 then begin
            ;  cgoplot, amp_fine_struct[0,day_i,freq_use,2,0,match_logic] - (obs_i*.2)
            amp_fine_struct_match[*,day_i,*,*,j,obs_i] = reform(amp_fine_struct[*,day_i,freq_use,*,j,match_logic])
          endif
        endfor
        
        zeros_amp = where(amp_fine_struct_match EQ 0, n_count)
        if n_count GT 0 then amp_fine_struct_match[zeros_amp] = !VALUES.F_NAN ;resistant_mean doesn't use NaN, force 0's to be NaN
        
        for tile_j=0,127 do begin
          for pol_j=0,1 do begin
            resistant_mean, reform(amp_fine_struct_match[pol_j,*,*,tile_j,j,obs_i]), 3, amp_mean, dimension=1
            ;cgoplot, amp_mean - (obs_i*.2), color='green'
            amp_mean_obs[freq_use, obs_i, j, pol_j, tile_j] = amp_mean ;mean for later, need weights though
          endfor
        endfor
      endfor
      
    ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
    endif
  endfor ;end pointing for
  
  ;phases=1
  if keyword_set(phases) then begin
    ;Create just fine structure for averaging purposes, phases only, using mean reference
    unwrapped_phase = phunwrap(atan(gain,/phase))
    new_phase = atan((exp(Complex(0,1)*unwrapped_phase)) / (exp(Complex(0,1)*atan(polygain,/phase))),/phase)
    zeros_phase = where(new_phase EQ 0, n_count)
    new_phase[zeros_phase]=!VALUES.F_NAN
    phase_mean =  reform(mean(mean(mean(new_phase,dimension=4,/NAN),dimension=4,/NAN),dimension=4,/NAN)) ;take the mean over the tiles, pointings, and obsids. perhaps day too? something to test
    for tile_i=0,127 do for pointing_i=0,4 do for obs_i=0,19 do phase_mean_array[*,*,*,tile_i,pointing_i,obs_i] = phase_mean
    
    phase_fine = atan((exp(Complex(0,1)*reform(atan(gain,/phase)))) / (exp(Complex(0,1)*(atan(polygain,/phase) + phase_mean_array))),/phase)
    zeros = where( atan(gain,/phase) EQ 0, n_count )
    if n_count GT 0 then phase_fine[zeros] = 0.
    phase_fine_match = FLTARR((size(phase_fine))[1:(size(phase_fine))[0]])
    
    for j=0,(size(pointing_num))[1]-1 do begin
      ;cgPS_Open,'/nfs/eor-00/h1/nbarry/obscheck.png',/quiet,/nomatch
    
    
      ;Create template for the plot
      obs_num=INTARR(day_num)
      for day_i=0, day_num-1 do obs_num[day_i] = N_elements(*obsid_day_pointing[j,day_i])
      ;cgplot, amp_fine_struct[0,0,freq_use,2,0,0], yrange=[-2.2,1.2], /nodata, $
      ;  title = 'Tile 2, all day averages sep by obsid count', xrange=[0,N_elements(freq_use)], $
      ;  xtitle='Unflagged freq index, no dig jump', ytitle='Norm fine struct'
      
      ;Loop over all days and plot fine structure with matched LSTs in a group. This explores
      ;how consistent obsids are and how they differ by LST. There should be differences due to
      ;differing modulation contamination.
      for obs_i=0, max(obs_num) - 1 do begin
      
        for day_i=0, day_num-1 do begin
          match_logic = where(*obsid_day_pointing_match[j,day_i] EQ obs_i,n_count)
          
          if n_count GT 0 then begin
            ;  cgoplot, amp_fine_struct[0,day_i,freq_use,2,0,match_logic] - (obs_i*.2)
            phase_fine_match[*,day_i,*,*,j,obs_i] = reform(phase_fine[*,day_i,*,*,j,match_logic])
          endif
        endfor
        
        zeros_amp = where(phase_fine_match EQ 0, n_count)
        if n_count GT 0 then phase_fine_match[zeros_amp] = !VALUES.F_NAN ;resistant_mean doesn't use NaN, force 0's to be NaN
        
        for tile_j=0,127 do begin
          for pol_j=0,1 do begin
            resistant_mean, reform(phase_fine_match[pol_j,*,*,tile_j,j,obs_i]), 3, phase_mean, dimension=1
            ;cgoplot, amp_mean - (obs_i*.2), color='green'
            phase_mean_obs[*, obs_i, j, pol_j, tile_j] = reform(phase_mean)
          endfor
        endfor
      endfor
      
    ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
    endfor
    
    freq_not_use = where((*obs.baseline_info).freq_use EQ 0)
    phase_mean_obs[freq_not_use,*,*,*,*]=0
    
    zeros = where(phase_mean_obs EQ 0,n_count)
    if n_count GT 0 then phase_mean_obs[zeros] = !VALUES.F_NAN
    phase_mean_pointing = mean(phase_mean_obs,dimension=2,/NAN)
    
  ;    for tile_j=0, 127 do begin
  ;      cgPS_Open,'/nfs/eor-00/h1/nbarry/longrun_ave_plots/phases/'+strtrim(string(tile_j, FORMAT='(I03)'),2)+'_allpointingXX.png',/quiet,/nomatch
  ;      cgplot, cal.freq[0:n_freq-1]/1E6, phase_mean_pointing[*,0,0,tile_j]+1., xtitle='Frequency (MHz)', ytitle='Norm phase, shifted', charsize=1, title='Averaged Tile '+ strtrim(string(tile_j),2) +', XX', yrange=[0,1.2]
  ;      for j=1,(size(pointing_num))[1]-1 do begin
  ;        cgoplot, cal.freq[0:n_freq-1]/1E6, phase_mean_pointing[*,j,0,tile_j]+1-(.2*j)
  ;      endfor
  ;      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  ;    endfor
    
  endif
  
  
  
  zeros_amp2 = where(amp_mean_obs EQ 0, n_count)
  amp_mean_obs2 = amp_mean_obs
  if n_count GT 0 then amp_mean_obs2[zeros_amp2] = !VALUES.F_NAN
  amp_mean_obs3 = mean(amp_mean_obs2, dimension=2,/NAN)
  ;freq_not_use = where((*obs.baseline_info).freq_use EQ 0)
  ;amp_mean_obs3[freq_not_use,*,*,*] = 0
  
  tau=1
  if keyword_set(tau) then begin
    freq_not_use = where((*obs.baseline_info).freq_use EQ 0)
    if (n_freq mod 2) EQ 0 then shift_num = n_freq/2-1 else shift_num = floor(n_freq/2.)
    T= (80000.)
    X = FINDGEN((n_freq - 1)/2) + 1
    x_axis = shift([0.0, X, n_freq/2, -n_freq/2 + X]/(n_freq*T)*1E9,shift_num)
    
    if keyword_set(amp) then begin
      infs = where(finite(amp_mean_obs3) EQ 0,n_count)
      if n_count GT 0 then amp_mean_obs3[infs] = 0
      amp_mean_obs3 = amp_mean_obs3 -1.
      amp_mean_obs3[freq_not_use,*,*,*] = 0
      fft_amp_mean = shift(fft(amp_mean_obs3,dimension=1),shift_num,0,0,0)
      
      amp_fine = FLTARR(2,24,n_freq,128,(size(pointing_num))[1],20)
      amp_fine[*,*,freq_use,*,*,*] = amp_fine_struct_match[*,*,freq_use,*,*,*]
      infs = where(finite(amp_fine) EQ 0, n_count)
      if n_count GT 0 then amp_fine[infs]=0
      
      amp_fine = amp_fine-1.
      amp_fine[*,*,freq_not_use,*,*,*]=0
      fft_amp = shift(fft(amp_fine,dimension=3),0,0,shift_num,0,0,0)
      
      if keyword_set(auto_cal) then begin
      
      endif
      
    endif
    
    if keyword_set(phases) then begin
      infs = where(finite(phase_mean_pointing) EQ 0)
      phase_mean_pointing[infs] = 0
      phase_mean_pointing[freq_not_use,*,*,*] = 0
      fft_phase_mean = shift(fft(phase_mean_pointing,dimension=1),shift_num,0,0,0)
      
      infs = where(finite(phase_fine) EQ 0, n_count)
      if n_count GT 0 then phase_fine[infs]=0
      
      phase_fine[*,*,freq_not_use,*,*,*]=0
      fft_phase = shift(fft(phase_fine,dimension=3),0,0,shift_num,0,0,0)
    endif
  endif
  
  
  
  hist=1
  if keyword_set(hist) then begin
  
    if keyword_set(phases) then begin
      infs = where(finite(phase_fine) EQ 0,n_count)
      if n_count GT 0 then phase_fine[infs]=0
      phase_fine_tot = total(phase_fine, 3)
      phase_fine_tot = reform(phase_fine_tot[0,*,*,*,*])
      zeros = where(phase_fine_tot EQ 0)
      dims = l_getdim(phase_fine_tot,zeros, mindim=1)
      for i=0, (size(dims))[2] - 1 do fft_phase[0,reform(dims[0,i]),*,reform(dims[1,i]),reform(dims[2,i]),reform(dims[3,i])]=!VALUES.F_NAN + Complex(0,1)*!VALUES.F_NAN
    endif
    if keyword_set(amp) then begin
      infs = where(finite(amp_fine) EQ 0,n_count)
      if n_count GT 0 then amp_fine[infs]=0
      amp_fine_tot = total(amp_fine, 3)
      amp_fine_tot = reform(amp_fine_tot[0,*,*,*,*])
      zeros = where(amp_fine_tot EQ -224., n_count)
      if n_count GT 0 then begin
        dims = l_getdim(amp_fine_tot,zeros, mindim=1)
        for i=0, (size(dims))[2] - 1 do fft_amp[0,reform(dims[0,i]),*,reform(dims[1,i]),reform(dims[2,i]),reform(dims[3,i])]=!VALUES.F_NAN + Complex(0,1)*!VALUES.F_NAN
      endif
    endif
    
    x_axis= FLTARR(401)
    for i=0,400 do x_axis[i]=float(i)*.01 - 2.
    x_axis2= FLTARR(401)
    for i=0,400 do x_axis2[i]=float(i)*.01
    x_axis3= FLTARR(601)
    for i=0,600 do x_axis3[i]=float(i)*.01 - 3.
    y_axis = shift([0.0, X, n_freq/2, -n_freq/2 + X]/(n_freq*T)*1E9,shift_num)
    
    for tile_i=0, 127 do begin
      result = LONG(FLTARR(401,256))
      result2 = LONG(FLTARR(81,256))
      error_full = FLTARR(401)
      sigma_full = FLTARR(401)
      for freq_i=0, 255 do begin
        result[*,freq_i] = histogram(reform(real_part(256.*fft_amp[0,*,freq_i,tile_i,*,*])),locations=locations,omax=omax,binsize=.01,/NAN,max=2.,min=-2.)
        result2[*,freq_i] = histogram(reform(real_part(256.*fft_amp[0,*,freq_i,tile_i,*,*])),locations=locations,omax=omax,binsize=.05,/NAN,max=2.,min=-2.)
        fit= gaussfit(locations, result2[*,freq_i],fit_params, sigma=temp, nterms=3)
        error_full[freq_i] = temp[2]
        sigma_full[freq_i] = fit_params[2]
      endfor
      cgPS_Open,'/nfs/mwa-00/h1/nbarry/longrun_ave_plots/fft/hist/4panel/amp/4panel_tile_'+strtrim(string(tile_i),2)+'_XX.png',/quiet,/nomatch
      quick_image, result[*,128:255], x_axis, y_axis[128:255], data_range=[0,50], data_aspect=1, savefile='/nfs/mwa-00/h1/nbarry/longrun_ave_plots/fft/hist/4panel/amp/4panel_tile_'+strtrim(string(tile_i),2)+'_XX', /png,$
        xtitle='Real part of calibration sol', ytitle='Time', title='Tile ' + strtrim(string(tile_i),2) + ', XX amp solution', $
        charsize=1, cb_title='counts', no_ps_close=1, multi_pos = [0,.5,.5,1]
      for poi_i=0, 4 do cgoplot, 256.*real_part(fft_amp_mean[128:255,poi_i,0,tile_i]),y_axis[128:255]
      ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
      ;result = LONG(FLTARR(401,256))
      ;for freq_i=0, 255 do result[*,freq_i] = histogram(reform(abs(256.*fft_phase[0,*,freq_i,tile_i,*,*])),locations=locations,omax=omax,binsize=.01,/NAN,max=4.,min=0)
      ;cgPS_Open,'/nfs/mwa-00/h1/nbarry/longrun_ave_plots/fft/hist/tile_'+strtrim(string(tile_i),2)+'_XX.png',/quiet,/nomatch
      ;quick_image, result[*,128:255], x_axis2, y_axis[128:255], data_range=[0,50], data_aspect=1,$
      ;  xtitle='Amplitude of calibration sol', ytitle='Time', title='Tile ' + strtrim(string(tile_i),2) + ', XX phase solution', $
      ;  charsize=1, cb_title='counts', no_ps_close=1, multi_pos = [0,0,.5,.5],/noerase
      ;for poi_i=0, 4 do cgoplot, 256.*abs(fft_phase_mean[128:255,poi_i,0,tile_i]),y_axis[128:255]
      cgplot, sigma_full[128:255],  y_axis[128:255], position = [0,.1,.5,.45],/noerase, xtitle = 'Sigma of real part distribution', ytitle='Time', title = 'Tile ' + strtrim(string(tile_i),2) + ', XX phase sigma',$
        aspect=1, charsize=.7,  yrange=[0,y_axis[255]]
      ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
        
      result = LONG(FLTARR(401,256))
      result2 = LONG(FLTARR(81,256))
      error_full = FLTARR(401)
      sigma_full = FLTARR(401)
      for freq_i=0, 255 do begin
        result[*,freq_i] = histogram(reform(imaginary(256.*fft_amp[0,*,freq_i,tile_i,*,*])),locations=locations,omax=omax,binsize=.01,/NAN,max=2.,min=-2.)
        result2[*,freq_i] = histogram(reform(imaginary(256.*fft_amp[0,*,freq_i,tile_i,*,*])),locations=locations,omax=omax,binsize=.05,/NAN,max=2.,min=-2.)
        fit= gaussfit(locations, result2[*,freq_i],fit_params, sigma=temp,nterms=3)
        error_full[freq_i] = temp[2]
        sigma_full[freq_i] = fit_params[2]
      endfor
      ;cgPS_Open,'/nfs/mwa-00/h1/nbarry/longrun_ave_plots/fft/hist/tile_'+strtrim(string(tile_i),2)+'_XX.png',/quiet,/nomatch
      quick_image, result[*,128:255], x_axis, y_axis[128:255], data_range=[0,50], data_aspect=1,$
        xtitle='Imaginary part of calibration sol', ytitle='Time', title='Tile ' + strtrim(string(tile_i),2) + ', XX amp solution', $
        charsize=1, cb_title='counts', no_ps_close=1, multi_pos = [.5,.5,1,1],/noerase
      for poi_i=0, 4 do cgoplot, 256.*imaginary(fft_amp_mean[128:255,poi_i,0,tile_i]),y_axis[128:255]
      
      
      ;      result = LONG(FLTARR(601,256))
      ;for freq_i=0, 255 do result[*,freq_i] = histogram(reform(atan(256.*fft_phase[0,*,freq_i,tile_i,*,*],/phase)),locations=locations,omax=omax,binsize=.01,/NAN,max=3.,min=-3.)
      ;cgPS_Open,'/nfs/mwa-00/h1/nbarry/longrun_ave_plots/fft/hist/tile_'+strtrim(string(tile_i),2)+'_XX.png',/quiet,/nomatch
      ;quick_image, result[*,128:255], x_axis3, y_axis[128:255], data_range=[0,50], data_aspect=1,$
      ;  xtitle='Phase of calibration sol', ytitle='Time', title='Tile ' + strtrim(string(tile_i),2) + ', XX phase solution', $
      ;  charsize=1, cb_title='counts', no_ps_close=1, multi_pos = [.5,0,1,.5],/noerase
      ;for poi_i=0, 4 do cgoplot, 256.*atan(fft_phase_mean[128:255,poi_i,0,tile_i],/phase),y_axis[128:255]
      cgplot, sigma_full[128:255], y_axis[128:255], position = [.5,.1,1,.45],/noerase, xtitle = 'Sigma of imaginary part distribution', ytitle='Time', title = 'Tile ' + strtrim(string(tile_i),2) + ', XX phase sigma',$
        aspect=1, charsize=.7, yrange=[0,y_axis[255]]
        
      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    endfor
  endif
  stop
  
  for tile_j=0, 127 do begin
    cgPS_Open,'/nfs/eor-00/h1/nbarry/longrun_ave_plots/'+strtrim(string(tile_j, FORMAT='(I03)'),2)+'_allpointingXX.png',/quiet,/nomatch
    cgplot, cal.freq[0:n_freq-1]/1E6, amp_mean_obs3[*,0,0,tile_j], xtitle='Frequency (MHz)', ytitle='Norm amplitude, shifted', charsize=1, title='Averaged Tile '+ strtrim(string(tile_j),2) +', XX', yrange=[0,1.2]
    for j=1,(size(pointing_num))[1]-1 do begin
      cgoplot, cal.freq[0:n_freq-1]/1E6, amp_mean_obs3[*,j,0,tile_j]-(.2*j)
    endfor
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  endfor
  
  
  stop
  
  freq_use = where((*obs.baseline_info).freq_use)
  
  amp_fine_struct = abs(gain * weight_invert(polygain))
  zeros_amp = where(amp_fine_struct EQ 0, n_count)
  if n_count GT 0 then amp_fine_struct[zeros_amp] = !VALUES.F_NAN
  amp_fine_struct_match=FLTARR(day_num,N_elements(freq_use),1)
  
  for day_i=0, day_num-1 do obs_num = N_elements(*obsid_day_pointing[0,day_i])
  cgplot, amp_fine_struct[0,0,freq_use,2,0,0], yrange=[-1.8,1.2], /nodata, title = 'Tile 2, all day averages sep by obsid count'
  
  for obs_i=0, max(obs_num) - 1 do begin
    for day_i=0, day_num-1 do begin
      match_logic = where(*obsid_day_pointing_match[0,day_i] EQ obs_i,n_count)
      if n_count GT 0 then begin
        cgoplot, amp_fine_struct[0,day_i,freq_use,2,0,match_logic] - (obs_i*.2)
        amp_fine_struct_match[day_i,*,*] = [[[amp_fine_struct_match[day_i,*,*]]],[[reform(amp_fine_struct[0,day_i,freq_use,2,0,match_logic])]]]
      endif
    endfor
    resistant_mean, reform(amp_fine_struct_match), 3, amp_mean, dimension=1
    cgoplot, amp_mean - (obs_i*.2), color='green'
  endfor
  
  cgoplot, abs(gain[0,0,freq_use,2,0,1] * weight_invert(polygain[0,0,freq_use,2,0,1])) - .2
  for day_i=0, day_num-1 do cgoplot, abs(gain[0,day_i,freq_use,2,0,1] * weight_invert(polygain[0,day_i,freq_use,2,0,1])) - .2
  resistant_mean, reform(amp_fine_struct[0,*,freq_use,2,0,1]), 3, amp_mean, dimension=1
  cgoplot, amp_mean - .2, color='green'
  
  cgoplot, abs(gain[0,0,freq_use,2,0,2] * weight_invert(polygain[0,0,freq_use,2,0,2])) - .4
  for day_i=0, day_num-1 do cgoplot, abs(gain[0,day_i,freq_use,2,0,2] * weight_invert(polygain[0,day_i,freq_use,2,0,2])) - .4
  resistant_mean, reform(amp_fine_struct[0,*,freq_use,2,0,2]), 3, amp_mean, dimension=1
  cgoplot, amp_mean - .4, color='green'
  
  cgoplot, abs(gain[0,0,freq_use,2,0,3] * weight_invert(polygain[0,0,freq_use,2,0,3])) - .6
  for day_i=0, day_num-1 do cgoplot, abs(gain[0,day_i,freq_use,2,0,3] * weight_invert(polygain[0,day_i,freq_use,2,0,3])) - .6
  resistant_mean, reform(amp_fine_struct[0,*,freq_use,2,0,3]), 3, amp_mean, dimension=1
  cgoplot, amp_mean - .6, color='green'
  
  cgoplot, abs(gain[0,0,freq_use,2,0,4] * weight_invert(polygain[0,0,freq_use,2,0,4])) - .8
  for day_i=0, day_num-1 do cgoplot, abs(gain[0,day_i,freq_use,2,0,4] * weight_invert(polygain[0,day_i,freq_use,2,0,4])) - .8
  resistant_mean, reform(amp_fine_struct[0,*,freq_use,2,0,4]), 3, amp_mean, dimension=1
  cgoplot, amp_mean - .8, color='green'
  
  cgoplot, abs(gain[0,0,freq_use,2,0,5] * weight_invert(polygain[0,0,freq_use,2,0,5])) - 1.
  for day_i=0, day_num-1 do cgoplot, abs(gain[0,day_i,freq_use,2,0,5] * weight_invert(polygain[0,day_i,freq_use,2,0,5])) - 1.
  resistant_mean, reform(amp_fine_struct[0,*,freq_use,2,0,5]), 3, amp_mean, dimension=1
  cgoplot, amp_mean - 1., color='green'
  
  cgoplot, abs(gain[0,0,freq_use,2,0,6] * weight_invert(polygain[0,0,freq_use,2,0,6])) - 1.2
  for day_i=0, day_num-1 do cgoplot, abs(gain[0,day_i,freq_use,2,0,6] * weight_invert(polygain[0,day_i,freq_use,2,0,6])) - 1.2
  resistant_mean, reform(amp_fine_struct[0,*,freq_use,2,0,6]), 3, amp_mean, dimension=1
  cgoplot, amp_mean - 1.2, color='green'
  
  cgoplot, abs(gain[0,0,freq_use,2,0,7] * weight_invert(polygain[0,0,freq_use,2,0,7])) - 1.4
  for day_i=0, day_num-1 do cgoplot, abs(gain[0,day_i,freq_use,2,0,7] * weight_invert(polygain[0,day_i,freq_use,2,0,7])) - 1.4
  resistant_mean, reform(amp_fine_struct[0,*,freq_use,2,0,7]), 3, amp_mean, dimension=1
  cgoplot, amp_mean - 1.4, color='green'
  
  cgoplot, abs(gain[0,0,freq_use,2,0,8] * weight_invert(polygain[0,0,freq_use,2,0,8])) - 1.6
  for day_i=0, day_num-1 do cgoplot, abs(gain[0,day_i,freq_use,2,0,8] * weight_invert(polygain[0,day_i,freq_use,2,0,8])) - 1.6
  resistant_mean, reform(amp_fine_struct[0,*,freq_use,2,0,8]), 3, amp_mean, dimension=1
  cgoplot, amp_mean - 1.6, color='green'
  
  stop
  
  pol_i=0
  tile_i=10
  j=2
  day_i=0
  
  cgplot, abs(polygain[pol_i,day_i,*,tile_i,j]),/nodata, yrange=[1,2]
  
  undefine, color_array
  cgLoadCT, 25, clip=[0,190], Ncolors=day_num+1
  Device, Decomposed=0
  for day_i=0,day_num-1 do cgoplot, abs(polygain[pol_i,day_i,*,tile_i,j]), color=day_i+1
  
  for day_i=0,day_num-1 do if ~keyword_set(color_array) then color_array=day_i+1 else color_array=[color_array,day_i+1]
  
  cgLegend, Title=day, $
    Color=color_array, Location=[0.8,0.87],charsize=1.2,VSpace=1.4,thick=2, $
    /center_sym
    
  Device, Decomposed=1
  
  
  stop
  
end