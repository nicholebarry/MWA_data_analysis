pro RTS_cal_solutions

  ;************begin email
  ;Email from Pietro that included the file used in this comparison
  ;I've finally managed to modify a RTS calibration solution in the way we discussed.
  ;Attached you'll find 24 Bandpass files and 24 DIJM files. You won't need the last ones,
  ;as it's just a bunch of ones and zeros, while the bandpass files have this structure:

  ;- first line represents the steps in frequency, it's internal to the RTS and I guess you can ignore it;
  ;- for each non flagged antenna you get 8 entries, the proper solution points
  ;(one for each fine channels, minus two at each edges and the central one) and the fit on those points,
  ;that is what the RTS actually uses and you get these for each polarisation (XX,XY,YX,YY), so 4x2=8 rows for each antenna;
  ;- you have amplitude and phase for each fine channel listed.
  ;************end of email

  ;First data point on x array is the tile number indexed from 1
  ;First data column on y array is the delta frequency step for that course channel in MHz
  ;data[1:26,1] ;XX amp
  ;data[1:26,2] ;XX amp fit
  ;data[1:26,3] ;XY amp
  ;data[1:26,4] ;XY amp fit
  ;data[1:26,5] ;YX amp
  ;data[1:26,6] ;XY amp fit
  ;data[1:26,7] ;YY amp
  ;data[1:26,8] ;YY amp fit
  ;they are at 32 fine frequencies per course channel (so double us) but flag 6 channels
  ;************


  ;************start of bp correction setup
  ;read in the memo bandpass gains for an individual coarse band (used in cotter)
  filename = '/nfs/mwa-00/h1/nbarry/MWA/IDL_code/MWA_data_analysis/subbandpass.txt'
  readcol, filename, memo_bp_gains
  ;read in the bandpass gains for an individual coarse band from the RTS made by Mitch
  filename = '/nfs/mwa-00/h1/nbarry/RTScal/rts_pfb_gain.txt'
  readcol, filename, rts_bp_gains
  
  ;The memo gains have 128 fine frequency channels, so break the fine frequencies into sets of 4 (for the RTS freq resolution)
  range_memo_low = where((findgen(128) mod 4) EQ 0)
  range_memo_high = where((findgen(128) mod 4) EQ 3)
  ;Average the memo gains in sets of four to create RTS freq resolution
  memo_bp_gains_averaged = FLTARR(32)
  for band_i=0, 31 do memo_bp_gains_averaged[band_i] = mean(memo_bp_gains[range_memo_low[band_i]:range_memo_high[band_i]])
  ;Stack them to create the full band of 30.72 MHz with resolution of 40 kHz
  memo_bp_gains_averaged_fullband = memo_bp_gains_averaged
  for band_i=0,22 do memo_bp_gains_averaged_fullband = [memo_bp_gains_averaged_fullband,memo_bp_gains_averaged]
  ;Stack the rts bandpass gains to create the full band of 30.72 MHz with resolution of 40 kHz
  rts_bp_gains_fullband = rts_bp_gains
  for band_i=0,22 do rts_bp_gains_fullband = [rts_bp_gains_fullband,rts_bp_gains]
  ;************end of bp correction setup
  
  ;************begin rts calibration read in
  rts_unfitted_cal=complex(FLTARR(2,384*2,128))
  rts_unfitted_cal_80kHz=complex(FLTARR(2,384,128))
  rts_fitted_cal=complex(FLTARR(2,384*2,128))
  rts_fitted_cal_80kHz=complex(FLTARR(2,384,128))
  
  ;The bandpass files are per coarse channel, and start index 0 till 23. Index 0 is the highest Hz, Index 23 is the lowest Hz
  bp_file_path=Reverse(findfile('/nfs/mwa-00/h1/nbarry/RTScal/*Bandpass*dat'))
  
  ;loop over all coarse channels
  for file_i=0, 23 do begin
    ;Read in the file given by Pietro for RTS cal. Read_col and textfast fail to read in complex values correctly given setup of file
    data_temp=read_ascii(bp_file_path[file_i], data_start=1, delimiter=',')
    data = data_temp.field01 ;the automatic structure name given by read_ascii
    
    ;Initialize the flagged number. Bart warned me that flagged tiles might not be in the file altogether, and that index would indicate that.
    flagged_num=0
    ;Create the frequency index range in the file, skipping over flagged channels, to create a fullband solution from the file
    freq_range=(file_i*32)+[2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29] ;(flagged 2 to start, 1 in middle, 2 at end)
    
    for tile_i=0,127 do begin
    
      ;if the tile is not flagged, fill out the values
      if data[0,tile_i*8+1-flagged_num*8.] EQ tile_i+1 then begin
        rts_unfitted_cal[0,freq_range,tile_i] = data[where((indgen(54) mod 2) EQ 0)+1,tile_i*8-flagged_num*8.] $
          *exp(Complex(0,1)*data[where((indgen(54) mod 2) EQ 1)+1,tile_i*8-flagged_num*8.]) ;XX amp and phase
        rts_fitted_cal[0,freq_range,tile_i] = data[where((indgen(54) mod 2) EQ 0)+1,tile_i*8-flagged_num*8.+1] $
          *exp(Complex(0,1)*data[where((indgen(54) mod 2) EQ 1)+1,tile_i*8-flagged_num*8.+1]) ;XX amp and phase fit
        rts_unfitted_cal[1,freq_range,tile_i] = data[where((indgen(54) mod 2) EQ 0)+1,tile_i*8-flagged_num*8.+6] $
          *exp(Complex(0,1)*data[where((indgen(54) mod 2) EQ 1)+1,tile_i*8-flagged_num*8.+6]) ;YY amp and phase
        rts_fitted_cal[1,freq_range,tile_i] = data[where((indgen(54) mod 2) EQ 0)+1,tile_i*8-flagged_num*8.+7] $
          *exp(Complex(0,1)*data[where((indgen(54) mod 2) EQ 1)+1,tile_i*8-flagged_num*8.+7]) ;YY amp and phase fit
      endif else flagged_num = flagged_num + 1
    ;move on to 2
    endfor
    
  endfor
  ;************end rts calibration read in
  
  ;Correct for bp difference between the two bp's
  for pol_i=0,1 do begin
    for tile_i=0,127 do begin
      rts_unfitted_cal[pol_i,*,tile_i] = abs(rts_unfitted_cal[pol_i,*,tile_i]) * (rts_bp_gains_fullband/memo_bp_gains_averaged_fullband)* $
        exp(Complex(0,1)*atan(rts_unfitted_cal[pol_i,*,tile_i],/phase))
      rts_fitted_cal[pol_i,*,tile_i] = abs(rts_fitted_cal[pol_i,*,tile_i]) * (rts_bp_gains_fullband/memo_bp_gains_averaged_fullband)* $
        exp(Complex(0,1)*atan(rts_fitted_cal[pol_i,*,tile_i],/phase))
    endfor
  endfor
  
  ;Average down to FHD freq resolution
  for pol_i=0,1 do begin
    for tile_i=0,127 do begin
      rts_unfitted_cal_80kHz[pol_i,*,tile_i] = (rts_unfitted_cal[pol_i, where((indgen(2*384) mod 2) EQ 0), tile_i] + rts_unfitted_cal[pol_i, where((indgen(2*384) mod 2) EQ 1), tile_i])/2.
      rts_unfitted_cal_80kHz[pol_i,where(((indgen(384)-8) mod 16) EQ 0),tile_i] = rts_unfitted_cal[pol_i, where(((indgen(2*384)-16) mod 32) EQ 0)+1, tile_i]
      rts_fitted_cal_80kHz[pol_i,*,tile_i] = (rts_fitted_cal[pol_i, where((indgen(2*384) mod 2) EQ 0), tile_i] + rts_fitted_cal[pol_i, where((indgen(2*384) mod 2) EQ 1), tile_i])/2.
      rts_fitted_cal_80kHz[pol_i,where(((indgen(384)-8) mod 16) EQ 0),tile_i] = rts_fitted_cal[pol_i, where(((indgen(2*384)-16) mod 32) EQ 0)+1, tile_i]
    endfor
  endfor
  
  ;correcting for the definition of gain (RTS: multiply, FHD: divide) (rec refers to reciprocal)
  rts_unfitted_cal_80kHz_rec = (1./abs(rts_unfitted_cal_80kHz))*exp(-Complex(0,1)*atan(rts_unfitted_cal_80kHz,/phase))
  rts_fitted_cal_80kHz_rec = (1./abs(rts_fitted_cal_80kHz))*exp(-Complex(0,1)*atan(rts_fitted_cal_80kHz,/phase))
  rts_unfitted_cal_80kHz_rec[where(finite(rts_unfitted_cal_80kHz_rec,/nan) EQ 1)] = 0
  rts_fitted_cal_80kHz_rec[where(finite(rts_fitted_cal_80kHz_rec,/nan) EQ 1)] = 0
  rts_unfitted_cal_80kHz_rec[*,255:383,*] = rts_unfitted_cal_80kHz_rec[*,255:383,*]*2.
  rts_fitted_cal_80kHz_rec[*,255:383,*] = rts_fitted_cal_80kHz_rec[*,255:383,*]*2.
  
  ;*******begin unapplied correction based off of Jones entries and spectral slope
  ;RTS does not whiten their data, but FHD does. FHD forces the spectral slope to be flat across the band
  
  ;setup
  spectral_correction = FLTARR(2,384,128)
  linfit_rts = FLTARR(2,384,128)
  jones_first_entry = FLTARR(24)
  
  ;Read in the first entry from the Jones matrix for each course band
  for file_i=0, 23 do begin
    jones_file_path=Reverse(findfile('/nfs/mwa-00/h1/nbarry/RTScal/*Jones*dat'))
    data2=read_ascii(jones_file_path[file_i], delimiter='+')
    jones_first_entry[file_i] = data2.field1[0]
  endfor
  
  ;Fit a line to the Jones entries over the whole band
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316296_obs.sav'
  center_freqs = (*obs.baseline_info).freq[where((findgen(384) mod 16) EQ 8)]
  linfit_coeff = linfit(center_freqs,jones_first_entry)
  
  ;Calculate both correction factors
  for pol_i=0,1 do begin
    for tile_i=0,127 do begin
      ;(Freq of alpha correction / Middle of band) ^ alpha divided by its mean gives the spectral correction across the band without changing the amplitude
      spectral_correction[pol_i,*,tile_i] =(((180E6)/(findgen(384)*80000.+1.67155E8))^(0.74))/mean(((180E6)/(findgen(384)*80000.+1.67155E8))^(0.74))
      
      ;linear fit divided by the mean of the linear fit to get a spectral correction across the band without changing amplitude
      linfit_rts[pol_i,*,tile_i] = (linfit_coeff[1] * (*obs.baseline_info).freq + linfit_coeff[0])/mean(linfit_coeff[1] * (*obs.baseline_info).freq + linfit_coeff[0])
    endfor
  endfor
  
  full_correction = spectral_correction * linfit_rts
  ;*******end unapplied correction
  
  ;unapplied correction, uncomment to apply
  ;rts_unfitted_cal_80kHz_rec2 = abs(rts_unfitted_cal_80kHz_rec)*full_correction*exp(Complex(0,1)*atan(rts_unfitted_cal_80kHz_rec,/phase))
  ;rts_fitted_cal_80kHz_rec2 = abs(rts_fitted_cal_80kHz_rec)*full_correction*exp(Complex(0,1)*atan(rts_fitted_cal_80kHz_rec,/phase))
  
  ;*****cal fhd read-in
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/calibration/1061316296_cal.sav'
  
  raw_fhd = complex(FLTARR(2,384,128))
  raw_fhd[0, *, *] = (*cal.gain[0] + *cal.gain_residual[0])
  raw_fhd[1, *, *] = (*cal.gain[1] + *cal.gain_residual[1])
  cal_fhd = complex(FLTARR(2,384,128))
  cal_fhd[0, *, *] = (*cal.gain[0])
  cal_fhd[1, *, *] = (*cal.gain[1])
  cal_mode = FLTARR(2,128,3)
  for pol_i=0, 1 do for tile_i=0, 127 do if cal.mode_params[pol_i,tile_i] NE !NULL then cal_mode[pol_i,tile_i,*] = (*cal.mode_params[pol_i,tile_i])[*]
  ;*****end cal fhd read-in
  
  ;*****rebuild the cable reflection mode as a function of frequency
  cal_freq_mode = complex(FLTARR(2,128,384))
  for pol_i=0,1 do for tile_i=0, 127 do cal_freq_mode[pol_i,tile_i,*] = cal_mode[pol_i,tile_i,1]*exp(-Complex(0,1)*2.*!Pi*(cal_mode[pol_i,tile_i,0]*findgen(384)/384)+Complex(0,1)*cal_mode[pol_i,tile_i,2])
  ;*****end rebuild calbe reflecion
  
  linfit_unfit_rts = FLTARR(2,384,128)
  linfit_unfit_fhd = FLTARR(2,384,128)
  
  ;***generate a linear fit to both rts and fhd gains to normalize and whiten both ad hoc (since non ad hoc methods are not working well enough)
  for pol_i=0, 1 do begin
    for tile_i=0, 127 do begin
      notmissing_index=where(abs(rts_unfitted_cal_80kHz_rec[pol_i,*,tile_i]) GT 0.0001,n_count)
      if n_count NE 0 then begin
        linfit_rts_coeff = linfit((findgen(384))[notmissing_index],abs(rts_unfitted_cal_80kHz_rec[pol_i,notmissing_index,tile_i]))
        linfit_unfit_rts[pol_i,*,tile_i] = linfit_rts_coeff[1]*findgen(384)+ linfit_rts_coeff[0]
      endif
      
      notmissing_index=where(abs(raw_fhd[pol_i,*,tile_i]) GT 0.0001,n_count)
      if n_count NE 0 then begin
        linfit_fhd_coeff = linfit((findgen(384))[notmissing_index],abs(raw_fhd[pol_i,notmissing_index,tile_i]))
        linfit_unfit_fhd[pol_i,*,tile_i] = linfit_fhd_coeff[1]*findgen(384)+ linfit_fhd_coeff[0]
      endif
    endfor
  endfor
  
  rts_unfitted_cal_80kHz_rec_norm = abs(rts_unfitted_cal_80kHz_rec)/(linfit_unfit_rts)*exp(Complex(0,1)*atan(rts_unfitted_cal_80kHz_rec,/phase))
  rts_fitted_cal_80kHz_rec_norm = abs(rts_fitted_cal_80kHz_rec)/(linfit_unfit_rts)*exp(Complex(0,1)*atan(rts_fitted_cal_80kHz_rec,/phase))
  
  raw_fhd[where(abs(rts_unfitted_cal_80kHz_rec_norm) EQ 0)] = 0
  unfitted_fhd_norm = abs(raw_fhd)/(linfit_unfit_fhd)*exp(Complex(0,1)*atan(raw_fhd,/phase))
  cal_fhd[where(abs(rts_fitted_cal_80kHz_rec_norm) EQ 0)] = 0
  fitted_fhd_norm = abs(cal_fhd)/(linfit_unfit_fhd)*exp(Complex(0,1)*atan(cal_fhd,/phase))
  ;***
  
  ;***Take the fourier transform of the absolute values of the normalized fhd and rts solutions, multiply by a convention factor,
  ;and take the residual between fhd and rts. Do it for raw and gain solutions
  fourier_unfitted_fhd = 384.*fft(abs(raw_fhd)/(linfit_unfit_fhd),dimension=2,/center)
  fourier_unfitted_rts = 384.*fft(abs(rts_unfitted_cal_80kHz_rec)/(linfit_unfit_rts),dimension=2,/center)
  fourier_unfitted_residual = fourier_unfitted_rts - fourier_unfitted_fhd
  
  fourier_fitted_fhd = 384.*fft(abs(cal_fhd)/(linfit_unfit_fhd),dimension=2,/center)
  fourier_fitted_rts = 384.*fft(abs(rts_fitted_cal_80kHz_rec)/(linfit_unfit_rts),dimension=2,/center)
  fourier_fitted_residual = fourier_fitted_rts - fourier_fitted_fhd
  ;***
  
  ;fourier the modefit for plotting purposes
  fourier_modefitted_fhd = 384.*fft(abs(cal_freq_mode+1.),dimension=3,/center)
  ;find the max in the "positive" part of the band
  void = max(fourier_modefitted_fhd[200:384], max_mode_value)
  
  cgPS_Open,'/nfs/mwa-00/h1/nbarry/RTScal/fitted_fft.png',/quiet,/nomatch
  cgplot, abs(fourier_modefitted_fhd[0,2,*]), yrange=[0,2],xrange=[384/2,384], title='Fitted, Tile 3, XX
  cgoplot, abs(fourier_fitted_fhd[0,*,2]), yrange=[0,2],xrange=[384/2,384],title='Fitted, Tile 3, XX',color='red'
  cgoplot, abs(fourier_fitted_rts[0,*,2]), yrange=[0,2],xrange=[384/2,384],title='Fitted, Tile 3, XX',color='blue'
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  cgPS_Open,'/nfs/mwa-00/h1/nbarry/RTScal/unfitted_fft.png',/quiet,/nomatch
  cgplot, abs(fourier_modefitted_fhd[0,2,*]), yrange=[0,2],xrange=[384/2,384], title='Unfitted, Tile 3, XX
  cgoplot, abs(fourier_unfitted_fhd[0,*,2]), yrange=[0,2],xrange=[384/2,384],title='Fitted, Tile 3, XX',color='red'
  cgoplot, abs(fourier_unfitted_rts[0,*,2]), yrange=[0,2],xrange=[384/2,384],title='Fitted, Tile 3, XX',color='blue'
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  cgPS_Open,'/nfs/mwa-00/h1/nbarry/RTScal/residual_unfitted_fft.png',/quiet,/nomatch
  cgplot, abs(fourier_modefitted_fhd[0,2,*]), yrange=[0,2],xrange=[384/2,384], title='Unfitted, Tile 3, XX
  cgoplot, abs(fourier_unfitted_residual[0,*,2]), yrange=[0,2],xrange=[384/2,384],title='Fitted, Tile 3, XX',color='green'
  cgoplot, abs(fourier_fitted_residual[0,*,2]), yrange=[0,2],xrange=[384/2,384],title='Fitted, Tile 3, XX',color='purple'

  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  
  ;***Take the residual between fhd and rts, both raw and gain solutions
  unfitted_residual = abs(rts_unfitted_cal_80kHz_rec_norm)-abs(unfitted_fhd_norm)
  unfitted_residual[where(abs(rts_unfitted_cal_80kHz_rec_norm) EQ 0)] = 0
  
  fitted_residual = abs(rts_fitted_cal_80kHz_rec_norm)-abs(fitted_fhd_norm)
  fitted_residual[where(abs(rts_fitted_cal_80kHz_rec_norm) EQ 0)] = 0
  ;***
  stop
  ;ranges = where((findgen(384) mod 16) EQ 0)
  
  
  
  
  ;*Plotting included
  mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  
  outdir = '/nfs/mwa-00/h1/nbarry/RTScal/plots/'
  freq_arr = (*obs.baseline_info).freq
  for pol_i=0,1 do begin
    for tile_i=0,127 do begin
    
      ;Naming convention
      if pol_i EQ 0 then pol_name='xx'
      if pol_i EQ 1 then pol_name='yy'
      
      cgPS_Open,outdir+'comparisons_'+pol_name+'_'+strtrim(string(tile_i),2)+'.png',/quiet,/nomatch
      
      cgplot, freq_arr/1E6,abs(rts_unfitted_cal_80kHz_rec[pol_i,*,tile_i]),xrange=[(freq_arr[0])/1E6,(freq_arr[383])/1E6], yrange=[.5,2], $
        title='Tile ' + strtrim(string(tile_i),2) + $
        ' ('+strtrim(string(UINT(cable_len[tile_i])),2) +'m cable), '+pol_name, XTICKFORMAT="(A1)", ytitle='Gain',charsize=1,color='black', $
        Position=[0.10, 0.35, 0.9, 0.90]
      titletotal='Unfit gains, RTS'
      colortotal=60B
      
      cgoplot, freq_arr/1E6,abs(raw_fhd[pol_i,*,tile_i]), color='blue'
      titletotal=[titletotal,'Unfit gains, FHD']
      
      cgLegend, Title=titletotal, $
        Color=['black','blue'],Length=.03,charsize=.7,$;Psym=[2,2,2,2,2,2],Length=0.0 $
        Location=[0.56,0.87]
        
      bottom_included=1
      if keyword_set(bottom_included) then begin
      
        cgPlot,(freq_arr)/1E6, unfitted_residual[pol_i,*, tile_i], $
          Position=[0.1, 0.1, 0.9, 0.34], color='black',xrange=[freq_arr[0]/1E6,freq_arr[383]/1E6],$
          XTICKFORMAT="(A1)",/NoErase, charsize=.75, yrange=[-.25,.25],ytitle='Residual',YTICKFORMAT="(F0.2)",ytickinterval=.05
      ;for obs_i=1,N_elements(*obs_ptr[j])-1 do cgoplot,real_part( abs(normal_sols[freq_use,tile_i,obs_i,pol_i]) -  $
      ;  abs(normal_input[freq_use,tile_i,obs_i,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,obs_i,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,obs_i,pol_i]) ),color=10B
      ;If keyword_set(titleB) then cgoplot,(freq_arr)/1E6,real_part( abs(B_sols[freq_use,tile_i,0,pol_i]) - mean_B[freq_use,tile_i,pol_i] ),color=20B
      ;If keyword_set(titleC) then cgoplot,(freq_arr)/1E6,real_part( abs(C_sols[freq_use,tile_i,0,pol_i]) - mean_C[freq_use,tile_i,pol_i] ),color=30B
      ;If keyword_set(titleD) then cgoplot,(freq_arr)/1E6,real_part( abs(D_sols[freq_use,tile_i,0,pol_i]) - mean_D[freq_use,tile_i,pol_i] ),color=40B
      ;If keyword_set(titleE) then cgoplot,(freq_arr)/1E6,real_part( abs(E_sols[freq_use,tile_i,0,pol_i]) - mean_E[freq_use,tile_i,pol_i] ),color=50B
          
      ; unphased_im=abs( abs(normal_sols[freq_use,tile_i,*,pol_i]) - $
      ;   abs(normal_input[freq_use,tile_i,*,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,*,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,*,pol_i]))
          
      ;cgPlot, (freq_arr)/1E6,mean(reform(unphased_im),dimension=2), $
      ;  Position=[0.1, 0.15, 0.9, 0.24], color=10B,xrange=[freq_arr[0]/1E6,freq_arr[335]/1E6],$
      ;  ytickinterval=.05,xtitle='Unflagged frequency MHz',/NoErase, charsize=.75, yrange=[0,.1],ytitle='Magnitude!cunphased';,YTICKFORMAT="(A1)"
          
      ;for obs_i=1,N_elements(*obs_ptr[j])-1 do cgoplot, imaginary( abs(normal_sols[freq_use,tile_i,obs_i,pol_i]) - $
      ;  abs(normal_input[freq_use,tile_i,obs_i,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,obs_i,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,obs_i,pol_i]) ),color=10B
          
      ;If keyword_set(titleB) then cgoplot,(freq_arr)/1E6,abs( abs(B_sols[freq_use,tile_i,0,pol_i]) - mean_B[freq_use,tile_i,pol_i] ),color=20B
      ;If keyword_set(titleC) then cgoplot,(freq_arr)/1E6,abs( abs(C_sols[freq_use,tile_i,0,pol_i]) - mean_C[freq_use,tile_i,pol_i] ),color=30B
      ;If keyword_set(titleD) then cgoplot,(freq_arr)/1E6,abs( abs(D_sols[freq_use,tile_i,0,pol_i]) - mean_D[freq_use,tile_i,pol_i] ),color=40B
      ;If keyword_set(titleE) then cgoplot,(freq_arr)/1E6,abs( abs(E_sols[freq_use,tile_i,0,pol_i]) - mean_E[freq_use,tile_i,pol_i] ),color=50B
          
      endif
      
      ;Device, Decomposed=1
      
      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
    endfor
  endfor
  
  
  
  
end