Pro cal_one_num_per_antenna

  restore, '/nfs/mwa-03/r1/EoR2013//fhd_nb_devel_Aug2015/calibration/1061316296_cal.sav'
  abs_mean = FLTARR(2,128)
  phase_mean = FLTARR(2,128)
  slope_gain = FLTARR(2,128)
  slope_phase = FLTARR(2,128)
  tot_mean = FLTARR(2,128)
  tot_mean_expand = FLTARR(8,128)
  
  restore, '/nfs/mwa-03/r1/EoR2013//fhd_nb_devel_Aug2015/metadata/1061316296_obs.sav'
  freq_arr=(*obs.baseline_info).freq
  
  
  
  for pol_i=0, 1 do begin
    for tile_i=0, 127 do begin
      non_zero = where(abs((*cal.gain[pol_i])[*,tile_i]) NE 0 ,n_count)
      if n_count GT 0 then begin
        abs_mean[pol_i,tile_i] = mean(abs((*cal.gain[pol_i])[non_zero,tile_i]))
        phase_mean[pol_i,tile_i] = mean(phunwrap(atan((*cal.gain[pol_i])[non_zero,tile_i],/phase)))
        
        slope_gain[pol_i,tile_i] = (mean(abs((*cal.gain[pol_i])[1:14,tile_i]))-mean(abs((*cal.gain[pol_i])[369:382,tile_i])))/(7.-375.)
      slope_phase[pol_i,tile_i] = (*cal.phase_params[pol_i,tile_i])[1]
      endif
      
    endfor
  endfor
  column_preheader = 'Real and imaginary components of the calibration mean per tile. Slope calculated as a function of channel number (0:383).' 
  column_headers = 'XX Real, XX Imaginary, XX Gain Slope, XX Phase Slope, YY Real, YY Imaginary, YY Gain Slope, YY Phase Slope'
  tot_mean = abs_mean*exp(Complex(0,1)*phase_mean)
  tot_mean_expand[0,*] = real_part(tot_mean[0,*])
  tot_mean_expand[1,*] = imaginary(tot_mean[0,*])
  tot_mean_expand[2,*] = slope_gain[0,*]
  tot_mean_expand[3,*] = slope_phase[0,*]
  tot_mean_expand[4,*] = real_part(tot_mean[1,*])
  tot_mean_expand[5,*] = imaginary(tot_mean[1,*])
  tot_mean_expand[6,*] = slope_gain[1,*]
  tot_mean_expand[7,*] = slope_phase[1,*]
  stop
  
end