pro zeroed_cal

  ;lower_band_only=1
  ;double_mode=1
  find_modes=1
  hyperfine_dft=1
  
  n_tiles=128
  n_pol=2
  n_poi=5
  
  ;Using preexisting file to extract information about which tiles have which cable length
  mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  ;Taking tile information and cross-matching it with the nonflagged tiles array, resulting in nonflagged tile arrays
  ;grouped by cable length
  cable_length_ref=cable_len[Uniq(cable_len,Sort(cable_len))]
  n_cable=N_Elements(cable_length_ref)
  tile_use_arr=Ptrarr(n_cable)
  FOR cable_i=0,n_cable-1 DO tile_use_arr[cable_i]=Ptr_new(where((cable_len EQ cable_length_ref[cable_i]) EQ 1))
  
  ;Calculate the shifting amount to center the fft
  if keyword_set(lower_band_only) then n_freq=256 else n_freq=384.
  if (n_freq mod 2) EQ 0 then shift_num = n_freq/2-1 else shift_num = floor(n_freq/2.)
  
  ;Calculate the fft axis in ns
  T= (80000.)
  X = FINDGEN((n_freq - 1)/2) + 1
  x_axis = shift([0.0, X, n_freq/2, -n_freq/2 + X]/(n_freq*T)*1E9,shift_num)
  
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/norm_gain_plus_phase_dig_poi.sav'
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/longrun_gain_plus_phase_dig_poi.sav'
  
  denom = (80000.*16.)^2. + (80000.)^2.
  
  coarse_1 = PTRARR(8,/allocate)
  coarse_2 = PTRARR(8,/allocate)
  
  for harmonic_i=1,7 do begin
    high_line = ( ((T*16.*float(harmonic_i))/denom)-(T/denom) )*1E9
    low_line = ( ((T*16.*float(harmonic_i))/denom)+(T/denom) )*1E9
    if harmonic_i EQ 1 then high_line_all=high_line else high_line_all=[high_line_all,high_line]
    if harmonic_i EQ 1 then low_line_all=low_line else low_line_all=[low_line_all,low_line]
  ;cgoplot, [high_line, high_line],[-1,1],  linestyle=2
  ;cgoplot, [low_line, low_line],[-1,1],  linestyle=2
  ;coarse_band_wh = where((abs(x_axis) LT low_line) AND (abs(x_axis) GT high_line),nc)
  ;if harmonic_i EQ 1 then coarse_bands_wh = coarse_band_wh else coarse_bands_wh = [coarse_bands_wh,coarse_band_wh]
  ;*coarse[harmonic_i] = where((abs(x_axis) GT low_line) AND (abs(x_axis) GT high_line))
  endfor
  
  for harmonic_i=0,6 do begin
    if harmonic_i EQ 6 then $
      inds = where((abs(x_axis) GT low_line_all[harmonic_i])) $
    else $
      inds = where((abs(x_axis) GT low_line_all[harmonic_i]) AND (abs(x_axis) LT high_line_all[harmonic_i+1]))
    *coarse_1[harmonic_i] = inds[where( inds LE shift_num)]
    *coarse_2[harmonic_i] = inds[where( inds GT shift_num)]
    
  endfor
  
  obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/1061316296_obs.sav','obs')
  if keyword_set(lower_band_only) then begin
    freq_use_logic = (*obs.baseline_info).freq_use
    freq_use_logic=freq_use_logic[0:255]
  endif else freq_use_logic = (*obs.baseline_info).freq_use
  freq_use=where(freq_use_logic,nf_use)
  freq_not_use=where(freq_use_logic EQ 0)
  
  if keyword_set(lower_band_only) then longrun_gain = longrun_gain[*,*,0:255,*]
  longrun_gain[*,*,freq_not_use,*]=0
  fft_longrun = shift(fft(abs(longrun_gain),dimension=3),0,0,shift_num,0)
  
  ;find where the coarse band lines are using simple data
  flat = FLTARR(n_freq)
  flat[*]=1.
  flat[freq_not_use]=0.
  fft_flat = shift(fft(flat),shift_num)
  wh_fft_flat = where(abs(fft_flat) GT .0001,n_count)
  ;wh_fft_flat_ar = INTARR(n_pol,n_poi,N_elements(wh_fft_flat),n_tiles)
  ;wh_fft_flat_ar[0,*,*,*] = rebin(wh_fft_flat,n_poi,N_elements(wh_fft_flat),n_tiles)
  ;wh_fft_flat_ar[1,*,*,*] = rebin(wh_fft_flat,n_poi,N_elements(wh_fft_flat),n_tiles)
  fft_flat_ar=complex(FLTARR(n_pol,n_poi,n_freq,n_tiles))
  for pol_i=0, n_pol-1 do for poi_i=0, n_poi-1 do for tile_i=0, n_tiles-1 do fft_flat_ar[pol_i,poi_i,*,tile_i]=fft_flat
  
  ;n_coarse_bands = N_elements(wh_fft_flat)
  
  
  
  if ~keyword_set(double_mode) then n_mode = 3 else n_mode=5
  mode_index = INTARR(n_pol,n_poi,n_mode,n_tiles) ;3 modes for each cable reflection: neg, center, and pos ;5 modes for each cable reflection: neg, center, pos, next to neg, next to pos
  
  cable_len_to_fit = [90.,150.,230.,320.,400.,524.]
  c_light=299792458.
  mode_i_arr=Fltarr(n_pol,n_tiles)
  ; Set up to fit for the reflection amplitude and phase
  ; Get the nominal tile lengths and velocity factors:
  cable_filepath=filepath(obs.instrument+'_cable_length.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=cable_filepath,first_line=1
  tile_i_file=Reform(data_array[0,*])
  tile_name_file=Reform(data_array[1,*])
  cable_len=Reform(data_array[2,*])
  cable_vf=Reform(data_array[3,*])
  
  
  for cable_i=0,3 do begin
  
    ;***Calculate the mode
    tile_ref_flag=0>Reform(data_array[4,*])<1
    ; Choose which cable lengths to fit
    cable_cut_i=where(cable_len NE cable_len_to_fit[cable_i],n_cable_cut)
    IF n_cable_cut GT 0 THEN tile_ref_flag[cable_cut_i]=0
    
    reflect_time=2.*cable_len/(c_light*cable_vf) ; nominal reflect time
    bandwidth=(Max((*obs.baseline_info).freq)-Min((*obs.baseline_info).freq))*n_freq/(n_freq-1)
    
    FOR pol_i=0,n_pol-1 DO mode_i_arr[pol_i,where(tile_ref_flag)]=bandwidth*reflect_time[where(tile_ref_flag)]
    ;***
    
    if cable_i EQ 0 then freq_channels = (*coarse_1[0])+round(N_elements(*coarse_1[0])/2) $
    else freq_channels = *coarse_1[cable_i-1]
    
    temp = max( (abs(fft_longrun - fft_flat_ar))[*,*,freq_channels,*tile_use_arr[cable_i]] , max_subscript,dim=3)
    
    ;find the index of the 150m cable reflections (pos and neg in fft space)
    s = SIZE(abs(fft_longrun[*,*,freq_channels,*tile_use_arr[cable_i]]))
    ncol = s[1]
    nrow = s[2]
    nz = s[3]
    z = (max_subscript / ncol / nrow) mod nz
    
    mode_index[*,*,0,*tile_use_arr[cable_i]] = freq_channels[z]
    mode_index[*,*,1,*tile_use_arr[cable_i]] = shift_num
    mode_index[*,*,2,*tile_use_arr[cable_i]] = shift_num+(shift_num-freq_channels[z])
    
    if keyword_set(double_mode) then begin
      temp2 = (fft_longrun - fft_flat_ar)[*,*,freq_channels,*tile_use_arr[cable_i]]
      for pol_i=0, n_pol-1 do for poi_i=0, n_poi-1 do for tile_i=0, N_elements(*tile_use_arr[cable_i])-1 do $
        temp2[pol_i,poi_i,z[pol_i,poi_i,tile_i],tile_i] -= reform(fft_longrun[pol_i,poi_i,freq_channels[z[pol_i,poi_i,tile_i]],(*tile_use_arr[cable_i])[tile_i]])
      temp2 = max(abs(temp2),max_subscript,dim=3)
      
      s = SIZE(abs(fft_longrun[*,*,freq_channels,*tile_use_arr[cable_i]]))
      ncol = s[1]
      nrow = s[2]
      nz = s[3]
      z = (max_subscript / ncol / nrow) mod nz
      
      mode_index[*,*,3,*tile_use_arr[cable_i]] = freq_channels[z]
      mode_index[*,*,4,*tile_use_arr[cable_i]] = shift_num+(shift_num-freq_channels[z])
    endif
    
  endfor
  
  ones = FLTARR(n_freq)
  ones[*] = 1.
  gain_mode_fit = complex(FLTARR(n_pol,n_poi,n_freq,n_tiles))
  
  if keyword_set(hyperfine_dft) then begin
    for pol_i=0, n_pol-1 do begin
      for tile_i=0, n_tiles-1 do begin
        mode_i=mode_i_arr[pol_i,tile_i]
        IF mode_i EQ 0 THEN CONTINUE
        ; We are going to fit the actual mode to subtract.
        mode0=mode_i ; start with nominal cable length
        dmode=0.05 ; overresolve the FT used for the fit (normal resolution would be dmode=1)
        nmodes=101 ; range around the central mode to test
        modes=(dindgen(nmodes)-nmodes/2)*dmode+mode0 ; array of modes to try
        modes=rebin(modes,nmodes,nf_use) ; reshape for ease of computing
        
        for poi_i=0, n_poi-1 do begin
          gain_temp = rebin_complex(transpose(reform(abs(longrun_gain[pol_i,poi_i,freq_use,tile_i]))),nmodes,nf_use)
          flat_temp = rebin_complex(transpose(reform(flat[freq_use])),nmodes,nf_use)
          freq_mat=rebin(transpose(freq_use),nmodes,nf_use) ; freq_use matrix to multiply/collapse in fit
          test_fits=Total(exp(Complex(0,1)*2.*!Pi/n_freq*modes*freq_mat)*gain_temp,2) ; Perform DFT of gains to test modes
          flat_dft=Total(exp(Complex(0,1)*2.*!Pi/n_freq*modes*freq_mat)*flat_temp,2) ; Perform DFT of gains to test modes
          amp_use=max(abs(test_fits-flat_dft),mode_ind)/nf_use*2. ; Pick out highest amplitude fit (mode_ind gives the index of the mode)
          phase_use=atan((test_fits-flat_dft)[mode_ind],/phase) ; Phase of said fit
          mode_i=modes[mode_ind,0] ; And the actualy mode
          
          gain_mode_fit[pol_i,poi_i,*,tile_i]=amp_use*exp(-Complex(0,1)*2.*!Pi*(mode_i*findgen(n_freq)/n_freq)+Complex(0,1)*phase_use)+ones
        ;gain_arr_fit[*,tile_i]+=gain_mode_fit
        ;cal_return.mode_params[pol_i,tile_i]=Ptr_new([mode_i,amp_use,phase_use])
        endfor
      endfor
    endfor
  endif
  
  
  con_signal=complex(FLTARR(n_pol,n_poi,n_freq,n_tiles))
  
  for cable_i=0,3 do begin
    for pol_i=0,n_pol-1 do begin
      for poi_i=0,n_poi-1 do begin
        for tile_i=0, N_elements(*tile_use_arr[cable_i])-1 do begin
        
          tile_num = (*tile_use_arr[cable_i])[tile_i]
          
          ;mode1 = mode_index[pol_i,poi_i,tile_num]
          
          ;wh_fft_mode = [mode1,shift_num,shift_num+(shift_num-mode1)]
          con_loc_total = INTARR(FLOAT(n_mode) * N_elements(wh_fft_flat))
          for mode_i=0,n_mode-1 do begin
          
          
            con_locs = (mode_index[pol_i,poi_i,mode_i,tile_num] + wh_fft_flat ) - (shift_num)
            
            con_signal[pol_i,poi_i,mode_index[pol_i,poi_i,mode_i,tile_num],tile_num] = fft_longrun[pol_i,poi_i,mode_index[pol_i,poi_i,mode_i,tile_num],tile_num]
            
            con_vals = reform(fft_longrun[pol_i,poi_i,mode_index[pol_i,poi_i,mode_i,tile_num],tile_num] * fft_longrun[pol_i,poi_i,wh_fft_flat,tile_num])
            wh_flip = where((con_locs) GE n_freq, n_count)
            if n_count GT 0 then con_locs[wh_flip] =  (abs(con_locs[wh_flip])-n_freq)
            con_signal[pol_i,poi_i,abs(con_locs),tile_num] += (con_vals)
            
            con_loc_total[mode_i * N_elements(wh_fft_flat): (mode_i+1)*N_elements(wh_fft_flat)-1] = con_locs
          endfor
          wh_con = where(reform(abs(con_signal[pol_i,poi_i,abs(con_loc_total),tile_num]))/reform(abs(fft_longrun[pol_i,poi_i,abs(con_loc_total),tile_num])) GT 1., n_count)
          if n_count GT 0 then con_signal[pol_i,poi_i,abs(con_loc_total[wh_con]),tile_num]=fft_longrun[pol_i,poi_i,abs(con_loc_total[wh_con]),tile_num]
        endfor
      endfor
    endfor
  endfor
  
  if keyword_set(hyperfine_dft) then begin
    gain_mode_fit[*,*,freq_not_use,*]=0.
    con_signal = shift(abs(fft(abs(gain_mode_fit),dim=3)),0,0,shift_num,0)
  endif
  
  zero_longrun = fft_longrun
  for harmonic_i=6, 0, -1 do begin
    ;zero_longrun[*,*,*coarse_1[harmonic_i],*] = 0.
    zero_longrun[*,*,(*coarse_1[harmonic_i]),*] = con_signal[*,*,(*coarse_1[harmonic_i]),*]
    zero_longrun[*,*,(*coarse_2[harmonic_i]),*] = con_signal[*,*,(*coarse_2[harmonic_i]),*]
    
  ;zero_longrun[*,*,*coarse_2[harmonic_i],*] = 0.
  ;zero_longrun[*,*,(*coarse_2[harmonic_i])[8],where(*tile_use_arr[1] EQ 1)] = fft_longrun[*,*,(*coarse_2[harmonic_i])[8],where(*tile_use_arr[1] EQ 1)]
  endfor
  fft_zero_longrun = abs(fft(shift(zero_longrun,0,0,-shift_num,0),/inverse,dim=3))
  stop
  
  quick_image, reform(abs(fft_longrun[0,0,0:384/2,where(*tile_use_arr[5] EQ 1)])), x_axis[0:384/2], data_range=[0,.001], xrange=[0,6000], $
    xtitle='$\tau$-space', ytitle='524m tiles', title='XX, -2 pointing',/png, savefile='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/xx_524_-2'
  quick_image, reform(abs(fft_longrun[0,1,0:384/2,where(*tile_use_arr[5] EQ 1)])), x_axis[0:384/2], data_range=[0,.001], xrange=[0,6000], $
    xtitle='$\tau$-space', ytitle='524m tiles', title='XX, -1 pointing',/png, savefile='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/xx_524_-1'
  quick_image, reform(abs(fft_longrun[0,2,0:384/2,where(*tile_use_arr[5] EQ 1)])), x_axis[0:384/2], data_range=[0,.001], xrange=[0,6000], $
    xtitle='$\tau$-space', ytitle='524m tiles', title='XX, 0 pointing',/png, savefile='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/xx_524_0'
  quick_image, reform(abs(fft_longrun[0,3,0:384/2,where(*tile_use_arr[5] EQ 1)])), x_axis[0:384/2], data_range=[0,.001], xrange=[0,6000], $
    xtitle='$\tau$-space', ytitle='524m tiles', title='XX, 1 pointing',/png, savefile='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/xx_524_1'
  quick_image, reform(abs(fft_longrun[0,4,0:384/2,where(*tile_use_arr[5] EQ 1)])), x_axis[0:384/2], data_range=[0,.001], xrange=[0,6000], $
    xtitle='$\tau$-space', ytitle='524m tiles', title='XX, 2 pointing',/png, savefile='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/xx_524_2'
    
end