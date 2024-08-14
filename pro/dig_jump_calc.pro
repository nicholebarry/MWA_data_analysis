pro dig_jump_calc

  obs = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/metadata/1061316296_obs.sav','obs')
  freq_arr = (*obs.baseline_info).freq
  
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/longrun_gain_ave/norm_gain_plus_phase_dig_poi.sav' ;2,995,384,128
  freq_not_use=where((*obs.baseline_info).freq_use EQ 0,nf_use)
  gain[*,*,freq_not_use,*]=0
  fft_gain = 2*abs(fft(abs(gain),dim=3))
  
  N=384.
  T= (80000.)
  X = FINDGEN((N - 1)/2) + 1
  x_axis = [0.0, X, N/2, -N/2 + X]/(N*T)
  k_axis = x_axis[sort(x_axis)]
  
  pulse_width = freq_arr[255] - freq_arr[0]
  period = freq_arr[383] - freq_arr[0]
  
  fourier_result = complex(FLTARR(192))

  ;fourier_result[0]=2*pulse_width/period
  ;Just the positive portion of the band
  for k=1, 191 do $
    fourier_result[k] = sqrt((.01/(k*!Pi)*sin(2*!Pi*k*pulse_width/period))^2. + (.01/(k*!Pi)*(1-cos(2*!Pi*k*pulse_width/period)))^2.)
    ;fourier_result[k] =fft_gain[0,0,k,0]*k*!Pi/sqrt(sin(2*!Pi*k*pulse_width/period)^2. + (1-cos(2*!Pi*k*pulse_width/period))^2.)

  stop
  
   mode0=384./2. ; start with nominal cable length
   dmode=1 ; overresolve the FT used for the fit (normal resolution would be dmode=1)
   nmodes=384.*1. ; range around the central mode to test
   modes=(dindgen(nmodes)-nmodes/2)*dmode+mode0 ; array of modes to try
   modes=rebin(modes,nmodes,384) ; reshape for ease of computing
        ; These lines are silly. Need to rebin gain_arr, but can't do complex numbers directly.
   gainr=rebin(transpose(reform(abs(longrun_gain[0,0,*,3]))),nmodes,384) ; dimension manipulation, add dim for mode fitting
   ;gaini=rebin(transpose(reform(imaginary(longrun_gain[0,0,*,3]))),nmodes,384) ; dimension manipulation, add dim for mode fitting
   gain_temp = gainr;+Complex(0,1)*gaini
     freq_mat=rebin(transpose(Findgen(384)),nmodes,384) ; this too...
   test_fits=Total(exp(-Complex(0,1)*2.*!Pi/384*modes*freq_mat)*gain_temp,2)/384. ; Perform DFT of gains to test modes
              amp_use=max(abs(test_fits),mode_ind)/nf_use ; Pick out highest amplitude fit (mode_ind gives the index of the mode)
              phase_use=atan(test_fits[mode_ind],/phase) ; Phase of said fit
              mode_i=modes[mode_ind,0] ; And the actualy mode
  
  stop
  
  ;mode = .2*sin(2.*!Pi*(0.78125E-6)*freq_arr)+1.
  mode = .2*sin(2.*!Pi*(1.237E-6)*freq_arr)+1.
  fft_mode = shift(fft(mode),384/2)
  wh_fft_mode = where(abs(fft_mode) GT .0001,n_count)
  flat = FLTARR(384)
  flat[*]=1.
  freq_not_use=where((*obs.baseline_info).freq_use EQ 0,nf_use)
  flat[freq_not_use]=0.
  fft_flat = shift(fft(flat),384/2)
  wh_fft_flat = where(abs(fft_flat) GT .0001,n_count)
  
  mixed_signal = shift(fft(flat*mode),384/2)
  
  ;con_locs0 = wh_fft_mode[0] - wh_fft_flat
  ;con_vals0 = fft_mode[wh_fft_mode[0]] * fft_flat[wh_fft_flat]
  con_locs1 = ((wh_fft_mode[1] + wh_fft_flat))
  con_vals1 = fft_mode[wh_fft_mode[1]] * fft_flat[wh_fft_flat]
  con_locs2 = ((wh_fft_mode[2] + wh_fft_flat))
  con_vals2 = fft_mode[wh_fft_mode[2]] * fft_flat[wh_fft_flat]
  con_signal = FLTARR(384)
  ;con_signal[con_locs0] = con_signal[con_locs0] + con_vals0
  con_locs = [con_locs1 - 384/2, con_locs1 + 384/2];, con_locs2 - 384/2, con_locs2 + 384/2]
  ;con_locs = [con_locs1 , con_locs2 ]
  wh_flip = where((con_locs) GT 384, n_count)
  if n_count GT 0 then con_locs[wh_flip] =  (abs(con_locs[wh_flip])-384)
  ;if n_count GT 0 then con_locs[wh_flip] =  384-abs(con_locs[wh_flip]);-384)
  ;con_vals = [con_vals1, con_vals1];, con_vals2, con_vals2]
  con_vals = [con_vals1, con_vals2]
  con_signal[abs(con_locs)] = con_signal[abs(con_locs)] + abs(con_vals)
  ;con_signal[abs(con_locs1)] = con_signal[abs(con_locs1)] + abs(con_vals1)
;con_signal[abs(con_locs2)] = con_signal[abs(con_locs2)] + abs(con_vals2)
  
end