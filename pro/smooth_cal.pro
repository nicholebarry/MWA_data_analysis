pro smooth_cal

  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod_zenithpointing_notileflag/calibration/1061316176_cal.sav'
  n_pol=2
  n_tile=128
  n_freq=384
  
  perfect_phase=1
  
  FOR pol_i=0,n_pol-1 DO BEGIN
    FOR tile_i=0L,n_tile-1 DO BEGIN
      fit_params=*cal.amp_params[pol_i,tile_i]
      gain_fit=fltarr(n_freq)
      FOR di=0L,2 DO gain_fit+=fit_params[di]*findgen(n_freq)^di
      phase_params=*cal.phase_params[pol_i,tile_i]
      phase_fit=fltarr(n_freq)
      IF ~keyword_set(perfect_phase) then begin
        FOR di=0L,1 DO phase_fit+=phase_params[di]*findgen(n_freq)^di
      endif else phase_fit[*]=0
      (*cal.gain[pol_i])[*,tile_i]=gain_fit*Exp(Complex(0,1)*phase_fit)
    endfor
  endfor
  stop
end