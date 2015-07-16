pro fun_run

  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_poly_saved_run/calibration/1061311664_cal.sav'
  n_freq=384
  pol_i=0
  degree=1
  phase_degree=1
  
  FOR tile_i=0,99 DO BEGIN
    gain_fit=DBLARR(n_freq)
    phase_fit=DBLARR(n_freq)
    mode_params=DBLARR(3)
    amp_params=*(cal.amp_params[pol_i,tile_i])
    phase_params=*(cal.phase_params[pol_i,tile_i])
    FOR di=0L,degree DO gain_fit+=amp_params[di]*findgen(n_freq)^di
    FOR di=0L,phase_degree DO phase_fit+=phase_params[di]*findgen(n_freq)^di
    If (cal.mode_params[pol_i,tile_i]) NE !NULL THEN mode_params=*(cal.mode_params[pol_i,tile_i])
    
    gain=(gain_fit*Exp(Complex(0,1)*phase_fit)+mode_params[1]*Exp(-Complex(0,1)*2.*!Pi*mode_params[0]*findgen(n_freq)/n_freq+Complex(0,1)*(mode_params[2])))
    if tile_i EQ 0 then cgplot,abs(gain),yrange=[1,2] else cgoplot,abs(gain)
    
  ENDFOR
  
end