pro beam_findzeros

  n_freq=48
  n_baselines=1
  n_pol=2
  
  n_uv_res=101
  
  psf = getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013cal_largerdim20/beams/1061319472_beams.sav','psf')
  ;psf = getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013cal_largerdim44/beams/1061311664_beams.sav','psf')
  ;psf = getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/beams/1061316296_beams.sav','psf')
  
  for pol_i=1, n_pol-1 do begin
    base_array = complex(FLTARR(psf.dim*psf.dim))
    for base_i=0, n_baselines-1 do begin
      ;for freq_i=0, n_freq-1 do begin
      for freq_i=n_freq-1, n_freq-1 do begin
        for uv_i=0, n_uv_res-1 do begin
          for uv_j=0, n_uv_res-1 do begin
            array = (*(*(*psf.beam_ptr)[pol_i,freq_i,base_i])[uv_i,uv_j])
            base_array = max([[abs(base_array)],[abs(array)]],dim=2)
          endfor
        endfor
      endfor
    endfor
    psf_single=Reform(base_array,psf.dim,psf.dim)
    stop
  endfor
  
  
end