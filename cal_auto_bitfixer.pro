pro cal_auto_bitfixer



  cal_auto=getvar_savefile('/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_Aug2017_autocal1/calibration/1061316296_cal.sav','cal')
  
  ;obs=getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_Aug2017_std/metadata/1061316296_obs.sav','obs')
  
  center_freq = where((INDGEN(384) mod 16) EQ 8)
  end_freq = where((INDGEN(384) mod 16) EQ 14)
  beg_freq = where((INDGEN(384) mod 16) EQ 1)
  n_tile=128
  
  for coarse_i=0,N_elements(center_freq)-1 do $
    (*cal_auto.gain[0])[beg_freq[coarse_i]:end_freq[coarse_i],*] *= $
    ;rebin_complex(weight_invert((*cal_auto.gain[0])[center_freq[coarse_i],*] * 0.14), 14,n_tile)
    rebin_complex(weight_invert(0.14), 14,n_tile)
  
  
  
end