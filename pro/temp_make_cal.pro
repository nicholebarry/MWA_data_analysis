pro temp_make_cal

  cal = getvar_savefile('/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_Aug2017_zeroedbp/calibration/1061316296_cal.sav','cal')
  cal2 = getvar_savefile('/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_Aug2017_globalbp_w_cable_w_digjump/calibration/1061316296_cal.sav','cal')
  
  for tile_i=0, 127 do begin
    cgps_open, '/nfs/eor-00/h1/nbarry/thesis_images/zerocal/'+strtrim(tile_i,2)+'.png',/quiet,/nomatch
    cgplot, abs((*cal.gain[0])[*,tile_i]), title = strtrim(tile_i,2)
    cgoplot, abs((*cal2.gain[0])[*,tile_i]), color='blue'
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  endfor
  
end