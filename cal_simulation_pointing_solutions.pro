pro cal_simulation_pointing_solutions


  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23zenith.txt'
  ;filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/obs_id_6176.txt'
  
  readcol, filename, obs_id, format='A', /silent
  
 ;   filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23plusone.txt'
 ; readcol, filename, obs_id2, format='A', /silent
 ; obs_id=[obs_id,obs_id2]
  
  parsednumbers=N_elements(obs_id)
  
  cal_sols=complex(FLTARR(384,128,2,parsednumbers))
  
  for obs_i=0, parsednumbers-1 do begin
    ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_maxcalsources_nod_zenithpointing_notileflag/calibration/'+ obs_id[obs_i] +'_cal.sav'
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_overfit_cal_eor_ones_short_baselines_included_zenithpointing_2000/calibration/'+ obs_id[obs_i] +'_cal.sav'
    
    ;if obs_i EQ 0 then cgplot, atan((*cal.gain[0])[*,3],/phase),yrange=[-.01,.01] else cgoplot, atan((*cal.gain[0])[*,3],/phase),yrange=[-.01,.01]
    for pol_i=0,1 do begin
      cal_sols[*,*,pol_i,obs_i]=*cal.gain[pol_i]
    endfor
  endfor
  ;for pol_i=0,1 do *cal.gain[pol_i]=*cal.gain[pol_i]+*cal.gain_residual[pol_i]
  
  cal_final=complex(FLTARR(384,2))
  
  for pol_i=0,1 do begin
    for freq_i=0, 383 do begin
      for obs_i=0,parsednumbers-1 do begin
      
        if obs_i EQ 0 then cal_pertile_perpol=cal_sols[freq_i,*,pol_i,obs_i] else cal_pertile_perpol=[cal_pertile_perpol,cal_sols[freq_i,*,pol_i,obs_i]]
        
      endfor
      resistant_mean, real_part(cal_pertile_perpol), 2, res_mean_re
      resistant_mean, imaginary(cal_pertile_perpol), 2, res_mean_im
      cal_final[freq_i,pol_i]=res_mean_re+Complex(0,1)*res_mean_im
    endfor
  endfor
  
  stop
  
end