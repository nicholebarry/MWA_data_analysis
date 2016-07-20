PRO vis_noise_calc_simcheck,restore_all=restore_all,interleave=interleave,vis_arr=vis_arr,flag_arr=flag_arr,nt=nt
  ; Script to plot visibility sigmas given various time interleaving or 2 second time visibilities differenced at various time steps

  ;****Restore the necessary information from the standard run to run this script outside of FHD.
  if keyword_set(restore_all) then begin
    RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316176_obs.sav' ;restore obs structure
    vis_XX_restored = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod_zenithpointing_notileflag/vis_data/1061316176_vis_XX.sav', 'vis_ptr') ;restore array of calibrated visibilities
    vis_YY_restored = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod_zenithpointing_notileflag/vis_data/1061316176_vis_YY.sav', 'vis_ptr')
    
    ;Combine the calibrated visibilities in the correct format for the script
    vis_arr = [vis_XX_restored, vis_YY_restored]
    
    RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod_zenithpointing_notileflag/vis_data/1061316176_flags.sav'
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_maxcalsources_nod_zenithpointing_notileflag/calibration/1061316176_cal.sav'
    
  endif
  ;****End of restore the necessary information from the standard run to run this script outside of FHD.
  
  ;TEMP
  RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316176_obs.sav' ;restore obs structure
  
  ;****Setup
  n_pol=obs.n_pol
  n_freq=Long(obs.n_freq)
  noise_arr=fltarr(n_pol,n_freq)
  pol_i=0
  
  
  ;Use a modified script to get the time bins to use for the visibility differences. The return value is the
  ;visibilility bins the same size as the visibility array. Input time determines which bins to pick out,
  ;so time 1 gives even odd differences, time 3 groups bins as 0,6,12... and 3,9,15...
  bin_i=split_vis_flags(obs,flag_arr,bi_use=bi_use,/preserve_flags)
  
  FOR pol_i=0,n_pol-1 DO BEGIN
    data_diff =Imaginary( (*vis_arr[pol_i])[*,*bi_use[0]])-Imaginary((*vis_arr[pol_i])[*,*bi_use[1]]) ; only use imaginary part
    flag_diff = ((*flag_arr[pol_i])[*,*bi_use[0]]>0)*((*flag_arr[pol_i])[*,*bi_use[1]]>0)
    FOR fi=0L,n_freq-1 DO BEGIN
      ind_use=where(flag_diff[fi,*],n_use)
      IF n_use GT 1 THEN noise_arr[pol_i,fi]=Stddev(data_diff[fi,ind_use])/Sqrt(2.)
    ENDFOR
  ENDFOR
  
  stop
  
END