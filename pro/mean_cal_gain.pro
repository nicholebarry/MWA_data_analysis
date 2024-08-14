pro mean_cal_gain

  ;file_base = '/data5/mit_backup/fhd_apb_EoR0_high_sem1_1/calibration/'
  file_base = '/fred/oz048/MWA/CODE/FHD/fhd_nb_calstop_fullbeam_GLEAMtenth/calibration/'
  file_base_missing = '/fred/oz048/MWA/CODE/FHD/fhd_nb_calstop_fullbeam_GLEAMtenth_missing/calibration/'

  obs_list = '/home/nbarry/MWA/FHD/obs_list/beardsley_thesis_list.txt'
  ;obs_list = '/home/nbarry/MWA/FHD/obs_list/beardsley_thesis_list.txt'

  readcol, obs_list, obs_id, format='(A10)'

  restore, '/fred/oz048/MWA/CODE/FHD/fhd_nb_calstop_fullbeam_GLEAMtenth/metadata/1061316296_obs.sav'
  ;restore, '/data5/mit_backup/fhd_apb_EoR0_high_sem1_1/metadata/1061316296_obs.sav'
  freq_flag_inds = where((*obs.baseline_info).freq_use EQ 0)
  tile_flag_inds = where((*obs.baseline_info).tile_use EQ 0)

  for obs_i=0, 1028 do begin

    file_logic = file_test(file_base + obs_id[obs_i] + '_cal.sav')
    if file_logic then restore, file_base + obs_id[obs_i] + '_cal.sav' else restore, file_base_missing + obs_id[obs_i] + '_cal.sav'
    
    cal0 = (*cal.gain[0]); + (*cal.gain_residual[0])
    cal1 = (*cal.gain[1]); + (*cal.gain_residual[1])

    cal0[freq_flag_inds,*] = 0
cal0[0:16,*]=0
cal0[256:383,*]=0
    cal0[*,tile_flag_inds] = 0
    cal1[freq_flag_inds,*] = 0
cal1[0:16,*]=0
cal1[256:383,*]=0
    cal1[*,tile_flag_inds] = 0
    zero_inds = where(abs(cal0) EQ 0, n_count)
    if n_count GT 0 then cal0[zero_inds] = !Values.F_NAN
    zero_inds = where(abs(cal1) EQ 0, n_count)
    if n_count GT 0 then cal1[zero_inds] = !Values.F_NAN

    print, mean(abs(cal0[*,*]),/NAN), mean(abs(cal1[*,*]),/NAN)

  endfor

end
