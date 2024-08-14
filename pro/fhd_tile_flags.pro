pro fhd_tile_flags

   file_path = '/fred/oz048/MWA/CODE/FHD/fhd_nb_calstop_fullbeam_GLEAMtenth/metadata/'
   file_path_missing = '/fred/oz048/MWA/CODE/FHD/fhd_nb_calstop_fullbeam_GLEAMtenth_missing/metadata/'
   obs_file = '/home/nbarry/MWA/FHD/obs_list/btl_noalltv_noocc4_minustwo.txt'
   readcol, obs_file, obs_arr, format='A', /silent

   for obs_i=0,N_elements(obs_arr)-1 do begin
     restore, file_path + obs_arr[obs_i] + '_obs.sav'
     if file_test(file_path_missing + obs_arr[obs_i] + '_obs.sav') then begin
       restore, file_path_missing + obs_arr[obs_i] + '_obs.sav'
     endif
     tile_flags = where((*obs.baseline_info).tile_use EQ 0,n_count)
     if n_count GT 0 then print, obs_arr[obs_i], tile_flags else print, obs_arr[obs_i]
   endfor


end
