pro run_check

        filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/beardsley_thesis_list.txt'
        obs_array='empty'
        readcol, filename, obs_array, format='A', /silent
        
        filename = '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly/calibration/*cal.sav'
        output_obsids = file_basename(file_search(filename), '_cal.sav')
        
        match, obs_array, output_obsids, suba, subb
        logic_arr = INTARR(N_elements(obs_array))
        logic_arr[suba]=1
        logic_inds = where(logic_arr EQ 0,n_count)
        missing_obs = obs_array[logic_inds]
        
        if n_count GT 0 then print, transpose(missing_obs) else print, "No missing obs"

        

end