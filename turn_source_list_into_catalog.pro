pro turn_source_list_into_catalog

  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23.txt'
  readcol, filename, obs_ids, format='A', /silent
  
  for obs_i=0, N_elements(obs_ids)-1 do begin
    restore,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_March2016_small/output_data/'+obs_ids[obs_i]+'_source_array.sav'
    catalog = source_array
    save, catalog, filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_March2016_small/output_data/'+obs_ids[obs_i]+'_source_array2.sav'
  endfor
  
end
