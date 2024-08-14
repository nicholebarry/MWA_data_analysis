pro combine_vis_in_freq

  cadence=20.
  n_freq=384.
  ;obsid='1061316296'
  obsid='1163765528'
  file_path_root='/nfs/mwa-04/r1/EoRuvfits/analysis/calibration_sim/fhd_nb_hash_phaseII_'
  n_chunks=floor(n_freq/cadence)
  
  for chunk_i=0, n_chunks-1 do begin
  
    freq_chunk = ulong((chunk_i+1)*(cadence))
    freq_chunk_start = ulong((chunk_i)*cadence)+1
    freq_chunk_end = ulong((chunk_i+1)*(cadence))+1
    if freq_chunk GT n_freq then freq_chunk=ulong(n_freq)-1
    file_path=file_path_root+strtrim(freq_chunk,2)
    
    restore, file_path + '/vis_data/'+obsid+'_vis_XX.sav'
    if chunk_i EQ 0 then vis_total_XX = Pointer_copy(vis_ptr)
    (*vis_total_XX)[freq_chunk_start:freq_chunk_end-1,*] = (*vis_ptr)[freq_chunk_start:freq_chunk_end-1,*]
    
    restore, file_path + '/vis_data/'+obsid+'_vis_YY.sav'
    if chunk_i EQ 0 then vis_total_YY = Pointer_copy(vis_ptr)
    (*vis_total_YY)[freq_chunk_start:freq_chunk_end-1,*] = (*vis_ptr)[freq_chunk_start:freq_chunk_end-1,*]
    
  endfor
  
  stop
  save, vis_total_XX, filename='/nfs/mwa-04/r1/EoRuvfits/analysis/calibration_sim/fhd_nb_hash_phaseII/vis_data/'+obsid+'_vis_XX.sav'
  save, vis_total_YY, filename='/nfs/mwa-04/r1/EoRuvfits/analysis/calibration_sim/fhd_nb_hash_phaseII/vis_data/'+obsid+'_vis_YY.sav'
  stop
  
end