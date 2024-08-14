pro make_res_vis

  file_base = '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/'
  
  pol = ['XX','YY']
  textfast,obsid,/read,/string,file_path='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23_longrunstyle.txt'
  
  for obs_i=0, N_elements(obsid)-1 do begin
    for pol_i=0,1 do begin
    
      vis_ptr=PTR_NEW(/ALLOCATE_HEAP)
      vis_model_ptr = getvar_savefile('/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013zenith_calonly_autovis/cal_prerun/vis_data/'+obsid[obs_i]+'_vis_model_'+pol[pol_i]+'.sav', 'vis_model_ptr')
      vis_dirty_ptr = getvar_savefile(file_base+'vis_data/'+obsid[obs_i]+'_vis_'+pol[pol_i]+'.sav','vis_ptr')
      restore, file_base+'metadata/'+obsid[obs_i]+'_obs.sav'
      *vis_ptr = *vis_dirty_ptr - *vis_model_ptr
      zeroed_inds = where(*vis_dirty_ptr EQ 0, n_count)
      if n_count GT 0 then (*vis_ptr)[zeroed_inds]=0.
      save, vis_ptr, obs, pol_i,filename=file_base+'vis_res/'+obsid[obs_i]+'_vis_'+pol[pol_i]+'.sav'
      
    endfor
  endfor
  
end