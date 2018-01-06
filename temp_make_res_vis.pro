pro temp_make_res_vis

  file_base = '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_2013longrun_std/'
  
  pol = ['XX','YY']
  textfast,obsid,/read,/string,file_path='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23_longrunstyle.txt'
  
  for obs_i=0, N_elements(obsid)-1 do begin
    for pol_i=0,1 do begin
    
      vis_ptr = getvar_savefile(file_base+'vis_res/'+obsid[obs_i]+'_vis_'+pol[pol_i]+'.sav','vis_ptr')
      restore, file_base+'metadata/'+obsid[obs_i]+'_obs.sav'
      save, vis_ptr, obs, pol_i,filename=file_base+'vis_res/'+obsid[obs_i]+'_vis_'+pol[pol_i]+'.sav'
      
    endfor
  endfor
  
end