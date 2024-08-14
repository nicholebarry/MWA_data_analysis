pro cal_simulation_noise_sav


  filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23plusone.txt'
  ;filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/Aug23zenith.txt'
  ;filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/obs_id_6176.txt'
  
  readcol, filename, obs_id, format='A', /silent
  parsednumbers=N_elements(obs_id)
  
  vis_arr=PTRARR(2,/allocate)
  cal_sim_input = 'fhd_nb_sim_unflagged_nodiffuse_onebeam_zenithpointing_calvisflag_overfit'
  
  for obs_i=0, parsednumbers-1 do begin
  
    plusone=['1061317272','1061317400','1061317520','1061317640','1061317760','1061317888','1061318008','1061318128','1061318248', $
      '1061318376','1061318496','1061318616','1061318736','1061318864','1061318984']
      
    zenith = ['1061315448','1061315568','1061315688','1061315808','1061315936','1061316056','1061316176','1061316296','1061316424', $
      '1061316544','1061316664','1061316784','1061316912','1061317032','1061317152']
      
    match_index=where(STRMATCH(plusone, obs_id[obs_i]),n_count)
    If n_count GT 0 then obs_temp = zenith[match_index] else obs_temp=obs_id[obs_i]
    
    vis_XX_model = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/'+cal_sim_input+'/vis_data/'+obs_temp+'_vis_model_XX.sav', 'vis_model_ptr')
    vis_YY_model = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/'+cal_sim_input+'/vis_data/'+obs_temp+'_vis_model_YY.sav', 'vis_model_ptr')
    *vis_arr[0]=*vis_XX_model
    *vis_arr[1]=*vis_YY_model
    
      obs = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/'+cal_sim_input+'/metadata/'+obs_temp+'_obs.sav', 'obs')
    
    vis_add_noise_simulation, cal_sim_input, vis_arr, obs_id[obs_i],obs
  endfor
  
  
end