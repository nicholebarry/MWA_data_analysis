pro cal_plotter_wrapper
  
  ;input directory
  dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_cross_manual/'

  ;read obs list
  filename='/home/nbarry/MWA/pipeline_scripts/bash_scripts/ozstar/obs_list/btl_noalltv_noocc4_minusone.txt'
  readcol, filename, obs_list, format='A', /silent
  n_obs = N_elements(obs_list)

  for obs_i=0, n_obs-1 do begin

    ;read sav files 
    restore, dir+'metadata/'+obs_list[obs_i]+'_obs.sav'
    restore, dir+'calibration/'+obs_list[obs_i]+'_cal.sav'
    plot_cals,cal,obs,file_path_base=dir+'output_images/'+obs_list[obs_i]

  endfor



end
