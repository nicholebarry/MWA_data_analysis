pro iono_check

  ;obs, median offset in arcmin, PCA eigenvalue
  ;textfast,data_array,/read,file_path='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/MWA_data_analysis/cjordan_iono_metric_2013-2015.txt'
  textfast,data_array,/read,file_path='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/MWA_data_analysis/cjordan_iono_metric_bmetric.txt'
  textfast,obs_array,/read,/string,file_path='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/MWA_data_analysis/cjordan_iono_metric_bmetric.txt'
  
  
  ;CJordan's full metric
  full_metric = fltarr(N_elements(reform(data_array[0,*])))
  low_PCA = where(data_array[2,*] LT 0.6,n_count)
  if n_count GT 0 then full_metric[low_PCA] = 25.*data_array[1,low_PCA]
  high_PCA = where(data_array[2,*] GE 0.6,n_count)
  If n_count GT 0 then full_metric[high_PCA] = 25.*data_array[1,high_PCA] + 64.*(data_array[2,high_PCA]-0.6)
  
  bad_obs = where(full_metric GT 3.,n_count)
  for string_i=0, N_elements(data_array[0,*])-1 do obs_array[string_i] = (strsplit(obs_array[string_i],' ',/extract))[0]
  bad_obs_name = obs_array[bad_obs]
  print, n_count
  
  stop
  
end