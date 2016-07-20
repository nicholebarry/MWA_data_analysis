pro apparent_brightness

  catalog_path='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_eor_ones_short_baselines_included/output_data/1061316176_skymodel.sav'
  ;catalog_path='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_cal_noeor_ones_maxcalsources_nod_zenithpointing_notileflag_modelmake10/output_data/1061316176_skymodel.sav'
  
  skymodel = getvar_savefile(catalog_path, 'skymodel')
  
  flux_i = FLTARR(N_elements(skymodel.source_list))
  
  for source_i=0, N_elements(skymodel.source_list)-1 do begin
    flux_i[source_i]=skymodel.source_list[source_i].flux.I
  endfor
  
  print, mean(flux_i)
  print, median(flux_i)
  print, minmax(flux_i)
  result=histogram(flux_i,binsize=.1,locations=xbin)
  cgplot, xbin, result, /ylog, /xlog, yrange=[10^(-2.),10^6.]
  stop
  
  ;catalog_path='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/mwa_calibration_source_list_gleam_kgs_no_fornax.sav'
  catalog_path='/nfs/mwa-09/r1/djc/EoR2013/Aug23/calibration_sim/fhd_nb_sim_perfect_cal_noeor_ones_dimcalsources_nod_notileflag/output_data/1061316176_skymodel.sav'
  
  ;skymodel = getvar_savefile(catalog_path, 'catalog')
  ;skymodel = getvar_savefile(catalog_path, 'skymodel')
  restore, catalog_path, /relax
  ;skymodel = catalog
  
  flux_i = FLTARR(N_elements(skymodel.source_list))
  
  for source_i=0, N_elements(skymodel.source_list)-1 do begin
    flux_i[source_i]=skymodel.source_list[source_i].flux.I
  endfor
  
  print, mean(flux_i)
  print, median(flux_i)
  print, minmax(flux_i)
  result=histogram(flux_i,binsize=.1,locations=xbin)
  cgoplot, xbin, result, color='green'
  stop
  
end