pro EO_to_catalog_sav

  make_EO_catalog=1
  if keyword_set(make_EO_catalog) then begin
  
    ;Read the EO catalog from fits
    file_name = 'eor0_clustered_points'
    eo_fits = '/nfs/mwa-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/'+file_name+'.fits'
    eo_read = mrdfits(eo_fits, 1)
    n_source = (size(eo_read))[1]
    
    ;Initialize
    pre_source_id = FLTARR(n_source)
    pre_x = FLTARR(n_source)
    pre_y = FLTARR(n_source)
    pre_ra = FLTARR(n_source)
    pre_dec = FLTARR(n_source)
    pre_ston = FLTARR(n_source)
    pre_freq = FLTARR(n_source)
    pre_alpha = FLTARR(n_source)
    pre_gain = FLTARR(n_source)
    ;pre_flag = FLTARR(n_source)
    ;pre_xx = FLTARR(n_source)
    ;pre_yy = FLTARR(n_source)
    ;pre_yx = FLTARR(n_source)
    pre_I = FLTARR(n_source)
    ;pre_Q = FLTARR(n_source)
    ;pre_U = FLTARR(n_source)
    ;pre_V = FLTARR(n_source)
    
    for source_i=0,n_source-1 do begin
      ;Build the inputs for the catalog from all sources
      pre_source_id[source_i] = eo_read[source_i].all_sources_id
      pre_x[source_i] = eo_read[source_i].all_sources_x
      pre_y[source_i] = eo_read[source_i].all_sources_y
      pre_ra[source_i]=eo_read[source_i].all_sources_ra
      pre_dec[source_i]=eo_read[source_i].all_sources_dec
      pre_ston[source_i]=eo_read[source_i].all_sources_ston
      pre_alpha[source_i]=eo_read[source_i].all_sources_alpha
      pre_gain[source_i]=eo_read[source_i].all_sources_gain
      pre_freq[source_i]=eo_read[source_i].all_sources_freq
      pre_I[source_i]=eo_read[source_i].all_sources_I
    endfor
    
    ;Take catalog inputs and build a FHD catalog from the consistent sources only
    catalog = source_comp_init(n_sources = n_source, xvals = pre_x, yvals = pre_y, frequency = pre_freq, ra=pre_ra, dec=pre_dec, flux=pre_I, id=pre_source_id,$
      StoN=pre_ston, alpha=pre_alpha, gain_factor = pre_gain)
    save, catalog, filename='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/'+file_name+'.sav'
  endif
  
  
end