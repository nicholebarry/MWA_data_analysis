pro lobes_to_catalog_sav

make_lobes_catalog=1
  if keyword_set(make_lobes_catalog) then begin
  
    ;Read the lobes catalog from fits
    ;lobes_fits = '/home/nbarry/tertiary_gits/FHD/catalog_data/LoBES_EoR0_FULL_FINAL_SpecGaussModCat_03MAR2022.fits'
    lobes_fits = '/home/nbarry/tertiary_gits/FHD/catalog_data/GLMv2_LoBES03Mar2022_combine.fits'
    lobes_read = mrdfits(lobes_fits, 1)
    n_source = (size(lobes_read))[1]
print, n_source   
 
    ;Initialize
    pre_ra = FLTARR(n_source)
    pre_dec = FLTARR(n_source)
    pre_source_id = FLTARR(n_source)
    pre_Jy_sources = FLTARR(n_source)
    pre_err_Jy_sources = FLTARR(n_source)
    pre_major_dc = FLTARR(n_source)
    pre_minor_dc = FLTARR(n_source)
    pre_pa_dc = FLTARR(n_source)

    shape_struct={x:0.,y:0.,angle:0.}
    shape_struct_new=Replicate(shape_struct,n_source>1)
    
    for source_i=0,n_source-1 do begin
;      ;If the source flux is above the std of the flux, count it as a consistent source
;      if not ((gleam_read[source_i].int_flux_181 - gleam_read[source_i].err_int_flux_181) LT 0) and not (gleam_read[source_i].int_flux_181 LT 0) then begin
;        if keyword_set(consistent_sources) then consistent_sources = [consistent_sources, source_i] else consistent_sources=source_i
;      endif
      ;Build the inputs for the catalog from all sources
      pre_source_id[source_i] = source_i
      pre_ra[source_i]=lobes_read[source_i].RA
      pre_dec[source_i]=lobes_read[source_i].DEC
      pre_Jy_sources[source_i]=lobes_read[source_i].int_flx181
      pre_err_Jy_sources[source_i]=lobes_read[source_i].err_int_flx181
      ;shape_struct={x:0.,y:0.,angle:0.}
      pre_major_dc[source_i]=lobes_read[source_i].major_dc
      pre_minor_dc[source_i]=lobes_read[source_i].minor_dc
      pre_pa_dc[source_i]=lobes_read[source_i].pa_dc
    endfor
    freq=181
    shape_struct_new.x=pre_major_dc
    shape_struct_new.y=pre_minor_dc
    shape_struct_new.angle=pre_pa_dc
    ;Take catalog inputs and build a FHD catalog from the consistent sources only
    ;catalog = source_comp_init(id=pre_source_id[consistent_sources], frequency=freq, ra=pre_ra[consistent_sources], dec=pre_dec[consistent_sources], flux=pre_Jy_sources[consistent_sources])
    catalog = source_comp_init(id=pre_source_id, frequency=freq, ra=pre_ra, dec=pre_dec, flux=pre_Jy_sources, shape=shape_struct_new)
    save, catalog, filename='/home/nbarry/tertiary_gits/FHD/catalog_data/LoBES_GLEAM_EoR0_181MHz.sav'
  endif

stop
  
  ;Snippet to analysis resultant source lists and their differences
  ;;;;;;
  gleam_catalog = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_gleam_firstpass/output_data/1061316296_skymodel.sav','skymodel')
  patti_catalog = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_patti_catalog/output_data/1061316296_skymodel.sav','skymodel')
  
  gleam_catalog_source_count = (size(gleam_catalog.source_list))[1]
  patti_catalog_source_count = (size(patti_catalog.source_list))[1]
  
  gleam_RA = FLTARR(gleam_catalog_source_count)
  gleam_DEC= FLTARR(gleam_catalog_source_count)
  gleam_flux= FLTARR(gleam_catalog_source_count)
  
  patti_RA = FLTARR(patti_catalog_source_count)
  patti_DEC= FLTARR(patti_catalog_source_count)
  patti_flux= FLTARR(patti_catalog_source_count)
  
  for source_i=0,patti_catalog_source_count-1 do begin
    patti_RA[source_i] = patti_catalog.source_list[source_i].RA
    patti_DEC[source_i] = patti_catalog.source_list[source_i].DEC
    patti_flux[source_i] = patti_catalog.source_list[source_i].FLUX.I
  endfor
  
  for source_i=0,gleam_catalog_source_count-1 do begin
    gleam_RA[source_i] = gleam_catalog.source_list[source_i].RA
    gleam_DEC[source_i] = gleam_catalog.source_list[source_i].DEC
    gleam_flux[source_i] = gleam_catalog.source_list[source_i].FLUX.I
  endfor
  
  binsize=.01
  gleam_result=histogram(gleam_flux,binsize=binsize,locations=gleam_locations,omax=omax, /NAN,reverse_indices=ri)
  patti_result=histogram(patti_flux,binsize=binsize,locations=patti_locations,omax=omax, /NAN,reverse_indices=ri)

  x_arr_gleam=[gleam_locations[0],gleam_locations+binsize/2.,gleam_locations[N_elements(gleam_locations)-1]-binsize/2.]
  y_arr_gleam=[gleam_result[0],gleam_result,gleam_result[N_elements(gleam_result)-1]]
  
  x_arr_patti=[patti_locations[0],patti_locations+binsize/2.,patti_locations[N_elements(patti_locations)-1]-binsize/2.]
  y_arr_patti=[patti_result[0],patti_result,patti_result[N_elements(patti_result)-1]]
  
  cgplot, x_arr_gleam, y_arr_gleam, xrange=[0,1], xtitle="Integrated Flux Density (Jy)",ytitle='Count', title='Catalog comparison' 
  cgoplot, x_arr_patti, y_arr_patti, xrange=[0,2], color='blue'
  
  cgplot, x_arr_gleam, y_arr_gleam, xrange=[0,3], yrange=[10^0.,10^4.], xtitle="Integrated Flux Density (Jy)",ytitle='Count', title='Catalog comparison, Semi-log',/ylog
  cgoplot, x_arr_patti, y_arr_patti, xrange=[0,4], yrange=[10^0.,10^4.], /ylog,color='blue'
  ;;;;;;
  stop
  
end
