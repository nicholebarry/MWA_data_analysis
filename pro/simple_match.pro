pro simple_match

  ;catalog1_path = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_gleamcal/output_data/'
  ;catalog2_path = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_presidelobe/output_data/'

  catalog1_path = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_gleamcal_sidelobe_skip/output_data/'
  ;catalog2_path = '/nfs/mwa-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/'
  catalog2_path = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_southsidelobe/output_data/'
  
  ;catalog1 = getvar_savefile(catalog1_path + '1061316296_source_array.sav', 'source_array')
  ;catalog2 = getvar_savefile(catalog2_path + '1061316296_source_array.sav', 'source_array')
  
  catalog1 = getvar_savefile(catalog1_path + '1061316296_source_array.sav', 'source_array')
  ;catalog2 = getvar_savefile(catalog2_path + 'GLEAM_EGC_catalog_241brighter.sav', 'catalog')
  catalog2 = getvar_savefile(catalog2_path + '1131458856_source_array.sav', 'source_array')
  
  ra2 = FLTARR(N_elements(catalog2))
  dec2 = FLTARR(N_elements(catalog2))
  ra1 = FLTARR(N_elements(catalog1))
  dec1 = FLTARR(N_elements(catalog1))
  flux2 = FLTARR(N_elements(catalog2))
  flux1 = FLTARR(N_elements(catalog1))
  
  ;Northern
  ;subset_dec = [5.,15.]
  ;subset_ra = [-5.,5.]
  
  subset_dec = [-65.,-55.]
  subset_ra = [-5.,5.]
  
  deg_match = .05 ; .017
  
  ;Pull out arguements from catalog2 structure
  for source_j=0, N_elements(catalog2)-1 do begin
    if catalog2[source_j].ra GT 180. then ra2[source_j]=catalog2[source_j].ra - 360. else ra2[source_j]=catalog2[source_j].ra
    dec2[source_j] = catalog2[source_j].dec
    flux2[source_j] = catalog2[source_j].flux.I
  endfor
  
  ;Pull out arguements from catalog1 structure
  for source_i=0, N_elements(catalog1)-1 do begin
    if catalog1[source_i].ra GT 180. then ra1[source_i]=catalog1[source_i].ra - 360. else ra1[source_i] = catalog1[source_i].ra
    dec1[source_i] = catalog1[source_i].dec
    flux1[source_i] = catalog1[source_i].flux.I
  endfor
  
  ;Constrain sources to be within a ra range
  if keyword_set(subset_ra) then begin
    wh_subset_ra2 = where( ra2 GT subset_ra[0] AND ra2 LT subset_ra[1], ncount_sub2)
    if ncount_sub2 GT 0 then ra2 = ra2[wh_subset_ra2] else message, "No sources in RA range"
    flux2 = flux2[wh_subset_ra2]
    dec2 = dec2[wh_subset_ra2]
    
    wh_subset_ra1 = where( ra1 GT subset_ra[0] AND ra1 LT subset_ra[1], ncount_sub1)
    if ncount_sub1 GT 0 then ra1 = ra1[wh_subset_ra1] else message, "No sources in RA range"
    flux1 = flux1[wh_subset_ra1]
    dec1 = dec1[wh_subset_ra1]
  endif
  
  ;Constrain sources to be within a dec range
  if keyword_set(subset_dec) then begin
    wh_subset_dec2 = where( dec2 GT subset_dec[0] AND dec2 LT subset_dec[1], ncount_sub2)
    if ncount_sub2 GT 0 then dec2 = dec2[wh_subset_dec2] else message, "No sources in DEC range"
    flux2 = flux2[wh_subset_dec2]
    ra2 = ra2[wh_subset_dec2]
    
    wh_subset_dec1 = where( dec1 GT subset_dec[0] AND dec1 LT subset_dec[1], ncount_sub1)
    if ncount_sub1 GT 0 then dec1 = dec1[wh_subset_dec1] else message, "No sources in DEC range"
    flux1 = flux1[wh_subset_dec1]
    ra1 = ra1[wh_subset_dec1]
  endif
  
  
  for source_i=0, N_elements(ra1)-1 do begin
  
    ra_match = where(ra2 LT (ra1[source_i]+.017) AND ra2 GT (ra1[source_i]-.017), ncount1)
    
    if ncount1 GT 0 then $
      dec_match = where(abs(dec2[ra_match]) LT (abs(dec1[source_i])+deg_match) AND abs(dec2[ra_match]) GT (abs(dec1[source_i])-deg_match), ncount2)
      
      
    if ncount1 GT 0 then if N_elements(dec_match) EQ 1 then begin
      if N_elements(match2) EQ 0 then match2=ra_match[dec_match] else match2 = [match2, ra_match[dec_match]]
      if N_elements(match1) EQ 0 then match1=source_i else match1 = [match1, source_i]
    endif
  endfor
  
  bright = where((flux2[match2] GT .75),n_count)
  ratio = where(flux1[match1[bright]]/flux2[match2[bright]] GT .5 AND flux1[match1[bright]]/flux2[match2[bright]] LT 1.3)
  
  log_line = linfit(alog(flux1[match1[bright[ratio]]]),alog(flux2[match2[bright[ratio]]]), sigma=sigma,chisqr=chisqr)
  
  stop
  
end