pro component_plotter

  for day_i=0,0 do begin
    ;readcol, '/nfs/mwa-00/h1/nbarry/MWA/IDL_code/obs_list/Aug27.txt', obs_id, format='A', /silent
    if day_i EQ 0 then dir = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_presidelobe_Oct23_EoR1/deconvolution/'
    ;if day_i EQ 1 then dir = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_presidelobe_Aug27/deconvolution/'
    
    filenames = FILE_SEARCH(dir, '*_fhd.sav')
    obs_id = FILE_BASENAME(filenames, '_fhd.sav')
    
    if day_i EQ 0 then num_tot_obs = N_elements(obs_id) else num_tot_obs = num_tot_obs + N_elements(obs_id)
    
    ;fits_files=1
    if keyword_set(fits_files) then begin
      for obs_i=0, N_elements(obs_id) -1 do begin
        ;for obs_i=0, 1 do begin
      
        component_array = getvar_savefile(dir+obs_id[obs_i]+'_fhd.sav','component_array')
        
        dec = component_array.dec
        ra = component_array.ra
        id = component_array.id
        freq = component_array.freq
        flux_I = component_array.flux.I
        
        component_struct = {ID:id, Declination:dec, RightAscension:ra, FluxI:flux_I}
        
        ;    hdr=["FIELD    =  EoR0", $
        ;      "DATE     =  Aug23", $
        ;      "FREQ_MHZ =  " + strtrim(freq[0],2), $
        ;      "OBSID    =  " +obs_id[obs_i]]
        
        hdr=["COMMENT   = Field: EoR1, Date: Oct24, MHz: " + strtrim(freq[0],2)+", ObsID: " +obs_id[obs_i],'']
        
        mwrfits, component_struct,dir+obs_id[obs_i]+'_lowEoR1.fits',hdr,/create
      endfor
    endif
    
    ;cgplot, [11+45./60.,12+5./60.],[-25-10./60.,-25-25./60.],/NODATA
    ;ra_range_sculptor = [11.7,12.1]
    ;dec_range_sculptor = [-25.15,-25.45]
    
    ra_range_sculptor = [40,50]
    dec_range_sculptor =  [-10,-20]
    
    for obs_i=0, N_elements(obs_id) -1 do begin
      component_array = getvar_savefile(dir+obs_id[obs_i]+'_fhd.sav','component_array')
      dec = FLTARR(N_elements(component_array))
      ra = FLTARR(N_elements(component_array))
      flux_I = FLTARR(N_elements(component_array))
      
      for component_i=0, N_elements(component_array) -1 do begin
        if component_array[component_i].ra GT 240. then ra[component_i] = component_array[component_i].ra - 360. else ra[component_i] = component_array[component_i].ra
        dec[component_i] = component_array[component_i].dec
        flux_I[component_i] = component_array[component_i].flux.I
      endfor
      
      ra_inds_low = where(ra GT ra_range_sculptor[0],n_count_low)
      if n_count_low LT 0 then continue
      ra_inds_high = where(ra[ra_inds_low] LT ra_range_sculptor[1],n_count_high)
      if n_count_high LT 0 then continue
      
      dec_inds_low = where(dec[ra_inds_low[ra_inds_high]] GT dec_range_sculptor[1],n_count_low)
      if n_count_low LT 0 then continue
      dec_inds_high = where(dec[ra_inds_low[ra_inds_high[dec_inds_low]]] LT dec_range_sculptor[0],n_count_high)
      if n_count_high LT 0 then continue
      
      sculptor_inds = ra_inds_low[ra_inds_high[dec_inds_low[dec_inds_high]]]
      
      if (obs_i EQ 0) AND (day_i EQ 0) then begin
        sculptor_comp_ra = ra[sculptor_inds]
        sculptor_comp_dec = dec[sculptor_inds]
        sculptor_comp_flux_I = flux_I[sculptor_inds]
      endif else begin
        sculptor_comp_ra = [sculptor_comp_ra, ra[sculptor_inds]]
        sculptor_comp_dec = [sculptor_comp_dec, dec[sculptor_inds]]
        sculptor_comp_flux_I = [sculptor_comp_flux_I, flux_I[sculptor_inds]]
      endelse
      
    ;cgoplot, ra, dec, psym = 3
    endfor
  endfor
  
  ;n_bins = 200.
  n_bins = 1000.
  
  binsize = (ra_range_sculptor[1]-ra_range_sculptor[0])/n_bins
  n_bins_dec = CEIL((dec_range_sculptor[0]-dec_range_sculptor[1])/binsize)
  ra_result = histogram(sculptor_comp_ra, binsize=binsize, min=ra_range_sculptor[0], reverse_indices=ri_ra)
  
  pixels = FLTARR(n_bins,n_bins_dec)
  
  ;  n_bins_dec = Ceil((dec_range_sculptor[0]-dec_range_sculptor[1])/binsize)
  
  for bin_i=0,n_bins-2 do begin
    if ri_ra[bin_i] EQ ri_ra[bin_i+1] then continue
    ra_inds = ri_ra[ri_ra[bin_i]:ri_ra[bin_i+1]-1]
    dec_result = histogram(sculptor_comp_dec[ra_inds], binsize=binsize, min=dec_range_sculptor[1], reverse_indices=ri_dec)
    for bin_j=0,N_elements(dec_result)-2 do begin
      if ri_dec[bin_j] EQ ri_dec[bin_j+1] then continue
      dec_inds = ri_dec[ri_dec[bin_j]:ri_dec[bin_j+1]-1]
      pixels[bin_i,bin_j] = total(sculptor_comp_flux_I[ra_inds[dec_inds]])
    endfor
  endfor
  
  pixels = pixels / num_tot_obs ;To get into Jy
  pixel_centers_ra = FINDGEN(n_bins)*binsize+ra_range_sculptor[0]+binsize/2.
  pixel_centers_dec = FINDGEN(n_bins_dec)*binsize + dec_range_sculptor[1] + binsize/2.
  
  stop
  quick_image, pixels, pixel_centers_ra, pixel_centers_dec,xtitle='RA',ytitle='Dec',/log,cb_title='Log of Flux Density', savefile = '/nfs/mwa-00/h1/nbarry/random_section'
end