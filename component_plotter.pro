pro component_plotter

  for day_i=0,0 do begin
    ;readcol, '/nfs/mwa-00/h1/nbarry/MWA/IDL_code/obs_list/Aug27.txt', obs_id, format='A', /silent
    ;if day_i EQ 0 then dir = '/nfs/mwa-10/r1/EoRuvfits/analysis/fhd_nb_decon_gleamcal_sidelobe_skip/deconvolution/'
    if day_i EQ 0 then dir = '/nfs/mwa-08/d1/DiffuseSurvey2015/fhd_rlb_GLEAM_cal_decon_Nov2016/deconvolution/'
    
    filenames = FILE_SEARCH(dir, '*_fhd.sav')
    obs_id = FILE_BASENAME(filenames, '_fhd.sav')
    
    obs_id = ['1131454296', '1131462216', '1131473016', '1131478656', '1131478776', '1131534824', '1131557504', '1131709192']
    
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
    
   ;ra_range_sculptor = [-.5,0]
   ;dec_range_sculptor =  [-60.7,-61.1]
   
   ra_range_sculptor = [-3.2,-2.9]
   dec_range_sculptor =  [-28,-28.3]
    ;ra_range_sculptor = [-1.,1.]
    ;dec_range_sculptor =  [-40.2,-41.7]
    
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
  
  n_bins = 50.
  ;n_bins = 1000.
  
  binsize = (ra_range_sculptor[1]-ra_range_sculptor[0])/n_bins
  ;binsize=.0056
  ;n_bins = Ceil((ra_range_sculptor[1]-ra_range_sculptor[0])/binsize)
  
  n_bins_dec = CEIL((dec_range_sculptor[0]-dec_range_sculptor[1])/binsize)
  ra_result = histogram(sculptor_comp_ra, binsize=binsize, min=ra_range_sculptor[0], reverse_indices=ri_ra)
  n_bins = N_elements(ra_result)
  
  pixels = FLTARR(n_bins,n_bins_dec)
  pixel_weights = FLTARR(n_bins,n_bins_dec)
  
  ;  n_bins_dec = Ceil((dec_range_sculptor[0]-dec_range_sculptor[1])/binsize)
  
  for bin_i=0,n_bins-2 do begin
    if ri_ra[bin_i] EQ ri_ra[bin_i+1] then continue
    ra_inds = ri_ra[ri_ra[bin_i]:ri_ra[bin_i+1]-1]
    dec_result = histogram(sculptor_comp_dec[ra_inds], binsize=binsize, min=dec_range_sculptor[1], reverse_indices=ri_dec)
    for bin_j=0,N_elements(dec_result)-2 do begin
      if ri_dec[bin_j] EQ ri_dec[bin_j+1] then continue
      dec_inds = ri_dec[ri_dec[bin_j]:ri_dec[bin_j+1]-1]
      pixels[bin_i,bin_j] = total(sculptor_comp_flux_I[ra_inds[dec_inds]])
      pixel_weights[bin_i,bin_j] = N_elements(ra_inds[dec_inds])
    endfor
  endfor
  
  pixels1 = pixels / num_tot_obs ;To get into Jy
  pixels2 = pixels * weight_invert(pixel_weights) ;To get into Jy
  pixel_centers_ra = FINDGEN(n_bins)*binsize+ra_range_sculptor[0]+binsize/2.
  pixel_centers_dec = FINDGEN(n_bins_dec)*binsize + dec_range_sculptor[1] + binsize/2.
  
  stop
  ;quick_image, pixels2, pixel_centers_ra, pixel_centers_dec,xtitle='RA',ytitle='Dec',/log,cb_title='Log of Flux Density', savefile = '/nfs/mwa-00/h1/nbarry/random_section'
  
  img2 = GAUSS_SMOOTH(pixels1 ,1)
  maximg2 = max(img2,max_sub)
  s = SIZE(img2)
  ncol = s[1]
  max_col = max_sub MOD ncol
  max_row = max_sub / ncol
  
  pad_col = max_col - (ncol-max_col)
  pad_row = max_row - (s[2]-max_row)
  pad_col_left=0
  pad_col_right=0
  pad_row_left=0
  pad_row_right=0
  if pad_col LT 0 then pad_col_left=abs(pad_col) else pad_col_right=abs(pad_col)
  if pad_row LT 0 then pad_row_left=abs(pad_row) else pad_row_right=abs(pad_row)
  
  img3 = FLTARR(n_bins+abs(pad_col),n_bins_dec+abs(pad_row))
  pixel_centers_ra3 = FLTARR(n_bins+abs(pad_col))
  pixel_centers_dec3 = FLTARR(n_bins_dec+abs(pad_row))
  img3[pad_col_left:s[1]+pad_col_left-1,pad_row_left:s[2]+pad_row_left-1] = img2
  
  pixel_centers_dec3[pad_row_left:s[2]+pad_row_left-1] = pixel_centers_dec
  delta_dec = pixel_centers_dec[0] - pixel_centers_dec[1]
  if pad_row_left GT 0 then $
    for pad_i=0, pad_row_left-1 do pixel_centers_dec3[pad_row_left - 1 - pad_i] = pixel_centers_dec3[pad_row_left - pad_i] - abs(delta_dec)
  if pad_row_right GT 0 then $
    for pad_i=0, pad_row_right-1 do pixel_centers_dec3[s[2]+pad_row_left + pad_i] = pixel_centers_dec3[s[2]+pad_row_left + pad_i - 1] + abs(delta_dec)
  
  pixel_centers_ra3[pad_col_left:s[1]+pad_col_left-1] = pixel_centers_ra
    delta_ra = pixel_centers_ra[0] - pixel_centers_ra[1]
  if pad_col_left GT 0 then $
    for pad_i=0, pad_col_left-1 do pixel_centers_ra3[pad_col_left - 1 - pad_i] = pixel_centers_ra3[pad_col_left - pad_i] - abs(delta_ra)
  if pad_col_right GT 0 then $
    for pad_i=0, pad_col_right-1 do pixel_centers_ra3[s[1]+pad_col_left + pad_i] = pixel_centers_ra3[s[1]+pad_col_left + pad_i - 1] + abs(delta_ra)
  
  img3s = size(img3)
  
  
  
  img4 = FLTARR((img3s[1]/3.)-2,(img3s[2]/3.)-2)
  pixel_centers_dec4 = FLTARR((img3s[2]/3.)-2)
  pixel_centers_ra4 = FLTARR((img3s[1]/3.)-2)
  
  for bin_i=1, (img3s[1]/3.)-2 do $
    for bin_j=1, (img3s[2]/3.)-2 do $
    img4[bin_i-1,bin_j-1] = total(img3[bin_i*3.-1.:(bin_i+1.)*3.-2.,bin_j*3.-1.:(bin_j+1.)*3.-2.])
    
  for bin_j=1, (img3s[2]/3.)-2 do $
    pixel_centers_dec4[bin_j-1] = (pixel_centers_dec3[bin_j*3.-1.:(bin_j+1.)*3.-2.])[1]
;    
;  sub0 = where(pixel_centers_dec4 EQ 0,n_count)
;  dec_min = min((pixel_centers_dec4), sub_min)  
;  dec_max = max(abs(pixel_centers_dec4), sub_max)  
;  dec_delta = pixel_centers_dec4[sub_max] - pixel_centers_dec4[sub_max+1]
;  if n_count GT 0 then $
;    if pixel_centers_dec4[min(sub0)] GT dec_min then $
;    for i=0, N_elements(sub0)-1 do pixel_centers_dec4[sub0[N_elements(sub0)-1-i]] = pixel_centers_dec4[sub_min-i] + dec_delta
;    
  for bin_i=1, (img3s[1]/3.)-2 do $
    pixel_centers_ra4[bin_i-1] = (pixel_centers_ra3[bin_i*3.-1.:(bin_i+1.)*3.-2.])[1]
;    
;  sub0 = where(pixel_centers_ra4 EQ 0,n_count)
;  ra_min = min((pixel_centers_ra4), sub_min)  
;  ra_max = max(abs(pixel_centers_ra4), sub_max)  
;  ra_delta = pixel_centers_ra4[sub_max] - pixel_centers_ra4[sub_max+1]
;  if n_count GT 0 then begin
;    if pixel_centers_ra4[min(sub0)] GT ra_min then $
;    for i=0, N_elements(sub0)-1 do pixel_centers_ra4[sub0[N_elements(sub0)-1-i]] = pixel_centers_ra4[sub_max-i] + ra_delta    
;    if pixel_centers_ra4[min(sub0)] LT ra_min then $
;    for i=0, N_elements(sub0)-1 do pixel_centers_ra4[sub0[i]] = pixel_centers_ra4[sub_min-i] + ra_delta  
;  endif  
;    
  stop
  
  ;img3=rebin(img2[3:74,*],72/4,72/4)
  ;pixel_centers_dec3=rebin(pixel_centers_dec,72/4)
  ;pixel_centers_ra3=rebin(pixel_centers_ra[3:74],72/4)
  img5 = img4
  sub = where(img5 LT .1)
  img5[sub]=0
subs = where(img5 GT 0)
s = SIZE(img5)
ncol = s[1]
col = subs mod ncol
row = subs /ncol
print, img5[subs],pixel_centers_ra4[col], pixel_centers_dec4[row]

  
end