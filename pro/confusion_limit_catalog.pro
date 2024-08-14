pro confusion_limit_catalog

  id = 3
  ;schetcher=1
  
  ;Schetcher model using fits from Willott 2001, no evolution
  if keyword_set(schetcher) then begin
    str_to_m2 = 12.5664
    normdens_low = 10^(-7.12)
    alpha_low = .539
    Jystar_low = 10^(26.1)/str_to_m2*10^(-26.)
    normdens_high = 10^(-6.196)
    alpha_high = 2.27
    Jystar_high = 10^(26.95)/str_to_m2*10^(-26.)
    
    Jy_bins = FINDGEN(1000)*.01+.01
    density_low=normdens_low*(Jy_bins/Jystar_low)^(-alpha_low)*exp(-Jy_bins/Jystar_low)
    density_high=normdens_high*(Jy_bins/Jystar_high)^(-alpha_high)*exp(-Jystar_high/Jy_bins)
    density_tot=density_low+density_high
    density_norm=TOTAL(density_tot)
    
    num_sources=5000. ;Gives 20 sources between 9 and 10 Jy (we see 10 in real data)
    num_per_bin=density_tot/density_norm*num_sources
    num_per_bin=FIX(num_per_bin)
    num_sources=Total(num_per_bin) ;After making integers
    
    Jy_sources=FLTARR(Total(num_per_bin))
    
    Jy_sources[0:TOTAL(num_per_bin[0:1])-1]=Jy_bins[0]
    For source_bin_i=2, N_elements(num_per_bin)-2 do begin
      Jy_sources[TOTAL(num_per_bin[0:source_bin_i-1])-1:TOTAL(num_per_bin[0:source_bin_i])-1]=Jy_bins[source_bin_i-1]
    endfor
  endif else begin
  
    fita = exp(6.64084)
    fitb = -2.13177
    
    Jy_bins = FINDGEN(1000)*.01+.01
    density_tot=fita*Jy_bins^(fitb)
    density_norm=TOTAL(density_tot)
    
    ; num_sources=10000. ;Gives 20 sources between 9 and 10 Jy (we see 10 in real data)
    num_sources=5000. ;temp
    num_per_bin=density_tot/density_norm*num_sources
    num_per_bin=FIX(num_per_bin)
    num_sources=Total(num_per_bin) ;After making integers
    
    Jy_sources=FLTARR(Total(num_per_bin))
    
    Jy_sources[0:TOTAL(num_per_bin[0:1])-1]=Jy_bins[0]
    For source_bin_i=2, N_elements(num_per_bin)-2 do begin
      Jy_sources[TOTAL(num_per_bin[0:source_bin_i-1])-1:TOTAL(num_per_bin[0:source_bin_i])-1]=Jy_bins[source_bin_i-1]
    endfor
    
    
  endelse
  
  ;dec
  ;dec=-(65*RANDOMU(Seed,num_sources)-15) ;about 0
  ;dec=-(20*RANDOMU(Seed,num_sources)+15) ;about 0
  num_sources=8000
  Jy_sources=FLTARR(num_sources)
  Jy_sources[*]=1
  t = 2.*!Pi*RANDOMU(Seed,num_sources)
  u = (107.-4)/2.*sqrt(RANDOMU(Seed,num_sources));+14.*RANDOMU(Seed,num_sources)
  ;fold_index = where(u GT 7, n)
  r = u
  ;r[fold_index] = 14- u[fold_index]
  ra = r*cos(t)-.2
  dec = r*sin(t)-26.8
  
  ;ra -- approx two generators about 0 line with half the numbers
  ;ra=(50*RANDOMU(Seed,num_sources)-25)
  ; ra=(20*RANDOMU(Seed,num_sources)-10)
  ; ra2=(25*RANDOMU(Seed,num_sources/2)+335)
  ;ra=[ra1,ra2]
  
  
  freq = 182.475; MHz, the center of the band.
  
  ;q_frac = 0
  ;u_frac = 0
  ;v_frac = 0
  
  flux_i = Jy_sources
  ;flux_q = flux_i * q_frac
  ;flux_u = flux_i * u_frac
  ;flux_v = flux_i * v_frac
  
  ;catalog = source_comp_init(id=source_id, frequency=freq, ra=ra, dec=dec, flux=Jy_sources)
  
  uniform=1
  if keyword_set(uniform) then begin
  
    catalog = source_comp_init(id=source_id, frequency=freq, ra=ra, dec=dec, flux=Jy_sources)
    save_path=filepath('uniform_full.sav',root=rootdir('FHD'),subdir='catalog_data')
    stop
    save, catalog, filename=save_path
    annuli=1
    if keyword_set(annuli) then begin
      primary = where(r LE 20.)
      negative_primary = where(r GT 20.)
      nulls = where(r GT 20. AND r LE 26.)
      negative_nulls = where(r LE 20. OR r GT 26.)
      sidelobe = where(r GT 26.)
      negative_sidelobe = where( r LE 26.)
      stop
      catalog = source_comp_init(id=source_id, frequency=freq, ra=ra[primary], dec=dec[primary], flux=Jy_sources[primary])
      save_path=filepath('uniform_primary.sav',root=rootdir('FHD'),subdir='catalog_data')
      save, catalog, filename=save_path
      catalog = source_comp_init(id=source_id, frequency=freq, ra=ra[nulls], dec=dec[nulls], flux=Jy_sources[nulls])
      save_path=filepath('uniform_nulls.sav',root=rootdir('FHD'),subdir='catalog_data')
      save, catalog, filename=save_path
      catalog = source_comp_init(id=source_id, frequency=freq, ra=ra[sidelobe], dec=dec[sidelobe], flux=Jy_sources[sidelobe])
      save_path=filepath('uniform_sidelobe.sav',root=rootdir('FHD'),subdir='catalog_data')
      save, catalog, filename=save_path
      
      catalog = source_comp_init(id=source_id, frequency=freq, ra=ra[negative_primary], dec=dec[negative_primary], flux=Jy_sources[negative_primary])
      save_path=filepath('uniform_negative_primary.sav',root=rootdir('FHD'),subdir='catalog_data')
      save, catalog, filename=save_path
      catalog = source_comp_init(id=source_id, frequency=freq, ra=ra[negative_nulls], dec=dec[negative_nulls], flux=Jy_sources[negative_nulls])
      save_path=filepath('uniform_negative_nulls.sav',root=rootdir('FHD'),subdir='catalog_data')
      save, catalog, filename=save_path
      catalog = source_comp_init(id=source_id, frequency=freq, ra=ra[negative_sidelobe], dec=dec[negative_sidelobe], flux=Jy_sources[negative_sidelobe])
      save_path=filepath('uniform_negative_sidelobe.sav',root=rootdir('FHD'),subdir='catalog_data')
      save, catalog, filename=save_path
    endif
  endif
  
  ;double_lobe=1
  if keyword_set(double_lobe) then begin
    save_path=filepath('dlt_point_sources'+number_formatter(id)+'.sav',root=rootdir('FHD'),subdir='catalog_data')
    save, catalog, filename=save_path
    
    brightest = reverse(sort(Jy_sources))
    top_brightest = brightest[0:499]
    double_index = Round(500*RANDOMU(Seed,100))
    sig = .0055
    lobe_sep = sig * sqrt(-2. * alog(1. - RANDOMU(Seed,100))) ;to get a rayleigh distributed random number
    deg_angle = 2.*!Pi*RANDOMU(Seed,100)
    ra_sep = lobe_sep * sin(deg_angle)
    dec_sep = lobe_sep * cos(deg_angle)
    
    ra[top_brightest[double_index]] = ra[top_brightest[double_index]] - ra_sep/2.
    dec[top_brightest[double_index]] = dec[top_brightest[double_index]] - dec_sep/2.
    ra = [ra,ra[top_brightest[double_index]] + ra_sep/2.]
    dec = [dec,dec[top_brightest[double_index]] + dec_sep/2.]
    Jy_sources = [Jy_sources,Jy_sources[top_brightest[double_index]]/2.]
    Jy_sources[top_brightest[double_index]] = Jy_sources[top_brightest[double_index]]/2.
    
    catalog = source_comp_init(id=source_id, frequency=freq, ra=ra, dec=dec, flux=Jy_sources)
    save_path=filepath('dlt_double_lobed'+number_formatter(id)+'.sav',root=rootdir('FHD'),subdir='catalog_data')
    save, catalog, filename=save_path
    
  endif
  
  if keyword_set(completeness) then begin
    complete_sources = where(Jy_sources GE .1, n_count)
    completeness_catalog = source_comp_init(id=source_id, frequency=freq, ra=ra[complete_sources], dec=dec[complete_sources], flux=Jy_sources[complete_sources])
    
    stop
    save_path=filepath('confusion'+number_formatter(id)+'.sav',root=rootdir('FHD'),subdir='catalog_data')
    save, catalog, filename=save_path
    
    catalog=completeness_catalog
    save_path=filepath('confusion_completeness'+number_formatter(id)+'.sav',root=rootdir('FHD'),subdir='catalog_data')
    save, catalog, filename=save_path
  endif
  
end