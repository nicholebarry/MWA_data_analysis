pro confusion_limit_catalog

  id = 3
  schetcher=1
  
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
    
    num_sources=10000. ;Gives 20 sources between 9 and 10 Jy (we see 10 in real data)
    num_per_bin=density_tot/density_norm*num_sources
    num_per_bin=FIX(num_per_bin)
    num_sources=Total(num_per_bin) ;After making integers
    
    Jy_sources=FLTARR(Total(num_per_bin))
    
    Jy_sources[0:TOTAL(num_per_bin[0:1])-1]=Jy_bins[0]
    For source_bin_i=2, N_elements(num_per_bin)-2 do begin
      Jy_sources[TOTAL(num_per_bin[0:source_bin_i-1])-1:TOTAL(num_per_bin[0:source_bin_i])-1]=Jy_bins[source_bin_i-1]
    endfor
    
    stop
  endelse
  
  ;dec
  dec=-(65*RANDOMU(Seed,num_sources)-15) ;about 0
  
  ;ra -- approx two generators about 0 line with half the numbers
  ra=(50*RANDOMU(Seed,num_sources)-25)
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
  
  catalog = source_comp_init(id=source_id, frequency=freq, ra=ra, dec=dec, flux=Jy_sources)
  
  complete_sources = where(Jy_sources GE .1, n_count)
  
  completeness_catalog = source_comp_init(id=source_id, frequency=freq, ra=ra[complete_sources], dec=dec[complete_sources], flux=Jy_sources[complete_sources])
  
  stop
  save_path=filepath('confusion'+number_formatter(id)+'.sav',root=rootdir('FHD'),subdir='catalog_data')
  save, catalog, filename=save_path
  
  catalog=completeness_catalog
  save_path=filepath('confusion_completeness'+number_formatter(id)+'.sav',root=rootdir('FHD'),subdir='catalog_data')
  save, catalog, filename=save_path
  
end