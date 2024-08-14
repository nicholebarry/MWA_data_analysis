pro ac_yaml_to_catalog_sav

make_ac_catalog=1
  if keyword_set(make_ac_catalog) then begin
  
    ;Read the GLEAM catalog from fits
    yaml_name = '/fred/oz048/achokshi/hyper-rm/data/fhd_aee_srclist/RM+49.0_HPX_EXACT.yaml'

    n_freq = 96
    flux_arr = dblarr(n_freq,4)
    freq_arr = dblarr(n_freq)
    freq_i=0

    openr, lun, yaml_name, /GET_LUN
    line = ''
    content = ''
    while ~eof(lun) do begin
      readf, lun, line

      ;RM0000:
      while strmid(line,0,2) eq 'RM' do begin
        source_id = ulong(strmid(line,2,4))
        
        readf, lun, line
        ;  - ra: 0.7538212592762423
        ra = double(strmid(line,8,22))
        readf, lun, line
        ;    dec: -25.749350800573403
        dec = double(strmid(line,9,22))
        ;    comp_type: point
        ;    flux_type:
        ;      list:
        readf, lun, line
        readf, lun, line
        readf, lun, line

        readf, lun, line
        ;        - freq: 167075000.0
        while strmid(line,10,4) eq 'freq' do begin
          freq = double(strmid(line,16,15))

          ;          i: 7.939277977185608
          ;          q: -2.0343656818049145
          ;          u: -1.2386478129826162
          ;          v: 1.1341825681693725

          readf, lun, line
          flux_arr[freq_i,0] = double(strmid(line,13,20))
          readf, lun, line
          flux_arr[freq_i,1] = double(strmid(line,13,20))
          readf, lun, line
          flux_arr[freq_i,2] = double(strmid(line,13,20))
          readf, lun, line
          flux_arr[freq_i,3] = double(strmid(line,13,20))

          freq_arr[freq_i] = freq

          freq_i = freq_i + 1



        ;        - freq: 167075000.0
          if ~eof(lun) then readf, lun, line

        endwhile

        freq_i = 0

        if source_id eq 0 then begin
          full_id = source_id
          full_ra = ra
          full_dec = dec
          full_flux_arr = flux_arr
        endif else begin
          full_id = [full_id,source_id]
          full_ra = [full_ra,ra]
          full_dec = [full_dec,dec]
          full_flux_arr = [[[full_flux_arr]],[[flux_arr]]]
        endelse

      ;RM0001:
      endwhile


    endwhile
    free_lun, lun

    flux_struct = {xx:0.,yy:0.,xy:Complex(0.),yx:Complex(0.),I:0.,Q:0.,U:0.,V:0.}
    full_flux_struct = replicate(flux_struct, N_elements(source_id))
    
    for freq_i=0, n_freq-1 do begin
      for source_i=0, N_elements(source_id)-1 do begin
        full_flux_struct[source_i].I=full_flux_arr[freq_i,0,source_i]
        full_flux_struct[source_i].Q=full_flux_arr[freq_i,1,source_i]
        full_flux_struct[source_i].U=full_flux_arr[freq_i,2,source_i]
        full_flux_struct[source_i].V=full_flux_arr[freq_i,3,source_i]
      endfor
      catalog = source_comp_init(id=full_id, frequency=freq_arr[freq_i], ra=full_ra, dec=full_dec, flux=full_flux_struct)
      save, catalog, filename='/fred/oz048/achokshi/hyper-rm/data/fhd_aee_srclist/RM+49.0_HPX_EXACT_'+strtrim(ulong(freq_arr[freq_i] / 1000),2)+'.sav'
    endfor
    
  endif
  
end
