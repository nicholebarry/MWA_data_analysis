pro bin_check

  n_ch = 1
  ch_bands = PTRARR(n_ch,/allocate)
  ;*ch_bands[0] = [0,95]
  ;*ch_bands[1] = [10,105]
  ;*ch_bands[2] = [0,127]
  ;*ch_bands[3] = [10,127]
  ;*ch_bands[4] = [10,117]
  ;*ch_bands[0]=[0,127]
  *ch_bands[0]=[10,127]
  ;*ch_bands[0]=[8,127]
  ;*ch_bands[3]=[0,117]
  ;*ch_bands[4]=[5,117]
  ;*ch_bands[1]=[10,117]
  ;*ch_bands[2]=[5,120]
  
  n_wedge=1
  wedge_angles = [110]
  ;wedge_angles = [115,125]
  
  n_kperp = 1
  kperp_range_ptr = PTRARR(n_kperp,/allocate)
  ;*kperp_range_ptr[0] = [10,60]
  ;*kperp_range_ptr[1] = [10,80]
  ;*kperp_range_ptr[2] = [10,100]
  ;*kperp_range_ptr[3] = [15,80]
  ;*kperp_range_ptr[4] = [12,80]
  ;*kperp_range_ptr[0]=[10,75]
  ;*kperp_range_ptr[1]=[15,75]
  ;*kperp_range_ptr[2]=[20,75]
  ;*kperp_range_ptr[0]=[10,80]
  ;*kperp_range_ptr[1]=[15,80]
  ;*kperp_range_ptr[5]=[20,80]
  ;*kperp_range_ptr[6]=[10,85]
  ;*kperp_range_ptr[7]=[15,85]
  ;*kperp_range_ptr[8]=[20,85]
  *kperp_range_ptr[0]=[20,80]
  ;*kperp_range_ptr[1]=[20,95]
  ;*kperp_range_ptr[2]=[20,75]
  
  ;generate=1
  if keyword_set(generate) then begin
  
    for ch_i=0, n_ch-1 do begin
      for wedge_i=0,n_wedge-1 do begin
        for kperp_i=0,n_kperp-1 do begin
        
          ps_wrapper, '/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/', 'beardsley_thesis_list',$
            /exact_obsname, ps_foldername='ps_bh',image_window_name='Blackman-Harris', pol_inc='yy',plot_stdset=0,max_uv_lambda=100,$
            wedge_angles=wedge_angles[wedge_i] , freq_ch_range=*ch_bands[ch_i], kperp_range_lambda_1dave=*kperp_range_ptr[kperp_i], /refresh_binning, /refresh_info
            
        endfor
      endfor
    endfor
    
  endif
  
  cut = 'beardsley_thesis_list'
  win = 'bh_'
  spec_window = 'swbh_'
  
  for ch_i=0, n_ch-1 do begin
    for wedge_i=0,n_wedge-1 do begin
      for kperp_i=0,n_kperp-1 do begin
      
        basefile = ['/nfs/mwa-01/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal_pskspan100/ps_bh/Combined_obs_'+cut+'_cubeXX__even_odd_joint_maxuv100_'+win+'ch'+$
          strtrim((*ch_bands[ch_i])[0],2)+'-'+strtrim((*ch_bands[ch_i])[1],2)]
        endfile = '_averemove_'+spec_window+'dencorr_no_'+strtrim(wedge_angles[wedge_i],2)+'deg_wedge_kperplambda'+$)
          strtrim((*kperp_range_ptr[kperp_i])[0],2)+'-'+strtrim((*kperp_range_ptr[kperp_i])[1],2)+'_1dkpower.idlsave'
          
        limit_percent=0.9772 ;2 sigma
        hubble_param=0.71
        
        cubes=['res']
        pols=['yy']
        
        file=basefile+'_'+cubes+'_'+pols+endfile
        restore,file
        
        n_k=n_elements(k_edges)
        k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
        delta=power*(k^3.)/(2.*!pi^2.)
        dsigma=(k^3.)/(2.*!pi^2.)*weight_invert(sqrt(weights))
        zeros = where(dsigma EQ 0, n_count)
        if n_count GT 0 then dsigma[zeros] = !VALUES.F_INFINITY
        limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
        lim=min(limits,ind)
        header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('$
          +'ch_bands: '+strtrim((*ch_bands[ch_i])[0],2)+'-'+strtrim((*ch_bands[ch_i])[1],2)+' wedge: '+$
          strtrim(wedge_angles[wedge_i],2)+' kperp: '+$
          strtrim((*kperp_range_ptr[kperp_i])[0],2)+'-'+strtrim((*kperp_range_ptr[kperp_i])[1],2)+')'
        print, header
        
      endfor
    endfor
  endfor
  
end