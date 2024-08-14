pro delay_spectrum, dir, obsid=obsid,spec_window_type=spec_window_type, plotfile=plotfile
  ;;;; This is a script to generate delay spectra from visibilities

  ;Make delay spectra plots. Can turn off to run without generating images
  png = 1 
  
  ;;; Defaults
  ; dir is the input directory where the visibilities are, and spec_window_type is a window function for the spectral dimension
  if ~keyword_set(dir) then dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_redo_redo_CasA_N13_rescaled/'
  if ~keyword_set(spec_window_type) then spec_window_type = 'Blackman-Harris7'

  ;;; Parameters for the observation id files
  ; The observation id files are separated by day and pointing. 
  dir_obsids = '/fred/oz048/MWA/data/2014/van_vleck_corrected/coarse_corr_no_ao/obs_lists_extra/'
  day = ['10','11','13','14','15','16','18','19','21','22','23','24','25','26','28','29','30','31','32','33']
  n_days = N_elements(day)
  pointingnames=['_minustwo_ssins','_minusone_ssins','_zenith_ssins','_plusone_ssins','_plustwo_ssins']
  n_pointings = N_elements(pointingnames)
  obs_ptr=PTRARR(n_days,n_pointings,/allocate_heap)

  ; Get observation ids and put them in a pointing pointer
  max_n_obs=0
  FOR poi_i=0,n_pointings-1 DO BEGIN
    max_n=0
    FOR day_i=0,n_days-1 DO BEGIN
      undefine, obs_temp
      readcol, dir_obsids + day[day_i] + pointingnames[poi_i]+'.txt', obs_temp, format='A', /silent
      if N_elements(obs_temp) GT 0 then *obs_ptr[day_i,poi_i]=obs_temp
      max_n_obs = max([N_elements(obs_temp),max_n_obs])
    ENDFOR
  ENDFOR

  ;;; Constants 
  speed_of_light = 299792458. ;m/s
  i_comp = Complex(0,1)
  z0_freq = 1420.4057517667e6 ;; Hz

  ;;; Get a representative frequency array for all data to calculate cosmology params 
  ;;; and conversion factors outside of the main loops
  restore, dir+'/metadata/'+(*obs_ptr[0,0])[0]+'_obs.sav'
  ; frequency array
  freq_arr = (*obs.baseline_info).freq
  ; number of frequencies
  nfreq = n_elements(freq_arr)
  ; resolution of frequencies
  freq_res = obs.freq_res
  ; unflagged indices in frequency array
  freq_inds = where((*obs.baseline_info).freq_use)
  ; redshifts
  redshifts = z0_freq/freq_arr - 1d

  ;;; Calculate comoving line-of-sight distances to be able to get the mean comoving distance in Mpc
  ; Scripts are from eppsilon
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los, hubble_param = hubble_param,$
    Ez=Ez, wedge_factor=wedge_factor
  ; freq_arr input into z_mpc should be in MHz
  z_mpc_mean = z_mpc(freq_arr / 1e6, hubble_param = hubble_param, f_delta = f_delta, $
    comov_dist_los = comov_dist_los, z_mpc_delta = z_mpc_delta)


  ;;; Calculate conversion factors
  ; conversion from cMpc to lambda for kperp
  kperp_lambda_conv = z_mpc_mean / (2.*!pi) 
  ; conversion from Jy to mK sterradian = 10^-26 * c^2 *10^3 / (2*f^2*kb)
  Jy_to_mK_str_conv = speed_of_light^2. / (2. * (freq_arr)^2. * 1.38065) 
  ; conversion from sterradian to Mpc^2
  str_to_Mpc2_conv = comov_dist_los^2.

  ;;; Calculate window function for freq fft
  ; Script from fhdps_utils
  window = spectral_window(nfreq, type=spec_window_type,/periodic)
  norm_factor = sqrt(nfreq/total(window^2.))
  window = window * norm_factor

  for day_i=0, n_days-1 do begin

    ; Array to hold the final delay power spectrum for each observation
    full_delay_arr = FLTARR(101,96,2, 1, n_pointings, max_n_obs)
    for poi_i=0, n_pointings-1 do begin

      obs_ids = *obs_ptr[day_i,poi_i]
      n_obs = N_elements(obs_ids)

      for obs_i=0,n_obs-1 do begin
        obsid = obs_ids[obs_i]

        ;Restore weights and metadata, skip if unavailable
        if ~file_test(dir+'/metadata/'+obsid+'_params.sav') then continue
        restore, dir+'/metadata/'+obsid+'_params.sav'
        restore, dir+'/vis_data/'+obsid+'_flags.sav'
        restore, dir+'/metadata/'+obsid+'_obs.sav'

        ;Get metadata
        nbl = n_elements(params.uu)
        flags = fltarr(nfreq,nbl,2)

        ;Only keep unflagged baselines between pols, 
        ; discounting the expected flagging for the edge channels and the dc channels
        for poli=0,1 do flags[*,*,poli] = *vis_weights[poli]
        flag_test = Total(Total(flags>0,1),2)
        bi_use=where(flag_test ge (2*nfreq-48*2. - 24))
        nbl = n_elements(bi_use)
        uu = params.uu[bi_use]
        vv = params.vv[bi_use]
        ww = params.ww[bi_use]

        bb = sqrt(uu^2.+vv^2.+ww^2.)

        ;Restore vis and create residual visibilities
        ;Only subtract defined data with models in order to avoid over subtraction
        data = Complex(fltarr(nfreq,nbl,2)) ; Stack pols
        restore, dir+'/vis_data/'+obsid+'_vis_XX.sav'
        restore, dir+'/vis_data/'+obsid+'_vis_model_XX.sav'
        data[*,*,0] = (*vis_ptr)[*,bi_use] 
        ;data[*,*,0] = (*vis_model_ptr)[*,bi_use] 
        temp = (*vis_model_ptr)[*,bi_use]
        temp = temp[freq_inds,*]
        data[freq_inds,*,0] = data[freq_inds,*,0] - temp
        restore, dir+'/vis_data/'+obsid+'_vis_YY.sav'
        restore, dir+'/vis_data/'+obsid+'_vis_model_YY.sav'
        data[*,*,1] = (*vis_ptr)[*,bi_use] 
        ;data[*,*,1] = (*vis_model_ptr)[*,bi_use]
        temp = (*vis_model_ptr)[*,bi_use]
        temp = temp[freq_inds,*]
        data[freq_inds,*,1] = data[freq_inds,*,1] - temp

        ;Could potentially delay filter in the future
        ;dataptr = ptrarr(2)
        ;vis_delay_filter, dataptr, vis_weights, params, obs
        ;data[*,*,0] = *dataptr[0]
        ;data[*,*,1] = *dataptr[1]

        undefine_fhd, vis_ptr, vis_model_ptr, flags, bi_use

        ; Phase to zenith to add coherently 
        w_mat = freq_arr#ww ; This should now be in wavelengths
        for poli=0,1 do data[*,*,poli] *= exp(i_comp * 2. * !pi * w_mat)

        ;Calculate final conv factor, where 2 is to get to Stokes I from xx/yy pols
        conv_factor = 2. * rebin(reform(Jy_to_mK_str_conv * str_to_Mpc2_conv, nfreq,1,1),nfreq,nbl,2,/sample) 
        
        ;Get data in proper units for fft, where dz and nfreq are fft factors
        dz = (max(comov_dist_los)-min(comov_dist_los))/(nfreq-1)
        data = conv_factor * dz * nfreq * data ; dz and nfreq to prepare for fft

        ;Aeff = 21. ; m^2 from Aaron E-W memo
        ;window_int = (z_mpc_mean * speed_of_light/mean(freq_arr))^2. / Aeff * (max(comov_dist_los)-min(comov_dist_los))
        ;window_int = window_int / 4. ; Checking back-of-the-envelope against eppsilon calcuation, I was off by ~4x
         ;data = data / sqrt(window_int) ; Should divide by window int after fft, but doing all units here instead.
        
        ave_beam_int = DBLARR(nfreq,2)
        for pol_i=0,1 do ave_beam_int[*,pol_i] = (*obs.primary_beam_sq_area[pol_i]) ; * obs.nf_vis,/nan) / total(obs.nf_vis,/nan)
        ;; convert rad -> Mpc^2, multiply by depth in Mpc to get Mpc^3
        window_int_beam_obs = 2. * ave_beam_int * z_mpc_mean^2. * (z_mpc_delta * nfreq)
        sqrt_window_int_beam_obs_expand = rebin(reform(sqrt(window_int_beam_obs),nfreq,1,2), nfreq, nbl, 2, /sample)
        data = data / temporary(sqrt_window_int_beam_obs_expand) ; Should divide by window int after fft, but doing all units here instead.

        ;ave removal
        data_mean = mean(data, dimension=1)
        data = data - (rebin(reform(real_part(data_mean),1,nbl,2), nfreq, nbl,2, /sample) + $
          dcomplex(0,1)*rebin(reform(imaginary(data_mean),1,nbl,2), nfreq, nbl,2, /sample))

        ; Apply spectral window function
        window_expand = rebin(reform(window,nfreq,1,1), nfreq, nbl, 2, /sample)
        data = data * window_expand
        data_mean = data_mean / norm_factor

        ;Average in frequency to reduce contamination from flagged channels
        n_avg = 2
        nfreq_avg = nfreq / n_avg
        freq_arr_avg = DBLARR(nfreq_avg)
        shift_val = nfreq_avg / 2 
        data_avg = dcomplex(DBLARR(nfreq_avg, nbl, 2))

        FOR freq_i=0L,nfreq_avg-1 DO data_avg[freq_i,*,*]=Total(data[freq_i*n_avg:(freq_i+1)*n_avg-1,*,*],1)
        FOR freq_i=0L,nfreq_avg-1 DO freq_arr_avg[freq_i]=Mean(freq_arr[freq_i*n_avg:(freq_i+1)*n_avg-1])

        spectra = (shift(fft(data_avg,dim=1),shift_val,0,0)) ; Shift only in fft direction.
        spectra[nfreq_avg/2,*,*] += (data_mean)
        spectra = abs(spectra)^2.

        undefine_fhd,data, data_avg
        ; fold over
        ndelay = shift_val
        spectra[(ndelay+1):*,*,*] += spectra[(ndelay-1):1:-1,*,*]
        spectra = spectra[ndelay:*,*,*]


      

          ; Bin up
          umag = sqrt(abs(uu)^2 + abs(vv)^2) * mean(freq_arr_avg)
          umin = min(umag)
          umax = max(umag)
          nbins = 100
          ubin = (umax-umin)/nbins
          uhist = histogram(umag, binsize=ubin, min=umin, omax=umax, locations=u_locs, reverse_indices=u_ri)
          u_centers = u_locs + ubin/2d
          u_edges = [u_locs, max(u_locs) + ubin]
          nbins = n_elements(uhist)
          kperp_edges = u_edges / kperp_lambda_conv
          delay_delta = 1./(max(freq_arr_avg)-min(freq_arr_avg))
          delay_max = (ndelay-1) * delay_delta
          ; Note extra 1e3 for hubble constant to get it in m/s/Mpc
          kpar_delta = delay_delta * 2*!pi * 100. * hubble_param * 1e3 * z0_freq * mean(Ez) / (speed_of_light * (1+mean(redshifts))^2.)
          kpar_edges = kpar_delta*findgen(ndelay+1)-kpar_delta/2.
          ; Make 2D image
          delay2d = fltarr(nbins,ndelay,2)
          for i=0,nbins-1 do begin
            if uhist[i] gt 0 then begin
              delay2d[i,*,*] = mean(spectra[*,u_ri[u_ri[i]:u_ri[i+1]-1],*],dim=2)
            endif
          endfor

          full_delay_arr[*,*,*, day_i, poi_i, obs_i] = delay2d

        if keyword_set(png) then begin

          wedge_amp = mean(wedge_factor) * !dpi / 180d * [20d, 90d]
          plotfile =  dir+'delay_plots/' + strtrim(obsid,2) + '_xx'
           kpower_2d_plots, power=delay2d[*,*,0], kperp_edges=kperp_edges, kpar_edges=kpar_edges, $
                        kperp_lambda_conv=kperp_lambda_conv, delay_params=1e9*[delay_delta, delay_max], $
                        hubble_param=hubble_param, plotfile=plotfile, /png, /hinv, /delay_axis, $
                        /baseline_axis, /plot_wedge_line, wedge_amp=wedge_amp, $
                        kperp_plot_range=[.007, .3], kpar_plot_range=[.003,4], full_title=strtrim(obsid,2) + ' xx', data_range=[6e12,3e13]

           plotfile =  dir+'delay_plots/' + strtrim(obsid,2) + '_yy'
           kpower_2d_plots, power=delay2d[*,*,1], kperp_edges=kperp_edges, kpar_edges=kpar_edges, $
                        kperp_lambda_conv=kperp_lambda_conv, delay_params=1e9*[delay_delta, delay_max], $
                        hubble_param=hubble_param, plotfile=plotfile, /png, /hinv, /delay_axis, $
                        /baseline_axis, /plot_wedge_line, wedge_amp=wedge_amp, $
                        kperp_plot_range=[.007, .3], kpar_plot_range=[.003,4], full_title=strtrim(obsid,2) + ' yy', data_range=[6e12,3e13]

        endif

      endfor
    endfor
    save, full_delay_arr, filename = dir+'full_delay_arr_'+day[day_i]+'.sav', kperp_edges, kpar_edges, kperp_lambda_conv, delay_delta, delay_max


  endfor



end
