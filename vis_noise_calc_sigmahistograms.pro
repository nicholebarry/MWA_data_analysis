PRO vis_noise_calc_sigmahistograms,restore_all=restore_all,interleave=interleave
  ; Script to plot visibility sigmas given various time interleaving or 2 second time visibilities differenced at various time steps

  ;****Restore the necessary information from the standard run to run this script outside of FHD.
  if keyword_set(restore_all) then begin
    RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316296_obs.sav' ;restore obs structure
    vis_XX_restored = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/vis_data/1061316296_vis_XX.sav', 'vis_ptr') ;restore array of calibrated visibilities
    vis_YY_restored = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/vis_data/1061316296_vis_YY.sav', 'vis_ptr')
    
    ;Combine the calibrated visibilities in the correct format for the script
    vis_arr = [vis_XX_restored, vis_YY_restored]
    
    RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/vis_data/1061316296_flags.sav'
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/calibration/1061316296_cal.sav'
    
  endif
  ;****End of restore the necessary information from the standard run to run this script outside of FHD.
  
  
  ;****Setup
  n_pol=obs.n_pol
  n_freq=Long(obs.n_freq)
  noise_arr=fltarr(n_pol,n_freq)
  pol_i=0
  
  sigma_arr=FLTARR(2,6)
  err_arr=FLTARR(2,6)
  ;****End of setup
  
  
  ;*************Begin looping through specified frequency indices, since sigma calculations depend on location within the course band
  for freq_i=0,1 do begin
  
    ;Index 5 is 1/3rd of the way through the first course band
    if freq_i EQ 0 then freq=5
    ;Index 14 is the edge of the first course band
    if freq_i EQ 1 then freq=14
    
    ;*****Begin looping through time steps, which will determine how much time is taken between visibilty differences,
    ;or how many 2 second intervals will be summed.
    for time_i=1, 6 do begin
    
      ;Use a modified script to get the time bins to use for the visibility differences. The return value is the
      ;visibilility bins the same size as the visibility array. Input time determines which bins to pick out,
      ;so time 1 gives even odd differences, time 3 groups bins as 0,6,12... and 3,9,15...
      bin_i=split_vis_flags_vissigma_timestep(obs,flag_arr,bi_use=bi_use,/preserve_flags,time=time_i)
      
      ;Setting this keyword sums the visibilities. For time 3, this sums group 1 as 0+1+2, 6+7+8... and
      ;group 2 as 3+4+5, 9+10+11...
      if keyword_set(interleave) then begin
      
        ;Take the first bin and save it for adding upon in the next for loop
        ;Use only the imaginary part in the calculation of noise, since the real and the imaginary parts
        ;result in the same sigma to about the one hundredth decimal, and a normalization will be applied
        ;later
        base0=Imaginary( (*vis_arr[pol_i])[*,*bi_use[0]])
        base1=Imaginary( (*vis_arr[pol_i])[*,*bi_use[1]])
        
        ;***Loop through each time bin and add it to the base visibility
        For ti=1, time_i-1 do begin
        
          ;find the bin associated with the next time step to sum, given through a modulo
          *bi_use[0]=where(bin_i mod (2*time_i) EQ ti,n_even)
          *bi_use[1]=where(bin_i mod (2*time_i) EQ (time_i+ti),n_odd)
          
          ;Make sure there are no leftover bins, and leave out those bins if so
          IF n_even LT n_odd THEN *bi_use[1]=(*bi_use[1])[0:n_even-1]
          IF n_odd LT n_even THEN *bi_use[0]=(*bi_use[0])[0:n_odd-1]
          
          ;Sum the current time step visibility to the base
          base0=base0+Imaginary( (*vis_arr[pol_i])[*,(*bi_use[0])] )
          base1=base1+Imaginary( (*vis_arr[pol_i])[*,(*bi_use[1])] )
          
        endfor
        ;***End of loop through each time bin and add it to the base visibility
        
        ;Calculate the difference of the summed visibilities, and normalize them
        data_diff_pre=(base0-base1)/sqrt(time_i)
        
      endif else data_diff_pre =Imaginary( (*vis_arr[pol_i])[*,*bi_use[0]])-Imaginary((*vis_arr[pol_i])[*,*bi_use[1]])
      
      ;Take the flags from just the first groups bin. This is an assumption, and the flags should be taken
      ;from both groups and from all the interleaves.
      flag_diff = ((*flag_arr[pol_i])[*,*bi_use[0]]>0)
      
      ;calculate the baseline length for each visibilty (in light travel time at this stage)
      uv=sqrt(cal.UU^2.+cal.VV^2.)
      ;Narrow down the baseline length to the correct dimension. Group 1 bins can represent both groups, since the baseline length
      ;is the same
      uv=uv[*bi_use[0]]
      
      ;The binsize is arbitrary, but it does give a pretty plot. Could be increased. Must be even.
      binsize=10.
      ;Find indices of unflagged visibilities.
      ind_use=where(flag_diff[freq,*] GT 0,n_use)
      ;Apply indices from the unflagged visibilities to both the visibility difference and the corresponding baseline
      IF n_use GT 1 THEN data_diff=reform(data_diff_pre[freq,ind_use])
      uv=uv[ind_use]
      
      ;Histogram the uv baselines to the binsize provided above, and record that in the reverse indices for applying
      ;to the visibility difference later. Convert the uv from light speed time to wavelengths.
      result=histogram(uv*3.*10^8.,binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
      
      ;Calculate the visibility sigma (normalized due to real/imag) for each bin using reverse indices
      vis_sigma=FLTARR(N_elements(result))
      for i=0, N_elements(result)-1 do if result[i] GT 0 then vis_sigma[i]=stddev(data_diff[ri[ri[i]:ri[i+1]-1]])/sqrt(2.)
      
      ;Make the x input and y input pretty for the plot
      y_arr=[vis_sigma[0], vis_sigma, vis_sigma[N_elements(vis_sigma)-1]]
      x_arr=[locations[0],locations+binsize/2,omax]
      
      ;Make
      cgPS_Open,'/nfs/eor-00/h1/nbarry/freq'+strtrim(freq,2)+'vissigma_inter_time'+strtrim(time_i,2)+'.png',/quiet,/nomatch
      cgplot, x_arr, y_arr, psym=10, Ystyle=8, xrange=[0,2000], title='Frequency index '+strtrim(freq,2)+', visibility sigma vs wavelength, interleave '+strtrim(time_i,2), xtitle='baseline ($\tex\lambda$)', ytitle='visibility sigma', charsize=1, position=[.15,.15,.85,.85]
      cgaxis, yaxis=1, yrange=[0,9000], /save, title='vis # in bin', charsize=1, color='blue'
      cgtext, .65,.75,' $\tex\sigma$!I'+strtrim(freq,2)+'!N='+strtrim(stddev(data_diff[*])/sqrt(2.),2), /normal, charsize=1
      cgoplot, x_arr,result, color='blue'
      cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
      sigma_arr[freq_i,time_i-1]=stddev(data_diff[*])/sqrt(2.)
      err_arr[freq_i,time_i-1]=N_elements(data_diff)
      
    endfor
  ;*****End of looping through time steps, which will determine how much time is taken between visibilty differences,
  ;or how many 2 second intervals will be summed.
    
  endfor
  ;*************End of looping through specified frequency indices, since sigma calculations depend on location within the course band
  
  ;****Begin making overall vis sigma with various time interleaving/time steps calculated above
  cgPS_Open,'/nfs/eor-00/h1/nbarry/vissigma_deltat_err.png',/quiet,/nomatch
  
  If ~keyword_set(interleave) then cgplot, [1,2,3,4,5,6],sigma_arr[0,*],psym=4, yrange=[20,40],xrange=[0,7],xtitle='$\tex\Delta$t steps', $
    ytitle='visibility sigma', title='Visibility sigma for various $\tex\Delta$t',charsize=1,color='red', $
    ERR_YLow=sigma_arr[0,*]/sqrt(err_arr[0,*]),ERR_YHigh=sigma_arr[0,*]/sqrt(err_arr[0,*]), ERR_color='red' else $
    cgplot, [1,2,3,4,5,6],sigma_arr[0,*],psym=4, yrange=[20,40],xrange=[0,7],xtitle='$\tex\Delta$t steps', $
    ytitle='visibility sigma', title='Visibility sigma for various $\tex\Delta$t',charsize=1,color='red'
    
  If ~keyword_set(interleave) then cgoplot, [1,2,3,4,5,6],sigma_arr[1,*],psym=5, color='blue',ERR_YLow=sigma_arr[0,*]/sqrt(err_arr[0,*]),ERR_YHigh=sigma_arr[0,*]/sqrt(err_arr[0,*]), ERR_color='blue' else $
    cgoplot, [1,2,3,4,5,6],sigma_arr[1,*],psym=5, color='blue'
    
  cglegend, title=['Index 5','Index 14'], color=['red','blue'],psym=[4,5], length=0,location=[.75,.25],charsize=1
  
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  ;****End making overall vis sigma with various time interleaving/time steps calculated above
  
END