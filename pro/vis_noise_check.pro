pro vis_noise_check

  vis_arr_XX = getvar_savefile('/fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/vis_data/1061316296_vis_model_XX.sav', 'vis_model_ptr')
  vis_arr_YY = getvar_savefile('/fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/vis_data/1061316296_vis_model_YY.sav', 'vis_model_ptr')
  vis_arr = PTRARR(2,/allocate)
  vis_arr = [Pointer_copy(vis_arr_XX),Pointer_copy(vis_arr_YY)]
  undefine_fhd, vis_arr_XX, vis_arr_YY
  obs = getvar_savefile('/fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/vis_data/1061316296_vis_XX.sav', 'obs')
  restore, '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/vis_data/1061316296_flags.sav'
  restore, '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_BH2grid_BH2degrid_GLEAMtenth_Z/metadata/1061316296_params.sav'

;  short_inds = where( sqrt(3.0E16*(params.uu^2. + params.vv^2)) LT 83.0, n_count)
 ; long_inds = where( sqrt(3.0E16*(params.uu^2. + params.vv^2)) GT 83.0, n_count)
;print, n_count
  t1 = SYSTIME(/SECONDS)
  vis_weights_use=split_vis_weights(obs,vis_weights,bi_use=bi_use,/preserve_weights,_Extra=extra)


  split_time =  SYSTIME(/SECONDS) - t1
  n_pol = obs.n_pol
  n_freq = obs.n_freq
  n_baselines = (size((*vis_arr[0])))[2]
  
  noise_arr_real = FLTARR(n_pol,n_freq)
  noise_arr_imag = FLTARR(n_pol,n_freq)
  noise_arr_abs = FLTARR(n_pol,n_freq)
  noise_arr_full = FLTARR(n_pol,n_freq)

  inds_imag = FLTARR(n_pol,n_freq)
  inds_real = FLTARR(n_pol,n_freq)
  inds_abs = FLTARR(n_pol,n_freq)
  inds_full = FLTARR(n_pol,n_freq)

if keyword_set(long_inds) then for pol_i=0, n_pol-1 do (*vis_weights_use[pol_i])[*,long_inds]=0

;  t1 = SYSTIME(/SECONDS)
;
;  FOR pol_i=0,n_pol-1 DO BEGIN
;    data_diff =REAL_part( (*vis_arr[pol_i])[*,*bi_use[0]])-REAL_PART((*vis_arr[pol_i])[*,*bi_use[1]]) ; only use imaginary part
;    vis_weight_diff = ((*vis_weights_use[pol_i])[*,*bi_use[0]]>0)*((*vis_weights_use[pol_i])[*,*bi_use[1]]>0)
;    FOR fi=0L,n_freq-1 DO BEGIN
;      ind_use=where(vis_weight_diff[fi,*],n_use)
;      IF n_use GT 1 THEN noise_arr_orig[pol_i,fi]=Stddev(data_diff[fi,ind_use])
;    ENDFOR
;  ENDFOR
  FOR pol_i=0,n_pol-1 DO BEGIN
    data_diff =(real_part( (*vis_arr[pol_i])[*,*bi_use[0]])-real_part((*vis_arr[pol_i])[*,*bi_use[1]]))^2.
    vis_weight_diff = ((*vis_weights_use[pol_i])[*,*bi_use[0]]>0)*((*vis_weights_use[pol_i])[*,*bi_use[1]]>0)
    FOR fi=0L,n_freq-1 DO BEGIN
      ind_use=where(vis_weight_diff[fi,*],n_use)
      inds_real[pol_i,fi] = n_use
      IF n_use GT 1 THEN noise_arr_real[pol_i,fi]=sqrt( (1./(N_elements(ind_use)-1.)) * total(data_diff[fi,ind_use]) )
      ;IF n_use GT 1 THEN noise_arr_intended[pol_i,fi]=sqrt( total(data_diff[fi,ind_use]) )
    ENDFOR
  ENDFOR
 
;  t_orig = SYSTIME(/SECONDS) - t1
  t1 = SYSTIME(/SECONDS)

  FOR pol_i=0,n_pol-1 DO BEGIN
    data_diff =(Imaginary( (*vis_arr[pol_i])[*,*bi_use[0]])-Imaginary((*vis_arr[pol_i])[*,*bi_use[1]]))
    vis_weight_diff = ((*vis_weights_use[pol_i])[*,*bi_use[0]]>0)*((*vis_weights_use[pol_i])[*,*bi_use[1]]>0)
    FOR fi=0L,n_freq-1 DO BEGIN
      ind_use=where(vis_weight_diff[fi,*],n_use)
      inds_imag[pol_i,fi] = n_use
      IF n_use GT 1 THEN noise_arr_imag[pol_i,fi]=sqrt( (1./(N_elements(ind_use)-1.)) * total(data_diff[fi,ind_use]^2.) )
        if (fi EQ 5) AND (pol_i EQ 0) then vis_slice_imag = data_diff[fi,ind_use]
      ;IF n_use GT 1 THEN noise_arr_intended[pol_i,fi]=sqrt( total(data_diff[fi,ind_use]) )
    ENDFOR
  ENDFOR

  t_intended = SYSTIME(/SECONDS) - t1
  t1 = SYSTIME(/SECONDS)

  FOR pol_i=0,n_pol-1 DO BEGIN
    data_diff =(abs( (*vis_arr[pol_i])[*,*bi_use[0]]) - abs((*vis_arr[pol_i])[*,*bi_use[1]]))^2.
    vis_weight_diff = ((*vis_weights_use[pol_i])[*,*bi_use[0]]>0)*((*vis_weights_use[pol_i])[*,*bi_use[1]]>0)
    FOR fi=0L,n_freq-1 DO BEGIN
      ind_use=where(vis_weight_diff[fi,*],n_use)
      inds_abs[pol_i,fi] = n_use
      IF n_use GT 1 THEN noise_arr_abs[pol_i,fi]=sqrt( (1./(N_elements(ind_use)-1.)) * total(data_diff[fi,ind_use]) )
      ;IF n_use GT 1 THEN noise_arr_abs[pol_i,fi]=sqrt( total(data_diff[fi,ind_use]) )
    ENDFOR
  ENDFOR

  t_abs = SYSTIME(/SECONDS) - t1

  print, "Number of defined baselines in the even/odd group"
  ind_use = where(vis_weight_diff,n_use)
  print, n_use
  undefine_fhd, vis_weight_diff, vis_weights, data_diff

  t1 = SYSTIME(/SECONDS)

  FOR pol_i=0,n_pol-1 DO BEGIN
      ;visibilities which are only defined for even-odd set
      vis_weight_full = fltarr(n_freq,n_baselines)

      vis_weight_full[*,*bi_use[0]] = ((*vis_weights_use[pol_i])[*,*bi_use[0]]>0)*((*vis_weights_use[pol_i])[*,*bi_use[1]]>0)
      vis_weight_full[*,*bi_use[1]] = ((*vis_weights_use[pol_i])[*,*bi_use[0]]>0)*((*vis_weights_use[pol_i])[*,*bi_use[1]]>0)
     

      ;ind_use = uniq(where(vis_weight_full GT 0)/n_freq)
      ;mean_vis = mean(abs(((*vis_arr[pol_i])[*,ind_use])),dim=2)
    FOR fi=0L,n_freq-1 DO BEGIN
      ind_use = where(vis_weight_full[fi,*] GT 0,n_use)
      ind_use_even = where(vis_weight_full[fi,*bi_use[0]] GT 0,n_use)
      ind_use_odd = where(vis_weight_full[fi,*bi_use[1]] GT 0,n_use)
      if n_use GT 0 then begin
	mean_vis = mean(((*vis_arr[pol_i])[fi,ind_use]))
        noise_arr_full[pol_i,fi]=sqrt( (1./(N_elements(ind_use)-1.)) * (total((real_part((*vis_arr[pol_i])[fi,ind_use]) - real_part(mean_vis))^2.) + $
          total((imaginary((*vis_arr[pol_i])[fi,ind_use]) - imaginary(mean_vis))^2.)) )
        if (fi EQ 5) AND (pol_i EQ 0) then vis_slice = (*vis_arr[pol_i])[fi,ind_use]
        if (fi EQ 5) AND (pol_i EQ 0) then vis_slice_even = (*vis_arr[pol_i])[fi,(*bi_use[0])[ind_use_even]]
        if (fi EQ 5) AND (pol_i EQ 0) then vis_slice_odd = (*vis_arr[pol_i])[fi,(*bi_use[1])[ind_use_odd]]
      endif
      inds_full[pol_i,fi] = n_use
    ENDFOR
  ENDFOR

;cgScatter2D, real_part(vis_slice), imaginary(vis_slice), XTitle='Real', YTitle='Imaginary', title='Distribution of visibilities in a frequency slice'

xrange = [-150,150];[Min(real_part(vis_slice)), Max(real_part(vis_slice))]
   yrange = [-150,150];[Min(imaginary(vis_slice)), Max(imaginary(vis_slice))]
   xbinsize = 5.
   ybinsize = 5.
   density = Hist_2D(real_part(vis_slice), imaginary(vis_slice), Min1=xrange[0], Max1=xrange[1], Bin1=xbinsize, $
                           Min2=yrange[0], Max2=yrange[1], Bin2=ybinsize)
   density_even = Hist_2D(real_part(vis_slice_even), imaginary(vis_slice_even), Min1=xrange[0], Max1=xrange[1], Bin1=xbinsize, $
                           Min2=yrange[0], Max2=yrange[1], Bin2=ybinsize)
   density_odd = Hist_2D(real_part(vis_slice_odd), imaginary(vis_slice_odd), Min1=xrange[0], Max1=xrange[1], Bin1=xbinsize, $
                           Min2=yrange[0], Max2=yrange[1], Bin2=ybinsize)
   density_emo = Hist_2D(real_part(vis_slice_even - vis_slice_odd), imaginary(vis_slice_even - vis_slice_odd), Min1=xrange[0], Max1=xrange[1], Bin1=xbinsize, $
                           Min2=yrange[0], Max2=yrange[1], Bin2=ybinsize)


maxDensity = Ceil(Max(density)/1e2) * 1e2
   scaledDensity = BytScl(density, Min=0, Max=maxDensity)
cgLoadCT, 33
   TVLCT, cgColor('gray', /Triple), 0
   TVLCT, r, g, b, /Get
   palette = [ [r], [g], [b] ]

cgPS_Open,'contours.png',/quiet,/nomatch

   cgImage, scaledDensity, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
      XTitle='Real', YTitle='Imaginary',charsize=.9

   thick = (!D.Name EQ 'PS') ? 6 : 2
   cgContour, density, LEVELS=maxDensity*[0.25, 0.5, 0.75], /OnImage, $
       C_Colors=['Tan','Tan', 'Brown'], C_Annotation=['Low', 'Avg', 'High'], $
       C_Thick=thick, C_CharThick=thick
      
   ; Display a color bar.
   cgColorbar, Position=[0.125, 0.875, 0.9, 0.925], Title='Density', $
       Range=[0, maxDensity], NColors=254, Bottom=1, OOB_Low='gray', $
       TLocation='Top',charsize=.9
cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent


maxDensity = Ceil(Max(density_even)/1e2) * 1e2
   scaledDensity_even = BytScl(density_even, Min=0, Max=maxDensity)
cgPS_Open,'contours_even.png',/quiet,/nomatch

   cgImage, scaledDensity_even, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
      XTitle='Real', YTitle='Imaginary',charsize=.9

   thick = (!D.Name EQ 'PS') ? 6 : 2
   cgContour, density_even, LEVELS=maxDensity*[0.25, 0.5, 0.75], /OnImage, $
       C_Colors=['Tan','Tan', 'Brown'], C_Annotation=['Low', 'Avg', 'High'], $
       C_Thick=thick, C_CharThick=thick

   ; Display a color bar.
   cgColorbar, Position=[0.125, 0.875, 0.9, 0.925], Title='Density', $
       Range=[0, maxDensity], NColors=254, Bottom=1, OOB_Low='gray', $
       TLocation='Top',charsize=.9
cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent

   scaledDensity_odd = BytScl(density_odd, Min=0, Max=maxDensity)
cgPS_Open,'contours_odd.png',/quiet,/nomatch

   cgImage, scaledDensity_odd, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
      XTitle='Real', YTitle='Imaginary',charsize=.9

   thick = (!D.Name EQ 'PS') ? 6 : 2
   cgContour, density_odd, LEVELS=maxDensity*[0.25, 0.5, 0.75], /OnImage, $
       C_Colors=['Tan','Tan', 'Brown'], C_Annotation=['Low', 'Avg', 'High'], $
       C_Thick=thick, C_CharThick=thick

   ; Display a color bar.
   cgColorbar, Position=[0.125, 0.875, 0.9, 0.925], Title='Density', $
       Range=[0, maxDensity], NColors=254, Bottom=1, OOB_Low='gray', $
       TLocation='Top',charsize=.9
cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent

   scaledDensity_emo = BytScl(density_emo, Min=0, Max=maxDensity)
cgPS_Open,'contours_emo.png',/quiet,/nomatch

   cgImage, scaledDensity_emo, XRange=xrange, YRange=yrange, /Axes, Palette=palette, $
      XTitle='Real', YTitle='Imaginary',charsize=.9

   thick = (!D.Name EQ 'PS') ? 6 : 2
   cgContour, density_emo, LEVELS=maxDensity*[0.25, 0.5, 0.75], /OnImage, $
       C_Colors=['Tan','Tan', 'Brown'], C_Annotation=['Low', 'Avg', 'High'], $
       C_Thick=thick, C_CharThick=thick

   ; Display a color bar.
   cgColorbar, Position=[0.125, 0.875, 0.9, 0.925], Title='Density', $
       Range=[0, maxDensity], NColors=254, Bottom=1, OOB_Low='gray', $
       TLocation='Top',charsize=.9
cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent

;hist = histogram(vis_slice_imag,binsize=xbinsize,max=150,min=-150,/NAN)
;cgplot, FINDGEN(61)*xbinsize + xrange[0], hist

  t_full = SYSTIME(/SECONDS) - t1



  print, "intended, abs, full"
  print, t_intended+split_time, t_abs+split_time, t_full+split_time

inds_zero = where(noise_arr_imag EQ 0)
noise_arr_real[inds_zero]= !VALUES.F_NAN
noise_arr_full[inds_zero]= !VALUES.F_NAN
noise_arr_imag[inds_zero]= !VALUES.F_NAN

stop

cgPS_Open,'standard_dev.png',/quiet,/nomatch
cgplot, (*obs.baseline_info).freq/1e6, noise_arr_imag[0,*], xrange=[165,200],xtitle='Frequency (MHz)',ytitle='Standard deviation of one polarization (Jy)', yrange=[20,60], charsize=.9
cgoplot, (*obs.baseline_info).freq/1e6, sqrt(noise_arr_imag[0,*]^2. + noise_arr_real[0,*]^2.)/sqrt(2), color='blue'     
cgoplot, (*obs.baseline_info).freq/1e6, noise_arr_full[0,*], color='green' 
al_legend, ['Equation 4','Equation 3','Equation 2'],color=['black','blue','green'], pos=[.5,.87], charsize=.9, box=0, /normal,linestyle=0
cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent

  stop

end
