pro bubbles_brad_uvf, restore_image=restore_image
;; Code to read in theoretical EoR simulations in sky space (Mpc) over redshift and convert them into
;; a uvf cube for typical MWA Phase I high band frequencies and dimension, in units of Jy

  t0=Systime(1)
  start_mem = memory(/current)

  ;Read in a sample obs file to get the desired frequencies
  restore, '/fred/oz048/MWA/CODE/FHD/fhd_nb_cross_RTS/metadata/1061316544_obs.sav'
  obs_freq = (*obs.baseline_info).freq
  n_freq=N_elements(obs_freq)

  ;General cosmology
  h = 0.68  ; hubble parameter
  nu_emitted = 1420.4057517*10^6.  ; Hz
  speed_of_light = 299792458.  ; meters/sec
  kB = 1.38064852 ;Boltzmann's constant
  Jy_to_mKstr = speed_of_light^2/(2*obs_freq^2*kB)
  redshifts = ( nu_emitted - obs_freq ) / obs_freq
  ;Set h to 1 to not include little h in all the calculations
  cosmology_measures, redshifts, comoving_dist_los=comoving_dist_los, hubble_param=1.;hubble_param
  comoving_dist_los = float(comoving_dist_los) ;double accuracy not required

  redshift_paths = file_search('/fred/oz078/bgreig/ForCath/FirstBox/BrightnessTemperatureLCBox/*dat')
  redshift_files = file_basename(redshift_paths)
  available_redshifts = strmid(redshift_files,37,6) ;redshift string is always at the same location in the file name

  ;Read in all of the images and reform to 2D
  box_range_Mpc = 7500 ;Mpc
  box_dim = 6400
  box_res_Mpc = float(box_range_Mpc)/float(box_dim) ;Mpc

restore_image=1
  if ~keyword_set(restore_image) then begin
    bubble_image_arr = FLTARR(box_dim, box_dim, 420)
    for z_i=0, 419 do begin
      bubble_image_filename = '/fred/oz078/bgreig/ForCath/FirstBox/BrightnessTemperatureLCBox/BrightnessTemperatureSlice_L'$
        +strtrim(box_range_Mpc,2)+'Mpc_z'+available_redshifts[z_i]+'_SliceNumber_DIM'+strtrim(box_dim,2)+'_'+strtrim(z_i,2)+'.dat'
      bubble_image = read_binary(bubble_image_filename, DATA_TYPE=4) ;float32 data type
      bubble_image = reform(bubble_image,box_dim,box_dim)
      bubble_image_arr[*,*,z_i] = bubble_image
    endfor
    undefine, bubble_image
  endif
  undefine, redshift_files, redshift_paths

  ;Approximate images to the desired redshifts
  fractional_index = FLTARR(n_freq)
  for z_i=0, n_freq-1 do begin
    temp = min(abs(float(available_redshifts) - redshifts[z_i]),index) ;find the closest index

    ;find the next closest neighbor index for the interpolation
    diff = float(available_redshifts[index]) - redshifts[z_i]
    if diff GT 0 then neighbor_index = index - 1 else neighbor_index = index + 1

    ;Get the desired redshifts in terms of the fractional index of the available redshifts
    fractional_index[z_i] = abs(diff) / abs(float(available_redshifts[index]) - float(available_redshifts[neighbor_index]))
    if diff GT 0 then fractional_index[z_i] = index - fractional_index[z_i] else fractional_index[z_i] = index + fractional_index[z_i]
  endfor

  ;Image of 7500Mpcx7500Mpc EoR for every frequency in a high band observation with the MWA
  if keyword_set(restore_image) then restore, '/fred/oz048/MWA/CODE/FHD/Brad_sim/r_Mpc_cube.sav' else begin
    bubble_image_arr_interpolated = interpolate(temporary(bubble_image_arr),temporary(fractional_index))
    save, bubble_image_arr_interpolated, obs_freq, filename = '/fred/oz048/MWA/CODE/FHD/Brad_sim/r_Mpc_cube.sav'
  endelse

  print, "Image made for desired frequencies"
  run_report, start_mem, t0, silent=silent

  undefine, fractional_index, available_redshifts, redshifts

  ;;For the model uv plane, we do a 2048x2048 uv plane with .5 lambda spacing for the central frequency in general
  ;; However, that is essentially full sky (2*!pi rad) and the MWA FoV is only 20 deg to FWHM at 200MHz (.35 rad)
  ;; A 7500x7500Mpc image is ~.85rad, thus should fill out the primary beam without tiling.
  ;; Instead of tiling to fill out the full sky, we will pad the image to achieve the correct resolution, and only 
  ;; lose information outside of the primary beam. 
  
  ;Using Morales & Hewitt 2004 conversions to get desired grid on the sky in radians
  dimension=2048.
  uv_res_inverserad = 0.5 ; lambda (inverse rad)
  uv_range_inverserad = uv_res_inverserad*dimension
  ;theta_res_rad = (2.*!pi) / uv_range_inverserad
  ;theta_range_rad = (2.*!pi) / uv_res_inverserad
  theta_res_rad = 1 / uv_range_inverserad
  theta_range_rad = 1 / uv_res_inverserad

  ;The theory box coordinates
  box_vector_Mpc = (FINDGEN(box_dim) - box_dim/2)*box_res_Mpc
  ;box_grid_Mpc = box_vector_Mpc#box_vector_Mpc ;Mpc

  bubble_image_grid = FLTARR(dimension,dimension,n_freq)
  uvf_cube_mKMpc2 = complex(FLTARR(dimension,dimension,n_freq))
  r_res_Mpc_arr = FLTARR(n_freq)
  pix_num_filled = FLTARR(n_freq)
stop
  for z_i=0, n_freq-1 do begin
    ;
    ;Loop over every frequency to interpolate to the appropriate grid
    ;and take the spatial ffts

    ;Mpc required for sky grid is redshift dependent.
    r_range_Mpc = theta_range_rad * comoving_dist_los[z_i]
    r_res_Mpc = theta_res_rad * comoving_dist_los[z_i]
    r_res_Mpc_arr[z_i] = r_res_Mpc
    r_vector_Mpc = (FINDGEN(dimension) - dimension/2)*r_res_Mpc

    ;Run the interpolation on just the area to be filled by theory instead of the full padded area
    defined_index = where(abs(r_vector_Mpc) LE (box_range_Mpc/2),n_index)
    if n_index LT 0 then print, "n_index calculation went wrong at z = " + strtrim(z_i,2)
    small_r_vector_Mpc = r_vector_Mpc[defined_index]
    
    fractional_index=FLTARR(n_index)
    for pixel_i=0, n_index-1 do begin
      ;
      ;Loop over every pixel desired and calculate its fractional index compared to the input cube

      temp = min(abs(small_r_vector_Mpc[pixel_i] - box_vector_Mpc),closest_index) ;find the closest index

      ;find the next closest neighbor index for the interpolation
      ;diff = small_r_vector_Mpc[pixel_i] - box_vector_Mpc[closest_index]
      diff = box_vector_Mpc[closest_index] - small_r_vector_Mpc[pixel_i]
      if diff GT 0 then neighbor_index = closest_index - 1 else neighbor_index = closest_index + 1

      ;edge case: just use the previous pixel difference to determine the fractional index
      if neighbor_index EQ box_dim then neighbor_index = closest_index - 1

      ;Get the desired sky Mpc in terms of the fractional index of the sky Mpc from the theory cube
      fractional_index[pixel_i] = abs(diff) / abs(box_vector_Mpc[closest_index] - box_vector_Mpc[neighbor_index])
      if diff GT 0 then fractional_index[pixel_i] = closest_index - fractional_index[pixel_i] $
        else fractional_index[pixel_i] = closest_index + fractional_index[pixel_i]
    endfor

    ;Range of the defined pixels to be placed in the padded image
    if (n_index mod 2) EQ 0 then range = [dimension/2 - (n_index)/2,dimension/2 + (n_index/2)] $
      else range = [dimension/2 - (n_index-1)/2,dimension/2 + (n_index/2)]

    ;Sky in radian (theta) grid in units of mK
    bubble_image_grid[range[0]:range[1],range[0]:range[1],z_i] = $
      interpolate(bubble_image_arr_interpolated[*,*,z_i],fractional_index,fractional_index,/grid)

pix_num_filled[z_i] = range[1] - range[0]
stop
    ;; Sky in radians is the Fourier dual of u,v.
    ;; Sky in Mpc is the Fourier dual of kx,ky
    ;;take fourier transform with respect to rx[Mpc],ry[Mpc] to get a uvf cube in mK*Mpc^2 units
    fourier_convension = (dimension)^2.*r_res_Mpc^2 ; See FT_memo by Beardsley, factor of N is IDL specific
    uvf_cube_mKMpc2[*,*,z_i] = fourier_convension * $
      fft(fft(bubble_image_grid[*,*,z_i], /CENTER,DIMENSION=1), /CENTER,DIMENSION=2)

  endfor

  save, uvf_cube_mKMpc2, uv_res_inverserad, uv_range_inverserad, obs_freq, $
    filename = '/fred/oz048/MWA/CODE/FHD/Brad_sim/uvf_cube_mKMpc2.sav'

  ;save the bubble image and then undefine it to save space
  save, bubble_image_grid, theta_res_rad, r_res_Mpc_arr, obs_freq, $
    filename = '/fred/oz048/MWA/CODE/FHD/Brad_sim/thetaf_mKstr_cube.sav'
  undefine, bubble_image_grid, bubble_image_arr_interpolated

  ;now have space to calculate the Jy uvf cube
  ;to go from mK Mpc^2 to mK, divide by Dm[Mpc]^2. To go from mK to Jy/str, divide by the Rayleigh Jeans law.
  ;to go from Jy/str to Jy/pixel, divide by 1/rad squared (area of grid pixel)
  uvf_cube_Jy = complex(FLTARR(dimension,dimension,n_freq))
  ;for z_i=0,n_freq-1 do uvf_cube_Jy[*,*,z_i] = uvf_cube_mKMpc2[*,*,z_i] / (comoving_dist_los[z_i])^2 / Jy_to_mKstr[z_i] / theta_res_rad^2
  ;for z_i=0,n_freq-1 do uvf_cube_Jy[*,*,z_i] = uvf_cube_mKMpc2[*,*,z_i] / r_res_Mpc_arr[z_i]^2 / Jy_to_mKstr[z_i] * theta_res_rad^2 / (box_res_Mpc / comoving_dist_los[z_i]) 
  for z_i=0,n_freq-1 do uvf_cube_Jy[*,*,z_i] = uvf_cube_mKMpc2[*,*,z_i] / comoving_dist_los[z_i]^2 / Jy_to_mKstr[z_i] * (2*pi!)^2/uv_res_inverserad^2 * 2;(4*!pi)/uv_res_inverserad^2 
  save, uvf_cube_Jy, uv_res_inverserad, uv_range_inverserad, obs_freq, pix_num_filled, $
    filename = '/fred/oz048/MWA/CODE/FHD/Brad_sim/uvf_Jy_cube.sav'
  undefine, uvf_cube_Jy

  ;Calculate the resolution in the LoS direction in Mpc, take the frequency fft
  comov_los_diff = comoving_dist_los - shift(comoving_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comoving_dist_los)-2]
  rz_delta_mpc = mean(comov_los_diff)
  uvn_cube_mKMpc3 = (n_freq*rz_delta_mpc)*fft(uvf_cube_mKMpc2, /CENTER,DIMENSION=3)

  save, uvn_cube_mKMpc3, obs_freq, $
    filename = '/fred/oz048/MWA/CODE/FHD/Brad_sim/uvn_cube_mKMpc3.sav'

  ;Calculate the k grid axes
  n_kx = dimension
  kx_delta_invmpc = (2.*!dpi)*uv_res_inverserad / mean(comoving_dist_los)
  kx_invmpc = (dindgen(n_kx)-n_kx/2) * kx_delta_invmpc
  n_ky = dimension
  ky_delta_invmpc = kx_delta_invmpc
  ky_invmpc = (dindgen(n_ky)-n_ky/2) * kx_delta_invmpc
  n_kz = n_freq
  kz_delta_invmpc = (2.*!dpi) / (max(comoving_dist_los) - min(comoving_dist_los) + rz_delta_mpc)
  kz_invmpc = (dindgen(n_kz)-n_kz/2) * kz_delta_invmpc

  ;This is the theoretical window function (full sky) when there is no instrument
  theory_window_function_Mpc3 = (2.*!pi)^3./(kx_delta_invmpc*ky_delta_invmpc*kz_delta_invmpc) ;Mpc^3

  ;power in mK2Mpc3 with axes of kx_invMpc, ky_invMpc, kz_invMpc
  powercube_mK2Mpc3 = abs(uvn_cube_mKMpc3)^2 / theory_window_function_Mpc3
  save, powercube_mK2Mpc3, kx_delta_invmpc, ky_delta_invmpc,kz_delta_invmpc, $
    filename = '/fred/oz048/MWA/CODE/FHD/Brad_sim/powercube_mK2Mpc3.sav'
 
  temp = sqrt(rebin(kx_invmpc^2d, n_kx, n_ky, n_kz) + $
    rebin(reform(ky_invmpc^2d, 1, n_ky), n_kx, n_ky, n_kz) + $
    rebin(reform(kz_invmpc^2d, 1, 1, n_kz), n_kx, n_ky, n_kz))
  kbin=0.0180673
  k_hist = histogram(temp, binsize = kbin, min = .01, omax = max_p, $
    locations = lin_k_locs, reverse_indices = k_ri)
  k_centers_invmpc = lin_k_locs + kbin/2d
    k_edges_invmpc = [lin_k_locs, max(lin_k_locs) + kbin]
  n_bins = n_elements(k_hist)
  power_1d = dblarr(n_bins)
  for i=0, n_bins-1 do begin
    inds =  k_ri[k_ri[i] : k_ri[i+1]-1]
    power_1d[i] = total(powercube_mK2Mpc3[inds])/N_elements(inds)
  endfor

  save, power_1d, k_centers_invmpc, k_edges_invmpc, $
    filename = '/fred/oz048/MWA/CODE/FHD/Brad_sim/power_1d.sav'

  print, "uvf cube made for desired dimension/frequencies"
  run_report, start_mem, t0, silent=silent
end
