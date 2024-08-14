pro healpix_pixel_choose, sample_healpix_file, output_file=output_file

  if ~keyword_set(output_file) then output_file='inds.sav'

  ;;A healpix file from FHD for the field in question. If no healpix pixel files exist, set
  ;;  restrict_hpx_inds=0 and specify the desired nside to make a full healpix cube for the
  ;;  field.
  ;;Should contain the healpix cube, inds, nside, n_avg, and obs structure
  restore, sample_healpix_file

  ;Get cosmology information
  freq = (*obs.baseline_info).freq
  n_freq = obs.n_freq
  n_freqbins = n_freq / n_avg
  frequencies = dblarr(n_freqbins)
  for i=0, n_freqbins-1 do begin
    frequencies[i] = mean(freq[i*n_avg:i*n_avg+(n_avg-1)]) / 1e6 ;; in MHz
  endfor
  z0_freq = 1420.40d ;; MHz
  redshifts = z0_freq/frequencies - 1d
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los, hubble_param = hubble_param
  
  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  z_mpc_mean = float(mean(comov_dist_los))
  kperp_lambda_conv = z_mpc_mean / (2.*!pi)

  ;; get pixel vectors
  pix2vec_ring, nside, hpx_inds, pix_center_vec
  ;; find mid point (work in x/y because of possible jumps in phi)
  vec_mid = [mean(pix_center_vec[*,0]), mean(pix_center_vec[*,1]), mean(pix_center_vec[*,2])]
  theta0 = acos(vec_mid[2])
  phi0 = atan(vec_mid[1], vec_mid[0])

  ;; To go to flat sky, rotate patch to zenith and flatten.
  ;; To get to current location, need to first rotate around z by
  ;; phi, then around y by -theta, then around z by -phi
  ;; use inverse to rotate back to zenith
  rot_matrix = get_rot_matrix(theta0, phi0, /inverse)
  new_pix_vec = rot_matrix ## pix_center_vec

  ;; then rotate to make as rectangular as possible
  pred_angle = healpix_rot(new_pix_vec[*,0], new_pix_vec[*,1])

  consv_delta_kperp_rad = 4.5* mean(frequencies*1e6) * z_mpc_mean / (3e8 * kperp_lambda_conv) ;use 4.5m to be conservative
  consv_xy_len = 2*!pi/consv_delta_kperp_rad
  radius = consv_xy_len/2.*sqrt(2)*1.1
  query_disc, nside, vec_mid, radius, listpix, nlist, /inc
  pix2vec_ring, nside, listpix, list_center_vec
  new_list_vec = rot_matrix ## list_center_vec
  x_list_rot = new_list_vec[*,0] * cos(pred_angle) - new_list_vec[*,1] * sin(pred_angle)
  y_list_rot = new_list_vec[*,0] * sin(pred_angle) + new_list_vec[*,1] * cos(pred_angle)
  ;cgplot, x_list_rot, y_list_rot, psym=3
              
  consv_lims = [-1*consv_xy_len/2., -1*consv_xy_len/2., consv_xy_len/2., consv_xy_len/2.]
  ;cgpolygon, reform(rebin(consv_lims[[0,2]], 2,2),4), reform(rebin(reform(consv_lims[[1,3]],1,2), 2,2),4), color='aqua'

  wh_listpix_close = where(x_list_rot ge consv_lims[0] and x_list_rot le consv_lims[2] and $
    y_list_rot ge consv_lims[1] and y_list_rot le consv_lims[3], count_list_close)            
  
  hpx_inds = listpix[wh_listpix_close]
  stop            
  save, file=output_file, nside, hpx_inds
  print, 'Inds file saved to ' + output_file


end
