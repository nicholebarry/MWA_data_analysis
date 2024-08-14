function Ez_inv, redshift
  Common cosmological_params, h, Om, Ol, Ok
  return, 1d/sqrt(Om*(1d + redshift)^3d + Ok*(1d + redshift)^2d + Ol)
end

pro comoving_distance_integral, redshift, hubble_param = hubble_param, omega_matter = omega_matter, omega_lambda = omega_lambda, $
                        omega_curv = omega_curv, hubble_distance = hubble_distance, Ez = Ez, comoving_dist_los = comoving_dist_los, $
                        wedge_factor = wedge_factor

  ;; using WMAP7 values as defaults
  if not keyword_set(hubble_param) then hubble_param = 0.71d
  if not keyword_set(omega_matter) then omega_matter = 0.27d
  if not keyword_set(omega_lambda) then omega_lambda = 0.73d
  if not keyword_set(omega_curv) then omega_curv = 0.d

  Common cosmological_params
  h = hubble_param
  Om = omega_matter
  Ol = omega_lambda
  Ok = omega_curv

  hubble_distance = 3000d / hubble_param
  
 
intreg = DBLARR(1200)

  for redshift=0, 1200-1 do begin

    Ez = sqrt(Om*(1d + redshift)^3d + Ok*(1d + redshift)^2d + Ol)
    intreg[redshift] = qromb('Ez_inv', 0, redshift, /double)

  endfor

  intreg_new = DBLARR(1200)
  Ok = 0.01d
  for redshift=0, 1200-1 do begin

    Ez = sqrt(Om*(1d + redshift)^3d + Ok*(1d + redshift)^2d + Ol)
    intreg_new[redshift] = qromb('Ez_inv', 0, redshift, /double)

  endfor
stop
  ;comoving_dist_los = hubble_distance * intreg

  ;wedge_factor = comoving_dist_los * Ez / (hubble_distance * (redshift+1d))
end