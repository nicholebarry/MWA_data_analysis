;;******* Spherical harmonic transform
i = Complex(0,1)
psf_image_dim = antenna1.psf_image_dim
input_phi_arr = (*antenna1.az_arr)
input_theta_arr = (*antenna1.za_arr)
n_theta = 60 
n_phi = n_theta*4
pix_hor=obs.dimension/antenna1.psf_scale

;spher_theta_arr = (FLTARR(n_phi)+1.)#((FINDGEN(n_theta)-0.5)*(90./n_theta))
spher_theta_arr = (FLTARR(n_phi)+1.)#((FINDGEN(n_theta)+1)*(90./n_theta))
spher_phi_arr = (FINDGEN(n_phi)*(360./n_phi))#(FLTARR(n_theta)+1.)  
y = pix_hor/2*sin(spher_theta_arr*!DPi/180)*sin(spher_phi_arr*!Dpi/180)
x = pix_hor/2*sin(spher_theta_arr*!DPi/180)*cos(spher_phi_arr*!Dpi/180)
image_power_beam_za_az = Interpolate(real_part(image_power_beam),x+(antenna1.psf_image_dim/2.),y+(antenna1.psf_image_dim/2.),cubic=-0.5)
full_image_power_beam_za_az = ([[image_power_beam_za_az[ 0 : n_phi-1 , 0 : n_theta-1]],$
  [reverse(image_power_beam_za_az[ 0 : n_phi-1, 0 : n_theta-2],2)]])
full_spher_theta_arr = ([[spher_theta_arr[ 0 : n_phi-1, 0 : n_theta-1]],$
  [spher_theta_arr[ 0 : n_phi-1, 0 : n_theta-2]+90.]])
;full_image_power_beam_za_az = [[full_image_power_beam_za_az],[reverse(full_image_power_beam_za_az)]]  
;full_spher_theta_arr = [[full_spher_theta_arr],[full_spher_theta_arr+180.]] 

full_image_power_beam_za_az = image_power_beam_za_az
full_spher_theta_arr = spher_theta_arr
C_lm=spherical_transform((full_image_power_beam_za_az),cos(full_spher_theta_arr*!Dpi/180))
;C_lm_neg=spherical_transform(reverse(full_image_power_beam_za_az,1),cos(reverse(full_spher_theta_arr,1)*!Dpi/180))

; Healpix option
nside=1024
;ang2vec,obs.obsdec,obs.obsra,cen_coords,/astro
;Query_disc,nside,cen_coords,180,hpx_inds_180,ninds,/deg
npix=NSIDE2NPIX(nside)
hpx_inds=INDGEN(npix,/long)
;Query_disc,nside,cen_coords,90,hpx_inds_90,ninds,/deg
pix2vec_ring,nside,hpx_inds,pix_coords
vec2ang,pix_coords,pix_theta,pix_phi;,/astro
;Jdate_use = obs.Jd0
;Eq2Hor,pix_ra,pix_dec,Jdate_use,alt_arr1,az_arr1,lat=obs.lat,lon=obs.lon,alt=obs.alt,/precess,/nutate,/refract
;pix_theta = 90. - (alt_arr1)
;pix_phi = (az_arr1)
;y = pix_hor/2*sin(pix_theta*!Dpi/180)*sin(pix_phi*!Dpi/180)
;x = pix_hor/2*sin(pix_theta*!Dpi/180)*cos(pix_phi*!Dpi/180)
y = pix_hor/2*sin(pix_theta)*sin(pix_phi)
x = pix_hor/2*sin(pix_theta)*cos(pix_phi)
pix_beam = Interpolate(real_part(image_power_beam),x+(antenna1.psf_image_dim/2.),y+(antenna1.psf_image_dim/2.),cubic=-0.5)
;pix_beam_full = dblarr(N_elements(hpx_inds_180))
;pix_beam_full[hpx_inds_90] = pix_beam
;healpix_quickimage, pix_beam,hpx_inds_180,nside,ordering='ring'
;!healpix.path.bin.f90='/apps/skylake/software/compiler/gcc/6.4.0/healpix/3.50//binf95/' 
ianafast, sqrt(pix_beam),/ring,binpath='binf95/anafast',alm1_out=file_path_fhd + '_alm.fits',nlmax=150,won=0
lun = fxposit(file_path_fhd + '_alm.fits', 0,/readonly)
data_struct=mrdfits(lun,1,data_header,/silent)
max_l=sxpar(data_header,'max-lpol')
;i = l^2 + l + m + 1
;C_lm = FLTARR(max_l, 2*max_l)
C_lm = DBLARR(max_l, max_l)
for l=0,max_l-1 do for m=0,l do begin $
  ind = where(data_struct.index EQ  (l^2 + l + m + 1), n_count) & $
  C_lm[l,m] = data_struct[ind].real & $
endfor
;l_inds = where(mean(abs(C_lm),dimension=2) GT 1e-5,n_count_l)
;if n_count_l GT 0 then max_l = max(l_inds) else message, 'No spherical harmonic coefficients for the image beam'

mesh_res=140;500
x_mesh = (FINDGEN(mesh_res)-(mesh_res/2.))#(FLTARR(mesh_res)+1.)    
y_mesh = (FLTARR(mesh_res)+1.)#(FINDGEN(mesh_res)-(mesh_res/2.)) 
phi_mesh = atan(-x_mesh,-y_mesh) 
r_mesh = sqrt(x_mesh^2+y_mesh^2) 
y_interp = r_mesh*(pix_hor/mesh_res)/sqrt(tan(phi_mesh)^2+1)
x_interp = sqrt( (r_mesh*(pix_hor/mesh_res))^2 - y_interp^2)
inds_x_neg = where(phi_mesh GT 0,n_x_neg) 
if n_x_neg GT 0 then x_interp[inds_x_neg] *= -1.
inds_y_neg = where(abs(phi_mesh) LT !Pi/2., n_y_neg) 
if n_y_neg GT 0 then y_interp[inds_y_neg] *= -1 
za_mesh = interpolate(input_theta_arr, reform(x_interp+psf_image_dim/2.,N_elements(x_interp)), $
  reform(y_interp+psf_image_dim/2.,N_elements(y_interp)),cubic=-0.5)*!Dpi/180
;pix_use_interp = where(za_mesh LT !pi/2)   
pix_use_interp = where((r_mesh*(pix_hor/mesh_res)) LT (pix_hor/2.),n_count)


u_res = 1/(!pi)
u_range = ceil(1/(!pi/mesh_res))
u_mesh = (FINDGEN(u_range) - (u_range/2.))#(FLTARR(u_range)+1.)
v_mesh = (FLTARR(u_range)+1.)#(FINDGEN(u_range) - (u_range/2.))
uv_mesh = sqrt(u_mesh^2+v_mesh^2)

C_transform = complex(FLTARR(N_Elements(pix_use_interp))) 

bessel=PTRARR(max_l,/allocate)
spher=complex(FLTARR((max_l,max_l))
B = complex(FLTARR(u_range,u_range))
for l=0,max_l-1 do begin & $
  *bessel[l] = sqrt(!pi/(2*uv_mesh))*BESELJ(uv_mesh,l+.5) & $
  for m=0,l do begin & $
    spher[l,m] = total(conj(SPHER_HARM(za_mesh[pix_use_interp],phi_mesh[pix_use_interp],l,m))) & $
  endfor & $
  B = B + (i^l * (*bessel[l]) * total(C_lm[l,*]*spher[l,*])) & $
endfor

;for u=u_range/2.,u_range-1 do begin & $
;  for v=u_range/2.,u_range-1 do begin & $
;    for l=0,max_l-1 do begin & $ 
;        B[u,v] = total(i^l * (*bessel[l])[u,v] * m_total) & $
;      for m=0,l do begin & $
;        m_total[l] = total(C_lm[l,m] * conj(*spher[l,m])) & $
;      endfor & $
;    endfor & $
;  endfor & $
;endfor
;for l=0, (size(C_lm))[1]-1 do for u=u_range/2.,u_range-1 do for v=u_range/2.,u_range-1 do *bessel[l] = sqrt(!pi/(2*uv_mesh[u,v]*sqrt(cos(za_mesh)^2+cos(phi_mesh)^2)))*BESELJ(uv_mesh[u,v] * sqrt(cos(za_mesh)^2+cos(phi_mesh)^2),l+.5)

;Grid version
for l=0,(size(C_lm))[1]-1 do for m=-l,l do begin $
  spher = SPHER_HARM(za_mesh[pix_use_interp],phi_mesh[pix_use_interp],l,m) & $
  C_transform = C_transform + real_part(C_lm[l,m] *( (spher) )) & $;+ (-1.)^m * conj(spher))) & $
endfor

;Healpix version
;i = l^2 + l + m + 1
for l=0,max_l-1 do for m=0,l do begin $
  spher = SPHER_HARM(za_mesh[pix_use_interp],phi_mesh[pix_use_interp],l,m) & $
  C_transform = C_transform + i^l * C_lm[l,m] * conj(spher) & $
endfor
  ;C_transform = C_transform + (-1.)^(m+2) * C_lm[l,m] * real_part(spher) & $


uv_test=complex(FLTARR(mesh_res,mesh_res))
uv_test[pix_use_interp]=C_transform
quick_image, uv_test
  C_transform = C_transform + i^l * (C_lm[l,m] * conj(spher) + (-1.)^m * C_lm_neg[l,m] * spher) & $
  ;C_transform = C_transform + (C_lm[l,m]*spher+C_lm_neg[l,m]*((-1)^m)*conj(spher)) & $
;;**************


;; ********* Polar transform
psf_image_dim = antenna1.psf_image_dim
input_phi_arr = (*antenna1.az_arr)
input_theta_arr = (*antenna1.za_arr)
n_theta = 40
n_phi = n_theta*4

spher_theta_arr = (FLTARR(n_phi)+1.)#((FINDGEN(n_theta))*(90./n_theta))
spher_phi_arr = (FINDGEN(n_phi)*(360./n_phi))#(FLTARR(n_theta)+1.)
pix_hor=obs.dimension/antenna1.psf_scale
y = pix_hor/2*sin(spher_theta_arr*!DPi/180)*sin(spher_phi_arr*!Dpi/180)
x = pix_hor/2*sin(spher_theta_arr*!DPi/180)*cos(spher_phi_arr*!Dpi/180)

image_power_beam_za_az = Interpolate(image_power_beam,x+(antenna1.psf_image_dim/2.),y+(antenna1.psf_image_dim/2.),cubic=-0.5)
full_image_power_beam_za_az = ([[image_power_beam_za_az],[reverse(image_power_beam_za_az[ *, 0 : n_theta-2],2)]])
full_spher_theta_arr = ([[spher_theta_arr],[spher_theta_arr[ *, 0 : n_theta-2]+90.]])
full_image_power_beam_za_az = [[full_image_power_beam_za_az],[reverse(full_image_power_beam_za_az)]]  
full_spher_theta_arr = [[full_spher_theta_arr],[full_spher_theta_arr+180.]] 

image_fft=fft_shift(FFT(fft_shift(image_power_beam_za_az)))
;image_fft=FFT(image_power_beam_za_az)
uv_res_x = 1/(360.)
uv_res_y = 1/(180.)
uv_range_x = 1/(90./n_theta)
uv_range_y = 1/(360./n_phi)
x_axis = FINDGEN(uv_range_x/uv_res_x)*uv_res_x
y_axis = FINDGEN(uv_range_y/uv_res_y)*uv_res_y

mesh_res=uv_range_x/uv_res_x
x_mesh = (FINDGEN(mesh_res)-(mesh_res/2.))#(FLTARR(mesh_res)+1.)
y_mesh = (FLTARR(mesh_res)+1.)#(FINDGEN(mesh_res)-(mesh_res/2.))
phi_mesh = atan(-x_mesh,-y_mesh)
r_mesh = sqrt(x_mesh^2+y_mesh^2)
y_interp = r_mesh*(pix_hor/mesh_res)/sqrt(tan(phi_mesh)^2+1)
x_interp = sqrt( (r_mesh*(pix_hor/mesh_res))^2 - y_interp^2)
inds_x_neg = where(phi_mesh GT 0,n_x_neg)
if n_x_neg GT 0 then x_interp[inds_x_neg] *= -1.
inds_y_neg = where(abs(phi_mesh) LT !Pi/2., n_y_neg)
if n_y_neg GT 0 then y_interp[inds_y_neg] *= -1
za_mesh = interpolate(input_theta_arr, reform(x_interp+psf_image_dim/2.,N_elements(x_interp)), $
  reform(y_interp+psf_image_dim/2.,N_elements(y_interp)),cubic=-0.5)*!Dpi/180
pix_use_interp = where(za_mesh LT !pi/2)
image_power_beam_za_az = Interpolate(image_fft,x+(antenna1.psf_image_dim/2.),y+(antenna1.psf_image_dim/2.),cubic=-0.5)





