;put in polar, create flipped versions, transform back to uv
Jmat_use=Reform(*Jmat_interp[0,0,0],n_zen_slice,n_ang_slice)
Expand,Jmat_use,n_zen_ang,n_az_ang,Jmat_single
Expand,!Dpi/180*abs((reform(phi_arr,31,121))),n_zen_ang,n_az_ang,phi_single
;Expand,!Dpi/180*(reform(theta_arr,31,121)),n_zen_ang,n_az_ang,theta_single
Expand,(reform(theta_arr,31,121)),n_zen_ang,n_az_ang,theta_single
full_Jmat = [Jmat_single,reverse(Jmat_single)]
full_phi_single = [phi_single,reverse(phi_single)]
full_theta_single = [theta_single,reverse(theta_single)]

xvals_instrument2=full_theta_single*Sin(full_phi_single)
yvals_instrument2=full_theta_single*Cos(full_phi_single)
r_polar = sin(full_theta_single)
theta_polar = full_phi_single
u=r_polar*cos(theta_polar)
v=r_polar*sin(theta_polar)
fft_Jmat = (2*n_zen_ang*n_az_ang)*fft_shift(FFT(fft_shift(full_Jmat)))






;spherical transform on image beam
input_phi_arr = (*antenna1.az_arr)
input_theta_arr = (*antenna1.za_arr)
n_theta = 23
n_phi = n_theta*4
;; Avoid the pole by offsetting theta.
;; Spherical harmonic transform does not do well with multiple values at the pole.
spher_theta_arr = (FLTARR(n_phi)+1.)#((FINDGEN(n_theta)-0.5)*(90./n_theta))
spher_phi_arr = (FINDGEN(n_phi)*(360./n_phi))#(FLTARR(n_theta)+1.)
pix_hor=obs.dimension/antenna1.psf_scale
y = pix_hor/2*sin(spher_theta_arr*!DPi/180)*sin(spher_phi_arr*!Dpi/180)
x = pix_hor/2*sin(spher_theta_arr*!DPi/180)*cos(spher_phi_arr*!Dpi/180)
image_power_beam_za_az = Interpolate(image_power_beam,x+(antenna1.psf_image_dim/2.),y+(antenna1.psf_image_dim/2.),cubic=-0.5)

;; Remove the first theta value to avoid the pole
;; Build full sphere of beam as if the horizon wasn't present to fully define the basis
full_image_power_beam_za_az = ([[image_power_beam_za_az[ n_theta-1 : n_phi-1 , 1 : n_theta-1]],$
  [reverse(image_power_beam_za_az[ n_theta-1 : n_phi-1, 1 : n_theta-2],2)]])
full_spher_theta_arr = ([[spher_theta_arr[ n_theta-1 : n_phi-1, 1 : n_theta-1]],$
  [spher_theta_arr[ n_theta-1 : n_phi-1, 1 : n_theta-2]+90.]])
;; Repeat sphere to hyperresolve the spherical harmonic transform
;full_image_power_beam_za_az = [[full_image_power_beam_za_az,reverse(full_image_power_beam_za_az)],[full_image_power_beam_za_az,reverse(full_image_power_beam_za_az)]]
;full_spher_theta_arr = [[full_spher_theta_arr,full_spher_theta_arr],[full_spher_theta_arr+180.,full_spher_theta_arr+180.]]
full_image_power_beam_za_az = [[full_image_power_beam_za_az],[reverse(full_image_power_beam_za_az)]]
full_spher_theta_arr = [[full_spher_theta_arr],[full_spher_theta_arr+180.]]
;; Can I then take these extra flipped beam and convert it back to the xy grid to reduce horizon aliasing?
;full_spher_phi_arr = ([[spher_phi_arr[ n_theta-1 : n_phi-1, 1 : n_theta-1]],$
;  [spher_phi_arr[ n_theta-1 : n_phi-1, 1 : n_theta-2]]])
;full_spher_phi_arr = [[full_spher_phi_arr,full_spher_phi_arr],[full_spher_phi_arr,full_spher_phi_arr]]
full_spher_phi_arr = 
y = pix_hor/2*sin(full_spher_theta_arr*!DPi/180)*sin(full_spher_phi_arr*!Dpi/180)
x = pix_hor/2*sin(full_spher_theta_arr*!DPi/180)*cos(full_spher_phi_arr*!Dpi/180)

for rotate_i=1,3 do begin
inds = where(full_spher_theta_arr GE 90.*rotate_i,n_count)
if n_count GT 0 then begin
y[inds] = pix_hor/2*(rotate_i+abs(sin((full_spher_theta_arr[inds]-90.*rotate_i)*!DPi/180)))*sin(full_spher_phi_arr[inds]*!Dpi/180)
x[inds] = pix_hor/2*(rotate_i+abs(sin((full_spher_theta_arr[inds]-90.*rotate_i)*!DPi/180)))*cos(full_spher_phi_arr[inds]*!Dpi/180)
endif
endfor

psf_image_dim = antenna1.psf_image_dim
input_phi_arr = (*antenna1.az_arr)
input_theta_arr = (*antenna1.za_arr)

x_mesh = (FINDGEN(psf_image_dim)-(psf_image_dim/2.))#(FLTARR(psf_image_dim)+1.)
y_mesh = (FLTARR(psf_image_dim)+1.)#(FINDGEN(psf_image_dim)-(psf_image_dim/2.))
phi_mesh = atan(-x_mesh,-y_mesh) 

r_mesh = sqrt(x_mesh^2+y_mesh^2)
rotator=fix(r_mesh / (pix_hor/2))
r_mesh = r_mesh - (rotator * (pix_hor/2))
flip_inds = where((rotator mod 2),n_flip)
if n_flip GT 0 then r_mesh[flip_inds] = pix_hor/2. - r_mesh[flip_inds] 
y_flip = r_mesh/sqrt(tan(phi_mesh)^2+1.)
x_flip = sqrt(r_mesh^2 - y_flip^2)

;theta_mesh = reform(interpolate(input_theta_arr,reform(x_flip+psf_image_dim/2.,N_elements(x_flip)),$
;  reform(y_flip+psf_image_dim/2.,N_elements(y_flip)),cubic=-0.5),psf_image_dim,psf_image_dim)

inds_x_neg = where(phi_mesh GT 0,n_x_neg)
if n_x_neg GT 0 then x_flip[inds_x_neg] *= -1.
inds_y_neg = where(abs(phi_mesh) LT !Pi/2., n_y_neg)
if n_y_neg GT 0 then y_flip[inds_y_neg] *= -1

inf_image_power_beam = reform( interpolate(image_power_beam, reform(x_flip+psf_image_dim/2, N_elements(x_flip)),$
  reform(y_flip+psf_image_dim/2, N_elements(y_flip)), cubic=-0.5) ,psf_image_dim,psf_image_dim)



inf_image_power_beam = interpolate(image_power_beam, reform(theta_mesh,N_elements(theta_mesh))/obs.degpix,$
  reform(phi_mesh,N_elements(phi_mesh))/obs.degpix,cubic=-0.5)
inf_image_power_beam = reform(inf_image_power_beam, psf_image_dim, psf_image_dim)

;; Sphere harmonics code does not benefit from extra flipping/additions since cos(theta) repeats
full_image_power_beam_za_az = ([image_power_beam_za_az[ n_theta-1 : n_phi-1 , 1 : n_theta-1]])
full_spher_theta_arr = ([spher_theta_arr[ n_theta-1 : n_phi-1, 1 : n_theta-1]])

C_lm=spherical_transform((full_image_power_beam_za_az),cos(full_spher_theta_arr*!Dpi/180))

C_transform = complex(FLTARR(N_Elements(pix_use)))
for l=0,(size(C_lm))[1]-1 do for m=-l,l do C_transform = C_transform + C_lm[l,abs(m)]*SPHER_HARM((*antenna1.za_arr)[pix_use]*!Dpi/180,(*antenna1.az_arr)[pix_use]*!Dpi/180, l, m)

uv_test = complex(FLTARR(1400,1400))
uv_test[pix_use] = C_transform
uv_range = 1/(90./n_theta)
uv_delta = 1/(360.)
;70 pixels for 90 degrees range

;recreate the resolution of a 2D uv grid and fill in
uv_res = 1/(1400.*obs.degpix*!pi/180)
scale_factor = (!pi / uv_res / 2 )
n_pix_fine = round(1400*scale_factor)
expand, xvals_instrument, n_pix_fine,n_pix_fine, xvals_instrument_fine
expand, yvals_instrument, n_pix_fine,n_pix_fine, yvals_instrument_fine
zen_ang_inst_fine=Sqrt(xvals_instrument_fine^2+yvals_instrument_fine^2)
az_ang_inst_fine=270. - Atan(temporary(yvals_instrument_fine),temporary(xvals_instrument_fine))*!Radeg+180.

horizon_test_fine = where(zen_ang_inst_fine GE 90, n_horizon_test_fine,complement=pix_use_fine)
J_transform = complex(FLTARR(N_Elements(pix_use_fine)))
for l=0,(size(C_lm))[1]-1 do for m=-l,l do J_transform = J_transform + C_lm[l,abs(m)]*SPHER_HARM( zen_ang_inst_fine[pix_use_fine] * !Dpi/180, az_ang_inst_fine[pix_use_fine] * !Dpi/180, l, m)






;spherical transform for Jones
n_input_phi = N_elements(Uniq(phi_arr))
n_input_theta = N_elements(phi_arr) / n_input_phi
overres_J_factor = 2.

;overresolve the input Jones matrix and associated coords
J_input = transpose(reform((Jmat_arr[0,0,0,*]),n_input_theta,n_input_phi))
theta_input = !Dpi/180*transpose(reform(theta_arr,n_input_theta,n_input_phi))
expand, temporary(J_input), n_input_phi*overres_J_factor,n_input_theta*overres_J_factor, J_input_highres 
expand, temporary(theta_input), n_input_phi*overres_J_factor,n_input_theta*overres_J_factor, theta_input_highres 

;flip the beam so that above *and* below the horizon is defined to reduce aliasing
;remove the zenith point since it is a multi-defined pole in the input simulation 
full_J_input_highres = ([[J_input_highres[ n_input_theta*overres_J_factor-1 : n_input_phi*overres_J_factor-1 , 1 : n_input_theta*overres_J_factor-1]],$
  [reverse(J_input_highres[ n_input_theta*overres_J_factor-1 : n_input_phi*overres_J_factor-1, 1 : n_input_theta*overres_J_factor-2])]])
full_theta_input_highres = ([[theta_input_highres[ n_input_theta*overres_J_factor-1 : n_input_phi*overres_J_factor-1, 1 : n_input_theta*overres_J_factor-1]],$
  [reverse(theta_input_highres[ n_input_theta*overres_J_factor-1 : n_input_phi*overres_J_factor-1, 1 : n_input_theta*overres_J_factor-2]+ 90.*!Dpi/180)]])
C_lm=spherical_transform((full_J_input_highres),cos(full_theta_input_highres))

;recreate the resolution of a 2D uv grid and fill in
uv_res = 1/(1400.*obs.degpix*!pi/180)
scale_factor = (!pi / uv_res / 2 ) 
n_pix_fine = round(1400*scale_factor)
expand, xvals_instrument, n_pix_fine,n_pix_fine, xvals_instrument_fine
expand, yvals_instrument, n_pix_fine,n_pix_fine, yvals_instrument_fine
zen_ang_inst_fine=Sqrt(xvals_instrument_fine^2+yvals_instrument_fine^2)
az_ang_inst_fine=270. - Atan(temporary(yvals_instrument_fine),temporary(xvals_instrument_fine))*!Radeg+180.

horizon_test_fine = where(zen_ang_inst_fine GE 90, n_horizon_test_fine,complement=pix_use_fine)
J_transform = complex(FLTARR(N_Elements(pix_use_fine)))
for l=0,(size(C_lm))[1]-1 do for m=-l,l do J_transform = J_transform + C_lm[l,abs(m)]*SPHER_HARM( zen_ang_inst_fine[pix_use_fine] * !Dpi/180, az_ang_inst_fine[pix_use_fine] * !Dpi/180, l, m)


;***reg
;spherical transform on interp
n_input_phi = N_elements(Uniq(phi_arr))
n_input_theta = N_elements(phi_arr) / n_input_phi

;Input Jones matrix and associated coords
J_input = transpose(reform((Jmat_arr[0,0,0,*]),n_input_theta,n_input_phi))
theta_input = !Dpi/180*transpose(reform(theta_arr,n_input_theta,n_input_phi))

;flip the beam so that above *and* below the horizon is defined to reduce aliasing
;remove the zenith point since it is a multi-defined pole in the input simulation 
full_J_input = ([[J_input[ n_input_theta-1 : n_input_phi-1 , 1 : n_input_theta-1]],$
  [reverse(J_input[ n_input_theta-1 : n_input_phi-1, 1 : n_input_theta-2])]])
full_theta_input = ([[theta_input[ n_input_theta-1 : n_input_phi-1, 1 : n_input_theta-1]],$
  [reverse(theta_input[ n_input_theta-1 : n_input_phi-1, 1 : n_input_theta-2]+ 90.*!Dpi/180)]])
C_lm=spherical_transform((full_J_input),cos(full_theta_input))

;recreate the resolution of a 2D uv grid and fill in
uv_res = 1/(1400.*obs.degpix*!pi/180)
scale_factor = (!pi / uv_res / 2 )
n_pix_fine = round(1400*scale_factor)
expand, xvals_instrument, n_pix_fine,n_pix_fine, xvals_instrument_fine
expand, yvals_instrument, n_pix_fine,n_pix_fine, yvals_instrument_fine
zen_ang_inst_fine=Sqrt(xvals_instrument_fine^2+yvals_instrument_fine^2)
az_ang_inst_fine=270. - Atan(temporary(yvals_instrument_fine),temporary(xvals_instrument_fine))*!Radeg+180.

horizon_test_fine = where(zen_ang_inst_fine GE 90, n_horizon_test_fine,complement=pix_use_fine)
J_transform = complex(FLTARR(N_Elements(pix_use_fine)))
for l=0,(size(C_lm))[1]-1 do for m=-l,l do J_transform = J_transform + C_lm[l,abs(m)]*SPHER_HARM( zen_ang_inst_fine[pix_use_fine] * !Dpi/180, az_ang_inst_fine[pix_use_fine] * !Dpi/180, l, m)


;***
uv_test_fine = complex(FLTARR(n_pix_fine,n_pix_fine))
uv_test_fine[pix_use_fine]=J_transform
quick_image, abs(fft_test[499:899,499:899]),window=1
quick_image, (abs(uv_test_fine[803:2202,803:2202]))[499:899,499:899],window=1,/log






;spherical transform
J= transpose(reform((Jmat_arr[0,0,0,*]),31,121))  
theta = !Dpi/180*transpose(reform(theta_arr,31,121))
full_J = ([[J[30:120,1:30]],[reverse(J[30:120,1:29])]])
full_theta = ([[theta[30:120,1:30]],[reverse(theta[30:120,1:29]+ 90.*!Dpi/180)]])
C_lm=spherical_transform((full_J),cos(full_theta))   

zen_ang_inst=Sqrt(xvals_instrument[pix_use]^2+yvals_instrument[pix_use]^2) 
az_ang_inst=270. - Atan(yvals_instrument[pix_use],xvals_instrument[pix_use])*!Radeg+180.

uv_test = complex(FLTARR(1400,1400))  
result2 = complex(FLTARR(N_elements(pix_use)))
for l=0,39 do for m=-l,l do result2 = result2 + result[l,abs(m)]*SPHER_HARM( zen_ang_inst * !Dpi/180, az_ang_inst * !Dpi/180, l, m)
uv_test[pix_use]=result2


;compare
test = complex(FLTARR(1400,1400)) 
test[pix_use] = (*Jones_matrix[0,0,0]) 
fft_test = fft_shift(fft(fft_shift(test)))  
uv_res = 1/(1400.*obs.degpix*!pi/180)
scale_factor = (!pi / uv_res / 2 ) * overres_J_factor
n_pix_fine = round(1400*scale_factor)
expand, xvals_instrument, n_pix_fine,n_pix_fine, xvals_instrument_fine
expand, yvals_instrument, n_pix_fine,n_pix_fine, yvals_instrument_fine
zen_ang_inst_fine=Sqrt(xvals_instrument_fine^2+yvals_instrument_fine^2)
az_ang_inst_fine=270. - Atan(yvals_instrument_fine,xvals_instrument_fine)*!Radeg+180.
horizon_test_fine = where(zen_ang_inst_fine GE 89, n_horizon_test_fine,complement=pix_use_fine)
result2 = complex(FLTARR(N_Elements(pix_use_fine)))
for l=0,39 do for m=-l,l do result2 = result2 + result[l,abs(m)]*SPHER_HARM( zen_ang_inst_fine[pix_use_fine] * !Dpi/180, az_ang_inst_fine[pix_use_fine] * !Dpi/180, l, m)
uv_test_fine = complex(FLTARR(n_pix_fine,n_pix_fine))
uv_test_fine[pix_use_fine]=result2

quick_image, abs(fft_test[499:899,499:899]),window=1
quick_image, (abs(uv_test_fine[803:2202,803:2202]))[499:899,499:899],window=1,/log

;polar transform
phi = !Dpi/180*abs(transpose(reform(phi_arr,31,121))-270)
;full_phi = phi[30:120,1:30]
full_phi = phi[*,1:30]
full_theta = theta[*,1:30]
;full_phi = ([[phi[30:120,1:30]],[reverse(phi[30:120,1:29]+ 90.*!Dpi/180)]])
full_r = abs(cos(!Dpi/90 - full_theta))
J_new = J[*,1:30]

mesh=(FINDGEN(30)+1)*1/30.
J_interp = J_new
full_new_theta=full_theta
for phi_i=0,120 do J_interp[phi_i,*] = interpolate(J_new[phi_i,*],FINDGEN(30)*(full_r[phi_i,*]/reverse(mesh[*])),cubic=-0.5)
for phi_i=0,120 do full_new_theta[phi_i,*] = !Dpi/90 - acos(phi_i*full_r[phi_i,*]/reverse(mesh[*]))

J_interp_fft=fft_shift(FFT(fft_shift(J_interp))) 
quick_image, abs(fft_shift(fft(j_interp_fft,/inverse)))

norm_new_theta = reform(transpose(full_new_theta),121*30.)
norm_phi = reform(transpose(full_phi),121*30.)
;norm_J_interp = reform(transpose(J_interp),121*30.)

Expand,transpose(J_interp_fft),n_zen_ang,n_az_ang,Jmat_single2
Jones_matrix_polar = Interpolate(Jmat_single2,zen_ang_inst/interp_res,az_ang_inst/interp_res,cubic=-0.5)

Jmat_use=Reform(*Jmat_interp[0,0,0],n_zen_slice,n_ang_slice)
Expand,Jmat_use,n_zen_ang,n_az_ang,Jmat_single
Expand,!Dpi/180*abs((reform(phi_arr,31,121))),n_zen_ang,n_az_ang,phi_single
Expand,!Dpi/180*(reform(theta_arr,31,121)),n_zen_ang,n_az_ang,theta_single
r_single = abs(cos(!Dpi/90 - theta_single))
mesh=(FINDGEN(n_zen_ang))*1/n_zen_ang
J_interp = Jmat_single
for phi_i=0,n_az_ang-1 do J_interp[*,phi_i] = interpolate(Jmat_single[*,phi_i],FINDGEN(n_zen_ang)*(r_single[*,phi_i]/reverse(mesh[*])),cubic=-0.5)
;put in polar, create flipped versions, transform back to uv
