pro bubble_stats
  ;, bandwidth, dnu, grid_size_wavelength
  print, 'Starting...'
  
  ; parse command line args
  compile_opt strictarr
  args = Command_Line_Args(count=nargs)
  bandwidth = double(args[0])
  ;obs_id = '1061316296'
  dnu = double(args[1])
  ;output_directory = '/nfs/eor-00/h1/nbarry/'
  grid_size_wavelength = double(args[2])
  ;version = 'default'
  ;cmd_args={version:version}
  
  
  ; bubble_stats
  ; A script to compile some of the imaging/filtering I've been doing
  ; but also incorporate some stats ideas Matt McQuinn had
  
  ; Set up a few options to start
  z_0=6.7892 ; temporary for nichole
  xi=0.79
  ;grid_size=2.*0.8226  ; 0.8226=.5lambda for 1.82MHz
  grid_size=(299792458./(182.*10^6.))*grid_size_wavelength  ; 0.8226=.5lambda for 182MHz
  ;bandwidth = 30.72/16. * 10^6.  ; Hz
  ;bandwidth = 30.72/32. * 10^6.  ; Hz
  ;dnu = 80000.  ; Hz
  
  ; Read in data
  print, 'Reading Data...'
  t0=Systime(1)
  ng = 512. ;I have no fucking clue here
  data_type = 'double'    ;;;;;;;;;;;probably redundant??????
  
  case xi of
    ;0.68: filename1 = '/data4/MWA/data_for_nichole/21cm_sims/delta_21cm_1gpch_xi_0.68_z6.90.dat'
    ;0.79: filename1 = '/data4/MWA/data_for_nichole/21cm_sims/delta_21cm_1gpch_xi_0.79_z6.90.dat'
    ;0.89: filename1 = '/data4/MWA/data_for_nichole/21cm_sims/delta_21cm_1gpch_xi_0.89_z6.90.dat'
    0.68: filename1 = '/nfs/eor-00/h1/nbarry/MWA/development/bubbles/delta_21cm_1gpch_xi_0.68_z6.90.dat'
    0.79: filename1 = '/nfs/eor-00/h1/nbarry/MWA/development/bubbles/delta_21cm_1gpch_xi_0.79_z6.90.dat'
    0.89: filename1 = '/nfs/eor-00/h1/nbarry/MWA/development/bubbles/delta_21cm_1gpch_xi_0.89_z6.90.dat'
    else: begin
      print, 'xi not recognized'
      return
    end
  endcase
  
  ;OpenR, lun, filename1, /Get_Lun
  
  ;slice = zeros(ng,ng,ng);
  ;fread(fh1,4,'ubit8')  ;read file fh1 into a column vector with 4 elements using 8bit ints, but into what fricking variable? Maybe read it to get it out of the way
  ;readf, lun, read_temp  ;I think this part is unnecessary?
  
  ;sim_image = fread(fh1,ng^3,data_type);
  ;readf, lun, sim_image  ;I hope this is in the right form, needs to be checked
  sim_image=read_binary(filename1, DATA_START=4, DATA_TYPE=5)
  sim_image = reform(sim_image,[ng,ng,ng])
  
  ;Free_lun, lun
  
  ; set up axis
  h = 0.70  ; hubble parameter
  Lb = 1000./h  ; Mpc, image extent
  ;Lb = 1500./h ;Trying to get something here...
  ;Lb = 1000  ; h Mpc, keep hubble parameter in there
  dr = Lb/ng  ; Mpc, image resolution
  
  r_sim_axis_range=Lb-1/ng
  r_sim_axis_delta=Lb/ng
  r_sim_axis_start=(Lb-1/ng)/2
  r_sim_axis=(findgen(round(r_sim_axis_range/r_sim_axis_delta+1))*r_sim_axis_delta)-r_sim_axis_start
  
  ; cut out one row to make it an odd size cube
  
  image_dim = size(sim_image, /Dimensions)
  axis_dim = size(r_sim_axis, /Dimensions)
  sim_image = sim_image[Where(~Histogram([0],MIN=0,MAX=image_dim[1]-1),/NULL),Where(~Histogram([0],MIN=0,MAX=image_dim[1]-1),/NULL),*]
  r_sim_axis = r_sim_axis[Where(~Histogram([0],MIN=0,MAX=axis_dim-1),/NULL)]
  
  timing=Systime(1)-t0
  
  
  ; Set up axes, physics, and cosmology
  ; Physics and cosmology
  nu_emitted = 1420.4057517*10^6.  ; Hz
  speed_of_light = 3.*10^8.  ; meters/sec
  k_B = 1.3806503*10^(-23.)  ; Boltzmann constant in J/K
  lambda_emitted = speed_of_light/nu_emitted
  hubble_const = 70000. ; m / (Mpc s)
  hubble_param = hubble_const / 100000.
  OmegaM = 0.27
  Omegak = 0.
  OmegaL = 0.73
  
  
  rmax = 750.  ; meters what the heck is this????
  ground_size = round(rmax*1.2/grid_size)*2.-1.  ; pad it and make it odd
  ground0 = (ground_size+1.)/2.  ; index for (0,0) in the layout array
  umax = 2.*rmax  ; meters
  uv_size = (2.*round(1.2*umax/grid_size) - 1.)  ; number of cells
  
  uv_axis_range=(uv_size-1.)*grid_size
  uv_axis_delta=grid_size
  uv_axis_start=((uv_size-1.)*grid_size/2.)
  uv_axis=(findgen(round(uv_axis_range/uv_axis_delta+1))*uv_axis_delta)-uv_axis_start ; meters
  u0_ind = round(uv_size/2.)
  
  
  ; Observational parameters
  lambda = lambda_emitted * (1. + z_0)
  nu_0 = nu_emitted/(1.+z_0)
  
  band_range=bandwidth
  band_delta=dnu
  band_start=(nu_0-bandwidth/2.)
  band=(findgen(round(band_range/band_delta+1))*band_delta)+band_start ; array for all the frequency values
  
  band_dim=size(band)
  ;if ( max(band_dim[1:band_dim[0]]) mod 2 ) EQ 0 THEN BEGIN
  ;  ; even number, make it odd.
  ;  band=(findgen(round(band_range/band_delta+2))*band_delta)+band_start
  ;  ;band[N_elements(band)]=band[N_elements(band)-1.]+dnu
  ;  band=band-dnu/2.
  ;endif
  redshifts = nu_emitted/band - 1.
  
  comoving_dist_los_array=DBLARR(N_elements(redshifts))
  For i=0, N_elements(redshifts)-1 DO BEGIN
    cosmology_measures, redshifts[i], comoving_dist_los=comoving_dist_los2 ; array with distances to the observation (Mpc), redshifts is an array
    comoving_dist_los_array[i]=comoving_dist_los2
  ENDFOR
  comoving_dist_los=DBLARR(N_elements(redshifts))
  comoving_dist_los=comoving_dist_los_array
  
  Dm = mean(comoving_dist_los) ; mean distance to the observation (Mpc), I think this works
  Ez = sqrt(OmegaM*(1.+z_0)^3. + Omegak*(1.+z_0)^2. + OmegaL)  ; about 14 for z=8
  lambdas = lambda_emitted * (1.+redshifts)  ; array with wavelengths per frequency channel
  
  chi_lambda = comoving_dist_los * lambdas  ; A useful array for converting to wavenumber
  mean_chi_lambda=mean(chi_lambda[*])  ;Need to check means
  kx_axis = uv_axis * 2. * !Pi / mean_chi_lambda  ; kx values for each voxel
  ky_axis = kx_axis
  k_perp_bin = kx_axis[1]-kx_axis[0] ; the size of the kperp bins
  dchi = abs(comoving_dist_los[N_elements(comoving_dist_los)-1]-comoving_dist_los[0])/(N_elements(comoving_dist_los)-1)  ; spatial resolution in LoS direction (Mpc)
  maxchi = max(comoving_dist_los[*])-min(comoving_dist_los[*])   ; width of observation in LoS direction (Mpc)
  kz_max = 2.*!Pi / (2.*dchi)  ; factor of 2 for positive and negative
  dkz = 2.*!Pi / maxchi  ; resolution in k_parallel
  
  kz_axis_range=2.*kz_max
  kz_axis_delta=dkz
  kz_axis_start=(kz_max)
  kz_axis=(findgen(round(kz_axis_range/kz_axis_delta+1))*kz_axis_delta)-kz_axis_start ; kz (k_parallel) values for each voxel
  
  inst_perp_axis_range=2.*!Pi/(kx_axis[1]-kx_axis[0])+.01
  inst_perp_axis_delta=(2.*!Pi/(kx_axis[N_elements(kx_axis)-1]-kx_axis[0]))
  inst_perp_axis_start=!Pi/(kx_axis[1]-kx_axis[0])
  inst_perp_axis=(findgen(round(inst_perp_axis_range/inst_perp_axis_delta+1))*inst_perp_axis_delta)-inst_perp_axis_start
  
  ;inst_par_axis_range=(((!Pi+.01)/(kz_axis[1]-kz_axis[0]))+(!Pi/(kz_axis[1]-kz_axis[0]))) I think Adam made a mistake. I made the .01 look like the inst_perp_axis, and now I get the correct number of indices (384 instead of 385 for normal res)
  inst_par_axis_range=(((!Pi)/(kz_axis[1]-kz_axis[0])+.01)+(!Pi/(kz_axis[1]-kz_axis[0])))
  inst_par_axis_delta=(2.*!Pi/(kz_axis[N_elements(kz_axis)-1]-kz_axis[0]))
  inst_par_axis_start=(!Pi/(kz_axis[1]-kz_axis[0]))
  inst_par_axis=(findgen(round(inst_par_axis_range/inst_par_axis_delta+1))*inst_par_axis_delta)-inst_par_axis_start
  
  undefine,sample_cube,sample_cube_ideal,r_ground,n_k_ideal,ground_x,ground_y,ground_ideal,n_k
  
  ;toc
  
  t1=Systime(1)
  ; Code to manipulate cubes to get them in same units and such
  print, (t1-t0)/60., ' minutes'
  print, 'Manipulating cubes...'
  t0=Systime(1)
  ;tic
  
  ; cut and interp
  
  ; first cut in z-direction
  ;z_width = 2*round(max(inst_par_axis)/(r_sim_axis(2)-r_sim_axis(1)))+11;
  z_width = 2.*round(max(inst_par_axis)/(r_sim_axis[1]-r_sim_axis[0]))+2.
  
  ;z_width = 121; % roughly number of sim pixels to reach inst extent
  z_start = 1. ; the offset for slice to use. can change this for independent simulations
  z_image_cut = (findgen(round((z_width-1.)+1.)))+(z_start)
  
  ; try this...
  z_image_cut1 = z_image_cut
  sim_size = size(sim_image)
  sim_length = sim_size[1] ;Find the number of rows in sim_length, a 3D matrix, could be wrong here
  ntile=14. ;Gets to roughly the right size
  ;ntile=7.
  sim_image_scaled = fltarr(ntile*sim_length,ntile*sim_length,N_elements(z_image_cut))
  ;sim_image_scaled = make_array(ntile*sim_length,ntile*sim_length,N_elements(z_image_cut), VALUE=0)
  
  ; tiling the sim in a hacky way
  for i=1, ntile DO BEGIN
    for j=1, ntile DO BEGIN
      sim_image_scaled[(((i-1)*sim_length+1)-1):((i*sim_length)-1),(((j-1)*sim_length+1)-1):((j*sim_length)-1),*] = sim_image[*,*,z_image_cut-1]
    endfor
  endfor
  ; sim_image_scaled((sim_length+1):2*sim_length,1:sim_length,:) = sim_image(:,:,z_image_cut);
  ; sim_image_scaled(1:sim_length,(sim_length+1):2*sim_length,:) = sim_image(:,:,z_image_cut);
  ; sim_image_scaled((sim_length+1):2*sim_length,(sim_length+1):2*sim_length,:) = sim_image(:,:,z_image_cut);
  
  
  ;sim_image_scaled = sim_image(:,:,z_image_cut);
  
  z_sim_extent = z_width*(r_sim_axis[1]-r_sim_axis[0])
  z_sim_axis_range=(z_sim_extent-z_sim_extent/z_width)
  z_sim_axis_delta=(z_sim_extent/z_width)
  z_sim_axis_start=(z_sim_extent-z_sim_extent/z_width)/2.
  z_sim_axis=(findgen(round(z_sim_axis_range/z_sim_axis_delta+1))*z_sim_axis_delta)-z_sim_axis_start
  
  ; Create meshgrids
  r_sim_axis2_range=((ntile*max(r_sim_axis)+(ntile-1.)*dr)-(ntile*min(r_sim_axis)))
  r_sim_axis2_delta=dr
  r_sim_axis2_start=(ntile*min(r_sim_axis))
  r_sim_axis2=(findgen(round(r_sim_axis2_range/r_sim_axis2_delta+1))*r_sim_axis2_delta)+r_sim_axis2_start
  
  ;[x_sim_cube,y_sim_cube,z_sim_cube] = meshgrid(r_sim_axis2,r_sim_axis2,z_sim_axis) ;Let's just say this works out...
  t1=Systime(1)
  print, (t1-t0)/60., ' minutes'
  ;print, "First Mesh routine..."
  ;t0=Systime(1)
  
  ;x_sim_cube=dblarr(N_elements(r_sim_axis2),N_elements(r_sim_axis2),N_elements(z_sim_axis))
  ;y_sim_cube=dblarr(N_elements(r_sim_axis2),N_elements(r_sim_axis2),N_elements(z_sim_axis))
  ;z_sim_cube=dblarr(N_elements(r_sim_axis2),N_elements(r_sim_axis2),N_elements(z_sim_axis))
  ;x_sim_cube_xymesh = (meshgrid(N_elements(r_sim_axis2),N_elements(r_sim_axis2),1)*dr)+(ntile*min(r_sim_axis))
  ;y_sim_cube_xymesh = (meshgrid(N_elements(r_sim_axis2),N_elements(r_sim_axis2),2)*dr)+(ntile*min(r_sim_axis))
  ;z_sim_cube_xymesh = (meshgrid(N_elements(r_sim_axis2),N_elements(z_sim_axis),2)*(z_sim_extent/z_width))-(z_sim_extent-z_sim_extent/z_width)/2.
  
  ;FOR i=0, N_elements(z_sim_axis)-1 DO BEGIN
  ;  x_sim_cube[*,*,i] = x_sim_cube_xymesh
  ;  y_sim_cube[*,*,i] = y_sim_cube_xymesh
  ;ENDFOR
  ;FOR i=0, N_elements(r_sim_axis2)-1 DO z_sim_cube[i,*,*] = z_sim_cube_xymesh
  
  ;t1=Systime(1)
  ;print, (t1-t0)/60., ' minutes'
  ;print, "Second Mesh routine..."
  ;t0=Systime(1)
  
  ; We're forcing the z-axis to match the instrument, and the perp axis to
  ; match the simulation. for some reason...
  ;x_sim_cube_instz=dblarr(N_elements(r_sim_axis2),N_elements(r_sim_axis2),N_elements(inst_par_axis))
  ;y_sim_cube_instz=dblarr(N_elements(r_sim_axis2),N_elements(r_sim_axis2),N_elements(inst_par_axis))
  ;z_sim_cube_instz=dblarr(N_elements(r_sim_axis2),N_elements(r_sim_axis2),N_elements(inst_par_axis))
  ;x_sim_cube_instz_xymesh = (meshgrid(N_elements(r_sim_axis2),N_elements(r_sim_axis2),1)*dr)+(ntile*min(r_sim_axis))
  ;y_sim_cube_instz_xymesh = (meshgrid(N_elements(r_sim_axis2),N_elements(r_sim_axis2),2)*dr)+(ntile*min(r_sim_axis))
  ;z_sim_cube_instz_xymesh = (meshgrid(N_elements(r_sim_axis2),N_elements(inst_par_axis),2)*(2.*!Pi/(kz_axis[N_elements(kz_axis)-1]-kz_axis[0])))-(!Pi/(kz_axis[1]-kz_axis[0]))
  
  
  ;t1=Systime(1)
  ;print, (t1-t0)/60., ' minutes'
  ;print, "Filling instataneous sim cubes..."
  ;t0=Systime(1)
  
  ;FOR i=0, N_elements(inst_par_axis)-1 DO BEGIN
  ;  x_sim_cube_instz[*,*,i] = x_sim_cube_xymesh
  ;  y_sim_cube_instz[*,*,i] = y_sim_cube_xymesh
  ;ENDFOR
  ;FOR i=0, N_elements(r_sim_axis2)-1 DO BEGIN
  ;  z_sim_cube_instz[i,*,*] = z_sim_cube_instz_xymesh
  ;  z_sim_cube[i,*,*] = z_sim_cube_xymesh
  ;ENDFOR
  
  ; Inperpolate simulation using a smooth function
  ;sim_image_interp = grid3(x_sim_cube,y_sim_cube,z_sim_cube,sim_image_scaled,x_sim_cube_instz,y_sim_cube_instz,z_sim_cube_instz)
  
  ;t1=Systime(1)
  ;print, (t1-t0)/60., ' minutes'
  print, "Linear interpolation..."
  t0=Systime(1)
  
  xscale=indgen(n_elements(r_sim_axis2))
  zscale=(inst_par_axis-min(z_sim_axis))/(z_sim_axis[1]-z_sim_axis[0])
  sim_image_interp=interpolate(sim_image_scaled,float(zscale),/grid)
  
  ; Don't need all these cubes lying around
  undefine, x_sim_cube,y_sim_cube,z_sim_cube,sim_image_scaled,xi_image_scaled,x_sim_cube_instz,y_sim_cube_instz,z_sim_cube_instz
  
  ; Pad cubes
  drperp = r_sim_axis[1]-r_sim_axis[0]
  dkperp = kx_axis[1]-kx_axis[0]
  
  rperp_axis_range=((round(!Pi/dkperp/drperp)*drperp)-(round(-!Pi/dkperp/drperp)*drperp))
  rperp_axis_delta=drperp
  rperp_axis_start=(round(-!Pi/dkperp/drperp)*drperp)
  rperp_axis=(findgen(round(rperp_axis_range/rperp_axis_delta+1))*rperp_axis_delta)-rperp_axis_start
  
  kperp_axis_range=((round(!Pi/drperp/dkperp)*dkperp)-(round(-!Pi/drperp/dkperp)*dkperp))
  kperp_axis_delta=dkperp
  kperp_axis_start=(round(-!Pi/drperp/dkperp)*dkperp)
  kperp_axis=(findgen(round(kperp_axis_range/kperp_axis_delta+1))*kperp_axis_delta)-kperp_axis_start
  
  kperp_axis_end=kperp_axis[N_elements(kperp_axis)-1]
  kperp_axis_start=kperp_axis[0]
  
  ; These guys are going to kill my memory...Adam's was single precision typecasted, so that might need to be changed.
  
  t1=Systime(1)
  print, (t1-t0)/60., ' minutes'
  print, "Defining padded sim cube..."
  t0=Systime(1)
  
  inst_par_axis_size=N_elements(inst_par_axis)
  rperp_axis_size=N_elements(rperp_axis)
  undefine, sim_image,comoving_dist_los_array,comoving_dist_los,inst_perp_axis,inst_par_axis
  undefine, kz_axis,kx_axis,uv_axis,band
  undefine, x_sim_cube_instz_xymesh,y_sim_cube_instz_xymesh,z_sim_cube_instz_xymesh
  undefine, x_sim_cube_instz,y_sim_cube_instz,z_sim_cube_instz
  undefine, chi_lambda, kperp_axis, ky_axis, lambdas, z_image_cut, z_image_cut1
  undefine, rperp_axis,redshifts, r_sim_axis, r_sim_axis2, xscale, zscale, z_sim_axis
  
  sim_image_interp_size = size(sim_image_interp)
  for i=0, sim_image_interp_size[3]-1 DO BEGIN
    sim_image_padded = complex(fltarr(rperp_axis_size,rperp_axis_size,1))
    
    if rperp_axis_size GT sim_image_interp_size[1] THEN BEGIN
      ; Need to pad sim image
      sim_pad_size = round((rperp_axis_size - sim_image_interp_size[1])/2)
      sim_perp_ind_range=((sim_pad_size+sim_image_interp_size[1])-(sim_pad_size+1.))
      sim_perp_ind_start=(sim_pad_size+1.)
      sim_perp_ind=(findgen(round(sim_perp_ind_range+1)))+sim_perp_ind_start
      sim_image_padded[sim_perp_ind,sim_perp_ind,1] = sim_image_interp[*,*,i]
    endif else begin
      sim_pad_size = round((sim_image_interp_size[1]-rperp_axis_size)/2.)
      sim_perp_ind_range=((sim_pad_size+rperp_axis_size)-(sim_pad_size+1.))
      sim_perp_ind_start=(sim_pad_size+1.)
      sim_perp_ind=(findgen(round(sim_perp_ind_range+1)))+sim_perp_ind_start
      sim_image_padded=sim_image_interp[sim_perp_ind,sim_perp_ind,i]
      
    endelse
    
    
    jacobian = (2.*!Pi/(kperp_axis_end-kperp_axis_start))^2.
    ;fourier_image = jacobian * fftshift(fftn(sim_image));
    
    t1=Systime(1)
    print, (t1-t0)/60., ' minutes'
    print, "Store and fft..."
    t0=Systime(1)
    
    ; Finish things up
    ft_image=jacobian*fft(fft(sim_image_padded,/CENTER,DIMENSION=1),/CENTER,DIMENSION=2)
    ; for now just save the inner most 2048x2048
    ft_image_size=size(ft_image)
    ind0=round(ft_image_size[1]/2)
    ft_image=ft_image[(ind0-1023):(ind0+1024),(ind0-1023):(ind0+1024)]
    
    ; save data
    
    ;f=fopen('nicholes_real_cube.bin','wb');
    ;fwrite(f,real(ft_image),'single');
    ;fclose(f);
    ;f=fopen('nicholes_imag_cube.bin','wb');
    ;fwrite(f,imag(ft_image),'single');
    ;fclose(f);
    save, ft_image, FILENAME='/nfs/eor-00/h1/nbarry/MWA/development/bubbles/cubes/ft_image_half_', i, '.sav'
    
  endfor
  
  t1=Systime(1)
  print, (t1-t0)/60., ' minutes'
;return, ft_image
end
