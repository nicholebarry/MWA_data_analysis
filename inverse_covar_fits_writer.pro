pro inverse_covar_fits_writer

  home_dir='/nfs/mwa-00/h1/nbarry/'
  
  hdr=['COMMENT = Data sum vector as a function of z for kx 87, ky 6.','']
  mwrfits, reform(data_sum[87,6,*]),home_dir+'data_sum_vector_filledpixel.fits',hdr,/create
  
  hdr=['COMMENT = Covariance weighted data sum vector as a function of z for kx 87, ky 6.','']
  mwrfits, reform(wt_data_sum[87,6,*]),home_dir+'wt_data_sum_vector_filledpixel.fits',hdr,/create
  
  hdr=['COMMENT = Covariance matrix (just instrumental contributors) as a function of z for kx 87, ky 6.','']
  mwrfits, reform(covar_z),home_dir+'covar_z_filledpixel.fits',hdr,/create
  
  hdr=['COMMENT = Inverse covariance matrix (just instrumental contributors) as a function of z for kx 87, ky 6.','']
  mwrfits, reform(inv_covar_z),home_dir+'inv_covar_z_filledpixel.fits',hdr,/create
  
  hdr=['COMMENT = Inverse covariance matrix (just instrumental contributors) as a function of kz for kx 87, ky 6.','']
  mwrfits, reform(inv_covar_kz),home_dir+'inv_covar_kz_filledpixel.fits',hdr,/create
  
  hdr=['COMMENT = Covariance weighted sigma squared vector as a function of z for kx 87, ky 6.','']
  mwrfits, reform(wt_sum_sigma2[87,6,*]),home_dir+'wt_sum_sigma2_vector_filledpixel.fits',hdr,/create
  
  hdr=['COMMENT = Matrix with sigma squared as the diagonal as a function of kz for kx 87, ky 6.','']
  mwrfits, reform(var_z_inst),home_dir+'var_z_inst_filledpixel.fits',hdr,/create
  
  hdr=['COMMENT = Fisher matrix as a function of kz for kx 87, ky 6.','']
  mwrfits, reform(fisher_kz),home_dir+'fisher_kz_filledpixel.fits',hdr,/create
  
  hdr=['COMMENT = Square root of the fisher matrix as a function of kz for kx 87, ky 6.','']
  mwrfits, reform(sqrt_fisher_kz),home_dir+'sqrt_fisher_kz_filledpixel.fits',hdr,/create
  
  hdr=['COMMENT = M matrix as a function of kz,kz for kx 87, ky 6.','']
  mwrfits, reform(m_matrix),home_dir+'m_matrix_filledpixel.fits',hdr,/create
  
  hdr=['COMMENT = Weighted sigma squared normalization (diagonal of MFM matrix) as a function of kz,kz for kx 87, ky 6.','']
  mwrfits, diag_matrix(reform(wt_sum_sigma2_norm[87,6,*,*])),home_dir+'wt_sum_sigma2_norm_filledpixel.fits',hdr,/create
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/bandpower1_realpart_filledpixel.png',/quiet,/nomatch
  cgplot, abs(real_part(band_power_1[87,6,*])), title='Real band power term 1, kx kx=87, ky=6',/ylog, color='green', yrange=[10^0.,10^10.],charsize=1, ytitle='Unweighted power', xtitle='k par index'
  cgoplot, real_part(band_power_1[87,6,*])
  cglegend, title=['Positive','Negative'], color=['black','green'],location=[.65,.85], charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/bandpower1_imagpart_filledpixel.png',/quiet,/nomatch
  cgplot, abs(imaginary(band_power_1[87,6,*])), title='Imag band power term 1, kx kx=87, ky=6',/ylog, color='green', yrange=[10^0.,10^10.],charsize=1, ytitle='Unweighted power', xtitle='k par index'
  cgoplot, imaginary(band_power_1[87,6,*])
  cglegend, title=['Positive','Negative'], color=['black','green'],location=[.65,.85], charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/bandweights1_realpart_filledpixel.png',/quiet,/nomatch
  cgplot, abs(real_part(band_weights_1[87,6,*])), title='Real band weights term 1, kx kx=87, ky=6',/ylog, color='green', yrange=[10^0.,10^10.],charsize=1, ytitle='weights', xtitle='k par index'
  cgoplot, real_part(band_weights_1[87,6,*])
  cglegend, title=['Positive','Negative'], color=['black','green'],location=[.65,.85], charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/bandweights1_imagpart_filledpixel.png',/quiet,/nomatch
  cgplot, abs(imaginary(band_weights_1[87,6,*])), title='Imag band weights term 1, kx kx=87, ky=6',/ylog, color='green', yrange=[10^0.,10^10.],charsize=1, ytitle='weights', xtitle='k par index'
  cgoplot, imaginary(band_weights_1[87,6,*])
  cglegend, title=['Positive','Negative'], color=['black','green'],location=[.65,.85], charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/power_realpart_filledpixel.png',/quiet,/nomatch
  cgplot, abs(real_part(power_3d[87,6,*])), title='Real weighted power full, kx kx=87, ky=6',/ylog, color='green',charsize=1, ytitle='weights', xtitle='k par index'
  cgoplot, real_part(power_3d[87,6,*])
  cglegend, title=['Positive','Negative'], color=['black','green'],location=[.65,.85], charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
  cgPS_Open,'/nfs/eor-00/h1/nbarry/power_imagpart_filledpixel.png',/quiet,/nomatch
  cgplot, abs(imaginary(power_3d[87,6,*])), title='Imag weighted power full, kx kx=87, ky=6',/ylog, color='green',charsize=1, ytitle='weights', xtitle='k par index'
  cgoplot, imaginary(power_3d[87,6,*])
  cglegend, title=['Positive','Negative'], color=['black','green'],location=[.65,.85], charsize=1
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  
end