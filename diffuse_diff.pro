pro diffuse_diff, old_diffuse_filepath,new_diffuse_filepath,store_path,old_diffuse_filepath2=old_diffuse_filepath2,new_diffuse_filepath2=new_diffuse_filepath2
  ;diffuse_diff, 'EoR0_diffuse_2048_4096.sav','EoR0_diffuse_2048_2048.sav',old_diffuse_filepath2='EoR0_diffuse_1024_4096.sav',new_diffuse_filepath2='EoR0_diffuse_1024_2048.sav'


  If ~keyword_set(store_path) then begin
    store_path = 'diffuse_diff_July9.sav'
    print, 'Will save as diffuse_diff_July9.sav'
  endif
  
  old_diffuse_filepath='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/'+old_diffuse_filepath
  new_diffuse_filepath='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/'+new_diffuse_filepath
  
  If keyword_set(old_diffuse_filepath2) OR keyword_set(new_diffuse_filepath2) then begin
    old_diffuse_filepath2='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/'+old_diffuse_filepath2
    new_diffuse_filepath2='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/'+new_diffuse_filepath2
    old_diffuse2 = getvar_savefile(old_diffuse_filepath2, 'model_arr')
    new_diffuse2 = getvar_savefile(new_diffuse_filepath2, 'model_arr')
  endif
  
  
  old_diffuse = getvar_savefile(old_diffuse_filepath, 'model_arr')
  old_diffuse_nside = getvar_savefile(old_diffuse_filepath, 'nside')
  old_diffuse_hpx = getvar_savefile(old_diffuse_filepath, 'hpx_inds')
  
  new_diffuse = getvar_savefile(new_diffuse_filepath, 'model_arr')
  new_diffuse_nside = getvar_savefile(new_diffuse_filepath, 'nside')
  new_diffuse_hpx = getvar_savefile(new_diffuse_filepath, 'hpx_inds')
  
  If keyword_set(old_diffuse_filepath2) OR keyword_set(new_diffuse_filepath2) then begin
  
    binsize=.001
    
    ;result=histogram(*old_diffuse[0]/2.,binsize=binsize,locations=locations,omax=omax)
    result=histogram(*old_diffuse[0],binsize=binsize,locations=locations,omax=omax)
    y_arr_old=[0, result, 0]
    x_arr_old=[locations[0]-binsize/2,locations+binsize/2,locations[N_elements(locations)-1]+binsize]
    
    binCenters = locations + (binsize / 2.0)
    yfit = GaussFit(binCenters, result, coeff_old, NTERMS=3) ;max value, center, standard deviation
    coeff_old[1]=Round(coeff_old[1]*1000.)/1000.
    coeff_old[2]=Round(coeff_old[2]*100.)/100.
    
    result=histogram(*old_diffuse2[0],binsize=binsize,locations=locations,omax=omax)
    y_arr_old2=[0, result, 0]
    x_arr_old2=[locations[0]-binsize/2,locations+binsize/2,locations[N_elements(locations)-1]+binsize]
    
    binCenters = locations + (binsize / 2.0)
    yfit = GaussFit(binCenters, result, coeff_old2, NTERMS=3) ;max value, center, standard deviation
    coeff_old2[1]=Round(coeff_old2[1]*1000.)/1000.
    coeff_old2[2]=Round(coeff_old2[2]*100.)/100.
    
    ;result=histogram(*new_diffuse[0]/2./4.,binsize=binsize,locations=locations,omax=omax)
    result=histogram(*new_diffuse[0],binsize=binsize,locations=locations,omax=omax)
    y_arr_new=[0, result, 0]
    x_arr_new=[locations[0]-binsize/2,locations+binsize/2,locations[N_elements(locations)-1]+binsize]
    
    binCenters = locations + (binsize / 2.0)
    yfit = GaussFit(binCenters, result, coeff_new, NTERMS=3) ;max value, center, standard deviation
    coeff_new[1]=Round(coeff_new[1]*1000.)/1000.
    coeff_new[2]=Round(coeff_new[2]*100.)/100.
    
    ;result=histogram(*new_diffuse2[0]/4.,binsize=binsize,locations=locations,omax=omax)
    result=histogram(*new_diffuse2[0],binsize=binsize,locations=locations,omax=omax)
    y_arr_new2=[0, result, 0]
    x_arr_new2=[locations[0]-binsize/2,locations+binsize/2,locations[N_elements(locations)-1]+binsize]
    
    binCenters = locations + (binsize / 2.0)
    yfit = GaussFit(binCenters, result, coeff_new2, NTERMS=3) ;max value, center, standard deviation
    coeff_new2[1]=Round(coeff_new2[1]*1000.)/1000.
    coeff_new2[2]=Round(coeff_new2[2]*100.)/100.
    
    cgPS_Open,'/nfs/eor-00/h1/nbarry/pixel_distribution_psftransfer.png',/quiet,/nomatch
    cgplot, x_arr_old, y_arr_old, psym=10,title='Distribution of Stokes I pixel values',xtitle='Stokes I amplitude (Jy/pixel)', ytitle='Density', charsize=1, /ylog, xrange=[-.2,.2], yrange=[1,1E6]
    cgoplot,x_arr_old2, y_arr_old2, psym=10,color='blue'
    cgoplot,x_arr_new, y_arr_new, psym=10,color='green'
    ;cgoplot,x_arr_new2, y_arr_new2, psym=10,color='pink'
    ;cglegend, title=['dim 2048, nside 4096','dim 2048, nside 2048','dim 1024, nside 4096','dim 1024, nside 2048'], color=['black','blue','green','pink'],length=.05,location=[.155,.85],charsize=1
    cglegend, title=['dim 3072, psf 2048','dim 2048, psf 2048','dim 3072, psf 3072'], color=['black','blue','green'],length=.05,location=[.155,.85],charsize=1
    ;cglegend, title=['$\sigma$='+strtrim(string(coeff_old[2],Format='(F10.2)'),2)+', mean='+strtrim(string(coeff_old[1],Format='(F10.2)'),2), $
    ;'$\sigma$='+strtrim(string(coeff_old2[2],Format='(F10.2)'),2)+', mean='+strtrim(string(coeff_old2[1],Format='(F10.2)'),2) $
    ;,'$\sigma$='+strtrim(string(coeff_new[2],Format='(F10.2)'),2)+', mean='+strtrim(string(coeff_new[1],Format='(F10.2)'),2),$
    ;'$\sigma$='+strtrim(string(coeff_new2[2],Format='(F10.2)'),2)+', mean='+strtrim(string(coeff_new2[1],Format='(F10.2)'),2)], color=['black','blue','green','pink'] $
    ;,length=.05,location=[.62,.85],charsize=1
    cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    stop
  endif else begin
  
    IF old_diffuse_nside NE new_diffuse_nside then message, 'Different nsides'
    
    model_arr = PTRARR(2,/allocate)
    *model_arr[0] = *old_diffuse[0]-*new_diffuse[0]
    *model_arr[1] = *old_diffuse[1]-*new_diffuse[1]
    
    nside=old_diffuse_nside
    hpx_inds = old_diffuse_hpx
    stop
    save,model_arr,nside,hpx_inds, filename='/nfs/mwa-00/h1/nbarry/MWA/IDL_code/FHD/catalog_data/'+store_path
  endelse
  
end