pro adams_get_fourth_line_fits_lowband,dir=dir,fits=fits,obsids=obsids,deep2_inds=deep2_inds,cal_cut_inds=cal_cut_inds,dir2=dir2, overplot=overplot, phase=phase, $
    ninety_only=ninety_only, onefifty_only=onefifty_only

  dir = '/nfs/mwa-09/r1/abrahamn/128T/eor0low/fhd_arn_percablebandpass_patticatalog/calibration/'

  if keyword_set(dir2) then dir2 = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/'

  
  ;Things to set per run
  save_path=dir+'mode_params_stdcompare/'

  text_title=' tile, M!Ilow!N'
  
  ;if keyword_set(overplot) then text_title=' tile, M!Io,std!N and M!Ia,phasesep,extraremove,crosspoly!N'
  if keyword_set(overplot) then text_title=' tile, M!Ihigh!N and M!Ilow!N'
  
  cal_files_text='_cal.sav'
  if keyword_set(dir2) then cal_files2_text='_cal.sav'
  ;cal_dir_added=1
  
  
  cal_files=findfile(dir+'*_cal.sav') ; use this as a proxy to get the obsids
  obsids=file_basename(cal_files,'_cal.sav')
  
  obsids_name=obsids
  obsids=ulong(obsids)
  nobs=n_elements(obsids)
  ntile=128
  npol=2
  
  fits = fltarr(npol,nobs,ntile,3) ; last index for mode_i,amp,phase

  for obs_i=0,nobs-1 do begin
    print,obs_i
    cal_files=dir+obsids_name+cal_files_text
    cal = getvar_savefile(cal_files[obs_i],'cal')
    for pol_i=0,npol-1 do begin
      for tile_i=0,ntile-1 do begin
        if (cal.mode_params[pol_i,tile_i] ne !null) then begin
          fits[pol_i,obs_i,tile_i,*] = (*cal.mode_params[pol_i,tile_i])
        endif
      endfor
    endfor
  endfor
  file_mkdir, save_path, /NOEXPAND_PATH
  save,filename=save_path+'fits_save.sav',fits

  
  if keyword_set(dir2) then begin
    fits2 = fltarr(npol,nobs,ntile,3) ; last index for mode_i,amp,phase
    cal_files2=dir2+obsids_name+cal_files2_text
    If keyword_set(cal_dir_added) then cal_files2=dir2+'calibration/'+obsids_name+cal_files2_text
    for obs_i=0,nobs-1 do begin
      print,obs_i
      cal = getvar_savefile(cal_files2[obs_i],'cal')
      for pol_i=0,npol-1 do begin
        for tile_i=0,ntile-1 do begin
          if (cal.mode_params[pol_i,tile_i] ne !null) then begin
            fits2[pol_i,obs_i,tile_i,*] = (*cal.mode_params[pol_i,tile_i])
          endif
        endfor
      endfor
    endfor
    
  endif ;end dir2 if
  
  
  ; get tile info
  restore, '/nfs/eor-03/r1/EoR2013/fhd_nb_std_test_twopolyquad/metadata/1061316296_obs.sav'
  tile_names=(*obs.baseline_info).tile_names
  tile_inds=where(fits[0,0,*,0] ne 0)
  
  If keyword_set(ninety_only) then begin
    mode_filepath=filepath(obs.instrument+'_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
    textfast,data_array,/read,file_path=mode_filepath,first_line=1
    cable_len=Reform(data_array[2,*])
    tile_cable_inds=where(cable_len EQ 90)
    match, tile_inds, tile_cable_inds, suba
    tile_inds=tile_inds[suba]
  endif
  If keyword_set(onefifty_only) then begin
    mode_filepath=filepath(obs.instrument+'_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
    textfast,data_array,/read,file_path=mode_filepath,first_line=1
    cable_len=Reform(data_array[2,*])
    tile_cable_inds=where(cable_len EQ 150)
    match, tile_inds, tile_cable_inds, suba
    tile_inds=tile_inds[suba]
  endif
  
  ; organize the fits info I actually want
  ;fits=pol,obs,tile,fitparams
  fits_full=fits[*,*,tile_inds,1]*exp(Complex(0,1)*fits[*,*,tile_inds,2]); - Complex(0,1)*2.*!Pi*fits[*,*,tile_inds,0]*findgen(384)/384.)
  if keyword_set(dir2) then fits_full2=fits2[*,*,tile_inds,1]*exp(Complex(0,1)*fits2[*,*,tile_inds,2])
  
  
  nt=n_elements(tile_inds)
  print,'Making plots, saving to '+save_path
  if keyword_set(dir2) and ~keyword_set(overplot) then begin
    phase_fits_full=atan(fits_full2,/phase)/atan(fits_full,/phase)
    fits_full=fits_full2/fits_full
    fits=fits2/fits
  endif else begin
    phase_fits_full=atan(fits_full,/phase)
    if keyword_set(overplot) then phase_fits_full2=atan(fits_full2,/phase)
  endelse
  
  for tile_i=0,nt-1 do begin
    print,tile_i
    
    
    
    filename=save_path+'tile'+number_formatter(tile_names[tile_inds[tile_i]])+'.png'
    cgPS_Open,filename,scale_factor=2,/quiet,/nomatch
    if keyword_set(dir2) then begin
      maxabs2=FLTARR(2)
      maxabs2[0]=.04;max(abs(fits_full2))
      maxabs2[1]=.04;max(abs(fits_full))
      maxabs=max(maxabs2)
      maxphase2=FLTARR(2)
      maxphase2[0]=max(abs(phase_fits_full2))
      maxphase2[1]=max(abs(phase_fits_full))
      maxphase=max(maxphase2)
    endif else begin
      maxabs=.04;max(abs(fits_full))
      maxphase=max(abs(phase_fits_full))
    endelse
    
    cgplot,abs(fits_full[0,*,tile_i]),position=[.05,.65,.95,.90],color='red',linestyle=0,title='magnitude',charsize=.5,$
      thick=0,yrange=[0,maxabs],XTICKFORMAT="(A1)";,Psym=2
    cgplot,abs(fits_full[1,*,tile_i]),color='blue',linestyle=0,/overplot,thick=0;,Psym=2
    if keyword_set(overplot) then begin
      cgoplot,abs(fits_full2[0,*,tile_i]),position=[.05,.65,.95,.90],color='red',linestyle=6,title='magnitude',charsize=.5,$
        thick=0,Psym=2
      cgoplot,abs(fits_full2[1,*,tile_i]),color='blue',linestyle=6,thick=0,Psym=2
    endif
    

    ;cgLegend, Title=['xx','yy'], Color=['red','blue'],Length=0.03,Location=[0.85,0.87], charsize=.6
    for day=0,n_elements(date_inds_cal)-1 do begin
      cgplot,[date_inds_cal[day],date_inds_cal[day]],[0,maxabs],linestyle=1,/overplot,thick=0
    endfor
    
    if ~keyword_set(phase) then begin
      cgplot,real_part(fits_full[0,*,tile_i]),position=[.05,.35,.95,.60],/noerase,color='red',linestyle=0,charsize=.5,$
        thick=0,title='real and imaginary parts',yrange=[-maxabs,maxabs],XTICKFORMAT="(A1)";,Psym=2
      cgplot,imaginary(fits_full[0,*,tile_i]),color='blue',linestyle=0,/overplot,thick=0;,Psym=2
    endif else begin
      cgplot,phase_fits_full[0,*,tile_i],position=[.05,.35,.95,.60],/noerase,color='red',linestyle=0,charsize=.5,$
        thick=0,title='phase',yrange=[-maxphase,maxphase],XTICKFORMAT="(A1)";,Psym=2
      cgplot,phase_fits_full[1,*,tile_i],color='blue',linestyle=0,/overplot,thick=0;,Psym=2
    endelse
    
    if keyword_set(overplot) then begin
      if keyword_set(phase) then begin
        cgoplot,phase_fits_full2[0,*,tile_i],position=[.05,.65,.95,.90],color='red',linestyle=6,title='magnitude',charsize=.5,$
          thick=0,Psym=2
        cgoplot,phase_fits_full2[1,*,tile_i],color='blue',linestyle=6,thick=0,Psym=2
      endif else begin
        cgoplot,real_part(fits_full2[0,*,tile_i]),position=[.05,.35,.95,.60],/noerase,color='red',linestyle=6,charsize=.5,$
          thick=0,Psym=2
        cgoplot,imaginary(fits_full2[0,*,tile_i]),color='blue',linestyle=6,thick=0,Psym=2
      endelse
    endif
   
    ;if keyword_set(phase) then cgLegend, Title=['xx','yy'], Color=['red','blue'],Length=0.03,Location=[0.85,0.57], charsize=.6 $
    ;else cgLegend, Title=['real','imag'], Color=['red','blue'],Length=0.03,Location=[0.85,0.57], charsize=.6
    for day=0,n_elements(date_inds_cal)-1 do begin
      cgplot,[date_inds_cal[day],date_inds_cal[day]],[-maxabs,maxabs],linestyle=1,/overplot,thick=0
    endfor
    
    ind=where(fits[*,*,tile_inds,0] ne 0)
    if keyword_set(dir2) and ~keyword_set(overplot) then moderange=[.9,1.1] else moderange=[37,39]
    If keyword_set(ninety_only) then moderange=[21,24]
    
    cgplot,fits[0,*,tile_inds[tile_i],0],position=[.05,.05,.95,.3],/noerase,color='red',linestyle=0,charsize=.5,$
      thick=0,title='mode',yrange=moderange,XTICKFORMAT="(A1)";,Psym=2
    cgplot,fits[1,*,tile_inds[tile_i],0],color='blue',linestyle=0,/overplot,thick=0;,Psym=2
    
    if keyword_set(overplot) then begin
      cgoplot,fits2[0,*,tile_inds[tile_i],0],position=[.05,.65,.95,.90],color='red',linestyle=6,title='magnitude',charsize=.5,$
        thick=0,Psym=2
      cgoplot,fits2[1,*,tile_inds[tile_i],0],color='blue',linestyle=6,thick=0,Psym=2
    endif
 
    ;cgLegend, Title=['xx auto influenced','yy auto influenced'],$
    ;  Color=['red','blue'],Location=[0.65,0.1], charsize=.6, psym=[3,3],VSpace=.9,Length=0.03, $
    ;  linestyle=[0,0]
    ;If keyword_set(dir2) then cgLegend, Title=['xx original', 'yy original'],$
    ;  Color=['red','blue'],Location=[0.85,0.1], charsize=.6, psym=[2,2],VSpace=.9,Length=0.0, $
    ;  /Center_Sym
    
    cgLegend, Title=['xx','yy'],$
      Color=['red','blue'],Location=[0.65,0.1], charsize=.6, psym=[3,3],VSpace=.9,Length=0.03, $
      linestyle=[0,0]
    If keyword_set(dir2) then cgLegend, Title=['xx orig', 'yy orig'],$
      Color=['red','blue'],Location=[0.85,0.1], charsize=.6, psym=[2,2],VSpace=.9,Length=0.0, $
      /Center_Sym
      
    ;text_title=number_formatter(tile_names[tile_inds[tile_i]])+' tile, (M!Io!N/M!Ipost,p!N)'
    ;text_title=number_formatter(tile_names[tile_inds[tile_i]])+' tile, pre digital gain jump pointing modefit'
    ;cgText, 0.27,0.95,number_formatter(tile_names[tile_inds[tile_i]])+' tile, all frequencies, obs modefit', /Normal, charsize=1.2
    ;cgText, 0.24,0.95,number_formatter(tile_names[tile_inds[tile_i]])+' tile, pre digital gain jump pointing modefit', /Normal, charsize=1.2
    cgText, 0.5,0.95,number_formatter(tile_names[tile_inds[tile_i]])+text_title, /Normal, charsize=1.2, Alignment=0.5
    ;(M!Iauto,pre,p!N-M!o!N)
    
    cgPS_Close,/png,Density=75,Resize=100.,/allow_transparent,/nomessage
    
  endfor
end
