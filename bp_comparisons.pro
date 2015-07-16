pro bp_comparisons, day,diff=diff, line_plot_multiplier=line_plot_multiplier, cable_type=cable_type,longrun=longrun, lst_compare=lst_compare
  ;Program for making comparisons between B_c and B_cp. B_c is read in from txt files of a B_c run, and
  ;vice-a-versa for B_cp.  B_cp must not have the polyfit taken out beforehand for the comparison
  ;to made.  Currently, this program makes waterfalls of each obsid vs frequency vs difference (either
  ;grouped by pointing or by cable type), with
  ;an optional plotter for the bandpass plots on top of one another for verification.

  ;cable_type is a num 0 through 5 (0:90m,1:150m,etc)


  pol_name_array=['xx','yy']
  outdir='/nfs/eor-00/h1/nbarry/'+day+'_cable_bp/
  
  IF keyword_set(cable_type) then begin
    cable_beg=cable_type
    cable_end=cable_type
  endif else begin
    cable_beg=0
    cable_end=5
  endelse
  
  
  If keyword_set(longrun) then parsednames=[day+'_minustwo',day+'_minusone',day+'_zenith',day+'_plusone',day+'_plustwo',day+'_plusthree'] else  $
    parsednames=[day+'minustwo',day+'minusone',day+'zenith',day+'plusone',day+'plustwo',day+'plusthree']
  ;parsednames=[day+'_minusfour',day+'_minusthree',day+'_minustwo',day+'_minusone',day+'_zenith',day+'_plusone',day+'_plustwo',day+'_plusthree',day+'_plusfour']
    
  If keyword_set(longrun) then pointing_num=[-2,-1,0,1,2] else pointing_num=[-2,-1,0,1,2,3]
  ;pointing_num=[-4,-3,-2,-1,0,1,2,3,4]
  obs_num_per_pointing=INTARR((size(pointing_num))[1])
  If keyword_set(lst_compare) then obs_num_per_pointing_uncut=INTARR((size(pointing_num))[1])
  FOR j=0,(size(pointing_num))[1]-1 DO BEGIN
    If keyword_set(longrun) then filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/' + parsednames[j] + '.txt' else $
      filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/'+ parsednames[j] + '.txt'
      
    obs_temp='empty'
    readcol, filename, obs_temp, format='A'
    obs_num_per_pointing[j]=N_elements(obs_temp)
    If obs_temp[0] NE 'empty' Then begin
      If j NE 0 then obsid=[obsid,obs_temp] else obsid=obs_temp
    endif
    If keyword_set(lst_compare) then begin
      filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/longrun_queries/uncut/' + parsednames[j] + '.txt'
      obs_temp='empty'
      readcol, filename, obs_temp, format='A'
      obs_num_per_pointing[j]=N_elements(obs_temp)
      If obs_temp[0] NE 'empty' Then begin
        If j NE 0 then obsid_uncut=[obsid_uncut,obs_temp] else obsid_uncut=obs_temp
      endif
    endif
    
    
  ENDFOR
  obsid_count=N_elements(obsid)
  If keyword_set(lst_compare) then begin
    obsid_count=N_elements(obsid_uncut)
    for obs_i=0,obsid_count-1 do begin
      IF where(obsid_uncut[obs_i] EQ obsid) EQ -1 THEN obsid_uncut[obs_i]='empty'
    endfor
    obsid=obsid_uncut
  endif
  
  
  If keyword_set(longrun) then pointing_name=['-2','-1','0','1','2'] else pointing_name=['-2','-1','0','1','2','3']
  poi_num=N_elements(pointing_name)
  pointing_ptr=PTRARR(poi_num,2,6,/allocate)
  for poi_i=0,poi_num-1 do begin
    If keyword_set(longrun) then filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_gainrescor_plots/fhd_nb_longrun_predigjump/fhd_nb_longrun_'+day+'/'+pointing_name[poi_i]+'_bandpass.txt' else $
      filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_gainrescor_plots/fhd_nb_reg_'+day+'_nopoly_fromnocable/'+pointing_name[poi_i]+'_bandpass.txt'
    readcol, filename, freq, cable90xx, cable90yy, cable150xx, cable150yy, cable230xx, cable230yy, cable320xx, cable320yy, cable400xx, cable400yy, cable524xx, cable524yy
    *pointing_ptr[poi_i,0,0]=cable90xx
    *pointing_ptr[poi_i,1,0]=cable90yy
    *pointing_ptr[poi_i,0,1]=cable150xx
    *pointing_ptr[poi_i,1,1]=cable150yy
    *pointing_ptr[poi_i,0,2]=cable230xx
    *pointing_ptr[poi_i,1,2]=cable230yy
    *pointing_ptr[poi_i,0,3]=cable320xx
    *pointing_ptr[poi_i,1,3]=cable320yy
    *pointing_ptr[poi_i,0,4]=cable400xx
    *pointing_ptr[poi_i,1,4]=cable400yy
    *pointing_ptr[poi_i,0,5]=cable524xx
    *pointing_ptr[poi_i,1,5]=cable524yy
  endfor
  
  If ~keyword_set(obsid_count) then obsid_count=94
  cable_ptr=PTRARR(obsid_count,2,6,/allocate)
  for obs_i=0,obsid_count-1 do begin
    If obsid[obs_i] NE 'empty' then begin
      filename='/nfs/eor-00/h1/nbarry/'+day+'_cable_bp/calibration/'+strtrim(string(obsid[obs_i]),2)+'_bandpass.txt'
      readcol, filename, freq, cable90xx, cable90yy, cable150xx, cable150yy, cable230xx, cable230yy, cable320xx, cable320yy, cable400xx, cable400yy, cable524xx, cable524yy
    endif else begin
      cable90xx=FLTARR(384)
      cable90xx=FLTARR(384)
      cable90yy=FLTARR(384)
      cable150xx=FLTARR(384)
      cable150yy=FLTARR(384)
      cable230xx=FLTARR(384)
      cable230yy=FLTARR(384)
      cable320xx=FLTARR(384)
      cable320yy=FLTARR(384)
      cable400xx=FLTARR(384)
      cable400yy=FLTARR(384)
      cable524xx=FLTARR(384)
      cable524yy=FLTARR(384)
      cable90xx=FLTARR(384)
    endelse
    *cable_ptr[obs_i,0,0]=cable90xx
    *cable_ptr[obs_i,1,0]=cable90yy
    *cable_ptr[obs_i,0,1]=cable150xx
    *cable_ptr[obs_i,1,1]=cable150yy
    *cable_ptr[obs_i,0,2]=cable230xx
    *cable_ptr[obs_i,1,2]=cable230yy
    *cable_ptr[obs_i,0,3]=cable320xx
    *cable_ptr[obs_i,1,3]=cable320yy
    *cable_ptr[obs_i,0,4]=cable400xx
    *cable_ptr[obs_i,1,4]=cable400yy
    *cable_ptr[obs_i,0,5]=cable524xx
    *cable_ptr[obs_i,1,5]=cable524yy
  endfor
  
  n_cable=6
  
  mode_filepath='/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/instrument_config/mwa_cable_reflection_coefficients.txt'
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  ;Initialize arrays for excluding flagged tiles
  tile_use_cable=PTRARR(n_cable, /allocate)
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_std/metadata/1061316296_obs.sav'
  
  freq_use=where((*obs.baseline_info).freq_use)
  
  data_range_xx=PTRARR(cable_end+1,/allocate)
  for pol_i=0,1 do begin
  
    pointing_allcables=PTRARR(6,cable_end+1,/allocate)
    
    pol_name=pol_name_array[pol_i]
    pol=pol_i
    for cable_i=cable_beg,cable_end do begin
      cable_bp=FLTARR(obsid_count,384)
      pointing_bp=FLTARR(obsid_count,384)
      
      for obs_i=0,obsid_count-1 do begin
        If cable_i eq 0 then cable_name='90'
        If cable_i eq 1 then cable_name='150'
        If cable_i eq 2 then cable_name='230'
        If cable_i eq 3 then cable_name='320'
        If cable_i eq 4 then cable_name='400'
        If cable_i eq 5 then cable_name='524'
        
        cable_bp[obs_i,*]=*cable_ptr[obs_i,pol,cable_i]
        
        If obs_i LT obs_num_per_pointing[0] then pointing_bp[obs_i,*]=*pointing_ptr[0,pol,cable_i]
        If (obs_i GE obs_num_per_pointing[0]) AND (obs_i LT total(obs_num_per_pointing[0:1])) then pointing_bp[obs_i,*]=*pointing_ptr[1,pol,cable_i]
        If (obs_i GE total(obs_num_per_pointing[0:1])) AND (obs_i LT total(obs_num_per_pointing[0:2])) then pointing_bp[obs_i,*]=*pointing_ptr[2,pol,cable_i]
        If (obs_i GE total(obs_num_per_pointing[0:2])) AND (obs_i LT total(obs_num_per_pointing[0:3])) then pointing_bp[obs_i,*]=*pointing_ptr[3,pol,cable_i]
        If (obs_i GE total(obs_num_per_pointing[0:3])) AND (obs_i LT total(obs_num_per_pointing[0:4])) then pointing_bp[obs_i,*]=*pointing_ptr[4,pol,cable_i]
        If ~keyword_set(longrun) AND ((obs_i GE total(obs_num_per_pointing[0:4])) AND (obs_i LT total(obs_num_per_pointing[0:5]))) then pointing_bp[obs_i,*]=*pointing_ptr[5,pol,cable_i]
        
      endfor
      
      x_range=INDGEN((obsid_count+1))
      y_range=INDGEN(384+1)
      line_array=INTARR((obsid_count+1))
      line_array[*]=256
      
      ;pointings
      ;pointingnumbers=[16,15,15,15,15,18]
      
      pointing_line_y=[0,384]
      minustwo_line_x=[obs_num_per_pointing[0],obs_num_per_pointing[0]]
      minusone_line_x=[total(obs_num_per_pointing[0:1]),total(obs_num_per_pointing[0:1])]
      zenith_line_x=[total(obs_num_per_pointing[0:2]),total(obs_num_per_pointing[0:2])]
      plusone_line_x=[total(obs_num_per_pointing[0:3]),total(obs_num_per_pointing[0:3])]
      plustwo_line_x=[total(obs_num_per_pointing[0:4]),total(obs_num_per_pointing[0:4])]
      
      if keyword_set(longrun) then total_obs=total(obs_num_per_pointing[0:4]) else total_obs=total(obs_num_per_pointing[0:5])
      minustwo_num_x=(obs_num_per_pointing[0]/2)/(total_obs*1.9)+.2
      minusone_num_x=(obs_num_per_pointing[1]/2+obs_num_per_pointing[0])/(total_obs*1.9)+.2
      zenith_num_x=(obs_num_per_pointing[2]/2+total(obs_num_per_pointing[0:1]))/(total_obs*1.9)+.2
      plusone_num_x=(obs_num_per_pointing[3]/2+total(obs_num_per_pointing[0:2]))/(total_obs*1.9)+.2
      plustwo_num_x=(obs_num_per_pointing[4]/2+total(obs_num_per_pointing[0:3]))/(total_obs*1.9)+.2
      If ~keyword_set(longrun) then plusthree_num_x=(obs_num_per_pointing[5]/2+total(obs_num_per_pointing[0:4]))/(total_obs*1.9)+.2
      
      for freq_i=0,383 do If cable_bp[0,freq_i] EQ 0. then cable_bp[*,freq_i]=1.
      
      
      half_waterfall=1
      If keyword_set(half_waterfall) then begin
        ;If pol_i EQ 0 then *data_range_xx[cable_i]=round(100.*minmax((pointing_bp[*,freq_use]-cable_bp[*,freq_use])))/100.
        cgPS_Open,outdir+'plots/'+cable_name+'m_waterfall_difference_'+pol_name+'.png',/quiet,/nomatch
        quick_image, (pointing_bp-cable_bp),x_range,y_range, data_range=[-.03,.03], $
          xtitle='Pointing',ytitle='Frequency Index', title=cable_name+'m, (B!Icp!N-B!Ic!N), '+day+', '+pol_name,charsize=1
        cgoplot, line_array, Linestyle=2,thick=3
        cgoplot, minustwo_line_x, pointing_line_y,Linestyle=2
        cgoplot, minusone_line_x, pointing_line_y, Linestyle=2
        cgoplot, zenith_line_x, pointing_line_y, Linestyle=2
        cgoplot, plusone_line_x, pointing_line_y, Linestyle=2
        cgoplot, plustwo_line_x, pointing_line_y, Linestyle=2
        ;cgoplot, plusthree_line_x, pointing_line_y, Linestyle=2
        cgText, .74,.655,/Normal, 'Digital Gain Jump', color='black', charsize=1
        
        cgPolygon, [-5,-5,90,90],[-1,-15,-15,-1], color='white', /FILL, /window
        
        cgText, minustwo_num_x,.17,/Normal, '-2', color='black', charsize=1
        cgText, minusone_num_x,.17,/Normal, '-1', color='black', charsize=1
        cgText, zenith_num_x,.17,/Normal, '0', color='black', charsize=1
        cgText, plusone_num_x,.17,/Normal, '1', color='black', charsize=1
        cgText, plustwo_num_x,.17,/Normal, '2', color='black', charsize=1
        if ~keyword_set(longrun) then cgText, plusthree_num_x,.17,/Normal, '3', color='black', charsize=1
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      endif
      
      ;*pointing_allcables[0,cable_i]=pointing_bp[0:15,freq_use]-cable_bp[0:15,freq_use]
      ;*pointing_allcables[1,cable_i]=pointing_bp[16:30,freq_use]-cable_bp[16:30,freq_use]
      ;*pointing_allcables[2,cable_i]=pointing_bp[31:45,freq_use]-cable_bp[31:45,freq_use]
      ;*pointing_allcables[3,cable_i]=pointing_bp[46:60,freq_use]-cable_bp[46:60,freq_use]
      ;*pointing_allcables[4,cable_i]=pointing_bp[61:75,freq_use]-cable_bp[61:75,freq_use]
      ;*pointing_allcables[5,cable_i]=pointing_bp[76:93,freq_use]-cable_bp[76:93,freq_use]
      
      *pointing_allcables[0,cable_i]=pointing_bp[0:obs_num_per_pointing[0]-1,*]-cable_bp[0:obs_num_per_pointing[0]-1,*]
      *pointing_allcables[1,cable_i]=pointing_bp[obs_num_per_pointing[0]:total(obs_num_per_pointing[0:1])-1,*]-cable_bp[obs_num_per_pointing[0]:total(obs_num_per_pointing[0:1])-1,*]
      *pointing_allcables[2,cable_i]=pointing_bp[total(obs_num_per_pointing[0:1]):total(obs_num_per_pointing[0:2])-1,*]-cable_bp[total(obs_num_per_pointing[0:1]):total(obs_num_per_pointing[0:2])-1,*]
      *pointing_allcables[3,cable_i]=pointing_bp[total(obs_num_per_pointing[0:2]):total(obs_num_per_pointing[0:3])-1,*]-cable_bp[total(obs_num_per_pointing[0:2]):total(obs_num_per_pointing[0:3])-1,*]
      *pointing_allcables[4,cable_i]=pointing_bp[total(obs_num_per_pointing[0:3]):total(obs_num_per_pointing[0:4])-1,*]-cable_bp[total(obs_num_per_pointing[0:3]):total(obs_num_per_pointing[0:4])-1,*]
      if ~keyword_set(longrun) then $
        *pointing_allcables[5,cable_i]=pointing_bp[total(obs_num_per_pointing[0:4]):total(obs_num_per_pointing[0:5])-1,*]-cable_bp[total(obs_num_per_pointing[0:4]):total(obs_num_per_pointing[0:5])-1,*]
        
        
      line_plots=1
      poi_name=['-2','-1','0','1','2','3']
      If keyword_set(line_plots) then begin
        rgbcolors=[[69,190,207],[144,23,3],[240,87,249],[45,165,30],[42,25,77],[1,53,8],[167,138,28],[251,193,236],[52,86,193],[246,229,194],[148,33,145],[253,88,89],[239,233,110],$
          [14,128,63],[96,50,14],[176,247,108],[74,126,157],[47,14,36],[212,120,250],[202,29,168]]
          
        for poi_i=0,poi_num-1 do begin
        
          cgPS_Open,outdir+'lines/'+cable_name+'m_'+poi_name[poi_i]+'_line_'+pol_name+'.png',/quiet,/nomatch
          cgplot, (*pointing_allcables[poi_i,cable_i])[0,*],$
            xtitle='Frequency Index',ytitle='Gain', title=cable_name+'m, pointing '+poi_name[poi_i]+', B!Ic!N vs B!Icp!N, '+pol_name,charsize=1, psym=2, symsize=0.2, yrange=[.8,1.2]
          for obs_i=1,(size(*pointing_allcables[poi_i, cable_i]))[1]-1 do begin
            Device, Decomposed=0
            TVLCT, rgbcolors[0,obs_i], rgbcolors[1,obs_i], rgbcolors[2,obs_i],100
            cgoplot, (*pointing_allcables[poi_i,cable_i])[obs_i,*], color=100B,psym=2, symsize=0.2
            Device, Decomposed=1
          endfor
          cgoplot, *pointing_ptr[poi_i,pol_i,cable_i], color='black',psym=1;, symsize=0.2
          cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
          
        endfor
      ;cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      endif
      
      
    endfor
    
    if keyword_set(half_waterfall) then begin
      poi_name=['-2','-1','0','1','2','3']
      for poi_i=0,poi_num-1 do begin
        current_pointing=[*pointing_allcables[poi_i,0],*pointing_allcables[poi_i,1],*pointing_allcables[poi_i,2],*pointing_allcables[poi_i,3],*pointing_allcables[poi_i,4],*pointing_allcables[poi_i,5]]
        obs_in_pointing=N_elements(*pointing_allcables[poi_i,0])
        cable_line_y=[0,383]
        cable90_line_x=[obs_in_pointing,obs_in_pointing]
        cable150_line_x=[obs_in_pointing*2,obs_in_pointing*2]
        cable230_line_x=[obs_in_pointing*3,obs_in_pointing*3]
        cable320_line_x=[obs_in_pointing*4,obs_in_pointing*4]
        cable400_line_x=[obs_in_pointing*5,obs_in_pointing*5]
        
        
        cgPS_Open,outdir+'plots/waterfall_pointing'+poi_name[poi_i]+'_difference_'+pol_name+'.png',/quiet,/nomatch
        quick_image, current_pointing,x_range,y_range, data_range=[-.03,.03],$
          xtitle='Observations in a pointing, split by cable',ytitle='Frequency Index', title=poi_name[poi_i]+' pointing, (B!Icp!N-B!Ic!N), '+day+' '+pol_name,charsize=1
        cgoplot, line_array, Linestyle=2,thick=3
        cgoplot, cable90_line_x, cable_line_y,Linestyle=2
        cgoplot, cable150_line_x, cable_line_y, Linestyle=2
        cgoplot, cable230_line_x, cable_line_y, Linestyle=2
        cgoplot, cable320_line_x, cable_line_y, Linestyle=2
        cgoplot, cable400_line_x, cable_line_y, Linestyle=2
        ;cgoplot, plusthree_line_x, pointing_line_y, Linestyle=2
        cgText, .74,.655,/Normal, 'Digital Gain Jump', color='black', charsize=1
        ;cgPolygon, [.1,.7,.1,.7],[.15,.15,.17,.17], color='black', /FILL, /device,/window
        cgPolygon, [-5,-5,90,90],[-1,-15,-15,-1], color='white', /FILL, /window
        cgText, .235,.17,/Normal, '90m', color='black', charsize=1
        cgText, .325,.17,/Normal, '150m', color='black', charsize=1
        cgText, .415,.17,/Normal, '230m', color='black', charsize=1
        cgText, .5,.17,/Normal, '320m', color='black', charsize=1
        cgText, .59,.17,/Normal, '400m', color='black', charsize=1
        cgText, .68,.17,/Normal, '524m', color='black', charsize=1
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      endfor
    endif
    
  endfor
  
  
  
end