PRO longrun_bp_plotter, difference=difference, half_waterfall=half_waterfall, line_plots=line_plots

  pol=0
  pol_name='xx'
  ;outdir='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_gainrescor_plots/fhd_nb_longrun_predigjump/'
  ;outdir='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_gainrescor_plots/fhd_nb_longrun_polyremoved/'
  outdir='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_gainrescor_plots/'
  
  day_name_array=['Aug23','Aug27','Sep02','Sep04','Sep06','Sep09','Sep11','Sep13','Sep17','Sep19','Sep30','Oct02','Oct04','Oct08','Oct10','Oct15','Oct23','Oct25','Oct28','Oct29'];,'Oct31','Nov17','Nov18','Nov29','Nov30']
  day_num=N_Elements(day_name_array)
  ;day num, pol, cable
  ;cable_ptr=PTRARR(26,2,6,\allocate)
  cable_ptr=PTRARR(day_num,9,2,6,/allocate) ;day, pointing, pol, cable
  pointing_num=['-4','-3','-2','-1','0','1','2','3','4']
  poi_num =N_elements(pointing_num)
  
  restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_std/metadata/1061316296_obs.sav'
  freq_use=where((*obs.baseline_info).freq_use)
  
  for poi_i=0,poi_num-1 do begin
    for i_path=0,day_num-1 do begin
      dir_path=outdir+'fhd_nb_longrun_'+day_name_array[i_path]+'/'
      file_test_num=file_test(dir_path+pointing_num[poi_i]+'_bandpass.txt')
      If file_test_num eq 1 then begin
        readcol,dir_path+pointing_num[poi_i]+'_bandpass.txt',freq_arr,cable90xx,cable90yy,cable150xx,cable150yy,cable230xx,cable230yy,cable320xx,cable320yy,cable400xx,cable400yy,cable524xx,cable524yy
      endif else begin
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
        cable90xx[*]=10.
        cable90yy[*]=10.
        cable150xx[*]=10.
        cable150yy[*]=10.
        cable230xx[*]=10.
        cable230yy[*]=10.
        cable320xx[*]=10.
        cable320yy[*]=10.
        cable400xx[*]=10.
        cable400yy[*]=10.
        cable524xx[*]=10.
        cable524yy[*]=10.
      endelse
      *cable_ptr[i_path,poi_i,0,0]=cable90xx
      *cable_ptr[i_path,poi_i,1,0]=cable90yy
      *cable_ptr[i_path,poi_i,0,1]=cable150xx
      *cable_ptr[i_path,poi_i,1,1]=cable150yy
      *cable_ptr[i_path,poi_i,0,2]=cable230xx
      *cable_ptr[i_path,poi_i,1,2]=cable230yy
      *cable_ptr[i_path,poi_i,0,3]=cable320xx
      *cable_ptr[i_path,poi_i,1,3]=cable320yy
      *cable_ptr[i_path,poi_i,0,4]=cable400xx
      *cable_ptr[i_path,poi_i,1,4]=cable400yy
      *cable_ptr[i_path,poi_i,0,5]=cable524xx
      *cable_ptr[i_path,poi_i,1,5]=cable524yy
    endfor
  endfor
  
  mean_cable=FLTARR(poi_num,2,6,384)
  for cable_i=0, 5 do begin
    for poi_i=0,poi_num-1 do begin
      for pol_i=0,1 do begin
        cable_temp=FLTARR(day_num,384)
        for day_i=0,day_num-1 do begin
        
          cable_temp[day_i,*]=*cable_ptr[day_i,poi_i,pol_i,cable_i]
          
        endfor
        for freq_i=0,383 do begin
          cable_mean_temp=cable_temp[*,freq_i]
          mean_cable[poi_i,pol_i,cable_i,freq_i]=MEAN(cable_mean_temp[where(cable_mean_temp NE 10.)])
        ;IF mean_cable[poi_i,pol_i,cable_i,3] EQ 1. AND mean_cable[poi_i,pol_i,cable_i,380] EQ 1. THEN mean_cable[poi_i,pol_i,cable_i,*]=0
        endfor
      endfor
    endfor
  endfor
  
  

  If keyword_set(half_waterfall) then begin
  
    for poi_i=0,poi_num-1 do begin
      for cable_i=0,5 do begin
        cable_bp=FLTARR(day_num,384)
        ;pointing_bp=FLTARR(94,384)
        
        for day_i=0,19 do begin
          If cable_i eq 0 then cable_name='90'
          If cable_i eq 1 then cable_name='150'
          If cable_i eq 2 then cable_name='230'
          If cable_i eq 3 then cable_name='320'
          If cable_i eq 4 then cable_name='400'
          If cable_i eq 5 then cable_name='524'
          
          If keyword_set(difference) then cable_bp[day_i,*]=reform(mean_cable[poi_i,pol,cable_i,*])-*cable_ptr[day_i,poi_i,pol,cable_i] else $
            cable_bp[day_i,*]=*cable_ptr[day_i,poi_i,pol,cable_i]
          IF (*cable_ptr[day_i,poi_i,pol,cable_i])[0] EQ 10. then cable_bp[day_i,*]=0.
        ;IF (day_i EQ 5 OR day_i EQ 11) AND cable_i EQ 0 AND poi_i EQ 2 THEN stop
        ;IF cable_i EQ 0 AND poi_i EQ 2 THEN stop
          
        ;If obs_i LT 16 then pointing_bp[obs_i,*]=*pointing_ptr[0,1,cable_i]
        ;If (obs_i GE 16) AND (obs_i LT 31) then pointing_bp[obs_i,*]=*pointing_ptr[1,1,cable_i]
        ;If (obs_i GE 31) AND (obs_i LT 46) then pointing_bp[obs_i,*]=*pointing_ptr[2,1,cable_i]
        ;If (obs_i GE 46) AND (obs_i LT 61) then pointing_bp[obs_i,*]=*pointing_ptr[3,1,cable_i]
        ;If (obs_i GE 61) AND (obs_i LT 76) then pointing_bp[obs_i,*]=*pointing_ptr[4,1,cable_i]
        ;If obs_i GE 76 then pointing_bp[obs_i,*]=*pointing_ptr[5,1,cable_i]
        endfor
        
        x_range=INDGEN((day_num+1))
        y_range=INDGEN(N_elements(freq_use)+1)
        line_array=INTARR((day_num+1))
        line_array[*]=223
        
        ;pointings
        ;pointingnumbers=[16,15,15,15,15,18]
        ;pointing_line_y=[0,336]
        ;minustwo_line_x=[16,16]
        ;minusone_line_x=[31,31]
        ;zenith_line_x=[46,46]
        ;plusone_line_x=[61,61]
        ;plustwo_line_x=[76,76]
        ;plusthree_line_x=[94,94]
        
        
        If keyword_set(difference) then begin
          ytitle='Frequency Index (unflagged only)'
          data_range=[-.03,.03]
          title=cable_name+'m, pointing '+pointing_num[poi_i]+', B!Icp(average)!N-B!Icp!N, '+pol_name
          save_file='waterfalls/'+cable_name+'m_'+pointing_num[poi_i]+'_waterfall_difference_'+pol_name+'.png'
        endif else begin
          ytitle='Frequency Index (unflagged only)'
          data_range=[.8,1.2]
          title=cable_name+'m, pointing '+pointing_num[poi_i]+', B!Icp!N, '+pol_name
          save_file='waterfalls/'+cable_name+'m_'+pointing_num[poi_i]+'_waterfall_'+pol_name+'.png'
        endelse
        
        make_fft=1
        ;******************************Taking FFTs of bandpass fit block
        IF Keyword_set(make_fft) THEN BEGIN
          for day_i=0,19 do begin
            cable_bp[day_i,*] = shift(384.*ABS(FFT(cable_bp[day_i,*])),384/2)
            ;cable_bp[day_i,*] = 384.*ABS(FFT(cable_bp[day_i,*]))
          endfor
          ytitle='FFT of Frequency'
          data_range=[0,.1]
          title=cable_name+'m, pointing '+pointing_num[poi_i]+', FFT of B!Icp(average)!N-B!Icp!N, '+pol_name
          save_file='waterfalls/'+cable_name+'m_'+pointing_num[poi_i]+'_waterfall_difference_fft_'+pol_name+'.png'
        ENDIF
        ;******************************End of FFT of bandpass fit block
        
        
        cgPS_Open,outdir+save_file,/quiet,/nomatch
        quick_image, cable_bp[*,freq_use],x_range,y_range,data_range=data_range, $
          xtitle='Day Index',ytitle=ytitle, title=title,charsize=1
        If ~keyword_set(make_fft) then cgoplot, line_array, Linestyle=2
        ;cgoplot, minustwo_line_x, pointing_line_y,Linestyle=2
        ;cgoplot, minusone_line_x, pointing_line_y, Linestyle=2
        ;cgoplot, zenith_line_x, pointing_line_y, Linestyle=2
        ;cgoplot, plusone_line_x, pointing_line_y, Linestyle=2
        ;cgoplot, plustwo_line_x, pointing_line_y, Linestyle=2
        ;cgoplot, plusthree_line_x, pointing_line_y, Linestyle=2
        If ~keyword_set(make_fft) then cgText, .74,.655,/Normal, 'Digital Gain Jump', color='black', charsize=1
        ;cgPolygon, [-5,-5,90,90],[-1,-15,-15,-1], color='white', /FILL, /window
        ;cgText, .235,.17,/Normal, '-2', color='black', charsize=1
        ;cgText, .325,.17,/Normal, '-1', color='black', charsize=1
        ;cgText, .415,.17,/Normal, '0', color='black', charsize=1
        ;cgText, .5,.17,/Normal, '1', color='black', charsize=1
        ;cgText, .59,.17,/Normal, '2', color='black', charsize=1
        ;cgText, .68,.17,/Normal, '3', color='black', charsize=1
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
      endfor
    endfor
    
  endif
  
  

  If keyword_set(line_plots) then begin
  
    for poi_i=0,poi_num-1 do begin
      for cable_i=0,5 do begin
        cable_bp=FLTARR(day_num,384)
        ;pointing_bp=FLTARR(94,384)
        
        for day_i=0,19 do begin
          If cable_i eq 0 then cable_name='90'
          If cable_i eq 1 then cable_name='150'
          If cable_i eq 2 then cable_name='230'
          If cable_i eq 3 then cable_name='320'
          If cable_i eq 4 then cable_name='400'
          If cable_i eq 5 then cable_name='524'
          
          If keyword_set(difference) then cable_bp[day_i,*]=reform(mean_cable[poi_i,pol,cable_i,*])-*cable_ptr[day_i,poi_i,pol,cable_i] else $
            cable_bp[day_i,*]=reform(*cable_ptr[day_i,poi_i,pol,cable_i])
          ;cable_bp[day_i,*]=*cable_ptr[day_i,poi_i,pol,cable_i]
          IF (*cable_ptr[day_i,poi_i,pol,cable_i])[0] EQ 10. then cable_bp[day_i,*]=0.
          
        endfor
        
        x_range=INDGEN((day_num+1))
        y_range=INDGEN(N_elements(freq_use)+1)
        line_array=INTARR((day_num+1))
        line_array[*]=223
        
        rgbcolors=[[69,190,207],[144,23,3],[240,87,249],[45,165,30],[42,25,77],[1,53,8],[167,138,28],[251,193,236],[52,86,193],[246,229,194],[148,33,145],[253,88,89],[239,233,110],$
          [14,128,63],[96,50,14],[176,247,108],[74,126,157],[47,14,36],[212,120,250],[202,29,168]]
          
        If keyword_set(difference) then begin
          save_file='waterfalls/'+cable_name+'m_'+pointing_num[poi_i]+'_line_difference_'+pol_name+'.png'
          yrange=[-.05,.05]
          title=cable_name+'m, pointing '+pointing_num[poi_i]+', B!Icp(average)!N-B!Icp!N, '+pol_name
        endif  else begin
          save_file='waterfalls/'+cable_name+'m_'+pointing_num[poi_i]+'_line_'+pol_name+'.png'
          yrange=[.8,1.2]
          title=cable_name+'m, pointing '+pointing_num[poi_i]+', B!Icp!N, '+pol_name
        endelse
        
        
        cgPS_Open,outdir+save_file,/quiet,/nomatch
        cgplot, cable_bp[0,*],yrange=yrange,$
          xtitle='Frequency Index',ytitle='Gain', title=title,charsize=1, psym=2, symsize=0.2
        for day_i=1,day_num-1 do begin
          Device, Decomposed=0
          TVLCT, rgbcolors[0,day_i], rgbcolors[1,day_i], rgbcolors[2,day_i],100
          cgoplot, cable_bp[day_i,*], color=100B,psym=2, symsize=0.2
          Device, Decomposed=1
        endfor
        
        ;cgoplot, line_array, Linestyle=2
        ;cgoplot, minustwo_line_x, pointing_line_y,Linestyle=2
        ;cgoplot, minusone_line_x, pointing_line_y, Linestyle=2
        ;cgoplot, zenith_line_x, pointing_line_y, Linestyle=2
        ;cgoplot, plusone_line_x, pointing_line_y, Linestyle=2
        ;cgoplot, plustwo_line_x, pointing_line_y, Linestyle=2
        ;cgoplot, plusthree_line_x, pointing_line_y, Linestyle=2
        ;cgText, .74,.655,/Normal, 'Digital Gain Jump', color='black', charsize=1
        ;cgPolygon, [-5,-5,90,90],[-1,-15,-15,-1], color='white', /FILL, /window
        ;cgText, .235,.17,/Normal, '-2', color='black', charsize=1
        ;cgText, .325,.17,/Normal, '-1', color='black', charsize=1
        ;cgText, .415,.17,/Normal, '0', color='black', charsize=1
        ;cgText, .5,.17,/Normal, '1', color='black', charsize=1
        ;cgText, .59,.17,/Normal, '2', color='black', charsize=1
        ;cgText, .68,.17,/Normal, '3', color='black', charsize=1
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
      endfor
    endfor
    
    
  endif
  
  
  
END