pro plot_all_fits

  poi_name=['-2','-1','0','1','2','3']
  day='Aug23'
  parsednames=[day+'minustwo',day+'minusone',day+'zenith',day+'plusone',day+'plustwo',day+'plusthree']
  
  ;obs ptr is where obs of a pointing are stored. It can be made bigger than needed to account for the pointings of the longrun
  obs_ptr=PTRARR(6,/allocate)
  
  FOR j=0,(size(poi_name))[1]-1 DO BEGIN
  
    filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/' + parsednames[j] + '.txt'
    obs_temp='empty'
    readcol, filename, obs_temp, format='A', /silent
    If j NE 0 then parsednumbers=[parsednumbers,N_elements(obs_temp)] else parsednumbers=N_elements(obs_temp)
    *obs_ptr[j]=obs_temp
    
  ENDFOR
  
  restore,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/metadata/1061311664_obs.sav'
  freq_use=where((*obs.baseline_info).freq_use)
  freq_arr=(*obs.baseline_info).freq[where((*obs.baseline_info).freq_use)]
  tile_names=(*obs.baseline_info).tile_names
  
  mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  rgbcolors=[[128,187,255],[202,18,27],[226,139,176],[229,153,89],[12,111,188]]
  
  
  For j=0,N_elements(poi_name)-1 do begin
  
    normal_res=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    normal_sols=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    normal_input=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    unsplit_quad_res=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    unsplit_quad_sols=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    split_quad_res=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    split_quad_sols=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    split_nodig_quad_res=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    split_nodig_quad_sols=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    split_nodig_quad_fancy_res=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    split_nodig_quad_fancy_sols=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    
    
    bp_sol=complex(FLTARR(384,128,2))
    
    for obs_i=0,N_elements(*obs_ptr[j])-1 do begin
    
      obsid=(*obs_ptr[j])[obs_i]
      
      
      
      filename='/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_cable_cal_pointing_gainrescor_plots/fhd_nb_reg_Aug23_nopoly_fromnocable/'+poi_name[j]+'_bandpass.txt
      readcol, filename, freq_arr_input, cable90xx, cable90yy, cable150xx, cable150yy, cable230xx, cable230yy, cable320xx, cable320yy, cable400xx, cable400yy, cable524xx, cable524yy, /silent
      bandpass_saved_sol=FLTARR(13,384)
      bandpass_saved_sol[0,*]=freq_arr_input
      bandpass_saved_sol[1,*]=cable90xx
      bandpass_saved_sol[2,*]=cable90yy
      bandpass_saved_sol[3,*]=cable150xx
      bandpass_saved_sol[4,*]=cable150yy
      bandpass_saved_sol[5,*]=cable230xx
      bandpass_saved_sol[6,*]=cable230yy
      bandpass_saved_sol[7,*]=cable320xx
      bandpass_saved_sol[8,*]=cable320yy
      bandpass_saved_sol[9,*]=cable400xx
      bandpass_saved_sol[10,*]=cable400yy
      bandpass_saved_sol[11,*]=cable524xx
      bandpass_saved_sol[12,*]=cable524yy
      
      
      for tile_i=0,127 do begin
        if (cable_len[tile_i] EQ 90) then bp_sol[*,tile_i,0]=cable90xx & if (cable_len[tile_i] EQ 90) then bp_sol[*,tile_i,1]=cable90yy
        if (cable_len[tile_i] EQ 150) then bp_sol[*,tile_i,0]=cable150xx & if (cable_len[tile_i] EQ 150) then bp_sol[*,tile_i,1]=cable150yy
        if (cable_len[tile_i] EQ 230) then bp_sol[*,tile_i,0]=cable230xx & if (cable_len[tile_i] EQ 230) then bp_sol[*,tile_i,1]=cable230yy
        if (cable_len[tile_i] EQ 320) then bp_sol[*,tile_i,0]=cable320xx & if (cable_len[tile_i] EQ 320) then bp_sol[*,tile_i,1]=cable320yy
        if (cable_len[tile_i] EQ 400) then bp_sol[*,tile_i,0]=cable400xx & if (cable_len[tile_i] EQ 400) then bp_sol[*,tile_i,1]=cable400yy
        if (cable_len[tile_i] EQ 524) then bp_sol[*,tile_i,0]=cable524xx & if (cable_len[tile_i] EQ 524) then bp_sol[*,tile_i,1]=cable524yy
      endfor
      
      
      
      ;Normal by obs
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/calibration/'+obsid+'_cal.sav'
      normal_res[*,*,obs_i,0]=*cal.gain_residual[0]
      normal_sols[*,*,obs_i,0]=*cal.gain[0]
      normal_res[*,*,obs_i,1]=*cal.gain_residual[1]
      normal_sols[*,*,obs_i,1]=*cal.gain[1]
      normal_input[*,*,obs_i,0]=*cal.gain_residual[0]+*cal.gain[0]
      normal_input[*,*,obs_i,1]=*cal.gain_residual[1]+*cal.gain[1]
      
      ;One quad, polyscaled cross-added modefit (90 and 150) **MODDED to two quads, polyscaled crossadded modefit with extra flagging
      restore, '~/Aug23_onequad_polyscaled_90150/forinput/'+obsid+'_cal.sav'
      unsplit_quad_sols[*,*,obs_i,0]=(*cal.gain[0])*bp_sol[*,*,0]
      unsplit_quad_sols[*,*,obs_i,1]=(*cal.gain[1])*bp_sol[*,*,1]
      
      unsplit_quad_res[*,*,obs_i,0]=normal_input[*,*,obs_i,0]-unsplit_quad_sols[*,*,obs_i,0]
      unsplit_quad_res[*,*,obs_i,1]=normal_input[*,*,obs_i,1]-unsplit_quad_sols[*,*,obs_i,1]
      ;stop
      ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_poly_saved_run_v2_plusmodeunsplit/calibration/'+obsid+'_cal.sav'
      ;unsplit_quad_sols[*,*,obs_i,0]=(*cal.gain[0])
      ;unsplit_quad_sols[*,*,obs_i,1]=(*cal.gain[1])
      
      ;unsplit_quad_res[*,*,obs_i,0]=(*cal.gain_residual[0])
      ;unsplit_quad_res[*,*,obs_i,1]=(*cal.gain_residual[0])
      
      
      ;Two split quads, polyscaled cross-added modefit
      ;restore, '~/Aug23_polyratio_quadsplit/forinput/'+obsid+'_cal.sav'
      ;split_quad_sols[*,*,obs_i,0]=*cal.gain[0]*bp_sol[*,*,0]
      ;split_quad_sols[*,*,obs_i,1]=*cal.gain[1]*bp_sol[*,*,1]
      
      ;split_quad_res[*,*,obs_i,0]=normal_input[*,*,obs_i,0]-split_quad_sols[*,*,obs_i,0]
      ;split_quad_res[*,*,obs_i,1]=normal_input[*,*,obs_i,1]-split_quad_sols[*,*,obs_i,1]
      
      ;Two constrained split quads, polyscaled cross-added modefit
      ;restore, '~/Aug23_nodig_quad_polyscaled/forinput/'+obsid+'_cal.sav'
      ;split_nodig_quad_sols[*,*,obs_i,0]=*cal.gain[0]*bp_sol[*,*,0]
      ;split_nodig_quad_sols[*,*,obs_i,1]=*cal.gain[1]*bp_sol[*,*,1]
      
      ;split_nodig_quad_res[*,*,obs_i,0]=normal_input[*,*,obs_i,0]-split_nodig_quad_sols[*,*,obs_i,0]
      ;split_nodig_quad_res[*,*,obs_i,1]=normal_input[*,*,obs_i,1]-split_nodig_quad_sols[*,*,obs_i,1]
      
      ;Two split linears, polyscaled cross-added modefit
      ;restore, '~/Aug23_unscaledautos_phasetransfer_polyratio/forinput/'+obsid+'_cal.sav'
      ;split_nodig_linear_sols[*,*,obs_i,0]=*cal.gain[0]*bp_sol[*,*,0]
      ;split_nodig_linear_sols[*,*,obs_i,1]=*cal.gain[1]*bp_sol[*,*,1]
      
      ;split_nodig_linear_res[*,*,obs_i,0]=normal_input[*,*,obs_i,0]-split_nodig_linear_sols[*,*,obs_i,0]
      ;split_nodig_linear_res[*,*,obs_i,1]=normal_input[*,*,obs_i,1]-split_nodig_linear_sols[*,*,obs_i,1]
      
      ;*******One quad, constrained split, with regular modefit by obs
      restore, '~/Aug23_std_test_nodigtwopolyquad/'+obsid+'_cal.sav'
      split_nodig_quad_sols[*,*,obs_i,0]=*cal.gain[0]*bp_sol[*,*,0]
      split_nodig_quad_sols[*,*,obs_i,1]=*cal.gain[1]*bp_sol[*,*,1]
      
      split_nodig_quad_res[*,*,obs_i,0]=normal_input[*,*,obs_i,0]-split_nodig_quad_sols[*,*,obs_i,0]
      split_nodig_quad_res[*,*,obs_i,1]=normal_input[*,*,obs_i,1]-split_nodig_quad_sols[*,*,obs_i,1]
      
      ;*******One quad, constrained split, with cross-phase added polyscaled modefits by obs
      restore, '~/Aug23_std_test_nodigtwopolyquad_fancymodeobs/'+obsid+'_cal.sav'
      split_nodig_quad_fancy_sols[*,*,obs_i,0]=*cal.gain[0]*bp_sol[*,*,0]
      split_nodig_quad_fancy_sols[*,*,obs_i,1]=*cal.gain[1]*bp_sol[*,*,1]
      
      split_nodig_quad_fancy_res[*,*,obs_i,0]=normal_input[*,*,obs_i,0]-split_nodig_quad_fancy_sols[*,*,obs_i,0]
      split_nodig_quad_fancy_res[*,*,obs_i,1]=normal_input[*,*,obs_i,1]-split_nodig_quad_fancy_sols[*,*,obs_i,1]
      
    endfor
    
    stop
    
    for pol_i=0,1 do begin
      for tile_i=0,127 do begin
      
      
        if pol_i EQ 0 then pol_name='xx'
        if pol_i EQ 1 then pol_name='yy'
        
        cgPS_Open,'/nfs/eor-00/h1/nbarry//Aug23_std_test_nodigtwopolyquad_fancymodeobs/fit_plots/'+poi_name[j]+'_comparisons_'+pol_name+'_'+strtrim(string(tile_i),2)+'.png',/quiet,/nomatch
        
        cgplot, (freq_arr)/1E6,abs(normal_input[freq_use,tile_i,0,pol_i]),xrange=[(freq_arr[0])/1E6,(freq_arr[335])/1E6], yrange=[min(abs(normal_input[freq_use,tile_i,*,0]))-.1,max(abs(normal_input[freq_use,tile_i,*,0]))+.1], $
          title=poi_name[j] + ' pointing, Tile ' + strtrim(string(tile_names[tile_i]),2) + $
          ' ('+strtrim(string(UINT(cable_len[tile_i])),2) +'m cable), '+pol_name, XTICKFORMAT="(A1)", ytitle='Gain',charsize=1,color='light grey', $
          Position=[0.10, 0.35, 0.9, 0.90]
          
        For normal_i=1, N_elements(*obs_ptr[j])-1 do begin
          cgoplot, (freq_arr)/1E6,abs(normal_input[freq_use,tile_i,normal_i,pol_i]),color='light grey'
        endfor
        
        Device, Decomposed=0
        TVLCT, rgbcolors[0,0], rgbcolors[1,0], rgbcolors[2,0],10
        For normal_i=0, N_elements(*obs_ptr[j])-1 do begin
          cgoplot, (freq_arr)/1E6,abs(normal_sols[freq_use,tile_i,normal_i,pol_i]),color=10B
        endfor
        
        TVLCT, rgbcolors[0,1], rgbcolors[1,1], rgbcolors[2,1],20
        cgoplot, (freq_arr)/1E6,abs(unsplit_quad_sols[freq_use,tile_i,0,pol_i]),color=20B, thick=2
        
        ;TVLCT, rgbcolors[0,2], rgbcolors[1,2], rgbcolors[2,2],30
        ;cgoplot, (freq_arr)/1E6,abs(split_quad_sols[freq_use,tile_i,0,pol_i]),color=30B,thick=2
        
        TVLCT, rgbcolors[0,3], rgbcolors[1,3], rgbcolors[2,3],30
        cgoplot, (freq_arr)/1E6,abs(split_nodig_quad_sols[freq_use,tile_i,0,pol_i]),color=30B,thick=2
        
        TVLCT, rgbcolors[0,4], rgbcolors[1,4], rgbcolors[2,4],40
        cgoplot, (freq_arr)/1E6,abs(split_nodig_quad_fancy_sols[freq_use,tile_i,0,pol_i]),color=40B,thick=2
        
        ;for Abraham
        ;TVLCT, rgbcolors[0,4], rgbcolors[1,4], rgbcolors[2,4],50
        ;cgoplot, (freq_arr)/1E6,abs(mean(normal_input[freq_use,tile_i,*,pol_i],dimension=3)),color=50B,thick=2
        
        TVLCT, 205,201,201,60
        
        ;for Abraham
        cgLegend, Title=['Unfit gains','Cal fit obsid (+150m)','Split cal fit pointing (+90m*bp,+150m*bp)','Split continuous cal fit pointing (+150m bo)',$
          'Split continuous cal fit pointing (+90m*bo,+150m*bo)'], $
          ;cgLegend, Title=['Unfit gains','Cal fit obsid (+150m)','Cal fit pointing (+90m,+150m)','Split cal fit pointing (+150m)',$
          ;'Split continuous cal fit pointing (+150m)','Split cal linear fit pointing (+150m)'], $
          Color=[60B,10B,20B,30B,40B],Length=.03,charsize=.7,$;Psym=[2,2,2,2,2,2],Length=0.0 $
          Location=[0.56,0.87]
          
          
        if keyword_set(bottom_included) then begin
          normal_input_phase=Atan(normal_input,/phase)
          normal_phase=Atan(normal_sols,/phase)
          unsplit_quad_phase=Atan(unsplit_quad_sols,/phase)
          split_quad_phase=Atan(split_quad_sols,/phase)
          split_nodig_quad_phase=Atan(split_nodig_quad_sols,/phase)
          split_nodig_linear_phase=Atan(split_nodig_linear_sols,/phase)
          
          mean_normal=mean( abs(normal_input)*exp(Complex(0,1)*normal_input_phase)*exp(-Complex(0,1)*normal_phase), dimension=3)
          mean_unsplit_quad=mean( abs(normal_input)*exp(Complex(0,1)*normal_input_phase)*exp(-Complex(0,1)*unsplit_quad_phase), dimension=3)
          mean_split_quad=mean( abs(normal_input)*exp(Complex(0,1)*normal_input_phase)*exp(-Complex(0,1)*split_quad_phase), dimension=3)
          mean_split_nodig_quad=mean( abs(normal_input)*exp(Complex(0,1)*normal_input_phase)*exp(-Complex(0,1)*split_nodig_quad_phase), dimension=3)
          mean_split_nodig_linear=mean( abs(normal_input)*exp(Complex(0,1)*normal_input_phase)*exp(-Complex(0,1)*split_nodig_linear_phase), dimension=3)
          
          unphased_real=real_part( abs(normal_sols[freq_use,tile_i,*,pol_i]) - $
            abs(normal_input[freq_use,tile_i,*,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,*,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,*,pol_i]))
            
          cgPlot,(freq_arr)/1E6, mean(reform(unphased_real),dimension=2), $
            Position=[0.1, 0.25, 0.9, 0.34], color=10B,xrange=[freq_arr[0]/1E6,freq_arr[335]/1E6],$
            XTICKFORMAT="(A1)",/NoErase, charsize=.75, yrange=[-.05,.05],ytitle='Real!cunphased',YTICKFORMAT="(F0.2)",ytickinterval=.05
          ;for obs_i=1,N_elements(*obs_ptr[j])-1 do cgoplot,real_part( abs(normal_sols[freq_use,tile_i,obs_i,pol_i]) -  $
          ;  abs(normal_input[freq_use,tile_i,obs_i,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,obs_i,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,obs_i,pol_i]) ),color=10B
          cgoplot,(freq_arr)/1E6,real_part( abs(unsplit_quad_sols[freq_use,tile_i,0,pol_i]) - mean_unsplit_quad[freq_use,tile_i,pol_i] ),color=20B
          cgoplot,(freq_arr)/1E6,real_part( abs(split_quad_sols[freq_use,tile_i,0,pol_i]) - mean_split_quad[freq_use,tile_i,pol_i] ),color=30B
          cgoplot,(freq_arr)/1E6,real_part( abs(split_nodig_quad_sols[freq_use,tile_i,0,pol_i]) - mean_split_nodig_quad[freq_use,tile_i,pol_i] ),color=40B
          cgoplot,(freq_arr)/1E6,real_part( abs(split_nodig_linear_sols[freq_use,tile_i,0,pol_i]) - mean_split_nodig_linear[freq_use,tile_i,pol_i] ),color=50B
          
          unphased_im=abs( abs(normal_sols[freq_use,tile_i,*,pol_i]) - $
            abs(normal_input[freq_use,tile_i,*,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,*,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,*,pol_i]))
            
          cgPlot, (freq_arr)/1E6,mean(reform(unphased_im),dimension=2), $
            Position=[0.1, 0.15, 0.9, 0.24], color=10B,xrange=[freq_arr[0]/1E6,freq_arr[335]/1E6],$
            ytickinterval=.05,xtitle='Unflagged frequency MHz',/NoErase, charsize=.75, yrange=[0,.1],ytitle='Magnitude!cunphased';,YTICKFORMAT="(A1)"
            
          ;for obs_i=1,N_elements(*obs_ptr[j])-1 do cgoplot, imaginary( abs(normal_sols[freq_use,tile_i,obs_i,pol_i]) - $
          ;  abs(normal_input[freq_use,tile_i,obs_i,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,obs_i,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,obs_i,pol_i]) ),color=10B
            
          cgoplot,(freq_arr)/1E6,abs( abs(unsplit_quad_sols[freq_use,tile_i,0,pol_i]) - mean_unsplit_quad[freq_use,tile_i,pol_i] ),color=20B
          cgoplot,(freq_arr)/1E6,abs( abs(split_quad_sols[freq_use,tile_i,0,pol_i]) - mean_split_quad[freq_use,tile_i,pol_i] ),color=30B
          cgoplot,(freq_arr)/1E6,abs( abs(split_nodig_quad_sols[freq_use,tile_i,0,pol_i]) - mean_split_nodig_quad[freq_use,tile_i,pol_i] ),color=40B
          cgoplot,(freq_arr)/1E6,abs( abs(split_nodig_linear_sols[freq_use,tile_i,0,pol_i]) - mean_split_nodig_linear[freq_use,tile_i,pol_i] ),color=50B
          
        endif
        
        Device, Decomposed=1
        
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
      endfor
    endfor
    
  endfor
  
end