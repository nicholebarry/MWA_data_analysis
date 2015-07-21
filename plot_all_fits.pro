pro plot_all_fits, outdir, bottom_included=bottom_included

  ;This code plots a variety of full cal solutions over a light grey background of input gains.
  ;Which cal solutions to plot can be changed below in sections A through E. Definitions of
  ;titleA, titleB, titleC, titleD, titleE determine if they get plotted.
  ;bottom_included plots two small residual plots below the main cal fit plot,
  ;as unphased mag and real residual
  ;Outdir is the path to the locations of the outputs, and must be specified, and end in /


  ;****************Begin setup

  undefine, titleA, titleB, titleC, titleD, titleE
  
  ;Define names of pointings, day, and names of obs files
  poi_name=['-2','-1','0','1','2','3']
  day='Aug23'
  parsednames=[day+'minustwo',day+'minusone',day+'zenith',day+'plusone',day+'plustwo',day+'plusthree']
  
  ;obs ptr is where obs of a pointing are stored. It can be made bigger than needed to account for the pointings of the longrun
  obs_ptr=PTRARR(6,/allocate)
  
  ;Read in obsids of each pointing for the given day
  FOR j=0,(size(poi_name))[1]-1 DO BEGIN
  
    filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/obs_list/' + parsednames[j] + '.txt'
    obs_temp='empty'
    readcol, filename, obs_temp, format='A', /silent
    If j NE 0 then parsednumbers=[parsednumbers,N_elements(obs_temp)] else parsednumbers=N_elements(obs_temp)
    *obs_ptr[j]=obs_temp
    
  ENDFOR
  
  ;Restore a typical obs file to get unflagged frequencies and tile names
  restore,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/metadata/1061311664_obs.sav'
  freq_use=where((*obs.baseline_info).freq_use)
  freq_arr=(*obs.baseline_info).freq[where((*obs.baseline_info).freq_use)]
  tile_names=(*obs.baseline_info).tile_names
  
  ;Read in cable length data
  mode_filepath=filepath('mwa_cable_reflection_coefficients.txt',root=rootdir('FHD'),subdir='instrument_config')
  textfast,data_array,/read,file_path=mode_filepath,first_line=1
  cable_len=Reform(data_array[2,*])
  
  ;Define colors for plotting
  rgbcolors=[[128,187,255],[202,18,27],[226,139,176],[229,153,89],[12,111,188]]
  
  ;****************End setup
  
  
  ;****************Begin reading and defining cal solutions to plot
  
  ;Loop through each pointing
  For j=0,N_elements(poi_name)-1 do begin
  
    ;Define complex arrays of solutions and residuals, freq x tile x obs in current pointing
    normal_res=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    normal_sols=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    normal_input=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    B_res=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    B_sols=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    C_res=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    C_sols=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    D_res=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    D_sols=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    E_res=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    E_sols=complex(FLTARR(384,128,N_elements(*obs_ptr[j]),2))
    
    ;Define complex array of the typical read-in bandpass solution for special cases
    bp_sol=complex(FLTARR(384,128,2))
    
    
    ;Loop through each obs per current pointing
    for obs_i=0,N_elements(*obs_ptr[j])-1 do begin
    
      ;The current obsid given loop count and current pointing
      obsid=(*obs_ptr[j])[obs_i]
      
      ;Semi-perminant file location for bp solutions
      filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/instrument_config/'+poi_name[j]+'_bandpass.txt
      ;Read in and save bp solutions for each cable/pol type
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
      
      ;Sort the bandpass solutions into one array of freq x tile x pol
      for tile_i=0,127 do begin
        if (cable_len[tile_i] EQ 90) then bp_sol[*,tile_i,0]=cable90xx & if (cable_len[tile_i] EQ 90) then bp_sol[*,tile_i,1]=cable90yy
        if (cable_len[tile_i] EQ 150) then bp_sol[*,tile_i,0]=cable150xx & if (cable_len[tile_i] EQ 150) then bp_sol[*,tile_i,1]=cable150yy
        if (cable_len[tile_i] EQ 230) then bp_sol[*,tile_i,0]=cable230xx & if (cable_len[tile_i] EQ 230) then bp_sol[*,tile_i,1]=cable230yy
        if (cable_len[tile_i] EQ 320) then bp_sol[*,tile_i,0]=cable320xx & if (cable_len[tile_i] EQ 320) then bp_sol[*,tile_i,1]=cable320yy
        if (cable_len[tile_i] EQ 400) then bp_sol[*,tile_i,0]=cable400xx & if (cable_len[tile_i] EQ 400) then bp_sol[*,tile_i,1]=cable400yy
        if (cable_len[tile_i] EQ 524) then bp_sol[*,tile_i,0]=cable524xx & if (cable_len[tile_i] EQ 524) then bp_sol[*,tile_i,1]=cable524yy
      endfor
      
      
      
      ;A***Normal by obs (one quadratic polfit by obs, modefit using crosses for 150 by obs)
      titleA='Cal fit obsid (+150m)'
      restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/calibration/'+obsid+'_cal.sav'
      normal_res[*,*,obs_i,0]=*cal.gain_residual[0]
      normal_sols[*,*,obs_i,0]=*cal.gain[0]
      normal_res[*,*,obs_i,1]=*cal.gain_residual[1]
      normal_sols[*,*,obs_i,1]=*cal.gain[1]
      normal_input[*,*,obs_i,0]=*cal.gain_residual[0]+*cal.gain[0]
      normal_input[*,*,obs_i,1]=*cal.gain_residual[1]+*cal.gain[1]
      ;A***End of normal by obs
      
      
      ;B***Two quad, polyscaled cross-added modefit (90 and 150) with extra flagging (proven wrong for 90m, improvement for 150m)
      titleB='Two quad cal fit pointing (+150pf,+90pf)'
      restore, '/nfs/eor-00/h1/nbarry/Aug23_std_test_twopolyquad_fancymodeobs/'+obsid+'_cal.sav'
      B_sols[*,*,obs_i,0]=(*cal.gain[0])*bp_sol[*,*,0]
      B_sols[*,*,obs_i,1]=(*cal.gain[1])*bp_sol[*,*,1]
      
      B_res[*,*,obs_i,0]=normal_input[*,*,obs_i,0]-B_sols[*,*,obs_i,0]
      B_res[*,*,obs_i,1]=normal_input[*,*,obs_i,1]-B_sols[*,*,obs_i,1]
      ;B***End of Two quad, polyscaled cross-added modefit (90 and 150) with extra flagging
      
      
      ;C***Unknown
      ;titleC=
      ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_poly_saved_run_v2_plusmodeunsplit/calibration/'+obsid+'_cal.sav'
      ;C_sols[*,*,obs_i,0]=(*cal.gain[0])
      ;C_sols[*,*,obs_i,1]=(*cal.gain[1])
      
      ;C_res[*,*,obs_i,0]=(*cal.gain_residual[0])
      ;C_res[*,*,obs_i,1]=(*cal.gain_residual[0])
      ;C***Unknown
      
      
      ;D***One quad, constrained split, with regular modefit by obs
      titleD='Two quad cal fit pointing (+150e,+90e)'
      restore, '/nfs/eor-00/h1/nbarry/Aug23_std_test_towpolyquad_extrafancymodeobs/'+obsid+'_cal.sav'
      D_sols[*,*,obs_i,0]=*cal.gain[0]*bp_sol[*,*,0]
      D_sols[*,*,obs_i,1]=*cal.gain[1]*bp_sol[*,*,1]
      
      D_res[*,*,obs_i,0]=normal_input[*,*,obs_i,0]-D_sols[*,*,obs_i,0]
      D_res[*,*,obs_i,1]=normal_input[*,*,obs_i,1]-D_sols[*,*,obs_i,1]
    ;D***End of One quad, constrained split, with regular modefit by obs
      
      
    ;E***One quad, constrained split, with cross-phase added polyscaled modefits by obs
    ;titleE=
    ;restore, '~/Aug23_std_test_nodigtwopolyquad_fancymodeobs/'+obsid+'_cal.sav'
    ;E_sols[*,*,obs_i,0]=*cal.gain[0]*bp_sol[*,*,0]
    ;E_sols[*,*,obs_i,1]=*cal.gain[1]*bp_sol[*,*,1]
      
    ;E_res[*,*,obs_i,0]=normal_input[*,*,obs_i,0]-E_sols[*,*,obs_i,0]
    ;E_res[*,*,obs_i,1]=normal_input[*,*,obs_i,1]-E_sols[*,*,obs_i,1]
    ;E***End of One quad, constrained split, with cross-phase added polyscaled modefits by obs
      
    endfor
    
    ;****************End reading and defining cal solutions to plot
    
    
    for pol_i=0,1 do begin
      for tile_i=0,127 do begin
      
        ;Naming convention
        if pol_i EQ 0 then pol_name='xx'
        if pol_i EQ 1 then pol_name='yy'
        
        cgPS_Open,outdir+poi_name[j]+'_comparisons_'+pol_name+'_'+strtrim(string(tile_i),2)+'.png',/quiet,/nomatch
        
        cgplot, (freq_arr)/1E6,abs(normal_input[freq_use,tile_i,0,pol_i]),xrange=[(freq_arr[0])/1E6,(freq_arr[335])/1E6], yrange=[min(abs(normal_input[freq_use,tile_i,*,0]))-.1,max(abs(normal_input[freq_use,tile_i,*,0]))+.1], $
          title=poi_name[j] + ' pointing, Tile ' + strtrim(string(tile_names[tile_i]),2) + $
          ' ('+strtrim(string(UINT(cable_len[tile_i])),2) +'m cable), '+pol_name, XTICKFORMAT="(A1)", ytitle='Gain',charsize=1,color='light grey', $
          Position=[0.10, 0.35, 0.9, 0.90]
        titletotal='Unfit gains'
        colortotal=60B
        
        For normal_i=1, N_elements(*obs_ptr[j])-1 do begin
          cgoplot, (freq_arr)/1E6,abs(normal_input[freq_use,tile_i,normal_i,pol_i]),color='light grey'
        endfor
        
        Device, Decomposed=0
        
        If keyword_set(titleA) then begin
          TVLCT, rgbcolors[0,0], rgbcolors[1,0], rgbcolors[2,0],10
          For normal_i=0, N_elements(*obs_ptr[j])-1 do begin
            cgoplot, (freq_arr)/1E6,abs(normal_sols[freq_use,tile_i,normal_i,pol_i]),color=10B
          endfor
          titletotal=[titletotal, titleA]
          colortotal=[colortotal, 10B]
        endif
        
        If keyword_set(titleB) then begin
          TVLCT, rgbcolors[0,1], rgbcolors[1,1], rgbcolors[2,1],20
          cgoplot, (freq_arr)/1E6,abs(B_sols[freq_use,tile_i,0,pol_i]),color=20B, thick=2
          titletotal=[titletotal, titleB]
          colortotal=[colortotal, 20B]
        endif
        
        If keyword_set(titleC) then begin
          ;TVLCT, rgbcolors[0,2], rgbcolors[1,2], rgbcolors[2,2],30
          ;cgoplot, (freq_arr)/1E6,abs(C_sols[freq_use,tile_i,0,pol_i]),color=30B,thick=2
          titletotal=[titletotal, titleC]
          colortotal=[colortotal, 30B]
        endif
        
        If keyword_set(titleD) then begin
          TVLCT, rgbcolors[0,3], rgbcolors[1,3], rgbcolors[2,3],40
          cgoplot, (freq_arr)/1E6,abs(D_sols[freq_use,tile_i,0,pol_i]),color=40B,thick=2
          titletotal=[titletotal, titleD]
          colortotal=[colortotal, 40B]
        endif
        
        If keyword_set(titleE) then begin
          TVLCT, rgbcolors[0,4], rgbcolors[1,4], rgbcolors[2,4],50
          cgoplot, (freq_arr)/1E6,abs(E_sols[freq_use,tile_i,0,pol_i]),color=50B,thick=2
          titletotal=[titletotal, titleE]
          colortotal=[colortotal, 50B]
        endif
        
        TVLCT, 205,201,201,60
        
        
        cgLegend, Title=titletotal, $
          Color=colortotal,Length=.03,charsize=.7,$;Psym=[2,2,2,2,2,2],Length=0.0 $
          Location=[0.56,0.87]
          
          
        if keyword_set(bottom_included) then begin
          normal_input_phase=Atan(normal_input,/phase)
          normal_phase=Atan(normal_sols,/phase)
          B_phase=Atan(B_sols,/phase)
          C_phase=Atan(C_sols,/phase)
          D_phase=Atan(D_sols,/phase)
          E_phase=Atan(E_sols,/phase)
          
          mean_normal=mean( abs(normal_input)*exp(Complex(0,1)*normal_input_phase)*exp(-Complex(0,1)*normal_phase), dimension=3)
          mean_B=mean( abs(normal_input)*exp(Complex(0,1)*normal_input_phase)*exp(-Complex(0,1)*B_phase), dimension=3)
          mean_C=mean( abs(normal_input)*exp(Complex(0,1)*normal_input_phase)*exp(-Complex(0,1)*C_phase), dimension=3)
          mean_D=mean( abs(normal_input)*exp(Complex(0,1)*normal_input_phase)*exp(-Complex(0,1)*D_phase), dimension=3)
          mean_E=mean( abs(normal_input)*exp(Complex(0,1)*normal_input_phase)*exp(-Complex(0,1)*E_phase), dimension=3)
          
          unphased_real=real_part( abs(normal_sols[freq_use,tile_i,*,pol_i]) - $
            abs(normal_input[freq_use,tile_i,*,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,*,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,*,pol_i]))
            
          cgPlot,(freq_arr)/1E6, mean(reform(unphased_real),dimension=2), $
            Position=[0.1, 0.25, 0.9, 0.34], color=10B,xrange=[freq_arr[0]/1E6,freq_arr[335]/1E6],$
            XTICKFORMAT="(A1)",/NoErase, charsize=.75, yrange=[-.05,.05],ytitle='Real!cunphased',YTICKFORMAT="(F0.2)",ytickinterval=.05
          ;for obs_i=1,N_elements(*obs_ptr[j])-1 do cgoplot,real_part( abs(normal_sols[freq_use,tile_i,obs_i,pol_i]) -  $
          ;  abs(normal_input[freq_use,tile_i,obs_i,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,obs_i,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,obs_i,pol_i]) ),color=10B
          If keyword_set(titleB) then cgoplot,(freq_arr)/1E6,real_part( abs(B_sols[freq_use,tile_i,0,pol_i]) - mean_B[freq_use,tile_i,pol_i] ),color=20B
          If keyword_set(titleC) then cgoplot,(freq_arr)/1E6,real_part( abs(C_sols[freq_use,tile_i,0,pol_i]) - mean_C[freq_use,tile_i,pol_i] ),color=30B
          If keyword_set(titleD) then cgoplot,(freq_arr)/1E6,real_part( abs(D_sols[freq_use,tile_i,0,pol_i]) - mean_D[freq_use,tile_i,pol_i] ),color=40B
          If keyword_set(titleE) then cgoplot,(freq_arr)/1E6,real_part( abs(E_sols[freq_use,tile_i,0,pol_i]) - mean_E[freq_use,tile_i,pol_i] ),color=50B
          
          unphased_im=abs( abs(normal_sols[freq_use,tile_i,*,pol_i]) - $
            abs(normal_input[freq_use,tile_i,*,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,*,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,*,pol_i]))
            
          cgPlot, (freq_arr)/1E6,mean(reform(unphased_im),dimension=2), $
            Position=[0.1, 0.15, 0.9, 0.24], color=10B,xrange=[freq_arr[0]/1E6,freq_arr[335]/1E6],$
            ytickinterval=.05,xtitle='Unflagged frequency MHz',/NoErase, charsize=.75, yrange=[0,.1],ytitle='Magnitude!cunphased';,YTICKFORMAT="(A1)"
            
          ;for obs_i=1,N_elements(*obs_ptr[j])-1 do cgoplot, imaginary( abs(normal_sols[freq_use,tile_i,obs_i,pol_i]) - $
          ;  abs(normal_input[freq_use,tile_i,obs_i,pol_i])*exp(Complex(0,1)*normal_input_phase[freq_use,tile_i,obs_i,pol_i])*exp(-Complex(0,1)*normal_phase[freq_use,tile_i,obs_i,pol_i]) ),color=10B
            
          If keyword_set(titleB) then cgoplot,(freq_arr)/1E6,abs( abs(B_sols[freq_use,tile_i,0,pol_i]) - mean_B[freq_use,tile_i,pol_i] ),color=20B
          If keyword_set(titleC) then cgoplot,(freq_arr)/1E6,abs( abs(C_sols[freq_use,tile_i,0,pol_i]) - mean_C[freq_use,tile_i,pol_i] ),color=30B
          If keyword_set(titleD) then cgoplot,(freq_arr)/1E6,abs( abs(D_sols[freq_use,tile_i,0,pol_i]) - mean_D[freq_use,tile_i,pol_i] ),color=40B
          If keyword_set(titleE) then cgoplot,(freq_arr)/1E6,abs( abs(E_sols[freq_use,tile_i,0,pol_i]) - mean_E[freq_use,tile_i,pol_i] ),color=50B
          
        endif
        
        Device, Decomposed=1
        
        cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
        
      endfor
    endfor
    
  endfor
  
end