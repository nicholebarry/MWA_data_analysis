pro plot_limits

  n_papers=8
  papers_array = PTRARR(n_papers,/allocate)
  name_array = STRARR(n_papers)
  
  ;Organized by paper: mK value, center of k bin/range, and redshift
  ali_2015 = transpose([22.4^2.,.325,8.4])
  *papers_array[0] = ali_2015
  name_array[0] = 'Ali, 2015'
  
    paciga_2013 = transpose([6.15E4,.5,8.6])
  *papers_array[1] = paciga_2013
  name_array[1] = 'Paciga, 2013'
  
  parsons_2014 = transpose([41^2.,.27,7.7])
  *papers_array[2] = parsons_2014
  name_array[2] = 'Parsons, 2014'
  
  patil_2017 = FLTARR(15,3)
  patil_2017[*,0] = [131.5^2.,242.1^2.,220.9^2.,337.4^2.,407.7^2.,86.4^2.,144.2^2.,184.7^2.,$
    296.1^2.,342.0^2.,79.6^2.,108.8^2.,148.6^2.,224.0^2.,366.1^2.]
  patil_2017[*,1] = [.053,.067,.083,.103,.128,.053,.067,.083,.103,.128,.053,.067,.083,.103,.128]
  patil_2017[*,2] = [8.3,8.3,8.3,8.3,8.3,9.15,9.15,9.15,9.15,9.15,10.1,10.1,10.1,10.1,10.1]
    patil_2017_short = FLTARR(6,3)
  patil_2017_short[*,0] = [131.5^2.,220.9^2.,86.4^2.,184.7^2.,79.6^2.,148.6^2.]
  patil_2017_short[*,1] = [.053,.083,.053,.083,.053,.083]
  patil_2017_short[*,2] = [8.3,8.3,9.15,9.15,10.1,10.1]
  *papers_array[3] = patil_2017_short
  name_array[3] = 'Patil, 2017'
  
  jacobs_2015 = FLTARR(4,3)
  jacobs_2015[*,0] = [1.64E4,6.9E3,6.9E3,2.4E3]
  jacobs_2015[*,1] = [.2,.1,.2,.2]
  jacobs_2015[*,2] = [10.29,8.54,7.94,7.55]
  *papers_array[4] = jacobs_2015
  name_array[4] = 'Jacobs, 2015'
 
  dillon_2014 = FLTARR(11,3)
  dillon_2014[*,0] = [2.6E7,1.16E6,8.64E5,6.7E5,1.3E6,1.28E7,5.26E7,5.67E8,4.58E6,2.93E8,6.92E8]
  dillon_2014[*,1] = [.058,.06,.063,.065,.0678,.0712,.0715,.149,.15,.15,.089]
  dillon_2014[*,2] = [11.68,10.868,10.153,9.518,8.444,7.985,7.57,7.1896,6.84,6.52,6.23]
  *papers_array[5] = dillon_2014
  name_array[5] = 'Dillon, 2014' 
  
  dillon_2015 = FLTARR(3,3)
  dillon_2015[*,0] = [3.8E4,3.69E4,4.67E4]
  dillon_2015[*,1] = [.18,.18,.16]
  dillon_2015[*,2] = [6.4,6.8,7.25]
  *papers_array[6] = dillon_2015
  name_array[6] = 'Dillon, 2015'
  
    beardsley_2016 = FLTARR(9,3)
  beardsley_2016[*,0] = [3.67E4,2.70E4,3.56E4,3.02E4,4.70E4,3.22E4,3.2E4,2.6E4,2.5E4]
  beardsley_2016[*,1] = [.231,.27,.24,.24,.20,.24,.16,.14,.14]
  beardsley_2016[*,2] = [7.1,7.1,6.8,6.8,6.5,6.5,7.1,6.8,6.5]
  beardsley_2016_short = FLTARR(6,3)
  beardsley_2016_short[*,0] = [2.70E4,3.02E4,3.22E4,3.2E4,2.6E4,2.5E4]
  beardsley_2016_short[*,1] = [.27,.24,.24,.16,.14,.14]
  beardsley_2016_short[*,2] = [7.1,6.8,6.5,7.1,6.8,6.5]
  *papers_array[7] = beardsley_2016_short
  name_array[7] = 'Beardsley, 2016'
  
  fiducial_z = [6.,6.14,6.28,6.43,6.58,6.73,6.88,7.04,7.2 ,7.37,7.53,7.7,7.88,8.06,8.24,8.42,8.61,8.8,$
    9.,9.2,9.4,9.61,9.82,10.04,10.26,10.49,10.71,10.95,11.19,11.43,11.68,11.93,12.19,12.46,12.73,13.,13.28]
  fiducial_val = [0.2070598,0.4520246,0.8806551,1.628951,2.754617,4.320467,6.392936,9.135006,12.36413,$
    16.02066,19.54978,23.19749,26.69443,29.66309,31.89128,33.15714,33.30726,32.40067,30.45852, 27.84677,$
    24.84632,21.71832,18.96164,16.65482,15.0046,13.9677,13.60488,13.83113,14.61855,15.90612,17.75889,$
    20.20181,23.36801,27.49878,32.82986,39.85733, 49.35672]
  
  
  ;k_bin_arr = [reform(transpose(ali_2015[*,1])),reform(transpose(dillon_2015[*,1])),reform(transpose(beardsley_2016[*,1])),$
  ;  reform(transpose(parsons_2014[*,1])),reform(transpose(patil_2017[*,1])),reform(transpose(jacobs_2015[*,1])),reform(transpose(dillon_2014[*,1])),reform(transpose(paciga_2013[*,1]))]
  
  rgbcolors = [[200,85,182],$
    [140,166,58],$
    [121,103,214],$
    [95,162,113],$
    [200,89,121],$
    [89,149,203],$
    [199,111,62],$
    [161,116,188]]
  color_num=[10,20,30,40,50,60,70,80]
  color_byte = [10B,20B,30B,40B,50B,60B,70B,80B]
  sym_array = ['Filled Diamond','Filled Bowtie','Filled Up Triangle','Filled Lower Half Circle','Filled Standing Bar','Filled Circle','Filled Laying Bar','Filled Star']
  sym_num = [14,24,17,40,26,16,28,46]
  
  color_names = ['purple', 'teal', 'blue', 'firebrick','purple', 'teal', 'blue', 'firebrick']
  
  x_range= [6,11]
  y_range = [1E-1,1E8]
  lessthan = '!9' + string("74B) + '!X'
  greaterthan = '!9' + string("76B) + '!X'
  lessthanequal = '!9' + string("243B) + '!X'
  titles = ['0 ' + lessthan + ' k '+ lessthan +' 0.07', '0.07 ' + lessthanequal + ' k '+ lessthan +' 0.14', $
     '0.14 ' + lessthanequal + ' k '+ lessthan +' 0.21', 'k '+ lessthanequal +' 0.21']
  
  positions = [[0.10, 0.6, 0.45, 0.95],[0.55, 0.6, 0.90, 0.95],[0.10, 0.15, 0.45, 0.5],[0.55, 0.15, 0.90, 0.5]]
  locations = [[0.1,0.045],[0.31,0.045],[0.51,0.045],[0.72,0.045],[0.1,0.01],[.31,.01],[.51,.01],[.72,.01]]
  x_title=['','','redshift','redshift']
  y_title=['$\Delta$$\up2$ (mK)$\up2$','','$\Delta$$\up2$ (mK)$\up2$','']
  
  cgps_open,'/nfs/eor-00/h1/nbarry/test_limit_plot.pdf',/quiet,/nomatch
  
  for plot_i=0, 3 do begin
  
    cgplot, [x_range[0],x_range[1]],[y_range[0],y_range[1]] , xrange=x_range, yrange=y_range, /ylog, charsize=.9,$
      Position=positions[*,plot_i], /NoData, title=titles[plot_i],xtitle=x_title[plot_i],ytitle=y_title[plot_i],yGridStyle=1, yTicklen=.5,/Noerase
      
    cgplot, fiducial_z, fiducial_val, color='brown',$
          xrange=x_range, yrange=y_range,/ylog, charsize=1, Position=positions[*,plot_i],/NoErase,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
      
    if plot_i EQ 0 then cgtext, .23,.65, 'Fiducial 21cmFAST Model', color='brown',charsize=.8  ,/normal
      
    for papers_i=0, n_papers-1 do begin
    
      k_arr = *papers_array[papers_i]
      
      if plot_i EQ 0 then inds = where(k_arr[*,1] LT .07,n_count)
      if plot_i EQ 1 then inds = where((k_arr[*,1] LT .14) AND (k_arr[*,1] GE .07),n_count)
      if plot_i EQ 2 then inds = where((k_arr[*,1] LT .21) AND (k_arr[*,1] GE .14),n_count)
      if plot_i EQ 3 then inds = where(k_arr[*,1] GE .21,n_count)
      
      if n_count GT 0 then begin
      
        TVLCT, rgbcolors[0,papers_i], rgbcolors[1,papers_i], rgbcolors[2,papers_i], color_num[papers_i]
        
        cgplot, k_arr[inds,2], k_arr[inds,0], color=color_names[papers_i],psym=sym_array[papers_i],$
          xrange=x_range, yrange=y_range,/ylog, charsize=1, Position=positions[*,plot_i],/NoErase,XTICKFORMAT="(A1)",YTICKFORMAT="(A1)"
          
      endif
      
      if plot_i EQ 3 then begin
        cglegend, Title=name_array[papers_i], color=color_names[papers_i], location = locations[*,papers_i], psym=sym_num[papers_i], charsize=1,Length=0
        ;stop
      endif
      
    endfor
    
  endfor
  
  
  
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
  stop
  
end