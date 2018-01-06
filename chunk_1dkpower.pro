pro chunk_1dkpower

  n_cubes=10
  
  chunks=strtrim(INDGEN(n_cubes)+1,2)
  chunk_names=['1-100','101-200','201-300','301-400','401-500','501-600','601-700','701-800','801-900','901-1000']
  cut='beardsley_thesis_list_int_chunk'
  chans=''
  win=''
  
  basefile='/nfs/mwa-04/r1/EoRuvfits/analysis/fhd_nb_2013longrun_autocal/ps_chunks/Combined_obs_'+cut+chunks+'_cubeXX__even_odd_joint_'+win+chans
  
  cubes=['res']
  pols=['xx','yy']
  
  endfile='_averemove_bh_dencorr_no_horizon_wedge_kperplambda10-50_1dkpower.idlsave'
  
  limit_percent=0.9772 ;2 sigma
  hubble_param=0.71
  
  rgbcolors=[[151,86,194],$
    [103,176,67],$
    [203,126,193],$
    [85,162,112],$
    [195,79,124],$
    [158,152,61],$
    [106,127,203],$
    [201,132,67],$
    [69,176,207],$
    [203,79,66]]
    
  color_num = [10,12,14,16,18,20,22,24,26,28]
  
  for n_i=0, n_cubes-1 do $
    TVLCT, rgbcolors[0,n_i], rgbcolors[1,n_i], rgbcolors[2,n_i], color_num[n_i]
    
    colors = [10B,12B,14B,16B,18B,20B,22B,24B,26B,28B]
    
  for j=0,1 do begin
    for cube_i=0,n_cubes-1 do begin
    
      file=basefile[cube_i]+cubes+'_'+pols[j]+endfile
      restore,file
      
      ; find the best limit in this file
      n_k=n_elements(k_edges)
      k=(k_edges[1:(n_k-1)]+k_edges[0:(n_k-2)])/2.
      delta=power*(k^3.)/(2.*!pi^2.)
      dsigma=(k^3.)/(2.*!pi^2.)/sqrt(weights)
      ;limits=abs(delta)+2.*dsigma ; actual value plus two sigma
      limits=dsigma*(inverf(limit_percent-(1.-limit_percent)*erf(delta/dsigma/sqrt(2)))*sqrt(2))+delta
      lim=min(limits,ind)
      header='#Limit: '+number_formatter(lim)+' mK^2, at k = '+number_formatter(k[ind]/hubble_param)+' h Mpc^-1 ('+cut+' '+chans+' '+cubes+' '+pols[j]+')'
      print,header
      
      ;if (cube_i EQ 0) then cgPS_Open,'/nfs/mwa-00/h1/nbarry/'+pols[j]+'_supercut_limit_comparison'+chans+'.png',/quiet,/nomatch
      if cube_i EQ 0 then begin
        cgplot, k,abs(limits),/xlog,/ylog,psym=10, xrange=[.08,1.5], yrange=[10^3.,10^8.],ytitle='$\Delta$$\up2$ (mK$\up2$)', $
          xtitle='k (h Mpc$\up-1$)', title=pols[j] + ' Limit Comparisons, Std 100 Obs Chunks ', charsize =1.5, color=colors[cube_i]
      endif else begin
        cgoplot, k,abs(limits),/xlog,/ylog,psym=10,color=colors[cube_i]
      endelse
      cgoplot, k,dsigma,/xlog,/ylog,psym=10,linestyle=2,color=colors[cube_i]
      
    ;if (cube_i EQ 2) then cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
      
    endfor
    cglegend, title=chunk_names,color=colors, location=[.2,.85], charsize=1
    stop
  endfor
  stop
  
end