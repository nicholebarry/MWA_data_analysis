pro export_1dkpower

  ;coarse_harm_width = 3
  ;wedge_angles=120d
  n_cubes=2
  
  cut = STRARR(n_cubes)
  ;cut[1]='wedge_cut_plus_res_cut'
  ;cut[1] = 'beardsley_thesis_list'
  cut[1] = 'Aug23'
  ;cut='Beardsley2016_zenith'
  ;cut[0] = 'beardsley_thesis_list'
  cut[0] = 'Aug23'
  
  ;chans='ch0-95' & band='low'
  chans='ch96-191' & band='high'
  ;chans='ch48-143' & band='mid'
  
  basefile = STRARR(n_cubes)
  ;basefile[1]=['/nfs/mwa-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps_jan2016/Combined_obs_'+cut[1]+'_cubeXX__even_odd_joint_'+chans+'_']
  basefile[0]=['/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_notimeavg_delaymodel/ps/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+chans+'_']
  basefile[1]=['/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_notimeavg_delaymodel_quarter/ps/Combined_obs_'+cut[0]+'_cubeXX__even_odd_joint_'+chans+'_']
  
  cubes=['res']
  pols=['xx','yy']
  
  endfile = STRARR(n_cubes)
  ;endfile[1]='_bh_dencorr_no_120deg_wedge_cbw3_kperplambda10-60_kpar0.15-200_1dkpower.idlsave'
  endfile[1]='_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_kpar0.15-200_1dkpower.idlsave'
  endfile[0]='_averemove_bh_dencorr_no_120deg_wedge_kperplambda10-60_kpar0.15-200_1dkpower.idlsave'
  
  limit_percent=0.9772 ;2 sigma
  hubble_param=0.71
  
  for j=0,1 do begin
    for cube_i=0,n_cubes-1 do begin
    
      file=basefile[cube_i]+cubes+'_'+pols[j]+endfile[cube_i]
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
      
      power=[0.,power] ; first element a placeholder to match size of k_edges
      weights=[0.,weights]
      foo=transpose([[k_edges],[power],[weights]])
      outfile='~/'+band+'_'+cubes+'_'+pols[j]+'.txt'
      ;textfast,foo,[header,'',''],file_path=outfile,/write
      
      if (cube_i EQ 0) then cgPS_Open,'/nfs/mwa-00/h1/nbarry/'+pols[j]+'_delaymodel2_quarter_nb_limit_comparison'+chans+'.png',/quiet,/nomatch
      if cube_i EQ 0 then cgplot, k_edges,limits,/xlog,/ylog,psym=10, xrange=[.08,1.5], ytitle='$\Delta$$\up2$ (mK$\up2$)', $
        xtitle='k (h Mpc$\up-1$)', title=pols[j] + ' Limit Comparisons ' + chans, charsize =1.5
      if cube_i EQ 1 then cgoplot, k_edges,limits,/xlog,/ylog,psym=10, xrange=[.08,1.5],color='blue'
      if cube_i EQ 1 then cglegend, title=['no delay filter','delay filter'],color=['black','blue'], location=[.2,.85], charsize=1
      if (cube_i EQ 1) then cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
     ; stop
    ;textfast,foo,file_path=outfile,/write
    endfor
  endfor
  
end