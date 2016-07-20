filename='/nfs/eor-00/h1/nbarry/MWA/IDL_code/FHD/instrument_config/0_bandpass.txt'
textfast,bandpass_saved_sol,/read,file=filename

cgPS_Open,'/nfs/mwa-00/h1/nbarry/General_plots/savedcable_bp.png',/quiet,/nomatch

cgplot, bandpass_saved_sol[0,*]/1E6, bandpass_saved_sol[1,*], color='blue', xtitle='Frequency (MHz)',ytitle='Gain', psym = 2, symsize=0.2, $
  yrange=[0.8,1.2], xrange=[163,202] ;[bandpass_saved_sol[0,0],bandpass_saved_sol[0,383]]
cgoplot, bandpass_saved_sol[0,*]/1E6, bandpass_saved_sol[3,*], color='red', xtitle='Frequency [MHz]',ytitle='Gain', psym = 2, symsize=0.2
cgoplot, bandpass_saved_sol[0,*]/1E6, bandpass_saved_sol[5,*], color='green', xtitle='Frequency [MHz]',ytitle='Gain', psym = 2, symsize=0.2
cgoplot, bandpass_saved_sol[0,*]/1E6, bandpass_saved_sol[7,*], color='purple', xtitle='Frequency [MHz]',ytitle='Gain', psym = 2, symsize=0.2
cgoplot, bandpass_saved_sol[0,*]/1E6, bandpass_saved_sol[9,*], color='yellow', xtitle='Frequency [MHz]',ytitle='Gain', psym = 2, symsize=0.2
cgoplot, bandpass_saved_sol[0,*]/1E6, bandpass_saved_sol[11,*], color='black', xtitle='Frequency [MHz]',ytitle='Gain', psym = 2, symsize=0.2

n_cable=6
Title = STRARR(6)
cable_length=[90,150,230,320,400,524]
FOR cable_i=0,n_cable-1 DO Title[cable_i]=String(format='(A,"m cables")',Strn(Round(cable_length[cable_i])))
cgLegend, Title=Title, Color=['blue','red','green','purple','yellow','black'],Psym=Replicate(2,n_cable),Length=0.0,Location=[0.7,0.85]

cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage