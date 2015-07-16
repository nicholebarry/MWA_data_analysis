pro runtemp

  print, 'Normal, nonpointing'
  chi_squared_polyfit,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/calibration',5,unreduced=1,nopointing=1
  ;print, 'Normal, nonpointing, modefit 150'
  ;chi_squared_polyfit,'/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_pointing_May2015/calibration',5,unreduced=1,nopointing=1,mode_calc=1
  
  print, 'One quadratic'
  chi_squared_polyfit,'Aug23_onequad_polyscaled_90150/forinput',5,unreduced=1
  ;print, 'One quadratic 150'
  ;chi_squared_polyfit,'Aug23_onequad_polyscaled_90150/forinput',5,unreduced=1,mode_calc=1
  
  print, 'Two split linears'
  chi_squared_polyfit,'Aug23_pointing_nodigjump_v2_plusmodepointing_frombp',5,unreduced=1, split_jump=1
  ;print, 'Two split linears 150'
  ;chi_squared_polyfit,'Aug23_pointing_nodigjump_v2_plusmodepointing_frombp',5,unreduced=1, split_jump=1, mode_calc=1
  
  print, 'Two split quadratics'
  chi_squared_polyfit,'Aug23_cross_polyfit_quad_data',8,unreduced=1, split_jump=1
  ;print, 'Two split quadratics 150'
  ;chi_squared_polyfit,'Aug23_cross_polyfit_quad_data',8,unreduced=1, split_jump=1,mode_calc=1
  
  print, 'Two quads by pointing, dig removed'
  chi_squared_polyfit,'Aug23_cross_polyfit_quad_nodig_data',9,unreduced=1, split_jump=1
  ;print, 'Two quads by pointing, dig removed 150'
  ;chi_squared_polyfit,'Aug23_cross_polyfit_quad_nodig_data',9,unreduced=1, split_jump=1,mode_calc=1
  
end