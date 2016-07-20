pro source_list_reader

catalog_path='/nfs/mwa-09/r1/djc/EoR2013/Aug23/calibration_sim/fhd_nb_sim_perfect_cal_noeor_ones_dimcalsources_nod_notileflag/output_data/1061316176_source_list.txt'

readcol, catalog_path,id,x_loc,y_loc,RA,Dec,SN,radius,avg_beam,XX_apparent,YY_apparent,Stokes_I_fit,Stokes_Q_fit,Stokes_I_res,Stokes_Q_res,Extended

print, mean(avg_beam)
print, max(avg_beam)
print, min(avg_beam)

end