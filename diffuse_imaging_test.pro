pro diffuse_imaging_test

	restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_presidelobe/vis_data/1061316296_vis_XX.sav'
	restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_presidelobe/vis_data/1061316296_vis_model_XX.sav'
	vis_res = *vis_ptr - *vis_model_ptr
	restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_decon_July2016_presidelobe/vis_data/1061316296_flags.sav'
	
end