pro uv_cubes_fhd

  ;uv model from a 2048x2048 uv extent
  grid_uv_model_2048 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_firstpass_map_4096_maxbaseline/grid_data/1061316296_uv_model_XX.sav', 'grid_uv_model')
  ;grid_uv_model_1024 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_farextent1_maxbaseline_512_sourcein/grid_data/1061316176_uv_model_XX.sav', 'grid_uv_model')
  
  ;uv model from a 1024x1024 uv extent
  grid_uv_model_1024 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_firstpass_map_2048_maxbaseline/grid_data/1061316296_uv_model_XX.sav', 'grid_uv_model')
  ;grid_uv_model_2048 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_farextent2_maxbaseline_512_sourcein/grid_data/1061316176_uv_model_XX.sav', 'grid_uv_model')
  
  ;uv dirty from a 2048x2048 uv extent
  grid_uv_dirty_2048 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_firstpass_map_4096_maxbaseline/grid_data/1061316296_uv_XX.sav', 'grid_uv')
  ;grid_uv_dirty_1024 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_farextent1_maxbaseline_512_sourcein/grid_data/1061316176_uv_XX.sav', 'grid_uv')
  
  ;uv dirty from a 1024x1024 uv extent
  grid_uv_dirty_1024 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_firstpass_map_2048_maxbaseline/grid_data/1061316296_uv_XX.sav', 'grid_uv')
  ;grid_uv_dirty_2048 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_farextent2_maxbaseline_512_sourcein/grid_data/1061316176_uv_XX.sav', 'grid_uv')
  
  grid_uv_res_1024 = grid_uv_dirty_1024 - grid_uv_model_1024
  grid_uv_res_2048 = grid_uv_dirty_2048 - grid_uv_model_2048
  
  ;array of counts of visibilities for 1024x1024
  vis_count_1024 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_firstpass_map_2048_maxbaseline/grid_data/1061316296_vis_count.sav', 'vis_count')
 ;vis_count_1024 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_farextent1_maxbaseline_512_sourcein/grid_data/1061316176_vis_count.sav', 'vis_count')
  
  ;array of counts of visibilities for 2048x2048
  vis_count_2048 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_firstpass_map_4096_maxbaseline/grid_data/1061316296_vis_count.sav', 'vis_count')
  ;vis_count_2048 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_farextent2_maxbaseline_512_sourcein/grid_data/1061316176_vis_count.sav', 'vis_count')
  
  ;weights for 2048x2048
  weights_uv_2048 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_firstpass_map_4096_maxbaseline/grid_data/1061316296_uv_weights_XX.sav', 'weights_uv')
  ;weights_uv_2048 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_farextent2_maxbaseline_512_sourcein/grid_data/1061316176_uv_weights_XX.sav', 'weights_uv')
  
  ;weights for 1024x1024
  weights_uv_1024 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_firstpass_map_2048_maxbaseline/grid_data/1061316296_uv_weights_XX.sav', 'weights_uv')
  ;weights_uv_1024 = getvar_savefile('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_sim_perfect_farextent1_maxbaseline_512_sourcein/grid_data/1061316176_uv_weights_XX.sav', 'weights_uv')
  
  ;uniform?
  uniform_model_1024 = grid_uv_model_1024 / vis_count_1024
  uniform_model_2048 = grid_uv_model_2048 / vis_count_2048
  ;lots of nans
  uniform_model_1024[where(FINITE(uniform_model_1024,/NAN))]=0
  uniform_model_2048[where(FINITE(uniform_model_2048,/NAN))]=0
  
  ;uniform?
  uniform_dirty_1024 = grid_uv_dirty_1024 / vis_count_1024
  uniform_dirty_2048 = grid_uv_dirty_2048 / vis_count_2048
  ;lots of nans
  uniform_dirty_1024[where(FINITE(uniform_dirty_1024,/NAN))]=0
  uniform_dirty_2048[where(FINITE(uniform_dirty_2048,/NAN))]=0
  
  ;uniform?
  uniform_res_1024 = grid_uv_res_1024 / vis_count_1024
  uniform_res_2048 = grid_uv_res_2048 / vis_count_2048
  ;lots of nans
  uniform_res_1024[where(FINITE(uniform_res_1024,/NAN))]=0
  uniform_res_2048[where(FINITE(uniform_res_2048,/NAN))]=0
  
  ;Ian wanted this, it doesn't make sense to me
    uniform_weights_1024 = weights_uv_1024 / vis_count_1024
  uniform_weights_2048 = weights_uv_2048 / vis_count_2048
  ;lots of nans
  uniform_weights_1024[where(FINITE(uniform_weights_1024,/NAN))]=0
  uniform_weights_2048[where(FINITE(uniform_weights_2048,/NAN))]=0
  
  ;apparent?
  apparent_model_1024 = grid_uv_model_1024 / weights_uv_1024
  apparent_model_2048 = grid_uv_model_2048 / weights_uv_2048
  ;lots of nans
  apparent_model_1024[where(FINITE(apparent_model_1024,/NAN))]=0
  apparent_model_2048[where(FINITE(apparent_model_2048,/NAN))]=0
  
  ;apparent?
  apparent_dirty_1024 = grid_uv_dirty_1024 / weights_uv_1024
  apparent_dirty_2048 = grid_uv_dirty_2048 / weights_uv_2048
  ;lots of nans
  apparent_dirty_1024[where(FINITE(apparent_dirty_1024,/NAN))]=0
  apparent_dirty_2048[where(FINITE(apparent_dirty_2048,/NAN))]=0
  
  ;apparent?
  apparent_res_1024 = grid_uv_res_1024 / weights_uv_1024
  apparent_res_2048 = grid_uv_res_2048 / weights_uv_2048
  ;lots of nans
  apparent_res_1024[where(FINITE(apparent_res_1024,/NAN))]=0
  apparent_res_2048[where(FINITE(apparent_res_2048,/NAN))]=0
  
  ;apparent_model_2048_half = apparent_model_2048[512:1535,512:1535]
  ;uniform_model_2048_half = uniform_model_2048[512:1535,512:1535]
  ;apparent_dirty_2048_half = apparent_dirty_2048[512:1535,512:1535]
  ;uniform_dirty_2048_half = uniform_dirty_2048[512:1535,512:1535]
  ;apparent_res_2048_half = apparent_res_2048[512:1535,512:1535]
  ;uniform_res_2048_half = uniform_res_2048[512:1535,512:1535]
    apparent_model_2048_half = apparent_model_2048[1024:3071,1024:3071]
  uniform_model_2048_half = uniform_model_2048[1024:3071,1024:3071]
  apparent_dirty_2048_half = apparent_dirty_2048[1024:3071,1024:3071]
  uniform_dirty_2048_half = uniform_dirty_2048[1024:3071,1024:3071]
  apparent_res_2048_half = apparent_res_2048[1024:3071,1024:3071]
  uniform_res_2048_half = uniform_res_2048[1024:3071,1024:3071]
  
  weights_uv_2048_half = weights_uv_2048[1024:3071,1024:3071]
    vis_count_2048_half = vis_count_2048[1024:3071,1024:3071]
    uniform_weights_2048_half = uniform_weights_2048[1024:3071,1024:3071]

  
  ;quick_image, abs(apparent_model_1024) - abs(apparent_model_2048_half),title='Apparent model difference, small - large extent', xtitle= 'U', ytitle='V', window=1
  ;quick_image, abs(apparent_model_1024),title='Apparent model 1024x1024', xtitle= 'U', ytitle='V', window=2
  ;quick_image, abs(uniform_model_1024) - abs(uniform_model_2048_half), title='Uniform model difference, small - large extent', xtitle= 'U', ytitle='V',window=3
  ;quick_image, abs(uniform_model_1024),title='Uniform model 1024x1024', xtitle= 'U', ytitle='V', window=4
  
  ;quick_image, abs(apparent_dirty_1024) - abs(apparent_dirty_2048_half),title='Apparent dirty difference, small - large extent', xtitle= 'U', ytitle='V', window=5
  ;quick_image, abs(apparent_dirty_1024),title='Apparent dirty 1024x1024', xtitle= 'U', ytitle='V', window=6
  ;quick_image, abs(uniform_dirty_1024) - abs(uniform_dirty_2048_half), title='Uniform dirty difference, small - large extent', xtitle= 'U', ytitle='V',window=7
  ;quick_image, abs(uniform_dirty_1024),title='Uniform dirty 1024x1024', xtitle= 'U', ytitle='V', window=8
  
  quick_image, atan(apparent_res_1024,/phase) - atan(apparent_res_2048_half,/phase),(FINDGEN(2048)-1024)*.5,(FINDGEN(2048)-1024)*.5, $
    title='Phase, Apparent res difference, small - large extent', xtitle= 'U', ytitle='V', data_range=[-.005,.005],/png,savefile='/nfs/mwa-00/h1/nbarry/data_apparent_res_diff_512_phase'
  quick_image, atan(apparent_res_1024,/phase),(FINDGEN(2048)-1024)*.5,(FINDGEN(2048)-1024)*.5, $
    title='Phase, Apparent res 1024x1024', xtitle= 'U', ytitle='V',data_range=[-3.14,3.14],/png,savefile='/nfs/mwa-00/h1/nbarry/apparent_res_512_phase_data'
  quick_image, abs(uniform_res_1024) - abs(uniform_res_2048_half),(FINDGEN(2048)-1024)*.5,(FINDGEN(2048)-1024)*.5, $
    title='Phase, Uniform res difference, small - large extent', xtitle= 'U', ytitle='V',data_range=[-5E-3,5E-3],/png,savefile='/nfs/mwa-00/h1/nbarry/uniform_res_diff_512_sim_zoom2_phase'
  quick_image, abs(uniform_res_1024),(FINDGEN(2048)-1024)*.5,(FINDGEN(2048)-1024)*.5, $
    title='Real, Uniform res 1024x1024', xtitle= 'U', ytitle='V',/png,savefile='/nfs/mwa-00/h1/nbarry/uniform_res_512_data'
    
    ;weights
   quick_image, abs(weights_uv_1024)/abs(vis_count_1024) - abs(weights_uv_2048_half)/abs(vis_count_2048),(FINDGEN(2048)-1024)*.5,(FINDGEN(2048)-1024)*.5, $
    title='Imaginary, Weights difference, small - large extent', xtitle= 'U', ytitle='V', data_range=[-1E-11,1E-11],/png,savefile='/nfs/mwa-00/h1/nbarry/weights_diff_512_sim_zoom2_imag'
    quick_image, abs(uniform_weights_1024) - abs(uniform_weights_2048_half),(FINDGEN(2048)-1024)*.5,(FINDGEN(2048)-1024)*.5, $
    title='Abs, Natural weights difference, small - large extent', xtitle= 'U', ytitle='V', data_range=[-1E-13,1E-13],/png,savefile='/nfs/mwa-00/h1/nbarry/natural_weights_diff_512_sim_zoom2_abs'
  quick_image, imaginary(vis_count_1024) - imaginary(vis_count_2048_half),(FINDGEN(2048)-1024)*.5,(FINDGEN(2048)-1024)*.5, $
    title='Imaginary, Vis count difference, small - large extent', xtitle= 'U', ytitle='V', data_range=[-1,1],/png,savefile='/nfs/mwa-00/h1/nbarry/vis_count_diff_512_sim_zoom2_imag'
    
  stop
  
end