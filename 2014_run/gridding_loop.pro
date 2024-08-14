pro gridding_loop

  ;When the gridding_loop is called with an obs ID passed via a bash script, 
  ; parse command line args to get the obs ID. Otherwise, run with a hardcoded
  ; obs ID for testing purposes.
  compile_opt strictarr
  args = Command_Line_Args(count=nargs)
  IF keyword_set(args) then begin
    obs_id = args[0]
  endif else begin
    obs_id= '1088716176' ;A minustwo pointing obs ID
  endelse

  ;Read in the various save files required for gridding input, including
  ; visibility (vis_ptr, a pointer per pol with Nfreq x Nbaselines), 
  ; vis_weights (vis_weights, a pointer per pol with Nfreq x Nbaselines), 
  ; psf (psf, a struct with the beam_arr: pointer to a pointer array of Npol x Nfreq x Nbaselines, 
  ; which themselves are pointer arrays Nuvres x Nuvres), 
  ; params (params, struct of various Nbaseline metadata arrays)

  sav_dir = '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_gd_woden_pyfhd/'

  ;vis_weights pointer array
  restore, sav_dir + 'vis_data/' + obs_id + '_flags.sav'

  ;params struct
  restore, sav_dir + 'metadata/' + obs_id + '_params.sav'

  ;psf struct, hardcode for now
  restore, '/fred/oz048/MWA/CODE/FHD/fhd_nb_data_pointing_beam/beams/gauss_beam_pointing-2.sav'

  ;obs struct. In order to avoid being overwritten later by obs structs saved alongside
  ;visibilities, read it in under a different variable name
  obs_out = getvar_savefile(sav_dir + 'metadata/' + obs_id + '_obs.sav','obs')

  ;Various variables required for running gridding in the appropriate manner, including
  ; bi_use (pointer per pol of the even, odd indices)
  ; model_return (return the model vis)
  ; perserve_visibilities (don't overwrite the visibilities)
  ; variance_flag (grid the variances)
  ; weights_flag (grid the weights)
  restore, sav_dir + obs_id + '_variables.sav'

  n_pol = obs_out.n_pol
  pol_names = ['XX','YY']
  iter_names = ['even','odd']
  n_iter = 2 ;Hardcode to iterate over even,odd indicies
  n_freq = obs_out.n_freq
  freq_use=(*obs_out.baseline_info).freq_use ;Bool array for flagged/unflagged freqs
  n_avg = 2. ;Hardcode to grid two freqs together for 160kHz res uvf cubes, FHD/epp standard
  freq_bin_i2=Floor(lindgen(n_freq)/n_avg)
  nf=Max(freq_bin_i2)+1L

  FOR pol_i=0,n_pol-1 DO BEGIN

    ;vis_ptr for the specific pol
    restore, sav_dir + 'vis_data/' + obs_id + '_vis_' + pol_names[pol_i] + '.sav'

    ;vis_model_ptr for the specific pol
    restore, sav_dir + 'vis_data/' + obs_id + '_vis_model_' + pol_names[pol_i] + '.sav'

    FOR iter_i=0,n_iter-1 DO BEGIN

      FOR fi=0L,nf-1 DO BEGIN

        ;Get the two frequencies to be gridded together. Ignore flagged frequencies
        fi_use=where((freq_bin_i2 EQ fi) AND (freq_use GT 0),nf_use)

        variance_holo=variance_flag ;initialize
        weights_holo=weights_flag ;initialize

        ;Skip if no unflagged frequencies, otherwise grid
        IF nf_use EQ 0 THEN BEGIN
          n_vis=0 
        ENDIF ELSE BEGIN
          dirty_UV=visibility_grid(vis_ptr,vis_weights[pol_i],obs_out,0,psf,params,fi_use=fi_use,bi_use=*bi_use[iter_i],$
              polarization=pol_i,weights=weights_holo,variance=variance_holo,silent=1,mapfn_recalculate=0,$
              model_ptr=vis_model_ptr,n_vis=n_vis,/preserve_visibilities,model_return=model_return)
          save, dirty_UV, model_return, variance_holo, weights_holo, n_vis, filename = sav_dir + obs_id + '_output_' + string(fi,format='(I03)') + '_' + pol_names[pol_i] + '_' + iter_names[iter_i] + '.sav'
        ENDELSE
      ENDFOR
    ENDFOR
  ENDFOR

end