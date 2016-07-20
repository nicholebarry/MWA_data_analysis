PRO ucla_cal_format, cal_sav_file_path, time

  ;************************
  ;Takes a calibration save file and turns it into the UCLA format text file.
  ;
  ;Text columns are:
  ;1. antenna name
  ;2. antenna index
  ;3. frequency (in MHz)
  ;4. polarization (EE, EN, NE, & NN)
  ;5. time (in JD)
  ;6. the real part of the complex gain
  ;7. the imaginary part of the complex gain
  ;8. antenna flag (0=unflagged, 1=flagged)
  ;
  ;And the header will contain:
  ;1. the program used to produce the calibration solutions (e.g. RTS)
  ;2. the convention for applying the calibration solutions to date (e.g. divide uncalibrated data
  ;by the solutions to obtain calibrated data).
  ;3. The last line of the header will always be the column names, written as:
  ;# ANT NAME, ANT INDEX, FREQ (MHZ), POL, TIME (JD), RE(GAIN), IM(GAIN), FLAG
  ;
  ;Created by Nichole Barry, July 9 2016
  ;
  ;************************

  ;Reading in the cal file (with a shortcut for the MIT cluster)
  folder_test = file_test(cal_sav_file_path)
  if folder_test EQ 0 then begin
    cal_sav_file_path_full = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/'+cal_sav_file_path
    folder_test = file_test(cal_sav_file_path_full)
    if folder_test EQ 0 then begin
      message, 'Save file not found in '+cal_sav_file_path+' or '+cal_sav_file_path_full
    endif else cal_sav_file_path = cal_sav_file_path_full
  endif
  
  ;Get overall directory output structure, put in location of obs sav file for flags
  obsid = file_basename(cal_sav_file_path)
  obsid = strmid(obsid,0,10)
  obs_sav_file_path = file_dirname(file_dirname(cal_sav_file_path)) + '/metadata/' + obsid + '_obs.sav'
  
  cal = getvar_savefile(cal_sav_file_path,'cal')
  obs = getvar_savefile(obs_sav_file_path,'obs')
  
  n_tile = cal.n_tile
  n_freq = cal.n_freq
  
  ;Getting tile/freq flags. By default, flagged data is logical 0 and unflagged data is logical 1, and needs to be switched.
  tile_use=(*obs.baseline_info).tile_use
  freq_use=(*obs.baseline_info).freq_use
  tile_flags = where(tile_use EQ 0, n_count)
  tile_use = INTARR(n_tile)
  if n_count GT 0 then tile_use[tile_flags]=1
  freq_flags = where(freq_use EQ 0, n_count)
  freq_use = INTARR(n_freq)
  if n_count GT 0 then freq_use[freq_flags]=1
  
  ;Getting antenna names in sequential order
  ant_name=[]
  for tile_i=10,160,10 do for tile_j=1,8 do ant_name = [ant_name,tile_i+tile_j]
  
  ;Getting antenna index in sequential order
  ant_index = INDGEN(n_tile)
  
  ;Getting frequencies of cal file
  freq = cal.freq/10^6.
  
  if cal.n_pol EQ 2 then pol = ['EE','NN'] else pol=['EE','EN','NE','NN'] ;order of the 4 pol is unchecked
  
  ;initialize
  real_gain=[]
  imag_gain=[]
  ant_name_full_array=[]
  ant_index_full_array=[]
  freq_full_array=[]
  pol_full_array=[]
  tile_use_full_array=[]
  freq_use_full_array=[]
  
  
  
  for pol_i=0, cal.n_pol-1 do begin
  
    for freq_i=0,n_freq-1 do begin
      real_gain = [real_gain,reform(real_part((*cal.gain[pol_i])[freq_i,*]))]
      imag_gain = [imag_gain,reform(imaginary((*cal.gain[pol_i])[freq_i,*]))]
      ant_name_full_array = [ant_name, ant_name_full_array]
      ant_index_full_array = [ant_index, ant_index_full_array]
      tile_use_full_array = [tile_use_full_array, tile_use]
      
      for tile_i=0,n_tile-1 do begin
        freq_full_array = [freq_full_array,freq[freq_i]]
        pol_full_array = [pol_full_array,pol[pol_i]]
        freq_use_full_array = [freq_use_full_array, freq_use[freq_i]]
      endfor ;end tile loop
    ;freq_full_array[(pol_i+1)*128*384+(freq_i+1)*128:(pol_i+1)*128*384+(freq_i+1)*128+128]=freq[freq_i]
    endfor ;end freq loop
    
  endfor ;end pol loop
  time_full_array=DBLARR((size(real_gain))[1])
  time_full_array[*] = time
  
  flag_full_array = ULONG(LOGICAL_OR(tile_use_full_array, freq_use_full_array))
  
  wh_flags = where(flag_full_array EQ 1, n_count)
  if n_count GT 0 then begin
    real_gain[wh_flags] = 0
    imag_gain[wh_flags] = 0
  endif
  
  data = transpose([[strtrim(ant_name_full_array,2)],[strtrim(ant_index_full_array,2)],[strtrim(freq_full_array,2)],[pol_full_array],$
    [string(time_full_array,format='(f0.5)')],[strtrim(real_gain,2)],[strtrim(imag_gain,2)],[strtrim(flag_full_array,2)]])
    
  line1 = '# Program of origin: FHD'
  line2 = '# Convention: Divide uncalibrated data by these gains to obtain calibrated data.'
  line3 = '# Metadata: Observation '+obsid+', calibration/flagging same throughout the 56 time samples.'
  line4 = '# ANT NAME, ANT INDEX, FREQ (MHZ), POL, TIME (JD), RE(GAIN), IM(GAIN), FLAG'
  
  header = transpose([[line1],[line2],[line3],[line4]])
  
  file_path = file_dirname(cal_sav_file_path) + '/' + obsid + '_cal.txt'

  openw, lun, file_path, /get_lun
  printf, lun, header, data
  close, lun
  free_lun, lun
  
end
