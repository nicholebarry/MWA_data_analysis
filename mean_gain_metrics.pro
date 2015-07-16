pro mean_gain_metrics, obsid


  ;obsid='1061316296'
  ;obsid='1064761640'
  calibration_polyfit=2
  
  ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/metadata/' + obsid + '_obs.sav'
  ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/metadata/' + obsid + '_params.sav'
  ;restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_no_cable_cal_std/calibration/' + obsid + '_cal.sav'
  ;no_cable_cal_xx=*cal.gain_residual[0]
  ;no_cable_cal_yy=*cal.gain_residual[1]
  
  restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/' + obsid + '_obs.sav'
  restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/metadata/' + obsid + '_params.sav'
  restore, '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/calibration/' + obsid + '_cal.sav'
  cable_cal_xx=*cal.gain_residual[0]
  cable_cal_yy=*cal.gain_residual[1]
  
  cc_cal_res_avg=FLTARR(2)
  ;no_cal_res_avg=FLTARR(2)
  cc_cal_res_avg[0]=Mean(Abs(cable_cal_xx),/NAN)
  cc_cal_res_avg[1]=Mean(Abs(cable_cal_yy),/NAN)
  ;no_cal_res_avg[0]=Mean(Abs(no_cable_cal_xx),/NAN)
  ;no_cal_res_avg[1]=Mean(Abs(no_cable_cal_yy),/NAN)
  
  print, cc_cal_res_avg
  
  
end