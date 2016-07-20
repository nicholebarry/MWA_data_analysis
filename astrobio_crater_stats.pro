pro astrobio_crater_stats

  UA_upper=[160, 40, 0,0]
  MA_upper=[600,150,25,0]
  MA_lower=[160,40,0,0]
  LA_lower=[600,150,25,0]
  LA_upper=[1600,400,67,0]
  UH_lower=[1600,400,67,0]
  UH_upper=[3000,750,125,0]
  LH_lower=[3000,750,125,0]
  LH_upper=[4800,1200,200,25]
  UN_upper=[0,0,400,100]
  UN_lower=[0,0,200,25]
  MN_upper=[0,0,400,200]
  MN_lower=[0,0,0,100]
  LN_upper=[0,0,0,200]
  
  age_upper=LIST(LN_upper,MN_upper,UN_upper,LH_upper,UH_upper,LA_upper,MA_upper,UA_upper]
  age_lower=LIST(LN_lower,MN_lower,UN_lower,LH_lower,UH_lower,LA_lower,MA_lower,UA_lower)
  
  color_array=['green','blue','purple','brown','red','light blue','light purple','light green']
  
  for age_i=0,7 do begin
    largest_cummulative = Total(age_upper[age_i])
    largest_cummulative_lower = Total(age_lower[age_i])
    second_largest_cummulative = Total((age_upper[age_i])[0:2])
    second_largest_cummulative_lower = Total((age_lower[age_i])[0:2])
    third_largest_cummulative =Total((age_upper[age_i])[0:1])
    third_largest_cummulative_lower =Total((age_lower[age_i])[0:1])
    lowest_cummulative =age_upper[age_i])[0]
    lowest_cummulative_lower =age_lower[age_i])[0]
    
    if age_i EQ 0 then cgplot, [largest_cummulative, second_largest_cummulative, third_largest_cummulative, lowest_cummulative],[16,5,2,1], color=color_array[age_i] else $
      cgoplot, [largest_cummulative, second_largest_cummulative, third_largest_cummulative, lowest_cummulative],[16,5,2,1], color=color_array[age_i]
            
  endfor
  
end