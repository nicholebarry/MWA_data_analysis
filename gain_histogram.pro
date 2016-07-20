pro gain_histogram

allgains=abs((*cal.gain[0])[*,0])
for i=1,127 do allgains=[allgains,abs((*cal.gain[0])[*,i])]
result=histogram(allgains, binsize=binsize, reverse_indices=ri,/NAN,locations=locations,omax=omax)
y_arr=[result[0], result, result[N_elements(result)-1]]
x_arr=[locations[0],locations+binsize/2,omax]
cgplot, x_arr, y_arr, psym=10




end