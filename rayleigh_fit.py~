import idlsave
from scipy import optimize

f = idlsave.read('seti_all_3D_20_thesis_evenodd2_minustwo.sav')

popt, pcov = optimize.curve_fit(f, x[3:200], y[3:200]/100000000)
print popt


popt, pcov = optimize.curve_fit(f, x[3:150], y[3:150]/1000000)
print popt
