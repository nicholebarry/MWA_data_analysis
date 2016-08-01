#import idlsave
from scipy.io import readsav
from scipy import optimize
import numpy as n

def f(x,a,s):
	#s = n.argmax(y) + 3
	#s=39.
	#f = a * x/s**2 * n.exp([-xx**2 for xx in x]/(2*s**2))
	f = a * x/s**2 * n.exp(-x**2/(2*s**2))
	return f

#f = idlsave.read('seti_all_3D_20_thesis_evenodd2_minustwo.sav')
binned_y = readsav('/nfs/eor-00/h1/nbarry/vis_res/thesis/residual/seti_binned_diff_thesis_residual_total.sav')
y = binned_y['binned_diff_total']

x = n.array([i for i in range(len(y))])

popt, pcov = optimize.curve_fit(f, x[3:200], y[3:200]/10000000)
print popt


popt, pcov = optimize.curve_fit(f, x[3:150], y[3:150]/1000000)
print popt
