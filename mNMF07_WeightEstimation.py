from __future__ import division
import numpy as np

def weight_estimation(exp_mat, prot_peptcount, globalparam_list):
	# for param in globalparam_list:
	mm, mrange_min, mrange_max, mz_range, MZ_SCALE, \
	tt, gradient_starttime, gradient_endtime, gradient_time, TIME_SCALE, \
	window, shift =  globalparam_list[1]

	emat = exp_mat.todense()
	l1 = np.sum(emat, axis=1)
	l2 = np.sqrt(np.sum(np.square(emat), axis=1))
	del(emat)
	weight_est = (1/(mz_range**0.5))*np.sum(gradient_time**0.5 - (l1/(l2+1e-6)))/(gradient_time**0.5 -1)
	return np.around(weight_est, 5)