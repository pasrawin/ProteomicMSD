from __future__ import division
import numpy as np
from scipy import sparse
from scipy import spatial
from scipy.stats import norm

### FUNCTIONS ###
from mNMF01_Tools import time_index, mz_index, smoothingtime_mat
flatten = lambda fl: [item for sublist in fl for item in sublist]

def h_construction(outfile, WithRTPrediction_, WithCosineH_, initRT_tuple, initRT_width,\
					V_mat, W_mat, noise_number, globalparam_list, iso_maxnumber, gaussian_width, cos_cutoff):
	
	V_shape, W_shape = V_mat.shape, W_mat.shape
	Hpeak_mean = 1*W_shape[0]/W_shape[1] #proofed per row
	if WithRTPrediction_:
		if W_shape[1] - noise_number != len(initRT_tuple):
			print 'Warning W_shape is wrong'
			exit()
		H_mat, initRT_correct = h_prediction(initRT_tuple, initRT_width, \
											Hpeak_mean, noise_number, globalparam_list, gaussian_width)
		H_mat = H_mat.todense()
	else:
		H_mat = h_constant(V_shape, W_shape, Hpeak_mean, noise_number)
		initRT_correct = None

	if WithCosineH_:
		H_mat, H_cosmat, W_nonzerolist = h_cos(V_mat, W_mat, H_mat, iso_maxnumber, cos_cutoff)
	else:
		H_cosmat = None

	if noise_number > 0:
		noise_mat = np.zeros((noise_number,H_mat.shape[1]))
		noise_mat += Hpeak_mean/H_mat.shape[1]
		H_mat = np.vstack([H_mat, noise_mat])
		# change to sparse if better for memory
		if np.count_nonzero(H_mat)*3 < H_mat.size:
			H_mat = sparse.coo_matrix(H_mat)
	return H_mat, initRT_correct, H_cosmat

### subfunctions
def h_prediction(initRT_tuple, initRT_width, Hpeak_mean, noise_number, globalparam_list, gaussian_width):
	for param in globalparam_list:
		mm, mrange_min, mrange_max, mz_range, MZ_SCALE, \
        tt, gradient_starttime, gradient_endtime, gradient_time, TIME_SCALE, \
        window, shift =  globalparam_list[1]

	peakside_index = int(gaussian_width/2 * tt) #before smooth
	peakwidth_index = int((peakside_index * 2) + 1)
	peakSD_index = peakwidth_index / 4
	rt_index = time_index([x[2] for x in initRT_tuple], gradient_starttime, tt)
	keeprow, keepcol, keepdata = [], [], []
	for idx, rt in enumerate(rt_index):
		if rt < 0:
			rt = 0
		elif rt > TIME_SCALE-1:
			rt = TIME_SCALE-1
		H_row = np.zeros(TIME_SCALE)
		left, right = rt-peakside_index, rt+peakside_index
		peak_at = np.linspace(left, right, peakwidth_index)
		peak_int = norm.pdf(peak_at, rt, peakSD_index) * Hpeak_mean
		# delete the plot which is located out of range
		mask = [(peak_at>=0) & (peak_at<TIME_SCALE)]
		peak_at = peak_at[tuple(mask)]
		peak_int = peak_int[tuple(mask)]

		nonzero = np.arange(peak_at[0],peak_at[-1]+1)
		keeprow.append([idx]*len(nonzero))
		keepcol.append(peak_at)
		keepdata.append(peak_int)

	keeprow, keepcol, keepdata = flatten(keeprow), flatten(keepcol), flatten(keepdata)
	H_mat = sparse.coo_matrix((keepdata,(keeprow,keepcol)),shape=(len(initRT_tuple),TIME_SCALE))
	H_mat = smoothingtime_mat(H_mat, window, shift)

	initRT_correct = [x[2] for x in initRT_tuple]
	initRT_correct = [gradient_starttime if x < gradient_starttime else x for x in initRT_correct]
	initRT_correct = [gradient_endtime if x > gradient_endtime else x for x in initRT_correct]

	initRT_correct_keep = []
	for ind, (a, b, c) in enumerate(initRT_tuple):
		initRT_correct_keep.append((a,b,c,initRT_correct[ind]))
	del(initRT_correct)

	return H_mat, initRT_correct_keep

def h_constant(V_shape, W_shape, Hpeak_mean, noise_number):
	H_mat = np.zeros((W_shape[1]-noise_number, V_shape[1]))
	H_mat += Hpeak_mean/H_mat.shape[0]
	return H_mat

def h_cos(V_mat, W_mat, H_mat, iso_maxnumber, cos_cutoff):
	def cosine_calculation(V_cut, W_cut):
		np.warnings.filterwarnings('ignore')
		similarity = spatial.distance.cosine(V_cut, W_cut)
		similarity = 1 - similarity
		return similarity

	H_coslist, W_nonzerolist = [], []
	H_matrow, H_matcol = H_mat.shape[0], H_mat.shape[1]
	W_mat = W_mat.tocsc()
	V_mat = V_mat.tocsr()
	peaknumber_cut, extend_add = 4, [-2,-1,0.5]
	for k in np.arange(H_matrow): # use H no noise
		W_keep = np.sort(np.nonzero(W_mat[:,k])[0][:iso_maxnumber]) # sorted non zero indices
		W_nonzerolist.append(W_keep)
		W_keep = W_keep[:peaknumber_cut]
		V_keep = V_mat[W_keep,:].todense()[:peaknumber_cut]
		W_cut = W_mat[W_keep,k].todense()[:peaknumber_cut]

		W_keepextend = W_keep.copy()	
		for n in extend_add:
			front_extend = W_keep[0] + int(n*(W_keep[1] - W_keep[0]))
			if front_extend >= 0: # ignore if negative
				W_keepextend = np.sort(np.insert(W_keepextend, 0, front_extend))
				V_keepextend = V_mat[W_keepextend,:].todense()
				W_cutextend = W_mat[W_keepextend,k].todense()

		t_range = np.nonzero(H_mat[k,:])[1]
		t_min, t_max = t_range[0], t_range[-1]
		for t in np.arange(H_matcol):
			if t_min <= t <= t_max:
				V_cut = V_keep[:,t]
				V_cutextend = V_keepextend[:,t]
				similarity_list = []				
				similarity_list.append(cosine_calculation(V_cutextend, W_cutextend)) # 
				similarity_list.append(cosine_calculation(V_cut, W_cut)) # 3 + back
				arg_back, arg_forth = int(t-1), int(t+1)
				arg_back = 0 if arg_back < 0 else int(t-1)
				arg_forth = V_keep.shape[1]-1 if arg_forth > V_keep.shape[1]-1 else int(t+1)

				similarity_list.append(cosine_calculation(V_keep[:,arg_back], W_cut))
				similarity_list.append(cosine_calculation(V_keep[:,arg_forth], W_cut))
				similarity_list.append(cosine_calculation(np.sum(V_keep[:,arg_back:arg_forth+1],axis=1), W_cut))
				similarity_list = np.nan_to_num(similarity_list)
				H_coslist.append((np.median(similarity_list) + np.mean(similarity_list))/2)
			else:
				H_coslist.append(0)

	H_cosmat = np.asarray(H_coslist).reshape((H_matrow, H_matcol))
	# print ' chk H_mat non zeros', np.count_nonzero(H_mat)
	# print ' process dub cosine less than or equal:', cos_cutoff
	H_cosmatdub = H_cosmat.copy()
	H_cosmatdub[H_cosmatdub<cos_cutoff] = 0
	H_mat = np.multiply(H_mat, H_cosmatdub)
	# print ' chk H_cosmat non zeros', np.count_nonzero(H_cosmat), H_cosmat.shape
	del(H_cosmatdub)
	return H_mat, sparse.coo_matrix(H_cosmat), W_nonzerolist
			