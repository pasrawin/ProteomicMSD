from __future__ import division
import numpy as np
import math
from scipy import sparse
from scipy import spatial
from scipy.stats import norm
import timeit
import gc

### FUNCTIONS ###
from mNMF01_Tools import time_index, mz_index, smoothingtime_mat
flatten = lambda fl: [item for sublist in fl for item in sublist]

def h_construction(outfile, WithRTPrediction_, WithCosineH_, initRT_tuple, initRT_width,\
					V_mat, W_mat, noise_number, globalparam_list, iso_maxnumber, gaussian_width, cos_cutoff, output_file):
	gc.disable()
	V_shape, W_shape = V_mat.shape, W_mat.shape
	Hpeak_mean = 1*W_shape[0]/W_shape[1] #proofed per row
	if WithRTPrediction_:
		if W_shape[1] - noise_number != len(initRT_tuple):
			print 'Warning W_shape is wrong'
			exit()
		H_mat, initRT_correct = h_prediction(initRT_tuple, initRT_width, \
											Hpeak_mean, noise_number, globalparam_list, gaussian_width)
		# H_mat = H_mat.todense()
	else:
		H_mat = h_constant(V_shape, W_shape, Hpeak_mean, noise_number)
		initRT_correct = None

	if WithCosineH_:
		H_size_rowcut = 250
		V_mat = V_mat.tocsr()
		W_mat = W_mat.tocsc()
		H_mat = H_mat.tocsr()
		
		if H_mat.shape[0] > H_size_rowcut:
			H_mat_list, H_cosmat_list = [], []			
			for row in np.arange(0,H_mat.shape[0],H_size_rowcut):				
				if row + H_size_rowcut >= H_mat.shape[0]:
					row_end = H_mat.shape[0]
				else:
					row_end = row + H_size_rowcut
				H_mat_out, H_cosmatdub_out = h_cos(V_mat, W_mat[:,row:row_end], H_mat[row:row_end,:], iso_maxnumber, cos_cutoff)
				H_mat_list.append(H_mat_out)
				del(H_mat_out)
				H_cosmat_list.append(H_cosmatdub_out)
				del(H_cosmatdub_out)

			H_cosmat = sparse.vstack(H_cosmat_list)
			H_cosmat = H_cosmat.tocoo()
			del(H_cosmat_list)		
			H_mat = sparse.vstack(H_mat_list)
			del(H_mat_list)
		else:
			H_mat, H_cosmat = h_cos(V_mat, W_mat, H_mat, iso_maxnumber, cos_cutoff)	
	else:
		H_cosmat = None

	if noise_number > 0:
		noise_mat = np.zeros((noise_number,H_mat.shape[1]))
		noise_mat += Hpeak_mean/H_mat.shape[1]
		try:
			H_mat = sparse.vstack([H_mat, noise_mat])
			H_mat.tocoo()
		except:
			H_mat = np.vstack([H_mat, noise_mat])
			# change to sparse if better for memory
			if np.count_nonzero(H_mat)*3 < H_mat.size:
				H_mat = sparse.coo_matrix(H_mat)
	sparse.save_npz(str(output_file)+'_Hmat_init.npz', H_mat)
	return H_mat, initRT_correct, H_cosmat

### subfunctions
def h_prediction(initRT_tuple, initRT_width, Hpeak_mean, noise_number, globalparam_list, gaussian_width):
	gc.disable()
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
	gc.disable()
	def cosine_calculation(V_cut, W_cut):
		np.warnings.filterwarnings('ignore')
		similarity = spatial.distance.cosine(V_cut, W_cut)
		similarity = 1 - similarity
		return similarity
	
	H_coslist_k, H_coslist_cos, H_coslist_t = [], [], []
	H_matrow, H_matcol = H_mat.shape[0], H_mat.shape[1]
	peaknumber_cut, extend_add = 4, [-2,-1,0.5]
	for k in np.arange(H_matrow): # use H no noise
		# print k, timeit.default_timer()
		k_count, v_count = 0, 0
		W_keep = np.sort(np.nonzero(W_mat[:,k])[0][:peaknumber_cut]) # sorted non zero indices
		V_keep = V_mat[W_keep,:].todense()#[:peaknumber_cut]
		W_cut = W_mat[W_keep,k].todense()#[:peaknumber_cut]
		W_keepextend = W_keep.copy()	
		for n in extend_add:
			front_extend = W_keep[0] + int(n*(W_keep[1] - W_keep[0]))
			if front_extend >= 0: # ignore if negative
				W_keepextend = np.insert(W_keepextend, 0, front_extend)
				V_keepextend = V_mat[W_keepextend,:].todense()
				W_cutextend = W_mat[W_keepextend,k].todense()
		t_range = np.nonzero(H_mat[k,:])[1]
		t_min, t_max = t_range[0], t_range[-1]
		for t in np.arange(H_matcol):
			if t_min <= t <= t_max:
				V_cut = V_keep[:,t]
				# pass if V_cut is nothing or has only one peak
				if np.sum(V_cut) != 0 and np.count_nonzero(V_cut) >= 2 and V_cut[0] != 0:
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
					similarity_score = (np.median(similarity_list) + np.mean(similarity_list))/2
					del(similarity_list)
					if similarity_score > cos_cutoff:
						H_coslist_cos.append(similarity_score)
						H_coslist_t.append(t)
						k_count += 1
		H_coslist_k.extend([k]*k_count)
		# print '*', v_count
		
	H_cosmatdub = sparse.csr_matrix((H_coslist_cos, (H_coslist_k, H_coslist_t)), shape=(H_matrow, H_matcol))
	del(H_coslist_cos)
	del(H_coslist_t)
	del(H_coslist_k)
	try:
		# H_mat is sparse
		H_mat = H_cosmatdub.multiply(H_mat)
	except:
		# H_mat is dense
		H_mat = np.multiply(H_mat, H_cosmatdub.todense())
	return H_mat, H_cosmatdub
			
