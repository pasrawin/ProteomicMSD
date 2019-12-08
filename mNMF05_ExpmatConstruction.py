from __future__ import division
import numpy as np
import pandas as pd
from scipy import sparse

### Global functions ###
from mNMF01_Tools import time_index, mz_index, smoothingtime_mat, rescale_mat

def expmat_construction(exp_file, exp_paramlist, charge_list):
	mslev = 1
	for param in exp_paramlist:
		mm, mrange_min, mrange_max, mz_range, MZ_SCALE, \
        tt, gradient_starttime, gradient_endtime, gradient_time, TIME_SCALE, \
        window, shift =  exp_paramlist[mslev]    
	print ' from:', exp_file, ' using mass spectrogram from', gradient_starttime, 'to', gradient_endtime, 'minutes'
	exp_df = pd.read_pickle(exp_file)
	# transform exp_df_head from str to numeric
	exp_df_head = ['ind', 'mslev', 'bpmz', 'bpint', 'starttime'] 
	for each in exp_df_head:
		exp_df[each] = pd.to_numeric(exp_df[each])
	# drop out of range time
	exp_df = exp_df[exp_df['starttime'] >= gradient_starttime]
	exp_df = exp_df[exp_df['starttime'] < gradient_endtime]

	# combine array and its bp to list of float
	for bp, ar, combine in zip(['bpmz', 'bpint'], ['mzarray', 'intarray'], ['allmz', 'allint']):
		exp_df[combine] = exp_df[bp].apply(lambda x: [x]) + exp_df[ar]

	## Create index
	exp_df['starttime'] = time_index(exp_df['starttime'], gradient_starttime, tt)
	exp_df['allmz'] = mz_index(exp_df['allmz'].values, mrange_min, mrange_max, mm)
	exp_df = exp_df[['ind','starttime','allmz','allint']]
	time_col = []
	time_col_temp = []
	for index, row in exp_df.iterrows():
		# remove out of range m
		row['allint'] = [i for m, i in zip(row['allmz'], row['allint']) if m >= 0 and m < MZ_SCALE]
		row['allmz'] = [m for m in row['allmz'] if m >= 0 and m < MZ_SCALE]
		# use bincount to sum int at same mz_index to create time_index col with MZ_SCALE length
		timecol_array = np.bincount(row['allmz'],row['allint'], minlength=(MZ_SCALE))
		timecol_array[timecol_array < 1] = 0
		time_col_temp.append(timecol_array) # append each row, int sum
		if index % 500 == 0:
			time_col.extend(time_col_temp)
			time_col_temp = []
	# flush last
	time_col.extend(time_col_temp)
	exp_df['allint_overlap'] = time_col

	expdf_row = np.tile(np.arange(MZ_SCALE), exp_df.shape[0])
	expdf_col = np.repeat(exp_df['starttime'].values, MZ_SCALE)
	expdf_value = np.concatenate(exp_df['allint_overlap'].values)

	exp_mat = sparse.coo_matrix((expdf_value,\
								(expdf_row, expdf_col)), \
								shape=(MZ_SCALE, TIME_SCALE))

	exp_mat = smoothingtime_mat(exp_mat, window, shift)
	exp_mat, mat_mean = rescale_mat(exp_mat)
	return exp_mat, exp_paramlist, mat_mean

