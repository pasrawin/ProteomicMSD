from __future__ import division
import numpy as np
import pandas as pd
from scipy import sparse
from pyteomics import mzml

def text_reader(msconvert_file):
	## pattern
	msbasic_match = [
				'index: ',
				'cvParam: ms level, ',
				'cvParam: base peak m/z, ',
				'cvParam: base peak intensity, ',
				'cvParam: scan start time, ',			
				]
	msarray_match = [
				'cvParam: m/z array, m/z',
				'cvParam: intensity'
				]
	### end pattern
	ms1_head = msbasic_match + msarray_match
	df_ms1 = pd.DataFrame(columns=ms1_head)

	with open(msconvert_file) as infile:
		count, list_sub = 0, []
		for line in infile:			
			line = line.strip()
			for match in (msbasic_match):
				if match in line:
					list_sub.append(line)
			for matcharr in msarray_match:
				if matcharr in line:			
					list_sub.append(next(infile, '').strip())
					count +=1
				# after complete one index set, fill in df according to ms-level
				if count == 2:	
					if list_sub[1] == 'cvParam: ms level, 1':
						df = pd.DataFrame([list_sub,],columns=(ms1_head))
						df_ms1 = df_ms1.append(df)		
					count, list_sub = 0, []
				elif count > 2:
					print count
					print "Warning: count is not 2"
					exit()

		df_ms1 = df_ms1[ms1_head].reset_index(drop=True)

	toparray = '^binary: \[[0-9]+\] ' # binary: [any numbers]
	# Delete unnecessary words
	other_delete = [', m/z', ', number of detector counts', ', minute' ]
	for head_delete in msbasic_match + other_delete:
		df_ms1.replace(head_delete, '', inplace=True, regex=True)
	for head_delete in msarray_match:
		df_ms1.replace(toparray, '', inplace=True, regex=True)

	col_set = {'index: ':'ind',
				'cvParam: ms level, ':'mslev', 
				'cvParam: base peak m/z, ':'bpmz',
				'cvParam: base peak intensity, ':'bpint',
				'cvParam: scan start time, ':'starttime',
				'cvParam: m/z array, m/z':'mzarray',
				'cvParam: intensity':'intarray'
				}

	df_ms1.rename(columns=col_set, inplace=True)
	for ar in ['mzarray', 'intarray']:
		df_ms1[ar] = df_ms1[ar].apply(lambda x: np.asarray(x.strip().split(" ")).astype(float).tolist())
	return df_ms1

def mzml_reader(msconvert_file):
	ind,mslev,bpmz,bpint,starttime,mzarray,intarray = [],[],[],[],[],[],[]
	with mzml.read(msconvert_file) as reader:
		k_count = 0
		for each_dict in reader:
			# print each_dict
			if each_dict['ms level'] == 1:
				ind.append(each_dict['index'])
				bpmz.append(each_dict['base peak m/z'])
				bpint.append(each_dict['base peak intensity'])			
				mzarray.append(each_dict['m/z array'])
				intarray.append(each_dict['intensity array'])
				v_dict = each_dict['scanList'].values()[2][0]
				starttime.append(v_dict['scan start time'])
	
	mslev = [1]*len(ind)
	mzarray = [x.tolist() for x in mzarray]
	intarray = [x.tolist() for x in intarray]
	col_set = ['ind','mslev','bpmz','bpint','starttime','mzarray','intarray']
	df_ms1 = pd.DataFrame(list(zip(ind,mslev,bpmz,bpint,starttime,mzarray,intarray)), columns=col_set)
	return df_ms1

def msconvert_reader(msconvert_file):
	print ' from:', str(msconvert_file)
	if msconvert_file[-5:] == '.mzML':
		df_ms1 = mzml_reader(msconvert_file)
		dub = 5
	elif msconvert_file[-4:] == '.txt':
		df_ms1 = text_reader(msconvert_file)
		dub = 4
	store = df_ms1.to_pickle(msconvert_file[:-dub]+'.pkl')
	return str(msconvert_file[:-dub]+'.pkl')
	
def mz_index(data_array, mrange_min, mrange_max, mm):
	# print ' process mz_index'
	data_array = np.asarray([np.asarray(x) for x in data_array])
	mz_index_list = [np.around((mz - mrange_min)*mm, 0).astype(int) for mz in data_array] ## replace below loop
	return mz_index_list

def time_index(data_array, gradient_starttime, tt):
	# print ' process time_index'
	data_array = np.asarray([np.asarray(x) for x in data_array])
	time_index_list = [np.around((time - gradient_starttime)*tt,0).astype(int) for time in data_array]
	return time_index_list

def mz_scaling(mz_range, mm):
    return (mz_range*mm)+1

def time_scaling(gradient_time, tt):
    return (gradient_time*tt)+1

def smoothingtime_mat(coomat, window, shift):
	# print ' process smoothingtime_mat', timeit.default_timer()
	def smoothingtime_axis(length,window,shift):
		saxis = int((length/shift)-1)
		return saxis 
	coomat = coomat.tocsc()
	coomat_col = smoothingtime_axis(coomat.shape[1],window,shift)
	coomat_row = coomat.shape[0]
	n, keeprow, keepcol, keepdata = 0, [], [], []
	for col in range(coomat_col):
		eachnewcol = coomat[slice(None),n:n+window]
		eachnewcol_deno = (eachnewcol != 0).sum(axis=1)
		eachnewcol_deno[eachnewcol_deno == 0] = 1
		eachnewcol = eachnewcol.sum(axis=1)/eachnewcol_deno
		eachnewcol_nonzero = np.nonzero(eachnewcol)[0]
		keeprow.append(eachnewcol_nonzero)
		keepcol.append([col]*len(eachnewcol_nonzero))
		keepdata.append(np.ravel(eachnewcol[eachnewcol_nonzero]))
		n += shift

	keeprow = [item for sublist in keeprow for item in sublist]
	keepcol = [item for sublist in keepcol for item in sublist]
	keepdata = [item for sublist in keepdata for item in sublist]
	scoomat = sparse.coo_matrix((keepdata,(keeprow,keepcol)),shape=(coomat_row,coomat_col))
	return scoomat

def rescale_mat(mat):
	mat = sparse.csr_matrix(mat)
	### mean no zero
	mat_nozero = mat[mat>0]
	mat_mean = np.mean(mat_nozero)
	# normalizing mat by mean without zero
	mat = mat/np.mean(mat_nozero)
	if np.around(np.mean(mat[mat>0]),0) != 1:
		print 'Warning after rescale mean (must be 1) but: ', np.around(np.mean(mat[mat>0]),0)
		exit()	
	coomat = sparse.coo_matrix(mat)
	return coomat, mat_mean
