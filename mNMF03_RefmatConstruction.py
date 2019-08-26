from __future__ import division
import numpy as np
import pandas as pd
import timeit
import collections
from sklearn.preprocessing import normalize
from scipy import sparse

### Global functions ###
from mNMF01_Tools import time_index, mz_index

def refMS1_construction(refms1_df, M_header, iso_header, charge_list, iso_maxnumber, globalparam_list, eps):
	mslev = 1
	mm = globalparam_list[mslev][globalparam_list[0].index('mm')]
	mrange_min = globalparam_list[mslev][globalparam_list[0].index('mrange_min')]
	mrange_max = globalparam_list[mslev][globalparam_list[0].index('mrange_max')]
	MZ_SCALE = globalparam_list[mslev][globalparam_list[0].index('MZ_SCALE')]

	if np.all(refms1_df.prot.values != refms1_df.sort_values('prot').prot.values):
		print 'Warning: Prot is not alphabetically sorted'
		exit()
	# calc number of charged peptide ion per each protein
	prot_peptcount = refms1_df.reset_index().melt(['prot'], M_header).dropna()
	prot_peptcount = prot_peptcount['prot'].value_counts().sort_index().to_frame(name='ms1count')
	all_peptcount = prot_peptcount.values.sum()
	
	mziso_df = refms1_df.rename_axis('pept_id')
	# melt/pivot Mheader (all charge) into 'variable' col and its mz value into 'value' col
	## so every line is a singly charged peptide with this iso head abundance .melt([dfkeep], pivotthing)
	mziso_df = mziso_df.reset_index().melt(['prot', 'pept_id','pept','mod']+iso_header, M_header)
	mziso_df = mziso_df.rename(columns = {'value':'mz', 'variable':'charge'})
	mziso_df['charge'] = mziso_df['charge'].str[-1].astype(int)
	# drop NA, sort same pept_id (same prot) up from small charge first
	mziso_df = mziso_df.dropna().sort_values(['pept_id','charge'])
	# reset index after correct sort to use as col
	mziso_df = mziso_df.rename_axis('tempidx').reset_index().drop('tempidx',1)
	prot_peptcount = mziso_df['prot'].value_counts().sort_index().to_frame(name='ms1count')
	peptcount = prot_peptcount.values.sum()

	mzidx_header = []
	for idx, iso_head in enumerate(iso_header):
		mziso_df['mz_'+iso_head] = mziso_df['mz'].values+(idx/mziso_df['charge'].values)
		mziso_df['mzidx_'+iso_head] = mz_index(mziso_df['mz_'+iso_head].values, mrange_min, mrange_max, mm)
		mzidx_header.append('mzidx_'+iso_head)

	mziso_df['final_mzidx'] = mziso_df[mzidx_header].values.tolist()
	mziso_df['final_normisoab'] = mziso_df[iso_header].values.tolist()

	if iso_maxnumber > 1:
		mziso_df_col = np.repeat(mziso_df.index.values, iso_maxnumber)
		mziso_df_row = np.concatenate(mziso_df['final_mzidx'].values)
		mziso_df_value = np.concatenate(mziso_df['final_normisoab'].values)
	else: # if no isotope
		mziso_df_col = np.array(mziso_df.index.values)
		mziso_df_row = np.array(mziso_df['mzidx_isoab0'].values) 
		mziso_df_value = np.array(mziso_df['isoab0'].values) 

	if mziso_df.index.values.tolist() != range(len(mziso_df.index.values)):
		print 'Warning mziso_df reindex is wrong'
		exit()

	nonzero_idx = np.multiply([mziso_df_row>=0], [mziso_df_value>0]) # keep True True
	nonzero_idx = tuple(nonzero_idx)
	sreference = sparse.coo_matrix((mziso_df_value[nonzero_idx],\
									(mziso_df_row[nonzero_idx], mziso_df_col[nonzero_idx])), \
									shape=(MZ_SCALE, peptcount))
	
	sreference = normalize(sreference, norm='l1', axis=0)
	return sreference, mziso_df, prot_peptcount