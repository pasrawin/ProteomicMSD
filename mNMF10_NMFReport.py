from __future__ import division
import numpy as np
import pandas as pd
from scipy import sparse, stats, spatial, signal

### FUNCTIONS ###
from mNMF01_Tools import smoothingtime_mat

def nmf_savemat(V_mat, W_mat, Hmat, weight, outfile):
	# print ' process sparse nmf_savemat all weights: '
	if V_mat != None:
		sparse.save_npz(str(outfile)+'_Vmat_init.npz', V_mat)
	if W_mat != None:
		sparse.save_npz(str(outfile)+'_Wmat_init.npz', W_mat)
	sparse.save_npz(str(outfile)+'_Hmat_result'+str(weight)+'.npz', Hmat)
	# print ' report nmf_savemat done: ', timeit.default_timer()

def nmf_identification(outfile, V_mat, W_mat, refresh_H, finalMS1_df,\
						globalparam_list, prot_peptcount, noiseW_list, noise_mean, weight,\
						tracking_it, initRT, eachprot_sum, eachprot_mean, eachprot_timeprofile,\
						iso_maxnumber, cos_mat, peak_report):
	
	for param in globalparam_list:
		mm, mrange_min, mrange_max, mz_range, MZ_SCALE, \
		tt, gradient_starttime, gradient_endtime, gradient_time, TIME_SCALE, \
		window, shift =  globalparam_list[1]	
	protcount = prot_peptcount.index.values
	peptcount_cumsum = np.insert(np.cumsum(prot_peptcount.values), 0, 0)

	writer_ = pd.ExcelWriter(outfile+'.xlsx', engine='xlsxwriter')
	# pd.DataFrame({'noisemean': [noise_mean]}).to_excel(writer_, sheet_name='Noisemean')
	timelist_label = []
	for label in protcount:
		timelist_label.append(label)
	timelist_label.append('noise')
	
	df_eachprot_timeprofile = pd.DataFrame(np.asarray(eachprot_timeprofile).T, columns=[timelist_label])
	df_eachprot_timeprofile.to_excel(writer_, sheet_name='XIC')
	
	all_df = report_all(refresh_H, cos_mat, noiseW_list, noise_mean, shift, tt, gradient_starttime, finalMS1_df, peak_report)
	all_df.to_excel(writer_, sheet_name='Peakresult')
	
	writer_.save()
	# print ' report identification done:', timeit.default_timer()

def report_all(refresh_H, cos_mat, noiseW_list, noise_mean, shift, tt, gradient_starttime, finalMS1_df, peak_report):
	cos_mat = cos_mat.todok()
	peakall_list = [[] for n in range(peak_report)]
	peakbestiden_list = []
	eachpeptsum_list = []
	for k in np.arange(refresh_H.shape[0] - noiseW_list[0].shape[1]):
		eachpept_chrom = np.ravel(refresh_H[k,:].todense())
		eachpept_sum = np.sum(eachpept_chrom)
		eachpeptsum_list.append(eachpept_sum)
		argrelmax = signal.argrelmax(eachpept_chrom, axis=0, order=2)[0]
		relmax = eachpept_chrom[argrelmax]
		argrelmax_sort = argrelmax[np.argsort(-relmax)][:peak_report]
		timemax_sort = ((argrelmax_sort+1)*shift/tt) + gradient_starttime		
		relmax_sort = eachpept_chrom[argrelmax_sort]
		cosmax_sort = [cos_mat[k,x] for x in argrelmax_sort]
		for idx, time in enumerate(timemax_sort):
			peakall_list[idx].append((k, time, relmax_sort[idx], cosmax_sort[idx]))
	idenhead_df = finalMS1_df.loc[:,['prot','pept','mod','charge','mz','isoab0']]
	idenhead_df.loc[:,'k'] = np.arange(finalMS1_df.shape[0])
	idenhead_df.loc[:,'sum'] = eachpeptsum_list
	for idx in np.arange(peak_report):
		idex_label = idx+1
		each_df = pd.DataFrame(peakall_list[idx], \
			columns=['k','peaktime'+str(idex_label),'peakintensity'+str(idex_label), 'peakcosscore'+str(idex_label)])
		if idx == 0:
			idenall_df = idenhead_df.join(each_df.set_index('k'), on='k')
		else:
			idenall_df = idenall_df.join(each_df.set_index('k'), on='k')
		idenall_df.loc[:,'peakintensity'+str(idex_label)] \
		= idenall_df.loc[:,'peakintensity'+str(idex_label)]*idenall_df.loc[:,'isoab0']		
	idenall_df = idenall_df.drop(['k','isoab0'], axis=1)

	return idenall_df
