from __future__ import division
import numpy as np
import pandas as pd
from scipy import sparse
from scipy import special

### FUNCTIONS ###
from mNMF05_ExpmatConstruction import smoothingtime_mat
from mNMF10_NMFReport import nmf_identification, nmf_savemat

def tol_test(it, last_it, max_iter, tracking_it_all):
	tol_value = 1e-4		
	if it > 3: # error tolerance (run at least 4 + 1 it)
		tol = (tracking_it_all[-2] - tracking_it_all[-1])/tracking_it_all[-2]
		if tol <= tol_value:
			last_it = True #go to last it
			# print ' report convergence at it: ', it, tol, last_it
		# else:
			# print ' report tol: ', tol
	return it, last_it

def kl_calc(Vmat, Wmat, Hmat, Vmat_mean, peptcount_eachsqrt, eachprot_sum_it, weight, eps):
	if weight != 0:
		penaltyH_it = np.sum(np.multiply(peptcount_eachsqrt,np.log(np.asarray(eachprot_sum_it) + eps))) #log L1 group
		penaltyH_it = weight*penaltyH_it
		# print ' chk main penaltyH_it:', penaltyH_it
	else:
		penaltyH_it = 0

	V_approx = Wmat.dot(Hmat) + eps
	cost_it = special.kl_div(Vmat.todense(),V_approx)
	cost_it = np.sum(cost_it) + penaltyH_it
	del(V_approx)	
	return cost_it

def kl_nmf(Vmat, Wmat, Hmat, cos_mat, prot_peptcount,\
        weight, max_iter, eps, noise_number, noise_mean, Vmat_mean, globalparam_list, iso_maxnumber,\
        initRT, outfile, finalMS1_df, peak_report, WithStandardize_):
	print ' process NMF'
	for param in globalparam_list:
		mm, mrange_min, mrange_max, mz_range, MZ_SCALE, \
		tt, gradient_starttime, gradient_endtime, gradient_time, TIME_SCALE, \
		window, shift =  globalparam_list[1]

	peptcount_each = prot_peptcount.values
	if WithStandardize_:
		peptcount_eachsqrt = peptcount_each**0.5 # for inside standardize
	else:
		peptcount_eachsqrt = 1
	peptcount_cumsum = np.insert(np.cumsum(peptcount_each), 0, 0)
	protcount, peptcount = len(peptcount_each), np.sum(peptcount_each)

	it_all, tracking_it_all = [], []
	eachprot_timeprofile = []
	eachprot_sum, eachprot_mean, eachprot_std = [], [], []

	Wmat = Wmat.tocsr()
	Vmat = sparse.csr_matrix(Vmat)

	if sparse.issparse(Hmat):
		refresh_H = Hmat.copy().todense()
	else:
		refresh_H = Hmat.copy()

	last_it = False
	# dimension arrangement lkij
	for it in xrange(max_iter):	
		it_all.append(it)
		numer = Wmat.dot(refresh_H) + eps
		numer = Wmat.transpose().dot(Vmat/numer)
		refresh_H_copy = refresh_H.copy()
		refresh_H = np.multiply(refresh_H, numer)
		del(numer)
		eachprot_sum_it, eachprot_sumax_it = [], []
		for each in range(protcount):
			from_k = peptcount_cumsum[each]
			to_k = peptcount_cumsum[each+1]
			refresh_Hg = refresh_H[from_k:to_k,slice(None)].copy()
			if weight != 0:
				deno = 1 + (weight*peptcount_eachsqrt[each]/(np.sum(abs(refresh_Hg))+eps)) # log l1 group 
				refresh_H[from_k:to_k,:] = refresh_Hg/deno
				del(refresh_Hg)
				del(deno)
				eachprot_sum_it.append(np.sum(refresh_H[from_k:to_k,:]))
				eachprot_sumax_it.append(np.sum(refresh_H[from_k:to_k,:],axis=1))
				# print ' chk refresh_Hg: ', eachprot_sum_it[-1]
			# every prot at last it					
			if last_it:
				eachprottarget_ratio, eachpepttarget_ratio = [], []
				refresh_Hg_process = (refresh_H[from_k:to_k,slice(None)]*Vmat_mean)						
				eachprot_sum.append(np.sum(refresh_Hg_process))
				eachprot_mean.append(np.mean(refresh_Hg_process))
				eachprot_timeprofile.append(np.ravel(np.sum(refresh_Hg_process,axis=0)))
				if each == protcount-1:
					cost_it = kl_calc(Vmat, Wmat, refresh_H, Vmat_mean, peptcount_eachsqrt, eachprot_sum_it, weight, eps) # calc last before clean
					tracking_it_all.append(cost_it)
					eachprot_timeprofile.append(np.ravel(np.sum(refresh_H[-noise_number:,slice(None)],axis=0))*Vmat_mean) #append noise
					it, last_it = tol_test(it, last_it, max_iter, tracking_it_all)
					# print ' report monotonic tracking:', tracking_it_all, all(x>y for x, y in zip(tracking_it_all, tracking_it_all[1:]))
					
		cost_it = kl_calc(Vmat, Wmat, refresh_H, Vmat_mean, peptcount_eachsqrt, eachprot_sum_it, weight, eps)
		tracking_it_all.append(cost_it)
		if last_it: #reach tolerance
			tracking_it_all = tracking_it_all[:-1]
			break
		it, last_it = tol_test(it, last_it, max_iter, tracking_it_all)

	print ' process report'
	it_all = pd.DataFrame(np.column_stack([np.arange(len(tracking_it_all)),tracking_it_all]), columns=['it', 'tracking'])
	nmf_identification(outfile, Vmat*Vmat_mean, Wmat, sparse.csr_matrix(refresh_H*Vmat_mean), finalMS1_df,\
						globalparam_list, prot_peptcount, [Wmat.tocsc()[:,-noise_number:].todense()], noise_mean, weight,\
						it_all, initRT, eachprot_sum, eachprot_mean, eachprot_timeprofile,\
						iso_maxnumber, cos_mat, peak_report)
	nmf_savemat(Vmat*Vmat_mean, Wmat, sparse.csr_matrix(refresh_H*Vmat_mean), weight, outfile)

		




	