from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from sklearn.decomposition import NMF
from scipy import sparse

from mNMF01_Tools import mz_scaling, time_scaling
from mNMF05_ExpmatConstruction import expmat_construction

def noise_subspace(file, W_mat, noise_no, globalparam_list, charge_list, start, end):
    noiseparam_list = noise_selection(globalparam_list, start, end)
    noise_mat, noiseparam_list, noise_mean\
    = expmat_construction(file, noiseparam_list, charge_list)
    W_mat = noise_scikit(W_mat, noise_mat, noise_no, noiseparam_list)
    return W_mat, noise_mean

def noise_selection(globalparam_list, start, end):
    noiseparam_list = np.asarray(globalparam_list[1]) # preserve globalparam_list
    duration = end - start
    for idx, each in enumerate(globalparam_list[0]):
        if each == 'tt':
            tt_factor = noiseparam_list[idx]
        elif each == 'gradient_starttime':
            noiseparam_list[idx] = start
        elif each == 'gradient_endtime':
            noiseparam_list[idx] = end
        elif each == 'gradient_time':
            noiseparam_list[idx] = duration
        elif each == 'TIME_SCALE':
            noiseparam_list[idx] = time_scaling(duration,tt_factor)

    noiseparam_list = [globalparam_list[0], noiseparam_list.tolist()]
    if globalparam_list[0:2] == noiseparam_list:
        print ' Warning globalparam_list is wrong: ', globalparam_list
        exit()
    return noiseparam_list

def noise_scikit(W_mat, noise_mat, noise_no, noiseparam_list):
    noise_mat = noise_mat.todense()
    # print ' process Noise Extraction with shape: ', noise_mat.shape
    model = NMF(n_components=noise_no, beta_loss='kullback-leibler', solver='mu', init='random', max_iter=250)		
    W = model.fit_transform(noise_mat)
    del(noise_mat)
    W = W/np.sum(W, axis=0) # norm W list
    # H = model.components_
    W = sparse.coo_matrix(W)
    W_mat = sparse.hstack([W_mat, W]) 
    return W_mat