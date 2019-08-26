from __future__ import division
import numpy as np
import pandas as pd
from pyteomics import achrom

def rt_prediction(calib_file, finalMS1_df):	
	print ' from:', calib_file
	calib_df = pd.read_excel(calib_file)
	reten_seq = [str(x) for x in calib_df['Sequence']]
	reten_peak = calib_df['PepRtimePeak'].tolist()
	np.warnings.filterwarnings('ignore')
	RCs = achrom.get_RCs_vary_lcp(reten_seq, reten_peak)
	initRT_tuple = []
	for idx, row in finalMS1_df.iterrows():
		pept = str(row['pept'])
		rt = achrom.calculate_RT(pept, RCs)
		initRT_tuple.append((idx,pept,rt))
	 
	reten_start = calib_df['PepRtimeStart']
	reten_end = calib_df['PepRtimeEnd']
	initRT_width = np.mean((reten_end - reten_start).values)
	
	reten_validate = []
	for seq in reten_seq:
		rt_v = achrom.calculate_RT(seq, RCs)
		reten_validate.append(rt_v)
	
	predict_coef = np.corrcoef(reten_peak, reten_validate)
	return initRT_tuple, initRT_width



		