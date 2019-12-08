from __future__ import division
import numpy as np
import pandas as pd
from collections import Counter
import math
import itertools

### MAIN ###
def insilico_digest(file, dictionary_list, param_list, globalparam_list, iso_maxnumber):
    print ' from:', str(file)
    df = df_prep(file)
    aa_monomass, aa_composition, aa_fixmodmass, aa_varmodmass, other_mass, el_abundance = dictionary_list
    misclvge_from, misclvge_to, len_from, len_to, charge_from, charge_to = param_list
    protseq_dict = trypsin_cut(df, misclvge_from, misclvge_to, len_from, len_to)
    insilico_df = mz_calc(protseq_dict, aa_monomass, aa_fixmodmass, other_mass, globalparam_list,\
                            iso_maxnumber, charge_from, charge_to)
    insilico_df = comp_calc(insilico_df, aa_composition)
    insilico_df = mzvarmod_calc(insilico_df, aa_varmodmass, globalparam_list, charge_from, charge_to, iso_maxnumber)
    insilico_df = iso_calc(insilico_df, el_abundance, iso_maxnumber)
    mz_header, iso_header, charge_list, eachpeptcount_list = label_df(insilico_df, charge_from, charge_to)
    insilico_number = insilico_df[mz_header].count().sum()
    return insilico_df, mz_header, iso_header, charge_list

### FUNCTIONS ###
def df_prep(ref_file):
	if ref_file[-6:] == '.fasta':
		protein_info, protein_sequence = fasta_read(ref_file)
		df = pd.DataFrame(list(zip(protein_info, protein_sequence)), columns =['Protein', 'Sequence'])
		df = df.dropna(axis=0, how='any')
	else:
		df = pd.read_excel(ref_file)
		df = df.dropna(axis=0, how='any')
	return df

def fasta_read(fasta_file):
	sequence_templist = []
	protein_info, protein_sequence = [], []
	with open(fasta_file) as file:
	    for line in file:
	        line = line.strip()
	        # if nothing, new loop
	        if not line:
	           continue
	        if line.startswith(">"):
	            sequence_info = line[1:]
	            if sequence_info not in sequence_templist:
	            	# join previous protein sequence and clear sequence list
	                protein_sequence.append(''.join(sequence_templist))
	                protein_info.append(sequence_info)
	                sequence_templist = []
	            # new loop
	            continue
	        sequence_templist.append(line)
		# flush last sequence
	    if sequence_templist:
	        protein_sequence.append(''.join(sequence_templist))
	    # delete first blank
	    if protein_sequence[0] == '':
	    	protein_sequence = protein_sequence[1:]
	return protein_info, protein_sequence

def trypsin_cut(protseq_df, misclvge_from, misclvge_to, len_from, len_to):
	misclvge_to += 1
	len_to += 1
	peptides_report, peptides = {}, []
	for protseq in protseq_df.itertuples():
		prot_in, seq_in = protseq[1], protseq[2]		
		# define tryptic site
		cut_sites=[]
		for i in range(0,len(seq_in)-1):
		    if seq_in[i]=='K' and seq_in[i+1]!='P':
		        cut_sites.append(i+1)
		    elif seq_in[i]=='R' and seq_in[i+1]!='P':
		        cut_sites.append(i+1)
		if cut_sites != []:
			if (cut_sites[-1]!=len(seq_in)) and (cut_sites!=[]):
			    cut_sites.append(len(seq_in))
			if (cut_sites[0]!=0) and (cut_sites!=[]):
				cut_sites.insert(0,0)
	    # digest with misclvge
		if len(cut_sites)>1:
			for m in range(misclvge_from, misclvge_to):
				prot_m_key = prot_in + '##' + str(m)
				pept_value = []
				for j in range(0,len(cut_sites)-1-m):
					if len_from <= cut_sites[j+1+m] - cut_sites[j] <= len_to:
						pept_value.append(seq_in[cut_sites[j]:cut_sites[j+1+m]])
				# dict {key:unique digested peptide}
				peptides_report.update({prot_m_key:list(set(pept_value))})
		else: # no tryptic site in the protein sequence
		    peptides_report.update({prot_key+ '##0':[seq_in]})

	return peptides_report

def mz_calc(peptides_report, aa_monomass, aa_fixmodmass, other_mass, globalparam_list, iso_maxnumber, charge_from, charge_to):
	mrange_min = globalparam_list[1][globalparam_list[0].index('mrange_min')]
	mrange_max = globalparam_list[1][globalparam_list[0].index('mrange_max')]
	mrange_maxbound = int((mrange_max - (iso_maxnumber/charge_from))*charge_to)
	charge_to += 1
	# replace by fix mod
	if len(aa_fixmodmass) > 0:
		for f_k, f_v in aa_fixmodmass.iteritems():
			aa_monomass[f_k] = f_v
	# calc mass
	M_df = pd.DataFrame()
	for p_k, p_v in peptides_report.iteritems():
		for p_v_each in p_v:			
			prot, mcv = p_k.split('##',1)[0], p_k.split('##',1)[1]
			calc_m = sum(aa_monomass[aa] for aa in p_v_each) + other_mass["H20"]
			if mrange_min <= calc_m <= mrange_maxbound:
				mass_data = pd.DataFrame({'prot' : pd.Series(prot),
							'misclvge' : pd.Series(int(mcv)),
							'pept' : pd.Series(p_v_each),
							'[M]': pd.Series(calc_m)
							})
				M_df = M_df.append(mass_data, ignore_index=True)

	M_df = M_df[['prot', 'misclvge', 'pept', '[M]']] #***
	for charge in range(charge_from, charge_to):
		head_mz = 'Mcharge'+str(charge)
		calc_mz = (M_df[['[M]']] + (charge*other_mass["Proton"]))/charge
		M_df[head_mz] = calc_mz
		# remove every charge of M_df [each row, column:head_mz] = np.NaN
		M_df.loc[(M_df[head_mz] <= mrange_min, head_mz)] = np.NaN
		M_df.loc[M_df[head_mz] >= mrange_max-(iso_maxnumber/charge), head_mz] = np.NaN
		
	M_df = M_df.groupby('prot').apply(pd.DataFrame.sort_values, '[M]').reset_index(drop=True)
	charge = charge_to - charge_from
	# clean blank mz
	mz_header = []
	for charge in range(charge_from, charge_to):
		mz_header.append('Mcharge'+str(charge))
	M_df = M_df.dropna(subset = [mz_header], how='all')
	M_df = M_df.reset_index(drop=True)
	return M_df

def comp_calc(M_df, aa_composition):
	# apply each row function and set up a new column
	M_df['comp'] = M_df.apply(lambda row: sum((Counter(aa_composition[aa]) for aa in row['pept']), Counter()), axis=1)
	return M_df

def mzvarmod_calc(M_df, aa_varmodmass, globalparam_list, charge_from, charge_to, iso_maxnumber):
	M_df['mod'] = ""
	mrange_min = globalparam_list[1][globalparam_list[0].index('mrange_min')]
	mrange_max = globalparam_list[1][globalparam_list[0].index('mrange_max')]

	if len(aa_varmodmass) > 0: # if more than 1 mod
		for mod_info, v_tup in aa_varmodmass.iteritems():
			mass_change = v_tup[0] # mass change of mod
			comp_change_dict = v_tup[2] # comp change of mod
			mod_amino = v_tup[3]

			if v_tup[1] == 'startswith': # position of mod
				# define new df for this var mod
				modM_df = M_df[M_df['pept'].str.startswith(mod_amino)].copy()
				modM_df['mod'] += mod_info
				mod_pept = []
				for pept in modM_df['pept']:
					mod_pept.append('x'+str(pept))
				modM_df['pept'] = mod_pept

			elif v_tup[1] == 'all':
				mod_number = v_tup[4]	
				modM_df = M_df[M_df['pept'].str.contains(mod_amino)].copy()
				if mod_number == 1:
					pass
				elif mod_number == 2:
					mod_count_inpept = modM_df['pept'].str.count(mod_amino)					
					modM_df = modM_df[mod_count_inpept == mod_number].copy() # replace
					modM_df = modM_df[[str(1)+str(mod_info[1:]) not in each for each in modM_df['mod']]]
				else:
					print ' Warning mod_number is more than 2'
					exit()
				mod_pept, mod_index_inpept = [], []
				for pept in modM_df['pept']:
					each_index_inpept = [i+1 for i, aa in enumerate(pept) if aa == mod_amino]
					mod_index_inpept.append(each_index_inpept)
					if mod_number == 1:
						i = each_index_inpept[0]
						pept = str(pept[:i-1])+'x'+str(pept[i-1:])
					elif mod_number == 2:
						each_index_inpept.reverse() #instant mod from the end
						for i in each_index_inpept:
							pept = str(pept[:i-1])+'x'+str(pept[i-1:])						
						each_index_inpept.reverse()
					mod_pept.append(pept)
				modM_df['mod'] += [mod_info + str(each) for each in mod_index_inpept]
				modM_df['pept'] = mod_pept

			# modM_df['mod'] = mod_info
			modM_df['[M]'] += mass_change	
			for Mcharge in np.arange(charge_from, charge_to+1):
				modM_df['Mcharge'+str(Mcharge)] += mass_change/Mcharge
			# update composition by sum value of same key via Counter
			comp_update = []
			for idx, compdict_each in enumerate(modM_df['comp']):
				comp_update.append(dict(Counter(compdict_each) + Counter(comp_change_dict)))
			modM_df['comp'] = comp_update
			
			# insert mod under unmod by creating reindex list
			newindex_list = range(modM_df.shape[0]+M_df.shape[0])
			M_df_indlist, modM_df_indlist = [], []
			for idx in M_df.index.values:
				M_df_indlist.append(newindex_list.pop(0))
				if idx in modM_df.index.values:
					modM_df_indlist.append(newindex_list.pop(0))
				
			M_df.index = M_df_indlist
			modM_df.index = modM_df_indlist
			M_df = pd.concat([M_df, modM_df]).sort_index()
	
	for charge in range(charge_from, charge_to+1):
		head_mz = 'Mcharge'+str(charge)
		M_df.loc[(M_df[head_mz] <= mrange_min, head_mz)] = np.NaN
		M_df.loc[M_df[head_mz] >= mrange_max-(iso_maxnumber/charge), head_mz] = np.NaN

	return M_df

def iso_calc(M_df, el_abundance, iso_maxnumber):
	def dist_coef(el, sub_num):
		if comp_each[el] >= sub_num: #combinatorics of way to sub of each el
			coef = math.factorial(comp_each[el])/(math.factorial(comp_each[el]-sub_num)*math.factorial(sub_num))
			return coef			
		else:
			return 0
	def dist_calc(possible_el, space_tosub, iso_of_el_tosub):
		abundance = 0
		for el in possible_el: #combinatorics * new iso abundance / core iso abundance
			abundance += dist_coef(el,space_tosub)*(el_abundance[el].values()[iso_of_el_tosub]**space_tosub/el_abundance[el].values()[0]**space_tosub)
		return abundance
	def dist_ofisopeak(space_tosub):
		#define el and space able to sub
		el_space = [('S',4),('S',2),('S',1),('O',2),('O',1),('H',1),('C',1),('N',1),('x',0)]
		#combinatorics of el upto number of space_tosub. 0 is necessary
		combi_all = [p for p in itertools.combinations_with_replacement(el_space, space_tosub)]
		combi_dist_calc = 0
		for combi_each in combi_all:
			if np.sum([x[1] for x in combi_each]) == space_tosub: #iso of el to sub must == space_tosub
				combi_each = [i for i in combi_each if i[1] > 0] #delete (x,0)
				combi_each = Counter(combi_each) #find how many item per el
				component_dist_calc = 1
				for k, v in combi_each.iteritems():
					component_dist_calc *= dist_calc([k[0]],v,k[1]) #multiply each component
				combi_dist_calc += component_dist_calc
		return combi_dist_calc
	
	for idx, comp_each in M_df.loc[:,['comp']].itertuples():
		# add water
		comp_each['H'] += 2
		comp_each['O'] += 1
		# calc core as product(abundance power number in comp) of each comp
		core_abundance = 1
		for comp_k, comp_v in comp_each.iteritems():
			core_abundance *= el_abundance[comp_k].values()[0] ** comp_v	
		# locate core
		M_df.loc[idx,'isoab0'] = core_abundance
		# calc and locate iso
		for isopeak in range(iso_maxnumber-1):
			isopeak += 1 #start from 1, not calc core_abundance
			M_df.loc[idx,'isoab'+str(isopeak)] = dist_ofisopeak(isopeak) * core_abundance

	iso_header = [head for head in M_df.columns if 'isoab' in head]
	M_df['iso_sum'] = M_df[iso_header].sum(axis=1)
	M_df = M_df.drop('comp', axis=1).reset_index(drop=True)
	return M_df

def label_df(M_df, charge_from, charge_to):
	charge_to += 1
	mz_header, charge_list, eachpeptcount_list = [], [], []
	for charge in range(charge_from, charge_to):
		mz_header.append('Mcharge'+str(charge))
		charge_list.append(charge)
	
	for index, row in M_df[mz_header].iterrows():
		eachpeptcount_list.append(row.count())

	iso_header = [head for head in M_df.columns if 'isoab' in head]
	return mz_header, iso_header, charge_list, eachpeptcount_list
