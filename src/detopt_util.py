# detopt_util.py

import os, re, math, copy, random, sys
import pandas as pd
from ast import literal_eval as lit


################
# UTILITY FXNS #
################

#############
# NUMERICAL #
#############

def eq(x, y):
	return abs(x-y) <= 1e-6


def floatjoin(L, sep):
	return sep.join(list(map(str, L)))


#####################
# DATAFRAME QUERIES #
#####################

from pandas_util import *


def read_neoantigens(patient):
	df = pd.read_csv('/Users/wuchh/Desktop/DetOpt/data/neoantigens.tsv', 
					 sep='\t', header=0, index_col=False)
	neos = df[df['Patient'] == patient][['Gene', 'Product', 'Variant_key']]
	neos['Variant_key'] = neos['Variant_key'].map(lambda x: x.replace(' ', '_'))
	return neos


######################
# COPY NUMBER STATES #
######################

def has_aberrant_state(row):
	state_cols = [col for col in row.index if 'state' in col]
	states = set(list(row[state_cols]))
	return states != {'1|1'}


def sample_aberrant_muts(df):
	df = df[df.apply(has_aberrant_state, axis=1)]
	aberrants = df['mut_index'].unique().tolist()
	print(f"There are {len(aberrants)} aberrant muts.")
	return aberrants


##################
# TREE FUNCTIONS #
##################

def tree_df_to_vec(df):
	vec = {}
	for (i, row) in df.iterrows():
		node = row['NODE_ID']
		vec[node] = row['PARENT_ID']
	return vec


def vec_to_anc_matrix(vec):
	A = {}
	root_id = 'ROOT'
	for node_i in vec.keys():
		A[node_i] = {}
		for node_j in vec.keys():
			A[node_i][node_j] = 0
	for node_j in vec.keys():
		node_i = node_j
		while node_i != root_id:
			A[node_i][node_j] = 1
			node_i = vec[node_i]
	for node in vec.keys():
		A[node][node] = 0 # a node is NOT an ancestor of itself for our purposes
	# result: if A[i][j] = 1, node i is ancestor of node j
	return A


def children(node, vec):
	return [child for child in vec.keys() if vec[child] == node]


def descendants(node, A):
	# a node is not a descendant of itself
	return [desc for desc in A[node].keys() if A[node][desc] == 1 and node != desc]


def distance_from_root(node, vec):
	distance = 0
	if node == 'ROOT': return 0
	else:
		while node != 'ROOT':
			distance += 1
			node = vec[node]
		return distance


def difference_from_diploid(state):
	return abs(state[0] - 1) + abs(state[1] - 1)


def states_to_int(states):
	res = []
	for state in states:
		(a,b) = state
		res.append((int(a), int(b)))
	return res


def splitmuts(muts):
	if isinstance(muts, float): return [] # for NaN cases
	elif isinstance(muts, list): return muts # already split
	return muts.split(',')


def find_node(mut, tree, mut_col_name):
	tree2 = tree.copy(deep=True)
	tree2[mut_col_name] = tree2[mut_col_name].map(splitmuts)
	node = 'ROOT'
	for (i, row) in tree2.iterrows():
		if mut in row[mut_col_name]:
			node = row['NODE_ID']
	return node


def values_to_tuple(d):
	res = {}
	for k in d.keys():
		res[k] = tuple(d[k])
	return res


def collapse(L):
	res = []
	j = 0
	res.append(L[0])
	for i in range(1, len(L)):
		if L[i] != res[j]:
			res.append(L[i])
			j += 1
	return res


def asg_to_state_tree(mut, asg, tree, vec, A, B):
	'''
	input:
	- tree: tree data frame
	- asg: assignment of allelic copy number per node
	- vec: tree topology
	- mut: mutation ID to return a state tree for
	- A: maternal copies of mut per node
	- B: paternal copies of mut per node

	output: state tree formatted as chain of connected states
	'''
	mut_node = find_node(mut, tree)
	#print("MUT NODE:", mut_node)
	#print("STATE AT MUT NODE:",  (int(asg[mut_node][0]), int(asg[mut_node][1]), int(A[mut_node]+B[mut_node])))
	state_tree = []
	while mut_node != "ROOT":
		state_tree.insert(0, (int(asg[mut_node][0]), int(asg[mut_node][1]), int(A[mut_node]+B[mut_node])))
		mut_node = vec[mut_node]
	state_tree.insert(0, (1,1,0))

	state_tree = collapse(state_tree)
	return state_tree


def calculate_DCF(mut, asg, tree, vec, A, B, presence, node_freqs, samples, mut_col_name):
	'''
	input:
	- tree: tree data frame
	- asg: assignment of allelic copy number per node
	- vec: tree topology
	- mut: mutation ID to return a state tree
	- A: maternal copies of mut per node
	- B: paternal copies of mut per node

	output: DCF of the mutation
	'''
	dcfs = {}
	mut_node = find_node(mut, tree, mut_col_name)
	for sample in samples:
		dcfs[sample] = node_freqs[mut_node][sample] * presence[mut_node]
		for desc in descendants(mut_node, vec_to_anc_matrix(vec)):
			dcfs[sample] += node_freqs[desc][sample] * presence[desc]

		print(f"DCF of mut {mut} in sample {sample} = {dcfs[sample]}")
	return dcfs


def mean_quality(value, values):
	values2 = copy.deepcopy(list(values))
	for v in values2:
		if eq(v, 1e6): 
			print('removing', v)
			values2.remove(v)
	return sum([abs(value - v) for v in values2]) / len(values2)


def mut_isin(mut_id, row):
	if isinstance(row['CNA_IDS'], float): return False
	return mut_id in row['CNA_IDS'].split(',')


def correct_cna_placement(mutation_id, node_of_event, tree):	
	mutation_id = 'cna_' + mutation_id
	tree = tree[tree.apply(lambda row: mut_isin(mutation_id, row), axis=1)]
	#print(f"node of event:{node_of_event}; actual node: {tree['NODE_ID'].iloc[0]}")
	return (tree['NODE_ID'].iloc[0])


def avg_vaf(node, sample, muts_at_node, vafs):
	res = 0
	for mut in muts_at_node:
		res += float(vafs[mut][sample])
	return res / len(muts_at_node)


def avg_vaf_agg(node, sample, muts_at_node, decifer):
	# add all variant read counts
	# divide by all total read counts
	tot_var = 0
	tot_tot = 0
	for mut in muts_at_node:
		var = decifer[(decifer['mut_index'] == mut) & (decifer['sample'] == sample)].squeeze()['var_reads']
		ref = decifer[(decifer['mut_index'] == mut) & (decifer['sample'] == sample)].squeeze()['ref_reads']
		tot_var += int(var)
		tot_tot += (int(var) + int(ref))
	return tot_var / tot_tot


def get_cn_neutral_at_node(node, tree, decifer, mut_col_name):
	candidates = tree[tree['NODE_ID'] == node].squeeze()[mut_col_name].split(',')
	res = []
	for m in candidates:
		row = decifer[decifer['mut_index'] == m].iloc[0]
		if has_aberrant_state(row):
			continue
		else:
			res.append(m)
	return res


def calc_observed_prevalence(node, tree, decifer, vec, samples, vafs, mut_col_name):
	'''
	vafs : dict mapping mut id to dict mapping sample id to vaf

	output:
		prevalence : dict mapping sample ID to vafs
	'''
	# freq(node) = avg vaf(node) - sum_children(avg vaf(child))

	prevalence = {}
	muts_at_node = get_cn_neutral_at_node(node, tree, decifer, mut_col_name)
	#print("NEUTRAL MUTS: ", node, muts_at_node)

	for sample in samples:
		# get avg vaf over all mutations in node

		vaf_at_node = avg_vaf_agg(node, sample, muts_at_node, decifer)
		prevalence[sample] = vaf_at_node
		for desc in children(node, vec):
			muts_at_desc = get_cn_neutral_at_node(desc, tree, decifer, mut_col_name)
			prevalence[sample] -= avg_vaf_agg(desc, sample, muts_at_desc, decifer)

		prevalence[sample] = max(0., 2 * prevalence[sample])
		#print(f'prev[{sample}] = {prevalence[sample]}')
	return prevalence


def append_root_prevalence(freqs):
	samples = freqs['0'].keys()
	for node in freqs.keys():
		for sample in samples:
			freqs['ROOT'][sample] = 1.

	for node in freqs.keys():
		for sample in samples:
			freqs['ROOT'][sample] -= freqs[node][sample]
	return freqs


def freqs_from_true_tree(node, true_tree, nsamples):
	res = {}
	samples = true_tree.iloc[0]['SAMPLE_IDS'].split(',')
	# if topologies are different and this node has no match in groundtruth, mark them as 0.
	freqs = true_tree[true_tree['NODE_ID'] == node].squeeze()['SAMPLE_NODE_FREQUENCIES']
	if len(freqs) == 0: 
		freqs = [0.] * nsamples
	else:
		freqs = freqs.split(',')
	for i in range(nsamples):
		res[samples[i]] = float(freqs[i])
	return res


def parse_vafs(df, samples_chosen):
	'''
	parse vafs only from decifer df
	'''
	vafs = {}
	df = df[df['sample'] != 'MASTER_SAMPLE']
	df = df[df['sample'].isin(samples_chosen)]

	for mut_id in set(df['mut_index'].values):
		vafs[mut_id] = {}

	for (i, row) in df.iterrows():
		mut_id = row['mut_index']
		vaf = row['var_reads'] / (row['ref_reads'] + row['var_reads'])
		sample = row['sample']
		if sample in samples_chosen:
			vafs[mut_id][sample] = vaf
	return vafs


def parse_gt_tree(df, num_samples):
	'''
	parse groundtruth tree
	'''
	samples = df.iloc[0]['SAMPLE_IDS'].split(',')
	if 'MASTER_SAMPLE' in samples: samples.remove('MASTER_SAMPLE')	
	samples_chosen = samples[:num_samples]
	vec, freqs = {}, {}
	for (i, row) in df.iterrows():
		node = row['NODE_ID']
		freqs[node] = {}
		vec[node] = row['PARENT_ID']
		nodefreqs = row['SAMPLE_NODE_FREQUENCIES'].split(',')
		samples = row['SAMPLE_IDS'].split(',')
		for i in range(len(samples)):
			sample = samples[i]
			if sample in samples_chosen:
					freqs[node][sample] = float(nodefreqs[i])
	return vec, freqs, samples_chosen


def parse_pt_tree(df, num_samples, true_tree):
	'''
	parse pairtree tree, providing ground truth node frequencies & clusters
	'''
	samples_chosen = ['S' + str(i) for i in range(num_samples)]
	vec, freqs = {}, {'ROOT':{}}
	for (i, row) in df.iterrows():
		node = row['NODE_ID']
		vec[node] = row['PARENT_ID']

		samples = row['SAMPLE_IDS'].split(',')
		'''
		for i in range(len(samples)):
			sample = samples[i]
			if sample in samples_chosen:
					freqs[node][sample] = float(nodefreqs[i])
		'''
	for node in vec.keys():
		if node == 'ROOT': continue
		freqs[node] = freqs_from_true_tree(node, true_tree, num_samples)
	freqs = append_root_prevalence(freqs)

	return vec, freqs, samples_chosen

# NEW
def split_cell(x):
	try:
		return (x.split(","))
	except:
		return ([])

# NEW
def split_cell_tup(x):
	try:
		string_with_tuples = x.replace(" ", "")
		return lit("[" + string_with_tuples + "]")
	except:
		return ([])

def expl(true_tree):
	'''
	explode a ground truth tree by all possible columns
	'''
	#print("HERE")
	#print(true_tree['SNV_IDS'])
	true_tree['SNV_IDS'] = true_tree['SNV_IDS'].map(split_cell)
	true_tree['SNV_ALLELES'] = true_tree['SNV_ALLELES'].map(split_cell)

	res = true_tree.explode(['SNV_IDS', 'SNV_ALLELES'])
	res['CNA_IDS'] = res['CNA_IDS'].map(split_cell)
	res['CNA_CHANGES'] = res['CNA_CHANGES'].map(split_cell_tup)	
	res = res.explode(['CNA_IDS', 'CNA_CHANGES'])

	res['SAMPLE_IDS'] = res['SAMPLE_IDS'].map(split_cell)
	res['SAMPLE_NODE_FREQUENCIES'] = res['SAMPLE_NODE_FREQUENCIES'].map(split_cell)
	res = res.explode(['SAMPLE_IDS', 'SAMPLE_NODE_FREQUENCIES'])
	return res


def calc_mut_copies(exp, mut, node, node_of_mut, node_of_cna, anc):
	# to have any mutational copies, the current node must be the mutation's
	# node, or a descendant thereof
	if not (anc[node_of_mut][node] or node == node_of_mut):
		return 0
	elif not (anc[node_of_cna][node] or node == node_of_cna):
		return 1
	else:
		allele = exp[exp['SNV_IDS'] == mut]['SNV_ALLELES'].iloc[0]
		cna_change = exp[exp['CNA_IDS'] == 'cna_'+mut]['CNA_CHANGES'].iloc[0]
		if allele == 'M':
			cna = int(cna_change[0])
		else:
			assert(allele=='P')
			cna = int(cna_change[1])
		# we always start from 1 copy and have at most one cna change
		return 1 + cna


def calc_seg_copies(exp, mut, node, node_of_mut, node_of_cna, anc):
	# calculate #copies of genomic segments harboring a mut across the tree
	if not (anc[node_of_cna][node] or node == node_of_cna):
		# if current node is not impacted by cna, it's diploiad
		return 2
	cna_change = exp[exp['CNA_IDS'] == 'cna_'+mut]['CNA_CHANGES'].iloc[0]
	# we always start from 1 copy and have at most one cna change
	return 2 + int(cna_change[0]) + int(cna_change[1])


def parse_frac_cns(true_tree, exp, samples_chosen, mode='segmental'):
	'''
	parse fractional copy numbers only from ground truth tree
	parse them only for copy aberrant mutations
	return a dict that maps
	mutation -> sample -> mutational fractional copynumber
	mode can be 'mutational' or 'segmental'
	'''
	if mode not in ['mutational','segmental']: return {}
	frac_cns = {}
	vec = tree_df_to_vec(true_tree)
	anc = vec_to_anc_matrix(vec)
	cna_ids = [i for i in (set(exp['CNA_IDS'])) if not isinstance(i, float)]
	#print(cna_ids)
	impacted_muts = list(map(lambda x : x.replace('cna_', ''), cna_ids))

	for mut in impacted_muts:
		#print(f'looking up {mut}',list(set(list(exp['SNV_IDS'].values))))
		node_of_mut = exp[exp['SNV_IDS'] == mut]['NODE_ID'].iloc[0]
		node_of_cna = exp[exp['CNA_IDS'] == 'cna_'+mut]['NODE_ID'].iloc[0]
		frac_cns[mut] = {}
		for sample in samples_chosen:
			frac_cns[mut][sample] = 0.
			for node in vec.keys():
				if mode == 'mutational': 
					if node == 'ROOT': continue
					ncopies = calc_mut_copies(exp, mut, node, node_of_mut, node_of_cna, anc)
				else: # for segmental fractional CN, count the root
					ncopies = calc_seg_copies(exp, mut, node, node_of_mut, node_of_cna, anc)
				nodefreq = exp[((exp['NODE_ID'] == node) & 
								(exp['SAMPLE_IDS'] == sample))]['SAMPLE_NODE_FREQUENCIES']
				nodefreq = nodefreq.iloc[0]
				frac_cns[mut][sample] = frac_cns[mut][sample] + (ncopies * float(nodefreq))
	return frac_cns


def parse_frac_cns_patient(decifer, patient, mode='segmental'):
	state_cols = [col for col in decifer.columns if 'state' in col and col != 'normal_state'] # TODO remove unused var
	prop_cols = [col for col in decifer.columns if 'prop' in col and col != 'normal_prop']

	neoantigens = read_neoantigens(patient)
	frac_cns = {}
	for mut in decifer['mut_index'].values:
		#if mut not in neoantigens['Variant_key'].values: continue
		if (decifer[decifer['mut_index'] == mut].tumor1_state.unique() == ['1|1']) and (decifer[decifer['mut_index'] == mut].tumor2_state.unique() == ['1|1']): continue
		#if (decifer[decifer['mut_index'] == mut].tumor1_state.unique() == ['1|1']): continue
		frac_cns[mut] = {}

	for mut in frac_cns.keys():
		for sample in decifer[(decifer['mut_index'] == mut)]['sample'].values:
			row = decifer[(decifer['mut_index'] == mut) & (decifer['sample'] == sample)].squeeze()
			row = row.dropna()
			frac_cns[mut][sample] = 0.
			if mode == 'segmental':
				for prop_col in [col for col in prop_cols if col in row]:
					state_col = prop_col.replace('prop', 'state')
					tot_cn = sum(map(int, row[state_col].split('|'))) # BUG cannot split NaN copy number states
					frac_cns[mut][sample] = frac_cns[mut][sample] + (tot_cn * row[prop_col])
				frac_cns[mut][sample] = frac_cns[mut][sample] + (2. * row['normal_prop'])
			else:
				raise Exception("Not Yet Implemented")
	return frac_cns


def node_freqs_from_tree(node, tree):
	res = {}
	samples = tree['SAMPLE_IDS'].values[0].split(',')
	freqs = tree[tree['NODE_ID'] == node].squeeze()['SAMPLE_NODE_FREQUENCIES'].split(',')
	for i in range(len(samples)):
		sample = samples[i]
		if sample == 'MASTER_SAMPLE': continue
		freq = freqs[i]
		res[sample] = freq
	return res


##########################################
# Accuracy Scoring and Objective Helpers #
##########################################

def get_actual_node(mut, truedcfs, tree, true_tree, treetype):
	if treetype == 'ground_truth':
		return truedcfs[truedcfs['MUT_ID'] == mut]['NODE_ID'].iloc[0]
	else:
		assert(treetype == 'pairtree')
		return '0'


def update_accuracy(mut_correct, cna_correct, node_of_mut, node_of_cna,
					incorrect_muts, mut, truedcfs, tree, true_tree, tree_type, cna):
	#actual_node = truedcfs[truedcfs['MUT_ID'] == mut]['NODE_ID'].iloc[0]
	actual_node = get_actual_node(mut, truedcfs, tree, true_tree, 'ground_truth')
	if str(actual_node) == node_of_mut: 
		mut_correct += 1
	else:
		print(f"{mut} incorrectly assigned to Node {node_of_mut} instead of {actual_node}.")
		incorrect_muts.append(mut)

	if tree_type == 'ground_truth':
		if correct_cna_placement(mut, node_of_cna, tree):
			cna_correct += 1
		else:
			print(f"For {mut}: {cna} assigned incorrectly to {node_of_cna}.")
	else: # no method of scoring right now. so just add 1
		cna_correct += 1

	return mut_correct, cna_correct, incorrect_muts



##########################################
# Functions for Analyzing an optimal ILP #
##########################################

def grab_vars(model, nodes):
	A, B, presence, delta = {}, {}, {}, {}
	abs_var1, abs_var2 = {}, {}
	U_A, U_B, phi = {}, {}, {}
	for node in nodes:
		for v in model.getVars():
			if v.VarName == f'A_{node}':
				A[node] = v.x
			if v.VarName == f'B_{node}':
				B[node] = v.x
			if v.VarName == f'pres_{node}':
				presence[node] = v.x
			if v.VarName == f'delta_{node}':
				delta[node] = v.x
			if v.VarName == 'phi_A':
				phi['A'] = v.x
			if v.VarName == 'phi_B':
				phi['B'] = v.x
			if v.VarName == f'U_A_{node}':
				U_A[node] = v.x
			if v.VarName == f'U_B_{node}':
				U_B[node] = v.x
	return A, B, presence, delta, U_A, U_B, phi


def variant_content(presence, node_freqs, A, B, nodes, sample):
	return sum([presence[node]*node_freqs[node][sample]*(A[node]+B[node]) for node in nodes])


def total_content(asg, node_freqs, nodes, sample): 
	return sum([(asg[node][0] + asg[node][1]) * node_freqs[node][sample] for node in nodes])


def get_impliedvaf(mut, sample, model):
	A,B,presence,delta,U_A,U_B,phi = grab_vars(model, vec.keys())
	v = variant_content(presence,freqs,A,B,vec.keys(),sample)
	t = total_content(wrong_asg, freqs, vec.keys(), sample)
	return v/t


def all_impliedvafs(mut, model):
	return {s: get_impliedvaf(sample, )}


#################
# MISCELLANEOUS #
#################

def state2tup(state):
	return tuple(state.split('|'))


def node_of_cn(mut, df):
	node = df[df['CNA_IDS'].str.contains('cna_'+mut)].squeeze()['NODE_ID']
	return node


def node_of_mut(mut, df):
	df = df[df['NODE_ID'] != 'ROOT']
	node = df[df['SNV_IDS'].str.contains(mut)].squeeze()['NODE_ID']
	print(f'looking for {mut}:', node)
	if isinstance(node, pd.Series) and len(node) > 0: node = node.iloc[0]
	elif isinstance(node, pd.Series) and len(node) == 0: node = '0'
	else: node == '0'
	return node



