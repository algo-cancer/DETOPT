# detopt.py
# please note that this file uses some functions from detopt_dev.py

import os, json, argparse
import pandas as pd
import polars as pl
from ast import literal_eval as lit
from detopt_placements import optimise
from detopt_util import *
from itertools import combinations, permutations
import time
from copy import deepcopy

C_MAX = 100

parser = argparse.ArgumentParser(description='DETOPT (DETermining Optimal Placement in Tumor progression history)', add_help=False)
parser.add_argument('--help', action='help',
    help='show this help message and exit')
parser.add_argument('-s', '--snv_file', required=True, type=str,
                    help='*.input.tsv file containing information about read counts and allele-copy number calls of SNVs')
parser.add_argument('-t', '--tree_file', required=True, type=str,
                    help='*.tree file containing information for each subclone, the PARENT_ID, MUTATIONS_AT_NODE, SAMPLE_IDS, NODE_FREQUENCIES')
parser.add_argument('-h', '--samples', required=True, type=str,
                    help='number of samples')
parser.add_argument('-p', '--cna_weight', default=0.25, type=str,
                    help='regularization value/weight')
parser.add_argument('-d', '--data_dir', required=True, type=str,
                    help='directory containing required snv and tree files')
parser.add_argument('-o', '--out', required=True, type=str,
					help='output filename, optionally with filepath')
args = parser.parse_args()

decifer_file = args.snv_file
tree_file = args.tree_file
n_arg = int(args.samples)
DATA_DIR = args.data_dir

vaf_reg, cna_reg = 1, args.cna_weight

num_threads = os.environ.get('SLURM_CPUS_PER_TASK')
if not num_threads: 
	num_threads = 4
else: num_threads = int(num_threads)
print(f"Using {num_threads} threads.")


def parse_tree(df, num_samples, true_tree=[], is_true_tree=False):
	# subsample `num_samples` number of samples
	samples = df.iloc[0]['SAMPLE_IDS'].split(',')
	if 'MASTER_SAMPLE' in samples: samples.remove('MASTER_SAMPLE')	
	samples_chosen = samples[:num_samples]
	vec, freqs = {}, {}
	for (i, row) in df.iterrows():
		node = row['NODE_ID']
		freqs[node] = {}
		vec[node] = row['PARENT_ID']
		if not is_true_tree: # pairtree nodefreqs
			nodefreqs = row['NODE_FREQUENCIES'].split(',')
		else: # ground truth
			nodefreqs = row['SAMPLE_NODE_FREQUENCIES'].split(',')
		samples = row['SAMPLE_IDS'].split(',')
		
		#if is_true_tree:
		for i in range(len(samples)):
			sample = samples[i]
			if sample in samples_chosen:
					freqs[node][sample] = float(nodefreqs[i])
	"""if not is_true_tree:
		assert(true_tree != []), 'Must provide ground truth tree as well!'
		for node in vec.keys():
			if node == 'ROOT':
				continue
			freqs[node] = freqs_from_true_tree(node, true_tree) #calc_observed_prevalence(node, df, decifer_df, vec, samples_chosen, vafs, 'MUTATIONS_AT_NODE')
		freqs = append_root_prevalence(freqs)
	"""
	return vec, freqs, samples_chosen


def parse_data(df, mut):
	states, vafs = [], {}
	df = df[df['mut_index'] == mut]
	state_cols = [col for col in df.columns if 'state' in col and 'normal' not in col]
	#print(state_cols, df)

	for col in state_cols:
		#print(col)
		state = tuple(df.iloc[0][col].split('|')) # BUG cannot split NaN copy number states FIXED
		states.append(state)

	for (i, row) in df.iterrows():
		var = row['var_reads']
		tot = row['ref_reads'] + row['var_reads']
		if var > 0 and tot == 0:
			vaf = 1.
		elif var == 0 and tot == 0:
			vaf = 0.
		else:
			vaf = row['var_reads'] / (row['ref_reads'] + row['var_reads'])
		sample = row['sample']
		vafs[sample] = vaf
	return states, vafs

def parse_mutation_data(df, mut, samples_chosen):
	states, vafs, = [], {}
	df = df[df['sample'] != 'MASTER_SAMPLE']
	df = df[df['sample'].isin(samples_chosen)]
	df = df[df['mut_index'] == mut]
	state_cols = [col for col in df.columns if 'state' in col and 'normal' not in col]
	prop_cols = [col for col in df.columns if 'prop' in col and 'normal' not in col]

	assert(len(state_cols) == len(prop_cols))

	# populate `states`	
	for i in range(len(state_cols)):
		state = df.iloc[0][state_cols[i]]
		state = state2tup(state)
		states.append(state)

	# populate `vafs` per sample (each sample is a row in this df)
	for (i, row) in df.iterrows():
		# `vafs` first
		vaf = row['var_reads'] / (row['ref_reads'] + row['var_reads'])
		sample = row['sample']
		if sample in samples_chosen:
			vafs[sample] = vaf

	#print(states, vafs)
	return states, vafs

def propagate(assignment, anc, vec, nodes):
	'''
	fill in the copy numbers of all descendants of aberrant nodes as the same state
	fill in all other nodes as the diploid state

	propagate to descendants starting from the node closest to the root
	'''
	assignment_copy = copy.deepcopy(assignment)
	nodes_left = copy.deepcopy(nodes)
	nodes_left.remove('ROOT')

	nodes_ordered = sorted(assignment.keys(), key=lambda node: distance_from_root(node, vec))
	for aberrant_node in nodes_ordered:
		for node in descendants(aberrant_node, anc):
			assignment_copy[node] = assignment.get(node, assignment[aberrant_node]) #assignment[aberrant_node] # propagate the same state to descendants
			if node in nodes_left:
				nodes_left.remove(node)
	for node in nodes_left: # remaining nodes get diploid state
		assignment_copy[node] = (1,1)

	assignment_copy['ROOT'] = (1,1)
	return assignment_copy


def satisfies_isa(state_assignment, vec):
	# for CNAs, we use the weak ISA states that if all copies of a seg are lost
	# then they cannot be gained again. So no 0 -> 1 pattern
	for node in vec.keys():
		parent = vec[node]
		state_of_child = state_assignment[node]
		state_of_parent = state_assignment.get(parent, (1,1))
		if int(state_of_parent[0]) == 0 and int(state_of_child[0]) > 0:
			return False
		if int(state_of_parent[1]) == 0 and int(state_of_child[1]) > 0:
			return False
	return True


def assignCNAs_update(
		                mutation: str, 
	                    vec: dict, 
			            node_freqs: dict, 
			            states: list, 
			            vafs: dict,
						frac_cns: dict, 
						sample_node_frac: dict,
			            num_threads: int
			        ):
	"""given mutation and following, assign copy number state to each node in tree

	UPDATED BRUTEFORCE desc. TODO use `networkx` for easier handling of ances/desc relationships, etc.

	Parameters
	----------
	mutation: str
		identifier for mutation that is being placed
	vec: dict
		maps a node (subclone) to its ancestor (subclone), root has NaN ancestor
		key is node, value is ancestor
	node_freqs: dict
		frequency of each sample at the node
	states: list (tuple)
		copy number state assigned to mutation, of segment
	vafs: dict
		variant allele frequencies of mutation in each sample
	frac_cns: dict
		fractional copy number of each mutation in each sample {mut1: {sample1: cn1 , sample2: cn2}, mut2: {sample1: cn1 , sample2: cn2}}
	num_threads: int
		number of threads
	"""

	# assignment_scores = {} # maps assigned node to its score from our qilp	
	# anc = vec_to_anc_matrix(vec) # each key corresponds to a node, values are all other nodes with values 1 if is descendent of node, else 0

	# nodes = list(vec.keys())
	# non_root_nodes = [n for n in nodes if n != "ROOT"]

	states = [tuple(map(int, cns)) for cns in states]

	states = sorted(states, key=difference_from_diploid) # sorting by total number of copy number changes (either up/down) for both alleles

	states = list(set(states))

	objective, node_assigned, state_assignment = optimise(
												states, 
												vec,
												node_freqs, 
												vafs, 
												frac_cns[mutation],
												sample_node_frac,
												vaf_reg, 
												cna_reg
												)

	#print(f"\tOptimal assignment: Mutation {mutation} attaches to Node {node_assigned}")

	# bestasg = values_to_tuple(lit(bestasg)) # revert hashed assignment
	return objective, node_assigned, state_assignment


def get_impliedvaf(v, t, actual):
	if t == 0.0: return actual
	else: return v / t


def parse_frac_cns_polars_seg(true_tree, exp, lookup_pl, samples_chosen, mode='segmental'):
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
	impacted_muts = list(map(lambda x : x.replace('cna_', ''), cna_ids))

	for mut in impacted_muts:
		node_of_mut = exp[exp['SNV_IDS'] == mut]['NODE_ID'].iloc[0]
		node_of_cna = exp[exp['CNA_IDS'] == 'cna_'+mut]['NODE_ID'].iloc[0]
		frac_cns[mut] = {}
		for sample in samples_chosen:
			frac_cns[mut][sample] = 0.
			for node in vec.keys():
				if mode == 'mutational': 
					if node == 'ROOT': continue
					ncopies = calc_mut_copies(exp, mut, node, node_of_mut, node_of_cna, anc) # never goes here
				else: # for segmental fractional CN, count the root
					ncopies = calc_seg_copies(exp, mut, node, node_of_mut, node_of_cna, anc) 
				nodefreq = lookup_pl.filter((pl.col("NODE_ID") == node) & (pl.col("SAMPLE_IDS") == sample)).select("SAMPLE_NODE_FREQUENCIES").item()
				frac_cns[mut][sample] = frac_cns[mut][sample] + (ncopies * float(nodefreq))
	return frac_cns


def node_desc_presence(freqs, vec):
	# NEW oct 5, 2023 new constr. returns node sample prevalence summing frequency in node v and its descendents

	tree = deepcopy(vec)
	freqs_dict = deepcopy(freqs)

	node_frequency = pd.DataFrame.from_dict(freqs_dict).T
	node_frequency = node_frequency.drop("ROOT", axis=0)

	del tree['ROOT']

	while tree != {}:
		r = []
		for n in tree.keys():
			if n not in tree.values():
				r.append(n)
				if tree[n] == "ROOT": continue
				p = (tree[n])
				node_frequency.loc[p,:] = node_frequency.loc[p,:] + node_frequency.loc[n,:] 
		for s in r:        
			del tree[s]
			
	return node_frequency.T.to_dict()


def main():
	
	start_time = time.time()
	outstring = ""

	num_threads = os.environ.get('SLURM_CPUS_PER_TASK')
	if not num_threads: 
		num_threads = 4
	else: num_threads = int(num_threads)
	print(f"Using {num_threads} threads.")
	
    # ALL HARDCODED
	decifer_df = pd.read_csv(f'{DATA_DIR}/{decifer_file}', sep='\t', header=0, index_col=False)
	tree = pd.read_csv(f'{DATA_DIR}/{tree_file}', sep='\t', header=0, index_col=False)
	samples = tree.iloc[0]['SAMPLE_IDS'].split(',')
	if 'MASTER_SAMPLE' in samples: samples.remove('MASTER_SAMPLE')
	samples_chosen = samples[:n_arg] # consider shuffling
	#print(samples_chosen)
	muts_tested = sorted(sample_aberrant_muts(decifer_df))
	#all_mutations = list(set(decifer_df['mut_index'].values))
	print("MUTS TESTED:", sorted(muts_tested))
	#vafs = parse_vafs(decifer_df, samples_chosen)
	vec, freqs, samples_chosen = parse_tree(tree, n_arg, is_true_tree=False)
	sample_node_frac = node_desc_presence(freqs, vec)
	#exp = expl(tree)

	# Jul 25 speed up
	#quick_lookup = exp[['NODE_ID','SAMPLE_IDS','SAMPLE_NODE_FREQUENCIES']]
	#quick_lookup = quick_lookup.sort_values(['NODE_ID','SAMPLE_IDS']).drop_duplicates()
	#pl_quick = pl.from_pandas(quick_lookup)
	
	# frac_cns = parse_frac_cns(tree, exp, samples_chosen, 'segmental') # long running time, scales poorly with number of mutations and number of nodes
	#frac_cns = parse_frac_cns_polars_seg(tree, exp, pl_quick, samples_chosen, 'segmental')
	frac_cns = parse_frac_cns_patient(decifer_df, '111e6e61e1') # ??? placing only cn calls of neoantigens? ###CHANGE SAMPLE NAME HERE; fucks up without it # MAKE SURE TO ALSO CHANGE DETOPT UTIL.PY FILE
	muts_tested = frac_cns.keys()

	#print(frac_cns)
	print(vec)
	#print(sample_node_frac)
	
	#start_time = time.time()
	
	#for mut in muts_tested:
	for mut in ["11:70031751-70031751_G>T","11:6477707-6477707_G>A","9:4833182-4833182_C>G","16:88947822-88947822_C>T","9:5090566-5090566_G>T","6:152419920-152419920_T>A","3:178952085-178952085_A>T","16:89350644-89350644_G>C"]:
		print(mut)

		#print(mut) # 9:90252946-90252946 A>C not placed; reason is because there are more than 1 tumour states in decifer.input and it is looking at the one with 1|1
		
		# parse the data and solve the model for each mutation individually
		states, vafs_of_mut = parse_mutation_data(decifer_df, mut, samples_chosen)
		states = [s for s in states if s != ('1', '1')] # FIXME does this force mutational loci to only have 1 copy number state change event? How is it exploring assignments between?
		# ADD BACK print(states)
		if states == [('1', '1')]: continue # 9:90252946-90252946 A>C not placed (because of this)
		#states.remove(('1', '1'))
		# ADD BACK print(f"\nmutation: {mut}, abberant copy number states: {states}\n")

		#start_time_cn_segment = time.time()
		
		objective_score, node_of_mutation, state_assignment = assignCNAs_update(
			mut, vec, freqs, states,
			vafs_of_mut, frac_cns, 
			sample_node_frac,
			num_threads
		)
		
		#print(f"\tRuntime for segment {time.time() - start_time_cn_segment}")
		#for sample, vaf in vafs_of_mut.items():
		#	print(f"{sample: <45}{round(vaf, 5)}")

		# print("\n")
		 
		#node_of_mutation, bestasg, bestmodel = assignCNAs_BTFE(
	    #	mut, vec, freqs, states,
		#	vafs_of_mut, frac_cns, 
		#	num_threads
		#	)
		
		try:
			outstring += f"{mut}\t{objective_score}\t{node_of_mutation}\t{':'.join([f'{s[0][0]} {s[0][1]} {s[1]}' for s in state_assignment])}\n"
			print(f'\tMutation {mut} assigned to Node {node_of_mutation}')
		except TypeError:
			print(f"MUTATION: {mut} has no assignment satisfying constraints")
			
		### UNCOMMENT THIS print(f'\tMutation {mut} assigned to Node {node_of_mutation}') ###
		#print(f"vars {bestmodel}")
		#print(f"Best objective for mutation {mut} : {bestmodel.objVal}")

	
		"""mut_assignments[mut] = {}
		mut_assignments[mut]['Node'] = node_of_mutation
		mut_assignments[mut]['Objective'] = bestmodel.objVal
		A, B, delta, psi = grab_vars(bestmodel, bestasg.keys())
		for node in bestasg.keys():
			if node == 'ROOT': 
				bestasg[node] = (bestasg[node][0], bestasg[node][1], 0) # mutational CN = 0
			else: 
				bestasg[node] = (bestasg[node][0], bestasg[node][1], int(A[node] + B[node]))
		mut_assignments[mut]['Assignment'] = bestasg"""
        
	#print(f"\nRuntime: {time.time() - start_time}")

	with open(f"{args.out}.detopt.tsv", "w") as o:
		o.write(outstring)
	
		"""with open(f"{args.out}.detopt.tsv", "w") as o:
			o.write('Variant key\tNode\tObjective\tAssignment\n')
			for mut in mut_assignments:
				asg = mut_assignments[mut]['Assignment']
				node_assigned = mut_assignments[mut]['Node']
				line = [mut, node_assigned,
						mut_assignments[mut]['Objective'], asg]
				o.write(floatjoin(line, '\t') + '\n')"""

	"""print(f"Wrote {args.out}.detopt.tsv")"""
	print(f"Total Runtime: {time.time() - start_time}")

	#return A, B, delta, psi
	return 0
	


def obj(mut):
	# For testing: retrieve all objective values for a given mutation
	return {x : abs(vafs[mut][x] - impliedvafs[x]) for x in vafs[mut].keys()}


if __name__ == '__main__': 
	#A, B, delta, psi = 
	main()