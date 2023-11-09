""" detopt.py

`DETOPT` is a a combinatorial optimization method for DETermining 
Optimal Placement in Tumor progression history of single nucleotide 
variants (SNVs) from the genomic regions impacted by copy-number 
aberrations (CNAs) using multi-sample bulk DNA sequencing data
"""

__version__ = '0.1.0'
__author__ = 'Chih Hao Wu'

import argparse
from copy import deepcopy
import os
import re
import time
from tqdm import tqdm
import pandas as pd
# import polars as pl

from detopt_placements import optimise
from detopt_util import state2tup, difference_from_diploid


def set_n_threads() -> int:
    """Set the number of threads

    Set the number of threads to the number of `SLURM_CPUS_PER_TASK` allocated
    or 4, by default if `DETOPT` is not being run on a machine in a cluster.

    Returns:
        int: The number of threads to use
    """
    n_threads = os.environ.get('SLURM_CPUS_PER_TASK')
    if not n_threads: 
        n_threads = 4
    else: 
        n_threads = int(n_threads)
    print(f"Using {n_threads} threads")
    return n_threads


def parse_tree(df: pd.DataFrame, samples: list) -> tuple:
    """Parse `tree_file` to get tree topology and sample node frequencies

    <description TBD>

    Args:
        df (pd.DataFrame): a `pd.DataFrame` for `tree_file`
        samples (list): a (sub)set of samples in entries of df in columns `NODE_FREQUENCIES`

    Returns:
        tuple:
            vec (dict): A mapping from a node to its parent, given from the tree
            freqs (dict): The sample node frequencies for each sample, in each node
    """
    vec, freqs = {}, {}

    for _, row in df.iterrows():
        node = row.NODE_ID
        freqs[node] = {}
        vec[node] = row.PARENT_ID

        row_samples, row_nodefreqs = row.SAMPLE_IDS.split(','), row.NODE_FREQUENCIES.split(',')
	
        for i, sample in enumerate(row_samples):
            if sample in samples:
                freqs[node][sample] = float(row_nodefreqs[i])
    
    return (vec, freqs)


def node_desc_presence(freqs: dict, vec: dict) -> dict:
    """In a post-order traversal of the tree, calculate the descendant cell fraction of a
    subclone, respresented by nodes.

    The node sample prevalence, or the descendant cell fraction, can be obtained by summing 
    the frequency in node v and all its descendents.

    Args:
        freqs (dict): The sample node frequencies for each sample, in each node
        vec (dict): A mapping from a node to its parent, given from the tree

    Returns:
        node_prevalence (dict): The node sample prevalence of nodes in tree
    """
    tree_struc, freqs_dict = deepcopy(vec), deepcopy(freqs)
    
    node_frequency = pd.DataFrame.from_dict(freqs_dict).T
    node_frequency.drop("ROOT", axis=0, inplace=True)
    
    del tree_struc['ROOT']
    
    while tree_struc != {}:
        leaves = []
        for node in tree_struc.keys():
            if node not in tree_struc.values():
                leaves.append(node)
                if tree_struc[node] == "ROOT": continue
                parent = tree_struc[node]
                node_frequency.loc[parent,:] = (node_frequency.loc[parent,:] + 
                                                node_frequency.loc[node,:]
                                                )
        for l in leaves:        
            del tree_struc[l]
    
    return node_frequency.T.to_dict()       # node_prevalence


# TODO put back into detopt_util.py
def parse_frac_cns_patient(snv: pd.DataFrame, mode: str='segmental', 
                           *mutation_list: list
                           ) -> dict: 
    """<main>

    <desc>

    Args:
        snv (pd.DataFrame): A `pd.DataFrame` containing information about
        read counts and allele-copy number calls of SNVs
        mode (str): The mode to parse copy-number state. Currently, just 
        segmental fractional copy number state
        mutation_list (list): optional, a list which elements are a 
        subset of mutations

    Returns:
        frac_cns (dict): The fractional copy number state of each segment
        on which a mutation is on, for each sample

    Raises:
        ValueError: if `mode` == 'segmental'; other modes not yet implemented
    """
    prop_cols, state_cols = [], []
    frac_cns = {}

    for col in snv.columns:
        
        is_pattern = bool(re.match(r'.*(?:_prop|_state)', col))
        match is_pattern:
            case True:
                suffix = col.split('_')[-1]
                if suffix == 'prop' and 'normal' not in col:
                    prop_cols.append(col)
                elif 'normal' not in col:
                    state_cols.append(col)        
            case _:
                pass

    total_rows = len(snv)

    for i, row in snv.iterrows():

        #print(f"Progress: {i+1}/{total_rows}", end='\r')
        
        mut = row.mut_index
        if mutation_list and (mut not in mutation_list): continue       # if a list of SNVs is provided, then use only those SNVs
        if ({row[state] for state in state_cols} == {'1|1'}): continue  # skip copy-number neutral variants

        sample = row['sample']

        if mut in frac_cns.keys():  pass
        else:
            frac_cns[mut] = {}
 
        frac_cns[mut][sample] = 0.

        if mode == 'segmental':
            for prop_col, state_col in zip(prop_cols, state_cols):
                total_copy_number= sum(map(int, row[state_col].split('|')))
                frac_cns[mut][sample] = (frac_cns[mut][sample] + 
                                         (total_copy_number * row[prop_col])
                                         )
            frac_cns[mut][sample] = frac_cns[mut][sample] + (2. * row['normal_prop'])
        else:
            raise ValueError('Not yet implemented')

    #print(f"Progress (frac_cns): {i+1}/{total_rows}")
    return frac_cns


def parse_mutation_data(snv: pd.DataFrame, mut: str, samples: list) -> tuple:
    """<main>

    <desc>

    Args:
        snv (pd.DataFrame): A `pd.DataFrame` containing information about
        read counts and allele-copy number calls of SNVs
        mut (str): A uniquely identified SNV
        samples (list): A (sub)set of samples

    Returns:
        tuple:
            states (list): A list of all states (str, str) for the copy number 
            states of the region on which `mut` is present
            vafs (dict): The variant allele frequency (VAF) of `mut` in each sample 
    """
    snv = snv[(snv.mut_index == mut) & (snv['sample'].isin(samples))]

    # TODO consider placing as standalone function
    prop_cols, state_cols = [], []
    for col in snv.columns:
    
        is_pattern = bool(re.match(r'.*(?:_prop|_state)', col))
        match is_pattern:
            case True:
                suffix = col.split('_')[-1]
                if suffix == 'prop' and 'normal' not in col:
                    prop_cols.append(col)
                elif 'normal' not in col:
                    state_cols.append(col)        
            case _:
                pass

    assert len(state_cols) == len(prop_cols), "Number of state cols and prop cols should be the same"

    states, vafs = [0] * len(state_cols), {}

    # populate `states`
    for i, state_col in enumerate(state_cols):
        state = snv.iloc[0][state_col]
        states[i] = state2tup(state)

    # populate `vafs` per sample, where each sample is a row in `snv``
    for _, row in snv.iterrows():

        vaf = row.var_reads / (row.ref_reads + row.var_reads)
        sample = row['sample']

        if sample in samples:
            vafs[sample] = vaf

    return states, vafs


def assignCNAs_update(mutation: str, 
                      vec: dict, 
                      node_freqs: dict,
                      states: list, 
                      vafs: dict, 
                      frac_cns: dict, 
                      sample_node_frac: dict,
                      n_threads: int, 
                      cna_reg: float,
                      vaf_reg: float=1.
			        ):
    """<main>

    <desc>

    Args:
        mutation (str): A uniquely identified copy-number aberrant SNV
	    vec (dict): A mappping from each node (subclone) to its uniquely 
        identified parental node (subclone). The root has NaN ancestor, 
        for sake of non-empty parental node.
		node_freqs (dict): The sample node frequencies for each sample, 
        in each node
	    states (list): A list of all states (str, str) for the copy number 
        states of the region on which `mut` is present
        vafs (dict): The variant allele frequency (VAF) of `mut` in each sample
	    frac_cns (dict): The fractional copy number state of each segment
        on which a mutation is on, for each sample
        n_threads (int): The number of threads used by Gurobi (or eq. solver)

    Returns:
        tuple:
            objective (float): The objective score returned by `DETOPT`, 
            refer to the objective function in the manuscript
            node_assigned (int): The node (subclone) to which the mutational
            event was assigned for the `mut`
            state_assignment (int): The node (subclone) to which the copy
            number change event(s) was (were) assigned for the `mut`
    """
    states = [tuple(map(int, cns)) for cns in states]
    states = sorted(states, key=difference_from_diploid)        # sorting by total number of copy number changes for both alleles
    states = list(set(states))

    obj_val, model_vars = optimise(states,
                                   vec,
                                   node_freqs,
                                   vafs,
                                   frac_cns[mutation],
                                   sample_node_frac,
                                   vaf_reg,
                                   cna_reg
                                   )

    return obj_val, model_vars


def write_detopt_assignments(assignments: dict, max_n_states: int, filename_prefix: str):
    """<main>

    <desc>

    Args:

    Returns:

    """
    # write output for assignments 
    header_cols = ['mut_index','snv_assignment']
    for n in range(max_n_states):
        header_cols += [f'cna_state_{n}', f'cna_assignment{n}'] 
    header = '\t'.join(header_cols) + '\n'

    with open(f'{filename_prefix}_assignments.tsv', 'w') as file:
        file.write(header)

        for mut in assignments:
            cna_str = [f'{state}\t{asg}' for state, asg in assignments[mut]['cna'].items()]
            cna_str = '\t'.join(cna_str)
            file.write(f"{mut}\t{assignments[mut]['snv']}\t{cna_str}\n")

    return None


def write_detopt_copy_number(copies: dict, filename_prefix: str, mutant_copies: bool=False):
    """<main>

    <desc>

    Args:

    Returns:
    
    """
    # write output for allele-specific (mutant) copy number state assignments
    rows = [0] * len(copies)
    muts = copies.keys()

    for i, mut in enumerate(muts):

        alleles = [copies[mut]['a'], copies[mut]['b']]

        mut_copies_both_alleles = {}
        for node in copies[mut]['a'].keys():
            mut_copy_states = ','.join([str(int(allele[node])) for allele in alleles])
            mut_copies_both_alleles[node] = f"{mut_copy_states}"

        rows[i] = mut_copies_both_alleles

    mut_copies_df = pd.DataFrame(rows)
    mut_copies_df['mut_id'] = muts
    mut_copies_df.set_index('mut_id', inplace=True)

    if mutant_copies:
        mut_copies_df.to_csv(f'{filename_prefix}_mut_copy_numbers.tsv', sep='\t')
    else:
        mut_copies_df.to_csv(f'{filename_prefix}_tot_copy_numbers.tsv', sep='\t')
    
    return None


def main():
    parser = argparse.ArgumentParser(
                        prog='DETOPT',
                        description='(DETermining Optimal Placement in Tumor progression history)',
                        prefix_chars='-+',
                        add_help=False
                        )
    parser.add_argument('+h', '++help', action='help', 
                        help='show this help message and exit'
                        )
    parser.add_argument('-p', '--cna_reg', nargs='?', default=.25, 
                        type=float, help='regularization weight (default: %(default)s)'
                        )
    parser.add_argument('-h', '--n_samples', nargs='?', type=int,
                        help='number of samples'
                        )
    parser.add_argument('-d', '--data_dir', required=True, type=str,
                        help='directory containing required `snv_file` and `tree_file` files'
                        )
    parser.add_argument('-o', '--out', required=True, type=str,
                        help='output filename prefix, optionally with filepath'
                        )
    parser.add_argument('-s', '--snv_file', required=True, type=str,
                        help='`snv_file` file containing information about read counts\
                              and allele-copy number calls of SNVs'
                        )
    parser.add_argument('-t', '--tree_file', required=True, type=str,
                        help='`tree_file` file containing information about the base tree'
                        )
    args = parser.parse_args()

    start_time = time.time()

    snv_file, tree_file, data_dir  = args.snv_file, args.tree_file, args.data_dir
    n_samples = args.n_samples
    vaf_reg, cna_reg = 1, args.cna_reg

    n_threads = set_n_threads()

    if not data_dir.endswith('/'):
        data_dir = data_dir+'/'
    else: pass
    
    snv_df = pd.read_csv(f'{data_dir}{snv_file}', sep='\t', header=0, index_col=False)

    tree_df = pd.read_csv(f'{data_dir}{tree_file}', sep='\t', header=0, index_col=False)
    samples = tree_df.loc[0, 'SAMPLE_IDS'].split(',')
    if not n_samples:
        samples = samples[:n_samples]       # optional to use fewer samples
    else: pass

   #print(f"There are {len(samples)} samples")

    vec, freqs = parse_tree(tree_df, samples)
    sample_node_frac = node_desc_presence(freqs, vec)
    frac_cns = parse_frac_cns_patient(snv_df, mode='segmental')

    muts = frac_cns.keys()

    assignments = {}
    mut_copies, tot_copies = {}, {}
    max_n_states = 1

    for mut in tqdm(muts, colour='green'):
    
        # parse the data and solve the model for each mutation individually
        mut_states, mut_vafs = parse_mutation_data(snv_df, mut, samples)
        mut_states = [s for s in mut_states if s != ('1', '1')]

        n_unique_states = len(set(mut_states))
        if n_unique_states > max_n_states: max_n_states = n_unique_states

        obj_val, model_vars = assignCNAs_update(mutation=mut,
                                                vec=vec,
                                                node_freqs=freqs,
                                                states=mut_states,
                                                vafs=mut_vafs,
                                                frac_cns=frac_cns,
                                                sample_node_frac=sample_node_frac,
                                                n_threads=n_threads,
                                                vaf_reg=vaf_reg,
                                                cna_reg=cna_reg
                                                )

        if obj_val >= 0:
            node_assignment, cna_assignment, mut_a, mut_b, tot_a, tot_b = model_vars
            #print(f'{mut}:\tOBJ_VAL -> {obj_val}\tSNV -> {node_assignment}\tCN -> {cna_assignment}')

            assignments[mut] = {'snv': node_assignment, 
                                'cna': cna_assignment
                                }
            
            mut_copies[mut] = {'a': mut_a,
                               'b': mut_b
                               }
            
            tot_copies[mut] = {'a': tot_a,
                               'b': tot_b
                               }
        else:
            #print(f'{mut} has no assignments satisfying constraints')
            pass

    print(f"Total Runtime: {round(time.time() - start_time, 2)} seconds")     # TODO add progress bar

    # write outputs
    write_detopt_assignments(assignments, max_n_states, args.out)

    # write output for allele-specific mutant copy number state assignments 
    write_detopt_copy_number(mut_copies, args.out, mutant_copies=True)

    # write output for allele-specific total copy number state assignments 
    write_detopt_copy_number(tot_copies, args.out)

    print(f"Wrote outputs to files: {args.out}*") 

if __name__ == '__main__':
    main()
