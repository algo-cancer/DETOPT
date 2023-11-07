"""hatchetconvert.py

Conversion of (required) SEG-UNC file from HATCHet and a mutation tab-
separated file (.tsv) into a DeCiFer-like input file for DETOPT.

`pandas` is a required package within the Python environment to run this
script. TODO modify to use `polars` for faster passing through mutations

This file can also be impored as a standalone module and the following 
functions may be called independently from the command line (however, 
some functions are dependent on the values, objects returned from 
other functions):

    * cn_assignment :
    * main :

Changelog:
"""

# TODO: implement as oop, with class structure
# TODO: rename mutations, table or list

import os, re, argparse, itertools
import pandas as pd

MAX_COPYNUMBER = 6
WRITE_HEADER = out_cols = ['mut_index',
                           'sample',
                           'var_reads',
                           'ref_reads',
                           'normal_state',
                           'normal_prop'
                           ]

parser = argparse.ArgumentParser(
    prog="hatchetconvert.py",
    description="Conversion of HATCHet to DeCiFer format"
)

parser.add_argument("-v", "--verbose", action="store_true", 
                    help="verbosity control, default False")
parser.add_argument("-1", "--hatchet", type=str, required=True, 
                    help="HATCHet seg.ucn file")
parser.add_argument("-2", "--mutations", type=str, required=True, 
                    help="mutations, tab-separated file")
parser.add_argument("--neutral", action='store_true',
                    help="return copy number neutral mutations ONLY")
parser.add_argument("-o", "--writeout", type=str, required=True, 
                    help="write to file name prefix")

args = parser.parse_args()

# FIXME: replace later with makefile to run examples
#if os.getcwd() != "/Users/wuchh/Desktop/wuchh/data/detopt-updates":
#    raise Exception("run from parent-most directory instead. exiting")


def get_unique(dataframe: pd.DataFrame, colname: str) -> list:
    """Returns unique values from column in a `pandas` dataframe as list
    
    Parameters
    ----------
    dataframe : pd.DataFrame
        mutations in `pandas` dataframe
    colname : str
        name of column in dataframe
    """
    return dataframe[colname].unique().tolist()


def cn_assignment(mutations: list, samples: list, cn_status: pd.DataFrame) -> dict: 
    """Assign mutations copy number states using segment copy number state
    also, assign cell prevalences to each copy number segment
    
    Parameters
    ----------
    mutations :
    samples :
    cn_states :

    Returns
    -------
    muts_assigned_cn : dict
    """
    
    unmatched_mut = set()

    keys = itertools.chain(*(map(lambda s: [(i, s) for i in mutations], samples)))
    assign_cn = {k: {} for k in keys}


    for mut in mutations:
        
        for sample in samples:
            chr, start, end = re.split(':|-|_', mut)[:3]

            inside_seg = (
                (cn_status["SAMPLE"] == sample) & 
                (cn_status["#CHR"] == chr) &
                (cn_status["START"] <= int(start)) &
                #(cn_status["END"] >= int(end))
                (cn_status["END"] > int(end))
            )
            """ FIXME THIS IS PROBLEM, lies on boundary of two CN-segments
                #CHR     START       END                           SAMPLE cn_normal  u_normal cn_clone1  u_clone1 cn_clone2  u_clone2
                2619   11  55000000  56000000  4214-1Met_FrTu_November_15_2016       1|1  0.487131       1|1       0.0       1|1  0.512869
                2628   11  56000000  56500000  4214-1Met_FrTu_November_15_2016       1|1  0.487131       0|1       0.0       1|0  0.512869 
                11:56000000-56000000_G>T
                MAYBE, by default choose to put this in CN
            """

            if inside_seg.sum() == 0:
                unmatched_mut.add(mut)
            else:
                mut_CNs = cn_status[inside_seg].squeeze()

                clones = [i.lstrip("u_") for i in mut_CNs.index if i.startswith("u_")]

                for clone in clones:
                    cn1, cn2 = mut_CNs["cn_"+clone].split('|')
                    assign_cn[(mut, sample)][clone] = (cn1, cn2, mut_CNs["u_"+clone])
        
    print(f"{len(unmatched_mut)}/{len(mutations)} mutations are unmatched to copy number affected segments")
    #print(unmatched_mut)
    #print(assign_cn[('9:90746970-90746970_G>-', '4435-1Met_FrTu_March_03_2021')])
    return (assign_cn, unmatched_mut)


def index_column(df: pd.DataFrame, colname: str, sort=False) -> pd.Series:
    """From a column in a `pandas` dataframe, using a list of its unique values, 
    return a `pandas` series where each value is replaced by the index position of that value in the list
    
    Parameters
    ----------
    df : pd.DataFrame
        `pandas` dataframe with column to index
    colname : str
        name of column to index
    sort_arg : bool
        (optional) sort list by substring
        
    Return
    ------
    """

    uniq = get_unique(df, colname)
    if sort:
        uniq.sort(key=lambda val: val[:4])
    return df[colname].apply(lambda value: uniq.index(value))


def write_conversion(filename: str, mutations: pd.DataFrame, copynumber: dict) -> None:
    """write output to hatchetconvert.py TODOC
    
    Parameters
    ----------
    filename: str
        Name of file to write to, optionally includes parent directory
    mutations: pd.DataFrame
        Mutations with copy number segments
    copynumber: dict
        Copy number segments with cell prevalences
    """

    lines = []

    for _, row in mutations.iterrows():

        row = row.astype(str).to_dict() 
        line = list(row.values())

        mutation, sample = row["character_label"], row["sample_label"]

        cn_seg = copynumber[(mutation, sample)]

        if not bool(cn_seg):
            continue # mutation not in cna affected segment

        # ONLY CN-neutral TODO FIXME
        #if cn_seg['clone1'][:2] != ('1','1') or cn_seg['clone2'][:2] != ('1','1'): continue
        #if cn_seg['clone1'][:2] != ('1','1'): continue
        #else: pass

        cn_states = [cn_seg[clone] for clone in cn_seg.keys()]
        cn_by_allele = [(int(cn1), int(cn2)) for cn1, cn2, _ in cn_states]
    
        # ONLY CN-neutral 
        if args.neutral and len(set(cn_by_allele)) != 1:
            continue # only unique state is "1|1", replaces above; we use this to obtain cn_neutral mutations only

        # enforce maximum copy number
        for e, cn in enumerate(cn_by_allele):
            a, b = cn
    
            if a > MAX_COPYNUMBER or b > MAX_COPYNUMBER:
                print(f"terminated at ({a}, {b})")
                break
            
            if e == len(cn_by_allele) - 1: # TODO: review from here
                for a, b, frac in cn_states: 
                    line = line + [f'{a}|{b}', f'{round(frac, 4)}']
                lines.append("\t".join(line))
        
        else:
            continue

    header = WRITE_HEADER
    for i in range(1, len(cn_seg)):
        header.append(f'tumor{i}_state')
        header.append(f'tumor{i}_prop')
        
    header = '\t'.join(header)
    lines.insert(0, header)

    with open(f"{filename}.decifer.input", 'w') as stream:
        stream.write('\n'.join(lines))


def main():
    """Read and reformat HATCHet and mutation files.

    Variables
    ---------
    cn_states : pandas.DataFrame
        copy number states of segments from HATCHet
    mutations : pandas.DataFrame
        mutations
    """

    cn_states = pd.read_csv(args.hatchet, sep="\t")
    cn_states["#CHR"] = cn_states["#CHR"].astype(str)

    mutations = pd.read_csv(args.mutations, sep="\t")
    mutations["Variant key"] = mutations["Variant key"].apply(lambda label: label.replace(' ','_'))
    mutations.drop(["cDNA change", "AA change"], axis=1, inplace=True)

    mutations.rename(
        {
            'Variant key' : 'character_label',
            'nreads_reference' : 'ref',
            'nreads_variant' : 'var',
            'DNA_sample_name' : 'sample_label'
        }, 
        axis=1,
        inplace=True
    )

    # FIXME: unsure problem
    mutations = mutations.replace(
        {
            '7919980Met-2_FrTu_January_09_2020':'111e6e61e1Met-2_FrTu_January_09_2020',
            '7919980Met_FrTu_January_09_2020' : '111e6e61e1Met_FrTu_January_09_2020'
        }
    )

    # ???: what are these duplicates
    mutations = mutations \
        .sort_values(['character_label', 'sample_label']) \
        .drop_duplicates(keep="first") \
        .reset_index(drop=True) \
        .loc[:,["character_label", "sample_label", "ref", "var"]]

    samples = get_unique(mutations, "sample_label")
    mutation_id = get_unique(mutations, "character_label")

    print(f"Input has {len(samples)} samples, {len(mutation_id)} mutations.")

    # assign mutations to possible copy number states and prevalences of each CN state
    copynumbers, unmatched_mutations = cn_assignment(mutation_id, samples, cn_states)
    mutations = mutations[~mutations.character_label.isin(unmatched_mutations)] 

    # index the mutation table by sample and mutation and sort REQ: unsure why needed
    mutations.insert(loc=0, 
                     column="sample_index", 
                     value=index_column(mutations,"sample_label", sort=True))

    mutations.insert(loc=0, 
                     column="character_index", 
                     value=index_column(mutations,"character_label"))

    # reformatting
    mutations.sort_values(["sample_index","character_index"], inplace=True)
    mutations.drop(["sample_index","character_index"], axis=1, inplace=True)
    mutations = mutations.reindex(columns=['character_label', 'sample_label', 'var', 'ref'])
    mutations[["var","ref"]] = mutations[["var","ref"]].astype(int)

    write_conversion(args.writeout, mutations, copynumbers)
    

if __name__ == "__main__":
    main()