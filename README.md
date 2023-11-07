# DETOPT
This is the repository for `DETOPT` a combinatorial optimization method for **DET**ermining **O**ptimal **P**lacement in **T**umor progression history of single nucleotide variants (SNVs) from the genomic regions impacted by copy-number aberrations (CNAs) using multi-sample bulk DNA sequencing data.

# Contents

  1. [Installation](#install)
  2. [Usage](#usage)
     * [required input data](#input)
     * [output](#output)
  4. [Demos](#demos)
     * [**summary information**](real_data/demo/README.md#Multi-region-longitudinal-bulk-DNA-sequenced-samples-from-metastatic-breast-cancer-(mBrCa)-patient)
     * [**longitudinal consistency analysis**](real_data/demo/README.md#Longitudinal-consistency-of-base-tree-topology)
     * [**variant placements**](real_data/demo/README.md#Output-of-DETOPT-on-patient-4355-data)
  6. [Support](#support)
    
<a name="install"></a>
# Installation
Clone the repository locally 
```console
$ git clone https://github.com/algo-cancer/DETOPT.git
$ cd DETOPT
```

**DEPENDENCIES** </br>
* Python >= 3.10
* Required packages are managed by conda. If it has not been previously installed, follow installation directions for [`conda`](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Create the conda environment from the `.yaml` file, then activate it using the following commands.
```console
$ conda env create -f env/detopt.yaml
$ conda activate detopt
```

> [!NOTE]
>`DETOPT` requires an installation and a valid license for [`Gurobi`](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-) (freely available for academic users). Follow instructions for retrival and setting up of a license [here](https://support.gurobi.com/hc/en-us/articles/12872879801105).

<a name="usage"></a>
# Usage
```
(detopt) foo@bar:$ python detopt.py --help
usage: detopt.py [--help] -s SNV_FILE -t TREE_FILE -h SAMPLES -p CNA_WEIGHT -d DATA_DIR -o OUT

DETOPT (DETermining Optimal Placement in Tumor progression history)

options:
  --help                                    show this help message and exit
  -s SNV_FILE,    --snv_file    SNV_FILE    *.input.tsv file containing information about read counts and allele-copy number calls of SNVs
  -t TREE_FILE,   --tree_file   TREE_FILE   *.tree file containing information for each subclone, the PARENT_ID, MUTATIONS_AT_NODE, SAMPLE_IDS, NODE_FREQUENCIES
  -h SAMPLES,     --samples     SAMPLES     number of samples
  -p CNA_WEIGHT,  --cna_weight  CNA_WEIGHT  regularization value/weight
  -d DATA_DIR,    --data_dir    DATA_DIR    directory containing required decifer and tree files
  -o OUT,         --out         OUT         output filename, optionally with filepath
```

<a name="input"></a>
**PARAMETERS** </br> 
  **Flag**|**Argument**|**Symbol Used**|**Description**
  --------|------------------|---------------|---------------
  -s|`snv_file`| - | refer to extended description of required input files 
  -t|`tree_file`| - | refer to extended description of required input files
  -p|`cna_weight`| *&rho;* |regularization parameter, denoted as in the manuscript, introduced to balance the two objective terms. By default, *&rho;* is set to 0.25.
  -h|`num_samples`| *h* | the total number of samples from the multi-sample bulk DNA sequencing data 
  -d|`data_dir`| - | a relative path to the directory in which the input `.decifer.tsv` and `.tree` files are located; provide the path to the directory without the trailing '/'
  -o|`out`| - | the filename for which `DETOPT` will create a file containing the output in the current working directory, e.g., `<out>.detopt.tsv`. By prepending filename with a filepath, the file can be writting to another directory, e.g., `</path/to/out>.detopt.tsv`
  
<a name="input"></a>
<details>
  <summary> <b> extended description of required input files </b> </summary> </br>

  `DETOPT` requires two (2) input files. The first contains the read counts and copy-number states (as obtained from HATCHet[^2], or any allele- and clone-specific copy-number caller) for each mutation, including both copy-number neutral and aberrant SNVs, in each sample. The second describes the base tree topology, constructed by the use of only copy-number neutral SNVs, with inferred sample node (subclone) frequencies.
  
**SNV file**. `--snv_file <SNV_FILE>` This file contains information for each SNV, the read counts mapping to the variant and reference alleles and the allele- and clone-specific copy number calls. `hatchetconvert.py` is provided for the conversion of a `.seg.ucn` file from HATCHet ([example](real_data/demo/demo_data/best.seg.ucn)) and a mutation tab-separated `.tsv` file ([example](real_data/demo/demo_data/111e6e61e1_all_mutations.tsv)) into an input file for DETOPT with fields,

```
mut_index    sample    var_reads    ref_reads    normal_state    normal_prop    tumor1_state    tumor1_prop    tumor2_state    tumor2_prop
mut_1        sample_1  42           123          1|1             0.23           2|1             0.54           2|0             0.23  
```

* `mut_index`: unique mutation identifier 
* `sample`: unique sample identifier
* `var_reads`: number of reads mapping to the variant allele
* `ref_reads`: number of reads mapping to the reference allele
* `normal_state`: normal copy number state, '1|1'
* `normal_prop`: proportion of cells in the sample that have the `normal_state` copy number state
* `tumor_state`: aberrant copy number state, 'A|B' where A and B are number of copies of allele A and B, respectively. Values of A and B are not allele-specific.
* `tumor_prop`: proportion of cells in the sample that have the `tumor_state` copy number state

**Tree file**. `--tree_file <TREE_FILE>` This file contains information for each subclone, represented by a node in the base tree, the information for the following fields,
  
```
NODE_ID    PARENT_ID    MUTATIONS_AT_NODE    SAMPLE_IDS                 NODE_FREQUENCIES
0          1            mut_1, ..., mut_i    sample_1, ..., sample_j    f_1, ..., f_j 
```
     
* `NODE_ID`: `node` in the tree, representing a subclone
* `PARENT_ID`: parental node (also a subclone) of `node`
* `MUTATIONS_AT_NODE`: copy number neutral mutations assigned to `node`
* `SAMPLE_IDS`: list of samples in which `node` (subclone) is present
* `NODE_FREQUENCIES`: list of the inferred sample node (subclone) frequencies. For each sample in `SAMPLE_IDS`, the node (subclone) frequency is the fraction of all cells in a sample, including normal cells, that belong to that subclone.
</details>

<a name="output"></a>
**OUTPUTS**  </br>
> [!IMPORTANT]
>`DETOPT` returns a tab-separated `.detopt.tsv` file containing the inferred placements of copy number gain and loss events. If an SNV was impacted by more than one copy-number aberration event, each copy-number state and its inferred placement are designated by a unique index that has no particular signficance other than to distinguish these events.

<a name="demos"></a>
# Demos

**Example**|**Description**|**output**
-----------|---------------|-----------
[4355](real_data/demo/README.md) |Demo of `DETOPT` on metastatic breast cancer patient 4355[^1] with 18 samples|[here](real_data/demo/README.md#Output-of-DETOPT-on-patient-4355-data)

```console
(detopt) $ python src/detopt.py -s 4355.snvs.input -t 4355.tree -h 18 -d real_data/demo/demo_inputs -o real_data/demo/results
```

<a name="support"></a>
# Support
The software and corresponding documentations are maintained by the research group of Dr. S. Cenk Sahinalp. If you have encountered issues with the tool or have inquiries, please raise it on the [issue forum](https://github.com/algo-cancer/DETOPT/issues) or contact [Chih Hao Wu](mailto:chih.wu@nih.gov) and [Salem Malikić](mailto:salem.malikic@nih.gov). Please always refer to the GitHub repository of [DETOPT](https://github.com/algo-cancer/DETOPT) for the most updated version of the software and relevant documentation.

<!-- References -->
[^1]: Zacharakis, N., Huq, L.M., Seitter, S.J., Kim, S.P., Gartner, J.J., Sindiri, S., Hill, V.K., Li, Y.F., Paria, B.C., Ray, S., et al.: Breast cancers are immunogenic: immunologic analyses and a phase ii pilot clinical trial using mutation-reactive autologous lymphocytes. Journal of Clinical Oncology 40(16), 1741–1754 (2022) [https://doi.org/10.1200/jco.21.02170](https://doi.org/10.1200/jco.21.02170)
[^2]: https://github.com/raphael-group/hatchet
