"""detopt_placements.py

From a fixed tree constructed using `CITUP2` using copy number neutral 
mutations, assign the placement at a node where a copy number affected 
mutation and the node(s) at which the copy number state changes to a 
distinct state.

ChangeLog:
    11/21/2023: (Lines 536-537) Round gp.GRB.INTEGER values to nearest
    integer values
"""

import gurobipy as gp


ALLELES = ['A','B']


def optimise(states: list,
             parent: dict,
             node_freq: dict,
             vafs: dict,
             seg_frac_cns: dict,
             sample_node_frac: dict,
             vaf_reg: int=1,
             cna_reg: int=0.25,
             n_threads: int=4
             ) -> tuple:
    """<main>

    <desc>

    Args:
        states (list): A list of all states (str, str) for the copy number 
        states of the region on which `mut` is present
	    parent (dict): A mappping from each node (subclone) to its uniquely 
        identified parental node (subclone). The root has NaN ancestor, 
        for sake of non-empty parental node.
		node_freqs (dict): The sample node frequencies for each sample, 
        in each node
        vafs (dict): The variant allele frequency (VAF) of `mut` in each sample
	    seg_frac_cns (dict): The fractional copy number state of each segment
        on which a mutation is on, for each sample
        n_threads (int): The number of threads used by Gurobi (or eq. solver)
        sample_node_frac (dict): The fraction of all cells in a sample, 
        including normal cells, that belong to a node (subclone)
        vaf_reg (int): A regularisation term for weighting discrepancy of VAFs
        cna_reg (int): A regularisation term for weighting discrepancy of CNAs

    Returns:
        tuple:
            objective (float): The objective score returned by `DETOPT`, 
            refer to the objective function in the manuscript
            node_assigned (int): The node (subclone) to which the mutational
            event was assigned for the `mut`
            state_assignment (int): The node (subclone) to which the copy
            number change event(s) was (were) assigned for the `mut`
    """

    nodes = parent.keys()
    states = [(1, 1)] + [(int(i[0]), int(i[1])) for i in states]

    env_params = {'OutputFlag': 0,
                 'Threads': n_threads
                 }

    with gp.Env("gurobi.log", params=env_params) as env, gp.Model(env=env) as model:
        objective = 0.

        # VAR. delta_{iv}

        delta = {
            node: model.addVar(vtype=gp.GRB.BINARY, name=f'delta_i_{node}') 
            for node in nodes
            }

        # CONSTRS. delta_{iv}

        model.addConstr(delta['ROOT'] == 0)
        model.addConstr(
            gp.quicksum(list(delta.values())) == 1, 
            name='constraint 1'
            )

        # VAR. theta_{iv}

        theta = {
            state: {
                node: model.addVar(vtype=gp.GRB.BINARY, name=f'theta_{state}_{node}')
                for node in nodes
                } 
            for state in states
            }
        
        # CONSTRS. theta_{iv}

        model.addConstr(theta[(1, 1)]['ROOT'] == 1)
        model.addConstrs(
            (gp.quicksum(theta[state][node] for node in nodes) == 1 for state in states), 
            name='constraint 2'
            )
        model.addConstrs(
            (gp.quicksum(theta[state][node] for state in states) <= 1 for node in nodes), 
            name='constraint 3'
            )

        # VAR. psi_{iv}

        psi = {
            state:{
                node: model.addVar(vtype=gp.GRB.BINARY, name=f'psi_{state}_{node}')
                for node in nodes
                }
            for state in states
            }

        # CONSTRS. psi_{iv}

        for node in nodes:
            theta_sum = gp.quicksum(theta[state][node] for state in states)

            for state in states:
                model.addConstr(
                    psi[state][node] >= theta[state][node], 
                    name='constraint 4'
                    )

                if node != 'ROOT':
                    model.addConstr(
                        psi[state][node] >= psi[state][parent[node]]-theta_sum, 
                        name='constraint 5'
                        )

        model.addConstrs(
            (gp.quicksum(psi[state][node] for state in states) == 1 for node in nodes), 
            name='constraint 6'
            )

        # VAR. T_{*v}

        tot_copies = {
            allele : {
                node : model.addVar(vtype=gp.GRB.INTEGER, name=f'total_copies_{node}')
                for node in nodes
                }
            for allele in ['A', 'B']
            }
        
        # CONSTRS. T_{*v} 

        for node in nodes:
            for i, allele in enumerate(ALLELES):
                model.addConstr(
                    tot_copies[allele][node] == \
                        gp.quicksum(psi[state][node]*state[i] for state in states), 
                    name=f'constraint 7'
                    )
    
        for node in nodes:
            if node != 'ROOT':
                model.addConstrs(
                    (tot_copies[allele][node] <= \
                     1000*tot_copies[allele][parent[node]] for allele in ALLELES), 
                    name='constraint 8'
                    )

        # VAR. e_{*v}, l_{*v}, g_{*v}

        e = {
            node: {
                allele : model.addVar(vtype=gp.GRB.BINARY, name=f'e_{node}_{allele}')      
                for allele in ALLELES
                } 
            for node in nodes if node !='ROOT' 
            }
        
        l = {
            node: {
                allele : model.addVar(vtype=gp.GRB.BINARY, name=f'l_{node}_{allele}')
                for allele in ALLELES
                }
            for node in nodes if node != 'ROOT'
            }

        g = {
            node: {
                allele : model.addVar(vtype=gp.GRB.BINARY, name=f'g_{node}_{allele}')
                for allele in ALLELES
                }
            for node in nodes if node != 'ROOT'
            }

        # create auxillary variables for the differences in copy number state of an 
        # allele from the parent vertex to the child vertex, T_{*p} - T_{*v} 

        aux_diff_par_chd = {
            allele : {
                node: model.addVar(vtype=gp.GRB.INTEGER, name=f'diff_{allele}_{node}]', lb=-100)
                for node in nodes
                } 
            for allele in ALLELES
            }
        
        # create auxillary variables for the differences in copy number state of an 
        # allele from the child vertex to the parent vertex, T_{*v} - T_{*p} 

        aux_diff_chd_par = {
            allele : {
                node: model.addVar(vtype=gp.GRB.INTEGER, name=f'diff_{allele}_{node}]', lb=-100)
                for node in nodes
                } 
            for allele in ALLELES
            }
        
        # create auxillary variables for taking absolute value of the differences
        # for each allele, |T_{*v} - T_{*p}|

        aux_abs_diff_chd_par = {
            allele : {
                node: model.addVar(vtype=gp.GRB.INTEGER, name=f'abs_diff_{allele}_{node}]')
                for node in nodes
                } 
            for allele in ALLELES
            }
        
        for node in nodes:
            for allele in ALLELES:

                if node != 'ROOT':
                    v, p = node, parent[node]

                    # AUX CONSTRS.

                    model.addConstr(
                        aux_diff_par_chd[allele][v] == \
                            tot_copies[allele][p]-tot_copies[allele][v],
                        name=f'auxillary constraint 1'
                        ) 
                    model.addConstr(
                        aux_diff_chd_par[allele][v] == \
                            tot_copies[allele][v]-tot_copies[allele][p],
                        name=f'auxillary constraint 2'
                        )
                    model.addConstr(
                        aux_abs_diff_chd_par[allele][v] == \
                            gp.abs_(aux_diff_chd_par[allele][v]),
                        name='auxillary constraint 3'
                        )
                    
                    # CONSTRS. e_{*v}, l_{*v}, g_{*v}

                    model.addConstr(
                        e[v][allele] >= 1-aux_abs_diff_chd_par[allele][v],
                        name='constraint 9'
                        )
                    model.addConstr(
                        e[v][allele] <= 1-0.001*aux_abs_diff_chd_par[allele][v], 
                        name='constraint 10'
                        )   
                    
                    model.addConstr(
                        l[v][allele] >= 0.001*aux_diff_par_chd[allele][v], 
                        name='constraint 11'
                        )
                    model.addConstr(
                        l[v][allele] <= 1+0.001*aux_diff_par_chd[allele][v], 
                        name='constraint 12'
                        )
                    model.addConstr(
                        l[v][allele] <= 1-e[v][allele],
                        name='constraint 13'
                        )
                    
                    model.addConstr(
                        g[v][allele] == 1-e[v][allele]-l[v][allele], 
                        name='constraint 14'
                        )

        # VAR. mu_{*}

        mu = {
            allele: model.addVar(vtype=gp.GRB.BINARY, name=f'mu_{allele}')
            for allele in ALLELES
            }
        
        # VAR. *_{v}

        mut_iv = {
            allele: {
                node: model.addVar(vtype=gp.GRB.INTEGER, name=f'mu_{allele}_{node}]')
                for node in nodes
                }
            for allele in ALLELES
            }
        
        # CONSTRS. mu_{*}

        model.addConstr(
            mu['A'] + mu['B'] == 1, 
            name='constraint 15'
            )
        
        # CONSTRS. *_{v}

        model.addConstrs(
            (mut_iv[allele]['ROOT'] == 0 for allele in ALLELES), 
            name='constraint 16'
            )
     
        for node in nodes:
            if node != 'ROOT':

                for allele in ALLELES:
                    model.addConstr(
                        mut_iv[allele][node] <= tot_copies[allele][node]*mu[allele], 
                        name='constraint 17'
                        )
                    model.addConstr(
                        mut_iv[allele][node] >= delta[node]*mu[allele],
                        name='constraint 18'
                        )

        """Case 1
        """

        for node in nodes:
            if node != 'ROOT':
                v, p = node, parent[node]

                for allele in ALLELES:
                    model.addConstr(
                        (e[node][allele] == 1) >> \
                            (mut_iv[allele][v] >= mut_iv[allele][p]),
                        name='constraint 19'
                        )
                    model.addConstr(
                        (e[node][allele] == 1) >> \
                            (mut_iv[allele][v] <= mut_iv[allele][p]+delta[node]),
                        name='constraint 20'
                        )

        """Case 2
        """       

        for node in nodes:
            if node != 'ROOT':
                v, p = node, parent[node]

                for allele in ALLELES:
                    model.addConstr((l[node][allele] == 1) >> \
                                    (mut_iv[allele][v] <= mut_iv[allele][p]+delta[node]))
        
        """Case 3
        """

        # VAR. I_{v}

        mut_indic_iv = {node: model.addVar(vtype=gp.GRB.BINARY, name=f"node_{node}_has_mut") for node in nodes}
        
        # CONSTRS. I_{v}

        for node in nodes:
            if node != 'ROOT':

                model.addConstr(
                    mut_indic_iv[node] <= mut_iv['A'][node]+mut_iv['B'][node],
                    name='constraint 21'
                    )
                model.addConstr(
                    1000*mut_indic_iv[node] >= mut_iv['A'][node]+mut_iv['B'][node],
                    name='constraint 22'
                    )
                
        model.addConstr(mut_indic_iv["ROOT"] == 0) 
        
        # auxillary variables to break up quadratic expression

        aux_par_mut = {
            allele: {
                node: model.addVar(
                    vtype=gp.GRB.INTEGER, 
                    lb=-100, # TODO remove?
                    name=f'par_{node}_mut_{allele}'
                    )                   
                for node in nodes
                }
            for allele in ALLELES
            }
        
        aux_chd_mut = {
            allele: {
                node: model.addVar(
                    vtype=gp.GRB.INTEGER, 
                    lb=-100, # TODO remove?
                    name=f'{node}_mut_{allele}'
                    )                   
                for node in nodes
                }
            for allele in ALLELES
            }
        
        for node in nodes:
            if node != "ROOT":
                v, p = node, parent[node]

                for allele in ALLELES:
                    model.addConstr(
                        aux_par_mut[allele][v] == \
                            aux_diff_chd_par[allele][v]*mut_indic_iv[p],
                        name='constraint 23'
                        )
                    model.addConstr(
                        aux_chd_mut[allele][v] == \
                            aux_diff_chd_par[allele][v]*delta[v],
                        name='constraint 24'
                        )
                    model.addConstr(
                        (g[v][allele] == 1) >> \
                            (mut_iv[allele][v] >= mut_iv[allele][p]),
                        name='constraint 25'
                        )
                    model.addConstr(
                        (g[v][allele] == 1) >> \
                            (mut_iv[allele][v] <= mut_iv[allele][p]+aux_par_mut[allele][v]+aux_chd_mut[allele][v]),
                        name='constraint 26'
                        )

        # VAR. U_{*v}

        mut_indic_all_iv = {
            allele: {
                node: model.addVar(vtype=gp.GRB.BINARY, name=f'node_{node}_allele_{allele}_has_mut')
                for node in nodes
                }
            for allele in ALLELES
            }

        # CONSTRS. U_{*v}

        for node in nodes:
            if node != "ROOT":
                v, p = node, parent[node]

                for allele in ALLELES:
                    match allele:
                        case 'A':
                            other_allele = 'B'
                        case 'B':
                            other_allele = 'A'

                    model.addConstr(
                        1000*(1-mut_indic_all_iv[allele][v]) >= \
                            tot_copies[allele][v]-mut_iv[allele][v],
                        name='constraint 27'
                        )
                    model.addConstr(
                        mut_indic_all_iv[allele][v] <= mut_iv[allele][v],
                        name='constraint 28'
                    )
                    model.addConstr(
                        mut_indic_all_iv[allele][v] >= \
                            mut_iv[allele][v]-tot_copies[allele][v]+mut_indic_iv[v]-mu[other_allele],
                        name='constraint 29'
                    )
                    model.addConstr(
                        mut_iv[allele][v]+1000*(1-g[v][allele]) >= \
                            mut_iv[allele][p]+aux_diff_chd_par[allele][v]*mut_indic_all_iv[allele][v],
                        name='constraint 30'
                    )

        # Constrain delta_iv such that it satisfies presence-absence constraints 
        # see file /Users/wuchh/Downloads/20231004_164106.jpg

        for node in sample_node_frac.keys():
            if node != "ROOT":
                for s in sample_node_frac[node].keys():
                    pres_abs_node = int(sample_node_frac[node][s] > 0.03)
                    pres_abs_sample = int(vafs[s] >= 0.05)
                    pres_abs_node_only = int(node_freq[node][s] > 0.03) 
                    # TODO flexible to get more assignments; hard to satisfy, 
                    # especially for very clonal mutations

                    model.addConstr(pres_abs_node >= \
                                    delta[node]*pres_abs_sample, 
                                    name="constraint 31"
                                    )
                    model.addConstr(pres_abs_sample >= \
                                    mut_indic_iv[node]*pres_abs_node_only,
                                    name="constraint 32"
                                    )

        aux_vaf = {}
        aux_cna = {}

        for sample in vafs.keys():
            # define variant genomic content
            variant_genomic_content = \
                gp.quicksum([(mut_iv['A'][node]+mut_iv['B'][node])*node_freq[node][sample] for node in nodes])
            
            # define total genomic content
            total_genomic_content = \
                gp.quicksum([(tot_copies['A'][node]+tot_copies['B'][node])*node_freq[node][sample] for node in nodes])

            # define auxiliary variables to allow absolute value in objective
            aux_vaf[sample] = model.addVar(vtype=gp.GRB.CONTINUOUS, name=f'aux_vaf_{sample}')
            aux_cna[sample] = model.addVar(vtype=gp.GRB.CONTINUOUS, name=f'aux_cna_{sample}')

            # objective terms
            discrepancy_vaf = (variant_genomic_content/seg_frac_cns[sample])-vafs[sample] 
            discrepancy_cna = total_genomic_content-seg_frac_cns[sample]

            # absolute value of objective terms
            model.addConstr(discrepancy_vaf <= aux_vaf[sample])
            model.addConstr(-1*discrepancy_vaf <= aux_vaf[sample])
            model.addConstr(discrepancy_cna <= aux_cna[sample])
            model.addConstr(-1*discrepancy_cna <= aux_cna[sample])

            objective += (vaf_reg*aux_vaf[sample]+cna_reg*aux_cna[sample])

        model.setObjective(objective, gp.GRB.MINIMIZE)
        model.optimize()

        try:
            obj_val = model.objVal

            mut_a, mut_b, tot_a, tot_b = {}, {}, {}, {}
            cna_assignment = {}

            for node in nodes:
                if round(delta[node].x, ndigits=None) == 1: 
                    node_assignment = node
                
                for state in states:
                    if round(theta[state][node].x, ndigits=None) == 1 and state != (1, 1): 
                        cna_assignment[state] = node

                mut_a[node], mut_b[node] = round(mut_iv['A'][node].x, ndigits=None), round(mut_iv['B'][node].x, ndigits=None)
                tot_a[node], tot_b[node] = round(tot_copies['A'][node].x, ndigits=None), round(tot_copies['B'][node].x, ndigits=None)

            model_vars = node_assignment, cna_assignment, mut_a, mut_b, tot_a, tot_b 

            return obj_val, model_vars

        except AttributeError:
            return -1, -1       # TODO not a good way to exit, will be fixed
        