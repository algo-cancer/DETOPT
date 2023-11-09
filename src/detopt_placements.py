"""detopt_placements.py

From a fixed tree constructed using `CITUP2` using copy number neutral 
mutations, assign the placement at a node where a copy number affected 
mutation and the nodes at which the copy number state changes to a 
distinct state.
"""

import gurobipy as gp
from detopt_util import *


N_THREADS = 4 # for local version


def optimise(states: list,
             parent: dict,
             node_freq: dict,
             vafs: dict,
             seg_frac_cns: dict,
             sample_node_frac: dict,
             vaf_reg: int = 1,
             cna_reg: int = 0.25
             ) -> tuple:
    """Returns a tuple with the node to which the first instance of the 
    mutation is assigned,the nodes to which the first instance of each 
    distinct copy number (excl. 1,1 normal) have been assigned, and the 
    objective score of the tree given these assignments.

    Parameters
    ----------
    states : list
        list of the distinct copy number states to be assigned to nodes
    parent : dict
        maps a descendent (key, node) to its ancestor (value, node), root 
        has "NaN" ancestor
    node_freqs: dict
		frequency of each sample at the node; this is phi in the manuscript!
	vafs: dict
		variant allele frequencies of mutation in each sample
	frac_cns: dict
		fractional copy number of each mutation in each sample {mut1: 
        {sample1: cn1 , sample2: cn2}, mut2: {sample1: cn1 , sample2: cn2}}
    sample_node_frac:
    cna_reg
    """

    nodes = parent.keys()
    states = [(1,1)] + [(int(i[0]),int(i[1])) for i in states]

    with gp.Env(empty=True) as env:
        env.setParam('OutputFlag', 0)
        env.setParam('Threads', N_THREADS)
        env.start()

        with gp.Model("LP", env=env) as model:
            objective = 0.

            # TODO initialise delta variables
            # d_iv is 1 iff mutation i occurs for the first time at node v 
            # the root is a normal clone, so no mutations are present
            # however, additional constraints ensure the below is always true. For now, we keep it
            delta = {node: model.addVar(vtype=gp.GRB.BINARY, name=f'delta_{node}') for node in nodes}
            # delta["ROOT"] = 0 # to be removed
            model.addConstr(delta["ROOT"] == 0)
            # model.addConstr(delta["23"] == 1)
            # delta["1"] = 1 # BUG remove, used for testing # try mut at 1, cn gain on 8, on allele a
            # define constraint for delta values
            # first instance of a mutation can only be assigned to one node c1
            model.addConstr(gp.quicksum(list(delta.values())) == 1, name="c1")


            # TODO initialise theta variables; nested dictionary comprehension
            # t_iv is 1 iff copy number gain/loss to a new state occurs for the first time at node v
            # the root has a normal copy number state
            theta = {
                node: {
                    state: model.addVar(vtype=gp.GRB.BINARY, name=f"theta_{node}_{state}")
                    for state in states
                }
                for node in nodes
            }
            model.addConstr(theta["ROOT"][(1, 1)] == 1)
            #model.addConstr(theta["8"][(2, 1)] == 1)
            # theta["10"][(2, 1)] = 1 # BUG remove, used for testing
            # define constraint for theta values
            # each new copy number state only gets assigned once c2
            # each vertex can only be assigned one copy number state c3
            model.addConstrs((gp.quicksum(theta[node][state] for node in nodes) == 1 for state in states), name="c2")
            model.addConstrs((gp.quicksum(theta[node][state] for state in states) <= 1 for node in nodes), name="c3")


            # TODO initialise psi variables
            # p_iv is 1 iff copy number status of node v is Si, where each state in the state vector are the indexed set of states
            psi = {
                node: {
                    state: model.addVar(vtype=gp.GRB.BINARY, name=f"psi_{node}_{state}")
                    for state in states
                }
                for node in nodes
            }
            # define constraints for psi values
            for node in nodes:

                for state in states:
                    model.addConstr(psi[node][state] >= theta[node][state], name="c4")

                if node != "ROOT":

                    for state in states:   
                        model.addConstr(psi[node][state] >= psi[parent[node]][state] - sum(theta[node][s] for s in states), name="c5")
                        # sum -> this is on nodes at which there are no new copy number gains

            model.addConstrs((gp.quicksum(psi[node][state] for state in states) == 1 for node in nodes), name="c6")


            # TODO initialise T_av as total number of copies of allele a on node v, same for allele b
            copy_a = {node: model.addVar(vtype=gp.GRB.INTEGER, name=f"copies_A_{node}") for node in nodes}
            copy_b = {node: model.addVar(vtype=gp.GRB.INTEGER, name=f"copies_B_{node}") for node in nodes}

            # define these copy number counts
            for node in nodes:
                model.addConstr(copy_a[node] == sum([psi[node][state]*int(state[0]) for state in states])) # BUG check the definition of C
                model.addConstr(copy_b[node] == sum([psi[node][state]*int(state[1]) for state in states]))
            
            # define isa-like constraint for loh; num copies is 0 if 0 in ancestor node
            for node in nodes:
                if node != "ROOT":
                    model.addConstr(copy_a[node] <= 1000 * copy_a[parent[node]], name="c7a") # BUG have not checked this yet
                    model.addConstr(copy_b[node] <= 1000 * copy_b[parent[node]], name="c7b")
                  

            # TODO initialise e (equal), l (loss), and g (gain) variables, indicators
            e_a = {node: model.addVar(vtype=gp.GRB.BINARY, name=f"equal_A_{node}") for node in nodes if node != "ROOT"} # why 0.001 when this is binary # BUG fixed, exclude root nodes
            l_a = {node: model.addVar(vtype=gp.GRB.BINARY, name=f"loss_A_{node}") for node in nodes if node != "ROOT"}
            g_a = {node: model.addVar(vtype=gp.GRB.BINARY, name=f"gain_A_{node}") for node in nodes if node != "ROOT"}

            e_b = {node: model.addVar(vtype=gp.GRB.BINARY, name=f"equal_B_{node}") for node in nodes if node != "ROOT"}
            l_b = {node: model.addVar(vtype=gp.GRB.BINARY, name=f"loss_B_{node}") for node in nodes if node != "ROOT"}
            g_b = {node: model.addVar(vtype=gp.GRB.BINARY, name=f"gain_B_{node}") for node in nodes if node != "ROOT"}

            # define values of e, l, g;
            # e is 1 iff number of copies of allele A at node v and p(v) are equal 
            # l is 1 iff number of copies of allele A at node v is smaller than at p(v)
            # g is then, by default defined by having neither equal or loss

            # use auxillary variables to calculate difference between parent and current node cna
            # BUG most likely there is a cleaner way to implement this
            diff_par_a_aux = {node: model.addVar(vtype=gp.GRB.INTEGER, name=f"diff_A_{node}", lb=-100) for node in nodes}
            diff_par_b_aux = {node: model.addVar(vtype=gp.GRB.INTEGER, name=f"diff_B_{node}", lb=-100) for node in nodes}

            diff_child_a_aux = {node: model.addVar(vtype=gp.GRB.INTEGER, name=f"diff_cp_A_{node}", lb=-100) for node in nodes}
            diff_child_b_aux = {node: model.addVar(vtype=gp.GRB.INTEGER, name=f"diff_cp_B_{node}", lb=-100) for node in nodes}
           
            for node in nodes:
                if node != "ROOT":
                    model.addConstr(diff_par_a_aux[node] == copy_a[node] - copy_a[parent[node]], name="c8_aux_a")
                    model.addConstr(diff_par_b_aux[node] == copy_b[node] - copy_b[parent[node]], name="c8_aux_b")

                    model.addConstr(diff_child_a_aux[node] == copy_a[parent[node]] - copy_a[node])
                    model.addConstr(diff_child_b_aux[node] == copy_b[parent[node]] - copy_b[node])

            # auxillary for handling absolute values
            diff_abs_a_aux = {node: model.addVar(vtype=gp.GRB.INTEGER, name=f"diff_A_node") for node in nodes}
            diff_abs_b_aux = {node: model.addVar(vtype=gp.GRB.INTEGER, name=f"diff_B_node") for node in nodes}
            
            model.addConstrs((diff_abs_a_aux[node] == gp.abs_(diff_par_a_aux[node]) for node in nodes), name="absconstr_a")
            model.addConstrs((diff_abs_b_aux[node] == gp.abs_(diff_par_b_aux[node]) for node in nodes), name="absconstr_b")
            
            for node in nodes:
                if node != "ROOT":
                    model.addConstr(e_a[node] >= 1 - diff_abs_a_aux[node])
                    model.addConstr(e_a[node] <= 1 - 0.001 * diff_abs_a_aux[node])   
                    model.addConstr(e_b[node] >= 1 - diff_abs_b_aux[node])
                    model.addConstr(e_b[node] <= 1 - 0.001 * diff_abs_b_aux[node])

                    model.addConstr(l_a[node] >= 0.001 * diff_child_a_aux[node])
                    model.addConstr(l_a[node] <= 1 + 0.001 * diff_child_a_aux[node])
                    model.addConstr(l_a[node] <= 1 - e_a[node])
                    model.addConstr(l_b[node] >= 0.001 * diff_child_b_aux[node])
                    model.addConstr(l_b[node] <= 1 + 0.001 * diff_child_b_aux[node])
                    model.addConstr(l_b[node] <= 1 - e_b[node])

                    model.addConstr(g_a[node] == 1 - e_a[node] - l_a[node])
                    model.addConstr(g_b[node] == 1 - e_b[node] - l_b[node])

            # TODO initialise mu variables
            # a mutation i, mu_a is 1 iff the mutation i first appears in allele A, similarly for B
            mu_a = model.addVar(vtype=gp.GRB.BINARY, name="mu_A")
            mu_b = model.addVar(vtype=gp.GRB.BINARY, name="mu_B")
            # define constraint that mutation can first (and only) appear on one allele, either A or B
            model.addConstr(mu_a + mu_b == 1)
            #model.addConstr(mu_a == 1)
            #mu_a = 1 # TODO testing
            # TODO mut_a_iv is the number of copies of allele A that has mutation i at node v, similarly for B
            mut_a_iv = {node: model.addVar(vtype=gp.GRB.INTEGER, name=f"copies_mut_A_{node}") for node in nodes}
            mut_b_iv = {node: model.addVar(vtype=gp.GRB.INTEGER, name=f"copies_mut_B_{node}") for node in nodes}
            # since the root node has no mutations, the number of copies with mutation is zero
            #mut_a_iv["ROOT"] = 0
            model.addConstr(mut_a_iv["ROOT"] == 0)
            model.addConstr(mut_b_iv["ROOT"] == 0)
            #mut_a_iv["ROOT"] = 0
            #mut_b_iv["ROOT"] = 0
            
            #mut_a_iv["7"] = 2 # TODO testing
            # define constraint where the maximum number of copies of allele A with mutation is the total number of copies
            # at that node and whether the mutation has/does appear on allele A
            # additionally, the number of copies with mutation has to be at least one if the mutation first appeared on 
            # node v, it must either be on A or B
            for node in nodes:
                if node != "ROOT":
                    model.addConstr(mut_a_iv[node] <= copy_a[node] * mu_a)
                    model.addConstr(mut_a_iv[node] >= delta[node] * mu_a)

                    model.addConstr(mut_b_iv[node] <= copy_b[node] * mu_b)
                    model.addConstr(mut_b_iv[node] >= delta[node] * mu_b)

            # for a non-root node v and its parent p(v)
            # TODO case 1; copy number status is equal between v and p(v)  !!! have not checked any starting from here
            for node in nodes:
                if node != "ROOT":

                    model.addConstr((e_a[node] == 1) >> (mut_a_iv[node] >= mut_a_iv[parent[node]])) # BUG i think this works without exactly following implementation
                    model.addConstr((e_a[node] == 1) >> (mut_a_iv[node] <= mut_a_iv[parent[node]] + delta[node]))
                    model.addConstr((e_b[node] == 1) >> (mut_b_iv[node] >= mut_b_iv[parent[node]])) # BUG i think this works without exactly following implementation
                    model.addConstr((e_b[node] == 1) >> (mut_b_iv[node] <= mut_b_iv[parent[node]] + delta[node]))
            
            # TODO case 2; copy number status loss from p(v) to v
            for node in nodes:
                if node != "ROOT":

                    model.addConstr((l_a[node] == 1) >> (mut_a_iv[node] <= mut_a_iv[parent[node]] + delta[node]))
                    model.addConstr((l_b[node] == 1) >> (mut_b_iv[node] <= mut_b_iv[parent[node]] + delta[node]))
            
            
            # TODO introduce indicator variable for mutation if it is present in the parent node p(v)
            # it is 1 iff the mutation i is present in the parent node; this being 1, meaning delta_iv is 0 since it already exists
            # we initialise the variable here
            par_mut_indic = {node: model.addVar(vtype=gp.GRB.BINARY, name=f"node_{node}_has_mut") for node in nodes}
            # define constraints for I
            for node in nodes:
                if node != "ROOT":
                    model.addConstr(par_mut_indic[node] <= mut_a_iv[node] + mut_b_iv[node])
                    model.addConstr(1000 * par_mut_indic[node] >= mut_a_iv[node] + mut_b_iv[node])
            model.addConstr(par_mut_indic["ROOT"] == 0) # BUG WAS 1 UNLESS EXPLCITLY SET
            
            # add another set of auxillary variables to break up quadratic expression
            aux_contr_copy_a = {
                node: {
                    "parent": model.addVar(vtype=gp.GRB.INTEGER, name=f"contr_A_par", lb=-100), 
                    "child": model.addVar(vtype=gp.GRB.INTEGER, name=f"contr_A_child", lb=-100)
                    }
                for node in nodes}
            
            aux_contr_copy_b = {
                node: {
                    "parent": model.addVar(vtype=gp.GRB.INTEGER, name=f"contr_B_par", lb=-100), # bug fix objective function, lower bound allow negative
                    "child": model.addVar(vtype=gp.GRB.INTEGER, name=f"contr_B_child", lb=-100)
                    }
                for node in nodes}
            
            for node in nodes:
                if node != "ROOT":

                    model.addConstr(aux_contr_copy_a[node]["parent"] == diff_par_a_aux[node] * par_mut_indic[parent[node]])
                    model.addConstr(aux_contr_copy_a[node]["child"] == diff_par_a_aux[node] * delta[node])
                    model.addConstr(aux_contr_copy_b[node]["parent"] == diff_par_b_aux[node] * par_mut_indic[parent[node]])
                    model.addConstr(aux_contr_copy_b[node]["child"] == diff_par_b_aux[node] * delta[node])
            
            # TODO case 3; copy number gain from p(v) to v
            for node in nodes:
                if node != "ROOT":

                    model.addConstr((g_a[node] == 1) >> (mut_a_iv[node] >= mut_a_iv[parent[node]]))
                    model.addConstr((g_b[node] == 1) >> (mut_b_iv[node] >= mut_b_iv[parent[node]]))

                    model.addConstr((g_a[node] == 1) >> (mut_a_iv[node] <= mut_a_iv[parent[node]] + aux_contr_copy_a[node]["parent"] + aux_contr_copy_a[node]["child"]))
                    model.addConstr((g_b[node] == 1) >> (mut_b_iv[node] <= mut_b_iv[parent[node]] + aux_contr_copy_b[node]["parent"] + aux_contr_copy_b[node]["child"]))

            
            # TODO add indicator variable u_aiv; understand this
            # u_aiv is 1 iff at a node v if there is a nonzero number of copies of allele A and all these harbour a mutation i, similarly for B
            mut_on_allcopies_a = {node: model.addVar(vtype=gp.GRB.BINARY, name=f"all_copies_A_mut_{node}") for node in nodes}
            mut_on_allcopies_b = {node: model.addVar(vtype=gp.GRB.BINARY, name=f"all_copies_B_mut_{node}") for node in nodes}            
            # define constraints for u_aiv            
            for node in nodes:
                if node != "ROOT":

                    model.addConstr(1000 * (1 - mut_on_allcopies_a[node]) >= copy_a[node] - mut_a_iv[node])
                    model.addConstr(mut_on_allcopies_a[node] <= mut_a_iv[node])
                    model.addConstr(mut_on_allcopies_a[node] >= mut_a_iv[node] - copy_a[node] + par_mut_indic[node] - mu_b)
                    model.addConstr(mut_a_iv[node] + 1000 * (1 - g_a[node]) >= mut_a_iv[parent[node]] + diff_par_a_aux[node] * mut_on_allcopies_a[parent[node]]) # idk what this is for ?

                    model.addConstr(1000 * (1 - mut_on_allcopies_b[node]) >= copy_b[node] - mut_b_iv[node])
                    model.addConstr(mut_on_allcopies_b[node] <= mut_b_iv[node])
                    model.addConstr(mut_on_allcopies_b[node] >= mut_b_iv[node] - copy_b[node] + par_mut_indic[node] - mu_a)
                    model.addConstr(mut_b_iv[node] + 1000 * (1 - g_b[node]) >= mut_b_iv[parent[node]] + diff_par_b_aux[node] * mut_on_allcopies_b[parent[node]]) # idk what this is for ?

            # BUG test constraints up to here

            # Additional constraints Oct. 5, 2023
            # TODO Constrain variable delta_iv such that it satisfies presence-absence constraints # see file /Users/wuchh/Downloads/20231004_164106.jpg
            for v in sample_node_frac.keys():
                if v != "ROOT":
                    #print(v)
                    for s in sample_node_frac[node].keys():
                        pres_abs_node = int(sample_node_frac[v][s] > 0.03)
                        pres_abs_sample = int(vafs[s] >= 0.05)
                        pres_abs_node_only = int(node_freq[v][s] > 0.03) # TODO flexible to get more assignments; hard to satisfy, especially for very clonal mutations

                        #print(f"\t{pres_abs_node}({sample_node_frac[v][s]})\t{pres_abs_sample}({round(vafs[s], 3)})\t{pres_abs_node >= pres_abs_sample}")
                        model.addConstr(pres_abs_node >= delta[v] * pres_abs_sample, name="presence_absence_c1")

                        #print(f"\t{s}\t{pres_abs_sample}({round(vafs[s], 3)})\t{pres_abs_node_only}({node_freq[v][s]})\t{pres_abs_sample >= pres_abs_node_only}")
                        model.addConstr(pres_abs_sample >= par_mut_indic[v] * pres_abs_node_only, name="presence_absence_c2")

            aux_vaf = {}
            aux_cna = {}
            for sample in vafs.keys():
                # define variant genomic content
                variant_genomic_content = gp.quicksum([(mut_a_iv[node] + mut_b_iv[node]) * node_freq[node][sample] for node in nodes])
                # define total genomic content
                total_genomic_content = gp.quicksum([(copy_a[node] + copy_b[node]) * node_freq[node][sample] for node in nodes])

                # define auxiliary variables to allow absolute value in objective
                aux_vaf[sample] = model.addVar(vtype=gp.GRB.CONTINUOUS, name=f'aux_vaf_{sample}')
                aux_cna[sample] = model.addVar(vtype=gp.GRB.CONTINUOUS, name=f'aux_cna_{sample}')

                # objective terms
                discrepancy_vaf = (variant_genomic_content / seg_frac_cns[sample]) - vafs[sample] 
                discrepancy_cna = total_genomic_content - seg_frac_cns[sample]

                # enforce an absolute value objective
                model.addConstr(discrepancy_vaf <= aux_vaf[sample])
                model.addConstr(-1 * discrepancy_vaf <= aux_vaf[sample])
                model.addConstr(discrepancy_cna <= aux_cna[sample])
                model.addConstr(-1 * discrepancy_cna <= aux_cna[sample])

                objective += (vaf_reg * aux_vaf[sample] + cna_reg * aux_cna[sample])

            model.setObjective(objective, gp.GRB.MINIMIZE)
            model.optimize()

            try:
                obj_val = model.objVal

                mut_a, mut_b, tot_a, tot_b = {}, {}, {}, {}
                cna_assignment = {}
                for node in nodes:
                    if round(delta[node].x, ndigits=None) == 1: 
                        #print(f"SNV {node}")
                        node_assignment = node
                    
                    for state in states:
                        
                        if round(theta[node][state].x, ndigits=None) == 1 and state != (1, 1): 
                            #print(f"CNA {node}, {state}")
                            cna_assignment[state] = node

                    mut_a[node], mut_b[node] = mut_a_iv[node].x, mut_b_iv[node].x
                    tot_a[node], tot_b[node] = copy_a[node].x, copy_b[node].x

                model_vars = node_assignment, cna_assignment, mut_a, mut_b, tot_a, tot_b 

                return obj_val, model_vars

            except AttributeError:
                return -1, -1       # TODO not a good way to exit, will be fixed