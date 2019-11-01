# Created by woochanghwang at 2019-03-28
#######################################################################
# Recipies / wrapper functions using toolbox methods for disease,
# drug and network analysis from Emre
# e.g. 03/2019
#######################################################################
import toolbox.network_utilities as nu
import toolbox.mapping_utilities as mu
import csv, numpy, os, pickle
import random
import pandas as pd


def get_gene_list_from_file(gene_list_file):
    gene_df = pd.read_table(gene_list_file, sep='\t')
    gene_list = list(gene_df['Gene'])

    return gene_list

def main():
    print("LifeArc wrapper")
    network_file = "/Users/woochanghwang/PycharmProjects/LifeArc/General/src_drug/Data/human_protein_interactome.sif"
    disease_gene_file = "/Users/woochanghwang/PycharmProjects/LifeArc/General/src_drug/Data/disease_genes.tsv"
    drug_target_file = "/Users/woochanghwang/PycharmProjects/LifeArc/General/src_drug/Data/drug_target_interactions.txt"

    network = nu.create_network_from_sif_file(network_file, use_edge_data=False, delim=None,
                                                        include_unconnected=True)
    nodes = set(network.nodes())
    print("network lengths",len(nodes))
    drug_to_targets = get_drug_target_drugbank(drug_target_file, nodes=nodes)
    print(drug_to_targets)
    # disease_to_genes, disease_to_category = get_diseasome_genes(disease_gene_file, nodes=nodes)
    # gene_list_file =  "/Users/woochanghwang/PycharmProjects/LifeArc/ULK/result/GBM_ULK1_gene_score_by_RW_pvalue_FC_230119.tsv"

    disease_name = 'GBM'
    # disease_to_genes, disease_to_category = get_diseasome_genes_from_selectedGenes(gene_list_file, disease_name,
    #                                                                                disease_category=None, nodes=nodes)
    #
    # print("network edges:", network.edges())
    # output_file = "{}_drug_proximity.tsv".format(disease_name)
    # calculate_proximity_multiple(network, from_file=drug_target_file, to_file=gene_list_file, disease_mode=disease_name,
    #                              out_file=output_file)

    #######################################
    ## Temp for ULK1,ULK2
    #######################################
    gene_list_file_ulk1 = "/Users/woochanghwang/PycharmProjects/LifeArc/ULK/result/GBM_ULK1_gene_score_by_RW_pvalue_FC_230119.tsv"
    gene_list_file_ulk2 = "/Users/woochanghwang/PycharmProjects/LifeArc/ULK/result/GBM_ULK2_gene_score_by_RW_pvalue_FC_230119.tsv"
    gene_list_ulk1 = get_gene_list_from_file("/Users/woochanghwang/PycharmProjects/LifeArc/ULK/result/GBM_ULK1_gene_score_by_RW_pvalue_FC_230119.tsv")
    gene_list_ulk2 = get_gene_list_from_file("/Users/woochanghwang/PycharmProjects/LifeArc/ULK/result/GBM_ULK2_gene_score_by_RW_pvalue_FC_230119.tsv")

    gene_list = list(set(gene_list_ulk1).union(set(gene_list_ulk2)))
    disease_name = 'GBM'
    disease_to_genes , disease_to_category = get_diseasome_genes_from_selectedGenes(gene_list, disease_name, disease_category=None, nodes=nodes)
    print("disease_gene", disease_to_genes)

    print("network edges:", network.edges())
    output_file = "/Users/woochanghwang/PycharmProjects/LifeArc/ULK/result/drug/{}_drug_proximity_{}.tsv".format(disease_name,"ULK1_2")
    calculate_proximity_multiple(network,from_file=drug_target_file, to_file=gene_list ,disease_mode = disease_name, out_file=output_file)

    ##############################
    ## test for disease_gene



def calculate_closest_distance(network, nodes_from, nodes_to, lengths=None):
    values_outer = []
    if lengths is None:
        for node_from in nodes_from:
            values = []
            for node_to in nodes_to:
                # print("from - to", node_from, node_to)
                if not nu.check_has_path(network,node_from,node_to): continue
                val = nu.get_shortest_path_length_between(network, node_from, node_to)
                values.append(val)
            if len(values) == 0:    continue
            d = min(values)
            # print (d)
            values_outer.append(d)
    else:
        for node_from in nodes_from:
            values = []
            vals = lengths[node_from]
            for node_to in nodes_to:
                val = vals[node_to]
                values.append(val)
            d = min(values)
            values_outer.append(d)
    d = numpy.mean(values_outer)
    # print d
    return d


##### Proximity related #####

def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None,
                        n_random=1000, min_bin_size=100, seed=452456, lengths=None):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """

    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network
    nodes_to = set(nodes_to) & nodes_network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None  # At least one of the node group not in network
    d = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = nu.get_degree_binning(network, min_bin_size,
                                                    lengths)  # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bins=bins, n_random=n_random,
                                             min_bin_size=min_bin_size, seed=seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bins=bins, n_random=n_random, min_bin_size=min_bin_size,
                                           seed=seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(len(nodes_from_random))  # n_random
    for i, values_random in enumerate(random_values_list):
        nodes_from, nodes_to = values_random
        # values[i] = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
        values[i] = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    # pval = float(sum(values <= d)) / len(values) # needs high number of n_random
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, (m, s)  # (z, pval)


def calculate_proximity_multiple(network, from_file=None, to_file=None, disease_mode = 'whole', n_random=1000, min_bin_size=100, seed=452456,
                                 lengths=None, out_file="output.txt"):
    """
    Run proximity on each entries of from and to files in a pairwise manner
    output is saved in out_file (e.g., output.txt)
    disease_mode : whole - from disease_gene.tsv
    disease_mode : disease_name - from our analysis
    """
    nodes = set(network.nodes())
    drug_to_targets = get_drug_target_drugbank(from_file, nodes=nodes)
    # drug_to_targets = dict((drug, nodes & targets) for drug, targets in drug_to_targets.iteritems())
    if disease_mode == 'whole': # all disease
        disease_to_genes, disease_to_category = get_diseasome_genes(to_file, nodes=nodes)
    else :
        disease_to_genes, disease_to_category = get_diseasome_genes_from_selectedGenes(to_file,disease_mode, disease_category=None, nodes=nodes)

    # Calculate proximity values
    print(len(drug_to_targets), len(disease_to_genes))
    # Get degree binning
    bins = nu.get_degree_binning(network, min_bin_size)
    f = open(out_file, 'w')
    f.write("source\ttarget\tn.source\tn.target\td\tz\n")
    for drug, nodes_from in drug_to_targets.items():
        values = []
        for disease, nodes_to in disease_to_genes.items():
            print(drug, disease)
            d, z, (m, s) = calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None,
                                               nodes_to_random=None, bins=bins, n_random=n_random,
                                               min_bin_size=min_bin_size, seed=seed, lengths=lengths)
            values.append((drug, disease, z, len(nodes_from), len(nodes_to), d, m, s))
        # f.write("%s\t%s\t%f\t%f\t%f\t%f\n" % (drug, disease, z, d, m, s))
        values.sort(key=lambda x: x[2])
        for drug, disease, z, k, l, d, m, s in values:
            # f.write("%s\t%s\t%f\t%d\t%d\t%f\t%f\t%f\n" % (drug, disease, z, k, l, d, m, s))
            f.write("%s\t%s\t%d\t%d\t%f\t%f\n" % (drug, disease, k, l, d, z))
    f.close()
    return

def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=100, degree_aware=True, seed=None):
    if bins is None:
        # Get degree bins of the network
        bins = network_utilities.get_degree_binning(network, min_bin_size)
    nodes_random = network_utilities.pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware,
                                                                         seed=seed)
    return nodes_random


def get_diseasome_genes(diseasome_file, nodes=None, network=None):
    """
    If nodes is not None, keep only nodes in the network
    If network is not None, keep only LCC
    """
    disease_to_genes = {}
    disease_to_category = {}
    for line in open(diseasome_file):
        words = line.strip("\n").split("\t")
        disease = words[1].strip('"')
        category = words[0]
        genes = set(words[2:])
        if nodes is not None:
            genes &= nodes
            if len(genes) == 0:
                continue
        if network is not None:
            network_sub = network.subgraph(genes)
            genes = network_utilities.get_connected_components(network_sub, False)[0]
        disease_to_genes[disease] = genes
        disease_to_category[disease] = category
    return disease_to_genes, disease_to_category

def get_drug_target_drugbank(drug_target_file, nodes=None, network=None):
    """
    If nodes is not None, keep only nodes in the network
    If network is not None, keep only LCC
    drugbank file = [[drugid,entrez],[]...]
    groupby drugid
    return = drug = gene_set
    """

    drug_target_df = pd.read_table(drug_target_file,sep='\t')
    drug_target_df = drug_target_df.applymap(str)
    drug_targets_df = drug_target_df.groupby('DrugID').agg({'Drug_Target':list})
    print(drug_targets_df.head())
    drug_targets_dict = drug_targets_df.to_dict()['Drug_Target']

    drug_to_genes_dict = {}
    ## Check weather targets are in network

    for drug, genes in drug_targets_dict.items():
        genes = set(genes)
        if nodes is not None:
            genes &= nodes
            if len(genes) == 0:
                continue
        if network is not None:
            network_sub = network.subgraph(genes)
            genes = network_utilities.get_connected_components(network_sub, False)[0]
        drug_to_genes_dict[drug] = genes


    return drug_to_genes_dict

# def get_diseasome_genes_from_selectedGenes(gene_list_file, disease_name, disease_category=None, nodes=None, network=None ):
#     """
#     It is exactly same function as get_diseasome_genes
#     from gene list from my analysis
#     one disease type
#     :param gene_list_file: result from my anlysis (RWR)
#     :param disease_name: pre difined
#     :param disease_category: None
#     :param nodes: nodes in network
#     :param network:
#     :return: dict(disease_to_gene[disease]={genes})
#     """
#     gene_df = pd.read_table(gene_list_file, sep='\t')
#     gene_list = list(gene_df['Gene'])
#     geneID_df = mapping_utilities.mapping_genes_id_to(gene_list,id_from='symbol', id_to='entrez', species = 'human', as_dataframe = True)
#     geneID_list = list(geneID_df['_id'])
#
#     disease_to_genes = dict()
#     disease_to_category = dict()
#
#     # consider genes in only network
#     genes = set(geneID_list)
#     if nodes is not None:
#         genes &= nodes
#         # if len(genes) == 0:
#         #     continue
#     if network is not None:
#         network_sub = network.subgraph(genes)
#         genes = network_utilities.get_connected_components(network_sub, False)[0]
#     disease_to_genes[disease_name] = genes
#     disease_to_category[disease_name] = disease_category
#
#
#     return disease_to_genes , disease_to_category

def get_diseasome_genes_from_selectedGenes(gene_list, disease_name, disease_category=None, nodes=None,
                                               network=None):

    """
        It is exactly same function as get_diseasome_genes
        from gene list from my analysis
        one disease type
        :param gene_list_file: result from my anlysis (RWR) ==> gene_list
        :param disease_name: pre difined
        :param disease_category: None
        :param nodes: nodes in network
        :param network:
        :return: dict(disease_to_gene[disease]={genes})
    """
    # gene_df = pd.read_table(gene_list_file, sep='\t')
    # gene_list = list(gene_df['Gene'])

    geneID_df = mapping_utilities.mapping_genes_id_to(gene_list, id_from='symbol', id_to='entrez', species='human',
                                                      as_dataframe=True)
    geneID_list = list(geneID_df['_id'])

    disease_to_genes = dict()
    disease_to_category = dict()

    # consider genes in only network
    genes = set(geneID_list)
    if nodes is not None:
        genes &= nodes
        # if len(genes) == 0:
        #     continue
    if network is not None:
        network_sub = network.subgraph(genes)
        genes = network_utilities.get_connected_components(network_sub, False)[0]
    disease_to_genes[disease_name] = genes
    disease_to_category[disease_name] = disease_category


    return disease_to_genes , disease_to_category


if __name__ == '__main__':
    main()
