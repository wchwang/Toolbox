# Created by woochanghwang at 2019-03-28

import networkx, random, copy
import os, pickle, numpy
import pandas as pd
from pypathway import random_walk_with_restart

MAX_NUMBER_OF_TRIAL = 10
FINITE_INFINITY = 999999

def main():
    print("network utilities")
    return

def pairwise(iterable):
    import itertools
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def calc_network_similarity(G1, G2, method='VEO'):

    G1_nodes = list(G1.nodes)
    G2_nodes = list(G2.nodes)
    G1_edges = list(G1.edges)
    G2_edges = list(G2.edges)
    if method == 'VEO':
        # simiarity = (len(set(G1_nodes).intersection(G2_nodes))+len(set(G1_edges).intersection(set(G2_edges)))) / (len(set(G1_nodes).union(G2_nodes))+len(set(G1_edges).union(set(G2_edges))))
        # simiarity = (len(set(G1_nodes).intersection(G2_nodes)) + len(set(G1_edges).intersection(set(G2_edges)))) / (len(set(G1_nodes).union(G2_nodes)) + len(set(G1_edges).union(set(G2_edges))))
        similarity = (len(set(G1_nodes).intersection(set(G2_nodes))) + len(set(G1_edges).intersection(set(G2_edges)))) / (
                    len(set(G1_nodes).union(set(G2_nodes))) + len(set(G1_edges).union(set(G2_edges))))
    elif method == 'EO':
        similarity = len(set(G1_edges).intersection(set(G2_edges))) / len(set(G1_edges).union(set(G2_edges)))
    elif method == 'GED':
        d = len(G1_nodes)+len(G2_nodes)-2*len(set(G1_nodes).intersection(set(G2_nodes))) + len(G1_edges)+len(G2_edges)-2*len(set(G1_edges).intersection(set(G2_edges)))
        similarity = 1/1+d

    return similarity
def calc_RWR_top_path(top_path_addr,start_gene):

    # # init a graph
    # G = nx.Graph([[1, 2], [2, 3], [3, 5], [2, 5], [1, 4], [4, 5]])
    #
    # # heats
    # h = {1: 0, 2: 1, 3: 0, 4: 1, 5: 0}
    # dict(random_walk_with_restart(G, h, rp=0.7, n=-1).node)

    ### Make Graph

    with open(top_path_addr) as top_path_f:
        top_path_all = [x.strip().split('\t') for x in top_path_f.readlines()]

    top_G = networkx.Graph(top_path_all)

    # print(top_G.nodes())
    #######################
    ## gene_diff * RW
    ######################
    heats = dict()
    for node in top_G.nodes():
        if node == start_gene:
            heats[node] = 1
        else:
            heats[node] = 0

    # ###################
    # # RW(gene_diff)
    # ##################
    # heats = dict()
    # for node in top_G.nodes():
    #     heats[node] =

    RWR_result = dict(random_walk_with_restart(top_G, heats, rp=0.7, n=-1).node)

    # print(RWR_result[start_gene]['heat'])

    return RWR_result

def get_interactions_among_geneSet_from_PPI(gene_set, ppi_file_name,mode='both',must_have_genes=None):
    """

    :param gene_set: Gene set that we want in the context network
    :param ppi_file_name: whole interaction file
    :param mode:
    - both(defalut) : interaction[left, right] all in the gene set
    - one : interaction[left, right] one of them or both in the gene set
    :return:
    """
    with open(ppi_file_name) as string_f:
        string_interactions = [x.strip().split('\t') for x in string_f.readlines()]

    string_G = networkx.Graph(string_interactions[1:])
    core_string_interactions = []
    for interaction in string_interactions[1:]: # header == False
        if mode == 'both' :

            if len(set(interaction).intersection(set(gene_set))) == 2:
                core_string_interactions.append(interaction)

        if mode == 'one':
            if len(set(interaction).intersection(set(gene_set))) >= 1:
                core_string_interactions.append(interaction)

        if mode == "one_neighbor":
            if len(set(interaction).intersection(set(gene_set))) == 2:
                core_string_interactions.append(interaction)
            elif len(set(interaction).intersection(set(gene_set))) == 1:
                node_not_in_set = (list(set(interaction)-(set(gene_set))))
                if len(node_not_in_set) == 0 :  continue
                node_not_in_set = list(set(interaction)-(set(gene_set)))[0]
                node_in_set = set(interaction)-set([node_not_in_set])
                neighbors_of_node_not_in = string_G.neighbors(node_not_in_set)
                neighbors_of_node_not_in = list(set(neighbors_of_node_not_in) - node_in_set)
                # print("node not in set:", node_not_in_set)
                # print("neighbors set",neighbors_of_node_not_in)
                for neighbor in neighbors_of_node_not_in:
                    if neighbor in gene_set:
                        # print([node_not_in_set,neighbor])
                        core_string_interactions.append([node_not_in_set,neighbor])


    # print("all_core interactions:", len(core_string_interactions))

    core_G = networkx.Graph(core_string_interactions)

    print("all_core interactions, nodes:", len(core_G.edges()), len(core_G.nodes()))

    if must_have_genes != None:
        not_in_G = list(set(must_have_genes)-set(core_G.nodes()))
        print(not_in_G)
        if(len(not_in_G)) >0:
            for interaction in string_interactions[1:]:
                if len(set(interaction).intersection(set(not_in_G))) > 0 :

                    core_string_interactions.append(interaction)

    core_G_final = networkx.Graph(core_string_interactions)
    print("after all_core interactions, nodes:", len(core_G_final.edges()), len(core_G_final.nodes()))
    core_interactions = core_G_final.edges()

    return core_interactions

def get_directed_Graph_from_omnipath():

    '''
    1. make omnipath_network(using pandas, from_pandas_edgelist)
    :return:
    gene_symbol - gene_symbol
    edge_attr = ['is_stimulation','is_inhibition','dip_url','sources','references']
    '''


    df_omnipath = pd.read_table("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/omni_path_signal_4.tsv", sep='\t')

    print(list(df_omnipath))

    ## remove undirected edge(100 edges)
    df_omnipath_directed = df_omnipath.loc[df_omnipath['is_directed']==1]

    omni_G = networkx.from_pandas_edgelist(df_omnipath_directed, 'source_genesymbol', 'target_genesymbol', ['is_directed', 'is_stimulation','is_inhibition','dip_url','sources','references'],create_using=nx.DiGraph())

    print(len(omni_G.edges()))

    return omni_G

def make_directed_Graph_from_gene_list_using_Omnipath(geneSet,mode='All'):
    '''

    :param geneSet:
    :param mode: All- 'all directed' ,Action - 'only activation or inhibition'
    :return: omni path
    '''
    df_omnipath = pd.read_table("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/omni_path_signal_4.tsv",
                                sep='\t')

    print(list(df_omnipath))

    ## remove undirected edge(100 edges)
    df_omnipath_directed = df_omnipath.loc[df_omnipath['is_directed'] == 1]

    ## Select only edges in geneSet
    if mode == 'All':
        ## remove undirected edge(100 edges)
        df_omnipath_directed = df_omnipath.loc[df_omnipath['is_directed'] == 1]

    elif mode =='Action':
        df_omnipath_directed = df_omnipath.loc[df_omnipath['is_directed'] == 1]
        df_omnipath_directed = df_omnipath_directed.loc[(df_omnipath_directed['is_stimulation']==1)|
                                                        (df_omnipath_directed['is_inhibition']==1)]

    df_selected_omnipath_directed = df_omnipath_directed.loc[
        (df_omnipath_directed['source_genesymbol'].isin(geneSet)) &
        (df_omnipath_directed['target_genesymbol'].isin(geneSet))]

    omni_G = networkx.from_pandas_edgelist(df_selected_omnipath_directed , 'source_genesymbol', 'target_genesymbol',
                                           ['is_directed', 'is_stimulation', 'is_inhibition', 'dip_url', 'sources',
                                            'references'], create_using=networkx.DiGraph())

    return omni_G

def calc_RWR_gene_score( top_gene_RWR_result_dict):

    '''
    calculate gene score only using RW value
    :param cancer_type:
    :param top_gene_RWR_result_dict:
    :return:
    '''
    gene_diff_dict = dict()

    top_gene_score = []

    for gene, rw_value in top_gene_RWR_result_dict.items():
        gene_score = float(rw_value['heat'])

        top_gene_score.append([gene,gene_score])

    top_gene_score_sorted = sorted(top_gene_score, key=lambda gene: abs(float(gene[1])), reverse=True)

    # print(top_gene_score_sorted[:10])

    return top_gene_score_sorted

def normalize_RWR_geneScore(top_gene_score):

    top_gene_score_normalized = []
    top_gene_score_normalized.append(top_gene_score[0]) # insert start gene

    top_gene_except_start = top_gene_score[1:]  # except start gene
    max_score = max([x[1] for x in top_gene_except_start])
    min_score = min([x[1] for x in top_gene_except_start])

    for gene_score in top_gene_except_start:
        gene = gene_score[0]
        score = gene_score[1]
        normalized_score = (score-min_score)/(max_score-min_score)
        top_gene_score_normalized.append([gene,normalized_score])

    return top_gene_score_normalized


def calc_RWR_from_interaction(graph_path_addr, start_genes, normalized=False):
    '''

    :param graph_path_addr: interaction path
    :param start_genes: start genes[gene1, gene2,...]
    :param normalized:
    :return:
    '''

    with open(graph_path_addr) as graph_path_f:
        graph_path_all = [x.strip().split('\t') for x in graph_path_f.readlines()]

    graph_G = networkx.Graph(graph_path_all)

    heats = dict()
    for node in graph_G.nodes():
        if node in start_genes:
            heats[node] = 1
        else:
            heats[node] = 0

    print("Start RWR")
    RWR_result_dict = dict(random_walk_with_restart(graph_G, heats, rp=0.7, n=-1).node)

    # print(RWR_result[start_gene]['heat'])

    gene_rwr_score = []

    for gene, rw_value in RWR_result_dict.items():
        gene_score = float(rw_value['heat'])

        gene_rwr_score.append([gene, gene_score])

    print("Soring the result")
    gene_score_sorted = sorted(gene_rwr_score, key=lambda gene: abs(float(gene[1])), reverse=True)

    if normalized == True:
        print("Normalize Result")
        return normalize_RWR_geneScore(gene_score_sorted)
    if normalized == False:
        return gene_score_sorted

def calc_RWR_from_Graph(graph_G, start_genes, normalized=False, dictionary=False):
    '''

    :param graph_path_addr: interaction path
    :param start_gene: start genes[gene1, gene2,...]
    :param normalized:
    :return:
    '''

    # with open(graph_path_addr) as graph_path_f:
    #     graph_path_all = [x.strip().split('\t') for x in graph_path_f.readlines()]
    #
    # graph_G = networkx.Graph(graph_path_all)

    heats = dict()
    for node in graph_G.nodes():
        if node in start_genes:
            heats[node] = 1
        else:
            heats[node] = 0

    print("Start RWR")
    RWR_result_dict = dict(random_walk_with_restart(graph_G, heats, rp=0.7, n=-1).node)

    # print(RWR_result[start_gene]['heat'])

    gene_rwr_score = []

    for gene, rw_value in RWR_result_dict.items():
        gene_score = float(rw_value['heat'])

        gene_rwr_score.append([gene, gene_score])

    print("Soring the result")
    gene_score_sorted = sorted(gene_rwr_score, key=lambda gene: abs(float(gene[1])), reverse=True)

    if (normalized == True) & (dictionary==False):
        print("Normalize Result")
        return normalize_RWR_geneScore(gene_score_sorted)
    if (normalized == False) & (dictionary==False):
        return gene_score_sorted

    if dictionary == True:
        gene_score_dict = dict()
        for gene_score in gene_score_sorted:
            gene_score_dict[gene_score[0]] = gene_score[1]
        return gene_score_dict



def create_graph(directed=False):
    """
        Creates & returns a graph
    """
    if directed:
        g = networkx.DiGraph()
    else:
        g = networkx.Graph()
    return g

def check_has_path(G,node_from,node_to):
    return(networkx.has_path(G,node_from,node_to))

def get_subgraph(G, nodes):
    """
	NetworkX subgraph method wrapper
    """
    return G.subgraph(nodes)


def get_nodes_and_edges_from_sif_file(file_name, store_edge_type=False, delim=None, data_to_float=True):
    """
	Parse sif file into node and edge sets and dictionaries
	returns setNode, setEdge, dictNode, dictEdge
	store_edge_type: if True, dictEdge[(u,v)] = edge_value
	delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set()
    setEdge = set()
    dictNode = {}
    dictEdge = {}
    flag = False
    f = open(file_name)
    for line in f:
        if delim is None:
            words = line.rstrip("\n").split()
        else:
            words = line.rstrip("\n").split(delim)
        id1 = words[0]
        setNode.add(id1)
        if len(words) == 2:
            if data_to_float:
                score = float(words[1])
            else:
                score = words[1]
            dictNode[id1] = score
        elif len(words) >= 3:
            if len(words) > 3:
                flag = True
            id2 = words[2]
            setNode.add(id2)
            setEdge.add((id1, id2))
            if store_edge_type:
                if data_to_float:
                    dictEdge[(id1, id2)] = float(words[1])
                else:
                    dictEdge[(id1, id2)] = words[1]
    f.close()
    if len(setEdge) == 0:
        setEdge = None
    if len(dictNode) == 0:
        dictNode = None
    if len(dictEdge) == 0:
        dictEdge = None
    if flag:
        print("Warning: Ignored extra columns in the file!")
    return setNode, setEdge, dictNode, dictEdge

def create_network_from_sif_file(network_file_in_sif, use_edge_data=False, delim=None, include_unconnected=True):
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file(network_file_in_sif,
                                                                              store_edge_type=use_edge_data,
                                                                              delim=delim)
    g = create_graph()
    if include_unconnected:
        g.add_nodes_from(setNode)
    if use_edge_data:
        for e, w in dictEdge.items():
            u, v = e
            g.add_edge(u, v, w=w)  # ,{'w':w})
    else:
        g.add_edges_from(setEdge)
    return g

def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length(G, source_id, target_id)

def get_connected_components(G, return_as_graph_list=True):
    """
        Finds (strongly in the case of directed network) connected components of graph
        returnAsGraphList: returns list of graph objects corresponding to connected components (from larger to smaller)
        otherwise returns list of node list corresponding nodes in connected components
    """
    result_list = []

    if return_as_graph_list:
        result_list = networkx.connected_component_subgraphs(G)
    else:
        result_list = [c for c in sorted(networkx.connected_components(G), key=len, reverse=True)]

    return result_list


def get_degree_binning(g, bin_size, lengths=None):
    '''

    It tried to make number of gene list( with same degree) to bin size.
    If number of gene list with some degree , it combine with ohter genes with another degree to meet bin size.
    '''
    degree_to_nodes = {}
    for node, degree in g.degree():  # .items(): # iterator in networkx 2.0
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = degree_to_nodes.keys()
    # values.sort() #pytyon 2.x
    values = sorted(values)
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1
        # print i, low, high, len(val)
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins


def get_degree_equivalents(seeds, bins, g):
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bins:
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes

def pick_random_nodes_matching_selected(network, bins, nodes_selected, n_random, degree_aware=True, connected=False,
                                        seed=None):
    """
    Use get_degree_binning to get bins
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes()
    for i in range(n_random):
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network)
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
                # nodes_random.append(random.choice(equivalent_nodes))
                chosen = random.choice(equivalent_nodes)
                for k in range(20):  # Try to find a distinct node (at most 20 times)
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [random.choice(nodes)]
                k = 1
                while True:
                    if k == len(nodes_selected):
                        break
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(network.neighbors(node_random))
                    if node_selected in nodes_random:
                        continue
                    nodes_random.append(node_selected)
                    k += 1
            else:
                nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    return values


def randomize_graph(graph, randomization_type, allow_self_edges=False):
    """
    Creates a random network from given network as a networkx graph
    randomization_type:
        - "random": add same number of edges randomly between nodes of original graph
        - "preserve_topology": keep edges, shuffle nodes of original graph
        - "preserve_topology_and_node_degree": keep edges, shuffle nodes of original graph with the nodes of same degree
        - "preserve_degree_distribution": remove an edge between two random nodes with degrees k, l then add to two nodes with degrees k-1 & l-1, then shuffle nodes
        - "preserve_degree_distribution_and_node_degree": remove 2 random edges between a-b and c-d where degree(a)=degree(c) and degree(b)=degree(d) then add 2 edges between a-d and b-c, then shuffle nodes with the same degree
	- "erdos_renyi": creates a graph where edges are redistributed based on erdos renyi random model.
	- "barabasi_albert": creates a graph where edges are redistributed based on barabasi albert model (preferential attachment).
    """

    debug = False

    n_node = graph.number_of_nodes()
    n_edge = graph.number_of_edges()

    if randomization_type == "same_degree_sequence":
        # Takes ages to find a suitable conformation for large graphs
        sequence = graph.degree().values()
        new_graph = None
        while new_graph is None:
            new_graph = networkx.random_degree_sequence_graph(sequence)
        return new_graph

    if randomization_type == "graph_tool_correlated":
        try:
            import graph_tool
        except:
            raise ValueError("Graph tool package not installed")
            return
        new_graph = graph.copy()
        graph_tool.generation.random_rewire(new_graph, model='uncorrelated', n_iter=1, edge_sweep=True,
                                            parallel_edges=False, self_loops=False, vertex_corr=None,
                                            block_membership=None, alias=True, cache_probs=True, persist=False,
                                            ret_fail=False, verbose=False)
        return new_graph

    if randomization_type == "erdos_renyi":
        # raise Exception("Work in progress")
        p = float(2 * n_edge) / (n_node * n_node - 2 * n_node)
        # Chooses each of the possible [n(n-1)]/2 edges with probability p
        new_graph = networkx.erdos_renyi_graph(n_node, p)
        mapping = dict(zip(new_graph.nodes(), graph.nodes()))
        new_graph = networkx.relabel_nodes(new_graph, mapping)
        available_edges = graph.edges()

        # Map graph from random model to new graph
        for edge in new_graph.edges():
            if len(available_edges) > 0:
                edge_org = available_edges.pop()
                if debug:
                    print ("From random:", (edge[0], edge[1]))
                new_graph.add_edge(edge[0], edge[1], graph.get_edge_data(edge_org[0], edge_org[1]))
            # If the random model added too many edges
            else:
                if debug:
                    print ("Removing:", edge)
                new_graph.remove_edge(edge[0], edge[1])

        # If the random model failed to add enough edges
        nodes = new_graph.nodes()
        for edge_org in available_edges:
            source_id = random.choice(nodes)
            target_id = random.choice(nodes)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes)
            if debug:
                print ("Adding:", (source_id, target_id))
            new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge_org[0], edge_org[1]))
        return new_graph

    if randomization_type == "barabasi_albert":
        # raise Exception("Work in progress")
        if n_edge >= n_node:
            # A graph of n nodes is grown by attaching new nodes each with m edges that are preferentially attached to existing nodes with high degree
            new_graph = networkx.barabasi_albert_graph(n_node, n_edge / n_node)
            mapping = dict(zip(new_graph.nodes(), graph.nodes()))
            new_graph = networkx.relabel_nodes(new_graph, mapping)
        else:
            new_graph = networkx.create_empty_copy(graph)

        available_edges = graph.edges()
        degree_map = dict(networkx.degree(new_graph))
        nodes = new_graph.nodes()

        # Map graph from random model to new graph
        for edge in new_graph.edges():
            if len(available_edges) > 0:
                edge_org = available_edges.pop()
                if debug:
                    print ("From random:", (edge[0], edge[1]))
                new_graph.add_edge(edge[0], edge[1], graph.get_edge_data(edge_org[0], edge_org[1]))
            # If the random model added too many edges
            else:
                nodes_to_select = [id for id, d in degree_map.items() for j in range(d + 1)]
                source_id = random.choice(nodes())
                target_id = random.choice(nodes_to_select)
                if debug:
                    print ("Removing:", (source_id, target_id))
                new_graph.remove_edge(source_id, target_id)
                degree_map[source_id] -= 1
                degree_map[target_id] -= 1

            # If the random model failed to add enough edges
        for edge_org in available_edges:
            nodes_to_select = [id for id, d in degree_map.items() for j in range(d + 1)]
            source_id = random.choice(nodes)
            target_id = random.choice(nodes_to_select)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes_to_select)
            if debug:
                print ("Adding:", (source_id, target_id))
            new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge_org[0], edge_org[1]))
            degree_map[source_id] += 1
            degree_map[target_id] += 1

        return new_graph

    new_graph = networkx.create_empty_copy(graph)
    # new_graph.add_nodes_from(graph.nodes())

    if randomization_type == "random":
        nodes = new_graph.nodes()
        for edge in graph.edges():
            source_id = random.choice(nodes)
            target_id = random.choice(nodes)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes)
            new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge[0], edge[1]))

    elif randomization_type == "preserve_topology":  # shuffle_nodes
        nodes = graph.nodes()
        random_nodes = graph.nodes()
        random.shuffle(random_nodes)
        equivalences = dict([(nodes[i], random_nodes[i]) for i in range(len(nodes))])
        new_graph.add_edges_from([(equivalences[current_edge[0]], equivalences[current_edge[1]],
                                   graph.get_edge_data(current_edge[0], current_edge[1])) for current_edge in
                                  graph.edges()])

    elif randomization_type == "preserve_topology_and_node_degree":  # shuffle_nodes_within_same_degree
        nodes_by_degree = dict((degree, []) for u, degree in graph.degree())  # .values()
        graph_degree = dict(graph.degree())
        [nodes_by_degree[graph_degree[node]].append(node) for node in graph_degree]
        equivalences = {}
        for current_degree in nodes_by_degree.keys():
            nodes = nodes_by_degree[current_degree]
            random_nodes = list(nodes)
            random.shuffle(random_nodes)
            equivalences.update(dict([(nodes[i], random_nodes[i]) for i in range(len(nodes))]))
        new_graph.add_edges_from([(equivalences[current_edge[0]], equivalences[current_edge[1]],
                                   graph.get_edge_data(current_edge[0], current_edge[1])) for current_edge in
                                  graph.edges()])

    elif randomization_type == "preserve_degree_distribution":
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, graph.get_edge_data(current_node1, current_node2))
        max_degree = sorted(zip(*list(graph.degree()))[1])[-1]  # .values()
        nodes_by_degree = dict((degree, {}) for degree in range(max_degree + 1))
        graph_degree = dict(graph.degree())
        [nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree]
        n_perturbation = random.randint(2 * n_edge / 3, n_edge)  # Perturb at least 66% of the edges
        for i in range(n_perturbation):
            n_trial = 0
            while True:
                n_trial += 1
                if n_trial > MAX_NUMBER_OF_TRIAL:
                    if debug:
                        print("Warning: Max number of trials exceeded in perturbation ", i)
                    break
                source_id = random.choice(new_graph.nodes())
                source_degree = new_graph.degree(source_id)
                while source_degree < 1:
                    source_id = random.choice(new_graph.nodes())
                    source_degree = new_graph.degree(source_id)
                target_id = random.choice(new_graph.neighbors(source_id))
                target_degree = new_graph.degree(target_id)
                del nodes_by_degree[source_degree][source_id]
                nodes_by_degree[source_degree - 1].setdefault(source_id)
                if target_id == source_id:
                    target_degree -= 1
                del nodes_by_degree[target_degree][target_id]
                nodes_by_degree[target_degree - 1].setdefault(target_id)
                ## not very important to check for cases where new_source = source (v.v. for targets)
                new_target_id = random.choice(nodes_by_degree[target_degree - 1].keys())
                if source_id == target_id:
                    new_source_id = new_target_id
                else:
                    new_source_id = random.choice(nodes_by_degree[source_degree - 1].keys())
                if debug:
                    print(source_id, target_id, " / ", new_source_id, new_target_id)
                    print(source_degree, target_degree)
                ## check if going to add an existing edge or self edge
                if new_graph.has_edge(new_source_id, new_target_id) or (
                        not allow_self_edges and new_source_id == new_target_id):
                    del nodes_by_degree[target_degree - 1][target_id]
                    nodes_by_degree[target_degree].setdefault(target_id)
                    del nodes_by_degree[source_degree - 1][source_id]
                    nodes_by_degree[source_degree].setdefault(source_id)
                    continue
                if debug:
                    print("rm %s %s" % (source_id, target_id))
                edge_data = new_graph.get_edge_data(source_id, target_id)
                new_graph.remove_edge(source_id, target_id)
                if debug:
                    print("add %s %s" % (new_source_id, new_target_id))
                new_graph.add_edge(new_source_id, new_target_id, edge_data)
                del nodes_by_degree[target_degree - 1][new_target_id]
                nodes_by_degree[target_degree].setdefault(new_target_id)
                if new_source_id == new_target_id and source_id != target_id:
                    source_degree += 1
                del nodes_by_degree[source_degree - 1][new_source_id]
                nodes_by_degree[source_degree].setdefault(new_source_id)
                break
        randomize_graph(new_graph, "preserve_topology")

    elif randomization_type == "preserve_degree_distribution_and_node_degree":
        ## add edges as well
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, graph.get_edge_data(current_node1, current_node2))
        nodes_by_degree = dict((degree, {}) for u, degree in graph.degree())
        graph_degree = dict(graph.degree())
        [nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree]

        # if n_edge < MIN_NUMBER_OF_PERTURBATION:
        #    n_perturbation = random.randint(1, n_edge)
        # else:
        #    n_perturbation = random.randint(MIN_NUMBER_OF_PERTURBATION, n_edge)
        n_perturbation = random.randint(n_edge / 2, n_edge)
        for i in range(n_perturbation):
            source_id = random.choice(new_graph.nodes())
            source_degree = new_graph.degree(source_id)
            ## find a node for which another node with the same degree exists
            # available_neighbors = []
            n_trial = 0
            while True:  # (len(nodes_by_degree[source_degree]) < 2 or len(available_neighbors) < 1):
                n_trial += 1
                if n_trial > MAX_NUMBER_OF_TRIAL:
                    if debug:
                        print("Warning: Max number of trials exceeded in perturbation ", i)
                    break
                source_id = random.choice(new_graph.nodes())
                source_degree = new_graph.degree(source_id)
                if len(nodes_by_degree[source_degree]) < 2:
                    continue
                available_neighbors = []
                ## find a neighbor for which another node with the same degree exists
                for neighbor_id in new_graph.neighbors_iter(source_id):
                    if source_degree == new_graph.degree(neighbor_id):
                        if len(nodes_by_degree[new_graph.degree(neighbor_id)]) > 2:
                            available_neighbors.append(neighbor_id)
                    else:
                        if len(nodes_by_degree[new_graph.degree(neighbor_id)]) > 1:
                            available_neighbors.append(neighbor_id)
                if len(available_neighbors) < 1:
                    continue
                target_id = random.choice(available_neighbors)
                target_degree = new_graph.degree(target_id)
                ## select a new source node with different id
                n_trial2 = 0
                inner_break = False
                while True:
                    n_trial2 += 1
                    if n_trial2 > MAX_NUMBER_OF_TRIAL:
                        if debug:
                            print("Warning: Max number of trials exceeded in perturbation ", i)
                        inner_break = True
                        break
                    new_source_id = random.choice(nodes_by_degree[source_degree].keys())
                    while new_source_id == source_id:
                        new_source_id = random.choice(nodes_by_degree[source_degree].keys())
                    new_available_neighbors = []
                    ## find a neighbor as new target node for which id is different from target and has an id equivalent to target
                    for neighbor_id in new_graph.neighbors_iter(new_source_id):
                        if target_degree == new_graph.degree(neighbor_id):
                            new_available_neighbors.append(neighbor_id)
                    if len(new_available_neighbors) < 1:
                        continue
                    new_target_id = random.choice(new_available_neighbors)
                    if len(new_available_neighbors) > 1:
                        while new_target_id == target_id:
                            new_target_id = random.choice(new_available_neighbors)
                            # print new_available_neighbors, new_target_id
                    else:
                        new_target_id = new_available_neighbors[0]
                    break
                if inner_break:
                    break
                if debug:
                    print(source_id, target_id, " / ", new_source_id, new_target_id)
                if source_id == new_target_id or new_source_id == target_id:
                    continue
                if new_graph.has_edge(source_id, new_target_id) or new_graph.has_edge(new_source_id, target_id):
                    continue
                if debug:
                    print("rm %d %d" % (source_id, target_id))
                    print("rm %d %d" % (new_source_id, new_target_id))
                edge_data_1 = new_graph.get_edge_data(source_id, target_id)
                edge_data_2 = new_graph.get_edge_data(new_source_id, new_target_id)
                new_graph.remove_edge(source_id, target_id)
                new_graph.remove_edge(new_source_id, new_target_id)
                if debug:
                    print("add %d %d" % (source_id, new_target_id))
                    print("add %d %d" % (new_source_id, target_id))
                new_graph.add_edge(source_id, new_target_id, edge_data_1)
                new_graph.add_edge(new_source_id, target_id, edge_data_2)

    else:
        raise Exception("Unknown randomization type %s" % randomization_type)

    return new_graph


if __name__ == '__main__':
    main()