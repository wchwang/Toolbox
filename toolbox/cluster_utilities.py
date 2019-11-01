# Created by woochanghwang at 2019-04-09

import markov_clustering as mc
import networkx as nx

def mcl_parameter_QC(network,range_from=15,range_to=26):
    matrix = nx.to_scipy_sparse_matrix(network)
    # perform clustering using different inflation values from 1.5 and 2.5
    # for each clustering run, calculate the modularity
    for inflation in [i / 10 for i in range(range_from, range_to)]:
        result = mc.run_mcl(matrix, inflation=inflation)
        clusters = mc.get_clusters(result)
        Q = mc.modularity(matrix=result, clusters=clusters)
        print("inflation:", inflation, "modularity:", Q)

def mcl_clustering(network, inflation=2.0):
    matrix = nx.to_scipy_sparse_matrix(network,network.nodes())
    result = mc.run_mcl(matrix)  # run MCL with default parameters
    clusters = mc.get_clusters(result)  # get clusters

    return clusters, matrix

def draw_mcl_clustering(clusters,matrix,network):
    '''
    draw mcl clustering with node index
    need to convert node name to node index to find position
    :param clusters:
    :param matrix:
    :param network:
    :return:
    '''
    pos = nx.nx_pydot.graphviz_layout(network)
    pos_index =  dict()

    gene_list = list(network.nodes())
    for key, value in pos.items():
        key_id = gene_list.index(key)
        pos_index[key_id] = value

    mc.draw_graph(matrix, clusters, pos=pos_index, node_size=50, with_labels=False, edge_color="silver")

def get_cluster_classes(den, label='ivl'):
    from collections import defaultdict

    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        # print("C,pi", c,pi)
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l

    return cluster_classes

def draw_denrogram():
    from scipy.cluster.hierarchy import dendrogram, linkage


def main():
    pass


if __name__ == '__main__':
    main()