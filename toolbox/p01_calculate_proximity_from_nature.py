# Created by woochanghwang at 2019-03-27

import sys
sys.path.append('/home/wch23/Project/LifeArc/General')
import argparse, os
import toolbox.network_utilities as network_util
import toolbox.wrappers as wrappers

# import network_utilities, wrappers

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--network_file') #, required=True)
    parser.add_argument('-s', '--nodes_from') #, required=True)
    parser.add_argument('-t', '--nodes_to') #, required=True)
    parser.add_argument('-d', '--disease_mode')  # , required=True)
    parser.add_argument('-o', '--out_file') #, required=True)
    parser.add_argument('-n', '--n_random', type=int, default=1000)
    parser.add_argument('-m', '--min_bin_size', type=int, default=100)
    parser.add_argument('-x', '--n_seed', type=int, default=452456)
    parser.add_argument('-f', '--parameter_file', type=str, default=None)
    parser.add_argument('-p', '--parameter_file_prefix', type=str, default=None)
    parser.add_argument('-i', '--parameter_file_start_index', type=int, default=None)
    parser.add_argument('-j', '--parameter_file_end_index', type=int, default=None)
    args = parser.parse_args()
    # Run more than once for given input files

    network = network_util.create_network_from_sif_file(args.network_file, use_edge_data=False, delim=None,
                                                        include_unconnected=True)
    wrappers.calculate_proximity_multiple(network,from_file=args.nodes_from, to_file=args.nodes_to ,disease_mode =args.disease_mode, out_file=args.out_file)


    ###########################################





    # network_file = "../src_drug/Data/human_protein_interactome.sif"
    # nodes_from = "../src_drug/scratch/drug_target_interaction_temp_1.txt"
    # nodes_to = "../../ULK/result/GBM_ULK1_2_gene_score_by_RW_pvalue_FC_230119.tsv"
    # disease_name = "GBM"
    #
    # output_file = "../src_drug/Result/{}_drug_proximity_t_1_1.tsv".format(disease_name)
    #
    #
    # network = network_util.create_network_from_sif_file(network_file, use_edge_data=False, delim=None,
    #                                                     include_unconnected=True)
    # wrappers.calculate_proximity_multiple(network, from_file=nodes_from, to_file=nodes_to,
    #                                       disease_mode=disease_name, out_file=output_file)

    ###########################################
    # if args.parameter_file_prefix is not None:
    #     parameter_file_prefix = args.parameter_file_prefix
    #     i_start = args.parameter_file_start_index
    #     i_end = args.parameter_file_end_index
    #     calculate_proximity_multiple(parameter_file_prefix, i_start, i_end)
    #     return
    # # # Run from input parameter file
    # # elif args.parameter_file_prefix is not None:min_bin_size
    # #     network_file, nodes_from, nodes_to, out_file, , n_random, n_seed = get_parameters_from_file(
    # #         args.parameter_file_prefix + "%s.txt" % 'n')
    # # Run once with provided arguments
    # else:
    #     nodes_from = args.nodes_from.split(",")
    #     nodes_to = args.nodes_to.split(",")
    #     network_file = args.network_file
    #     n_random = args.n_random
    #     min_bin_size = args.min_bin_size
    #     n_seed = args.n_seed
    #     out_file = args.out_file
    # network = network_util.create_network_from_sif_file(network_file, use_edge_data=False, delim=None,
    #                                                              include_unconnected=True)
    # # print args
    # print(network_file, nodes_from, nodes_to, n_random, min_bin_size, n_seed, out_file)
    # values = wrappers.calculate_proximity(network, nodes_from=nodes_from, nodes_to=nodes_to, n_random=n_random,
    #                                       min_bin_size=min_bin_size, seed=n_seed)
    # if values is not None:  # not in network
    #     d, z, (m, s) = values
    #     # print z, d, (m, s)
    #     open(out_file, 'w').write("%f %f %f %f\n" % (z, d, m, s))
    return

def main_pree():

    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--network_file') #, required=True)
    parser.add_argument('-s', '--nodes_from') #, required=True)
    parser.add_argument('-t', '--nodes_to') #, required=True)
    parser.add_argument('-o', '--out_file') #, required=True)
    parser.add_argument('-n', '--n_random', type=int, default=1000)
    parser.add_argument('-m', '--min_bin_size', type=int, default=100)
    parser.add_argument('-x', '--n_seed', type=int, default=452456)
    parser.add_argument('-f', '--parameter_file', type=str, default=None)
    parser.add_argument('-p', '--parameter_file_prefix', type=str, default=None)
    parser.add_argument('-i', '--parameter_file_start_index', type=int, default=None)
    parser.add_argument('-j', '--parameter_file_end_index', type=int, default=None)
    args = parser.parse_args()
    # Run more than once for given input files
    if args.parameter_file_prefix is not None:
        parameter_file_prefix = args.parameter_file_prefix
        i_start = args.parameter_file_start_index
        i_end = args.parameter_file_end_index
        calculate_proximity_multiple(parameter_file_prefix, i_start, i_end)
        return
    # # Run from input parameter file
    # elif args.parameter_file_prefix is not None:
    #     network_file, nodes_from, nodes_to, out_file, min_bin_size, n_random, n_seed = get_parameters_from_file(
    #         args.parameter_file_prefix + "%s.txt" % 'n')
    # Run once with provided arguments
    else:
        nodes_from = args.nodes_from.split(",")
        nodes_to = args.nodes_to.split(",")
        network_file = args.network_file
        n_random = args.n_random
        min_bin_size = args.min_bin_size
        n_seed = args.n_seed
        out_file = args.out_file
        network = network_util.create_network_from_sif_file(network_file, use_edge_data=False, delim=None,
                                                                 include_unconnected=True)
    # print args
    print(network_file, nodes_from, nodes_to, n_random, min_bin_size, n_seed, out_file)
    values = wrappers.calculate_proximity(network, nodes_from=nodes_from, nodes_to=nodes_to, n_random=n_random,
                                          min_bin_size=min_bin_size, seed=n_seed)
    if values is not None:  # not in network
        d, z, (m, s) = values
        # print z, d, (m, s)
        open(out_file, 'w').write("%f %f %f %f\n" % (z, d, m, s))
    return


def b_calculate_proximity_multiple(parameter_file_prefix, i_start, i_end):
    network_file, nodes_from, nodes_to, out_file, min_bin_size, n_random, n_seed = get_parameters_from_file(
        parameter_file_prefix + "%s.txt" % i_start)
    network = network_util.create_network_from_sif_file(network_file, use_edge_data=False, delim=None,
                                                             include_unconnected=True)
    bins = network_util.get_degree_binning(network, min_bin_size, lengths=None)
    for i in range(i_start, i_end):
        if not os.path.exists(parameter_file_prefix + "%s.txt" % i):
            print("File does not exists for index (aborting):", i)
            break
        network_file, nodes_from, nodes_to, out_file, min_bin_size, n_random, n_seed = get_parameters_from_file(
            parameter_file_prefix + "%s.txt" % i)
        if os.path.exists(out_file):
            print("Skipping existing file for index:", i)
            continue
        print(network_file, nodes_from, nodes_to, n_random, min_bin_size, n_seed, out_file)
        values = wrappers.calculate_proximity(network, nodes_from=nodes_from, nodes_to=nodes_to, bins=bins,
                                              n_random=n_random, min_bin_size=min_bin_size, seed=n_seed)
        if values is not None:  # not in network
            d, z, (m, s) = values
            # print z, d, (m, s)
            open(out_file, 'w').write("%f %f %f %f\n" % (z, d, m, s))
    return

def calculate_proximity_multiple(parameter_file_prefix, i_start, i_end):
    network_file, nodes_from, nodes_to, out_file, min_bin_size, n_random, n_seed = get_parameters_from_file(
        parameter_file_prefix + "%s.txt" % i_start)
    network = network_util.create_network_from_sif_file(network_file, use_edge_data=False, delim=None,
                                                             include_unconnected=True)
    bins = network_util.get_degree_binning(network, min_bin_size, lengths=None)
    for i in range(i_start, i_end):
        if not os.path.exists(parameter_file_prefix + "%s.txt" % i):
            print("File does not exists for index (aborting):", i)
            break
        network_file, nodes_from, nodes_to, out_file, min_bin_size, n_random, n_seed = get_parameters_from_file(
            parameter_file_prefix + "%s.txt" % i)
        if os.path.exists(out_file):
            print("Skipping existing file for index:", i)
            continue
        print(network_file, nodes_from, nodes_to, n_random, min_bin_size, n_seed, out_file)
        values = wrappers.calculate_proximity(network, nodes_from=nodes_from, nodes_to=nodes_to, bins=bins,
                                              n_random=n_random, min_bin_size=min_bin_size, seed=n_seed)
        if values is not None:  # not in network
            d, z, (m, s) = values
            # print z, d, (m, s)
            open(out_file, 'w').write("%f %f %f %f\n" % (z, d, m, s))
    return


def get_parameters_from_file(file_name):
    network_file = None
    nodes_from = None
    nodes_to = None
    out_file = None
    min_bin_size = None
    n_random = None
    n_seed = None
    arguments = open(file_name).read().split()
    # print arguments
    i = 0
    if arguments[0] == "":
        i = 1
    while i < len(arguments):
        if arguments[i] == "-e":
            network_file = arguments[i + 1]
        if arguments[i] == "-s":
            nodes_from = arguments[i + 1].split(",")
        if arguments[i] == "-t":
            nodes_to = arguments[i + 1].split(",")
        if arguments[i] == "-o":
            out_file = arguments[i + 1]
        if arguments[i] == "-m":
            min_bin_size = int(arguments[i + 1])
        if arguments[i] == "-n":
            n_random = int(arguments[i + 1])
        if arguments[i] == "-x":
            n_seed = int(arguments[i + 1])
        # print arguments[i], arguments[i+1]
        i += 2
    return network_file, nodes_from, nodes_to, out_file, min_bin_size, n_random, n_seed

if __name__ == "__main__":
    main()

