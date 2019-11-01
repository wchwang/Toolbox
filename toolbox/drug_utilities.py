# Created by woochanghwang at 2019-06-26

from rdkit import Chem


###### Drug similarity related ######

def get_smiles_similarity(smiles1, smiles2, similarity="fingerprint"):
    from rdkit import DataStructs
    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit.Chem.AtomPairs import Pairs
    """
    fp_type: sim | sub
    metric: tanimoto | tversky
    """
    if len(smiles1) == 0 or len(smiles2) == 0:
        return None

    ms = [Chem.MolFromSmiles(smiles1), Chem.MolFromSmiles(smiles2)]

    if similarity == "fingerprint":
        fps = [FingerprintMols.FingerprintMol(x) for x in ms]

        d = DataStructs.FingerprintSimilarity(fps[0],fps[1],metric=DataStructs.TanimotoSimilarity)

    elif similarity == "atom":
        pairFps = [Pairs.GetAtomPairFingerprint(x) for x in ms]
        d = DataStructs.DiceSimilarity(pairFps[0], pairFps[1])
        # d = DataStructs.TanimotoSimilarity(pairFps[0],pairFps[1])
    # print(d)
    return d


# def get_target_similarity(targets1, targets2, target_to_occurrences=None):
#     """
#     Weighted jaccard, if target_to_occurrences (diseases / side effects) is not None
#     """
#     if len(targets1) == 0 or len(targets2) == 0:
#         return None
#     targets_common = targets1 & targets2
#     if target_to_occurrences is None:
#         # d = len(targets_common) / float(max(len(targets1), len(targets2))) # ~worse
#         d = len(targets_common) / float(len(targets1 | targets2))
#     else:
#         d = 0.0
#         for t in targets_common:
#             d += 1.0 / len(target_to_occurrences[t])
#         d /= len(targets1 | targets2)  # max(len(targets1), len(targets2))
#     return d
#
#
# def get_target_ppi_similarity(targets1, targets2, network):
#     if len(targets1) == 0 or len(targets2) == 0:
#         return None
#     vals = []
#     for target1 in targets1:
#         for target2 in targets2:
#             d = network_utilities.get_shortest_path_length_between(network, target1, target2)
#             vals.append(d)
#     d = numpy.exp(-numpy.mean(vals))
#     return d
#
#
# def get_drug_similarity(drug_to_values, method="target", network=None, dump_file=None):
#     """
#     values = targets or smiles
#     """
#     if dump_file is not None and os.path.exists(dump_file):
#         drug_to_drug_similarity = cPickle.load(open(dump_file))
#         return drug_to_drug_similarity
#     if network is None and method == "target-ppi":
#         raise ValueError("Network is required for target-ppi")
#     drug_to_drug_similarity = {}
#     drugs = drug_to_values.keys()
#     for i, drug1 in enumerate(drugs):
#         for j, drug2 in enumerate(drugs):
#             if i >= j:
#                 continue
#             # comb = tuple(sorted([drug1, drug2]))
#             val1 = drug_to_values[drug1]
#             val2 = drug_to_values[drug2]
#             d = None
#             if method == "target":
#                 d = get_target_similarity(val1, val2)
#             elif method == "target-ppi":
#                 d = get_target_ppi_similarity(val1, val2, network)
#             elif method == "chemical":
#                 d = get_smiles_similarity(val1, val2)
#             else:
#                 raise ValueError("Uknown method: %s" % method)
#             drug_to_drug_similarity.setdefault(drug1, {})[drug2] = d
#             drug_to_drug_similarity.setdefault(drug2, {})[drug1] = d
#     if dump_file is not None:
#         cPickle.dump(drug_to_drug_similarity, open(dump_file, 'w'))
#     return drug_to_drug_similarity
#
#
# ###### Chemical similarity based target prediction (~SEA) related ######
#
# def get_chemical_similarity_based_target_predictions(parameters, smiles_list, cutoff=0.9, method="fishers"):
#     """
#     method: any_smiles / at_least_one_above (inherit targets of any matching smiles w.r.t. cutoff) | majority_above (inherit the target majority of whose smiles is above cutoff) | all_above (inherit the target all of whose smiles are above cutoff)
#     """
#     parser = get_drugbank(parameters)
#     smiles_to_geneids, geneid_to_smiles_strings = get_drug_smiles_by_target(parameters, parser)
#     # print len(geneid_to_smiles_strings), geneid_to_smiles_strings.items()[:5]
#     all_smiles = reduce(lambda x, y: x | y, geneid_to_smiles_strings.values())
#     smiles_to_smiles_similarity = {}
#     for smiles1 in smiles_list:
#         smiles_to_smiles_similarity[smiles1] = {}
#         for smiles2 in all_smiles:
#             try:
#                 d = get_smiles_similarity(smiles1, smiles2, fp_type="sim", metric="tanimoto")
#             except:
#                 # print smiles1, smiles2 # chirality not possible
#                 continue
#             smiles_to_smiles_similarity[smiles1][smiles2] = d
#     smiles_to_geneids_predicted = {}
#     if method == "fishers":
#         smiles_to_query_smiles = {}
#         for smiles in all_smiles:
#             for smiles_query, smiles_to_value in smiles_to_smiles_similarity.iteritems():
#                 try:
#                     d = smiles_to_value[smiles]
#                 except:
#                     continue
#                 if d >= cutoff:
#                     smiles_to_query_smiles.setdefault(smiles, set()).add(smiles_query)
#         print
#         "Drugs with matching smiles:", len(smiles_to_query_smiles)
#         smiles_to_geneids_predicted = get_side_effect_targets_fishers(smiles_to_geneids, smiles_to_query_smiles,
#                                                                       cutoff=float(parameters.get("fdr_cutoff")),
#                                                                       correct_pvalues=True)
#     elif method == "any_smiles" or method == "at_least_one_above":
#         for smiles in smiles_list:
#             for smiles2, d in smiles_to_smiles_similarity[smiles].iteritems():
#                 if d >= cutoff:
#                     print
#                     smiles2, d
#                     geneids = smiles_to_geneids_predicted.setdefault(smiles, set())
#                     geneids |= smiles_to_geneids[smiles2]
#     elif method in ("majority_above", "all_above"):
#         for smiles in smiles_list:
#             for geneid, smiles_strings in geneid_to_smiles_strings.iteritems():
#                 values = []
#                 for smiles2 in smiles_strings:
#                     try:
#                         # print smiles2 # chirality not possible
#                         d = smiles_to_smiles_similarity[smiles][smiles2]
#                         values.append(d >= cutoff)
#                     except:
#                         continue
#                 if len(values) == 0:
#                     continue
#                 if method == "majority_above":
#                     if sum(values) / float(len(values)) > 0.5:  # any(values):
#
#                         smiles_to_geneids_predicted.setdefault(smiles, set()).add(geneid)
#                 elif method == "all_above":
#                     if all(values):
#                         smiles_to_geneids_predicted.setdefault(smiles, set()).add(geneid)
#             # If smiles is among known smiles, add known targets
#             if smiles in smiles_to_geneids:
#                 smiles_to_geneids_predicted[smiles] |= smiles_to_geneids[smiles]
#     else:
#         raise ValueError("Unknown method: %s!" % method)
#     return smiles_to_geneids_predicted


def main():
    smile1 = "CNC(=O)C1=NC=CC(OC2=CC3=C(C=C2)N=C(NC2=CC(=C(Cl)C=C2)C(F)(F)F)N3)=C1"
    smile2 = "[H][C@]1(O)N=C(O)C2=C3C4=CC=CC=C4N4C3=C3N(C5=CC=CC=C5C3=C12)[C@@]1(C)O[C@]4([H])C[C@@]([H])(NC)[C@@]1([H])OC"

    d = get_smiles_similarity(smile1,smile2)
    d2 = get_smiles_similarity(smile1,smile2,similarity="atom")



if __name__ == '__main__':
    main()