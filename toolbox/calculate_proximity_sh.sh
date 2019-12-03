#!/bin/bash

#network_file="../src_drug/Data/human_protein_interactome.sif"
network_file="../src_drug/Data/human_protein_interactome_with_STRING.sif"
#nodes_from_prefix="/home/wch23/Project/LifeArc/General/src_drug/scratch/drug_target_interaction"
nodes_from_prefix="/home/wch23/Project/LifeArc/General/data/Drug/Drugbank_drug_target_interaction_db06176.tsv"

nodes_to="/home/wch23/Project/LifeArc/SOX2/result/SOX2_STAT1/Sig.Genes/LUSC_SOX2_gene_score_by_RW_pvalue_0.001_FC_CHOL_DepMap_filter.tsv"
#nodes_to="/home/wch23/Project/LifeArc/ITGB5/result_2/CHC/Sig.Genes/CHOL_ITGB5_gene_score_by_RW_pvalue_0.0075_FC_CHOL_DepMap_filter.tsv"
#nodes_to="/home/wch23/Project/LifeArc/ULK/result/GBM_LGG_v2/Sig_Gene/GBM_LGG_ULK_conbined_gene_score_by_RW_pvalue_FC_GBM_DepMap_filter.tsv"
#nodes_to="/home/wch23/Project/LifeArc/ULK/result/GBM_LGG_v2/Sig_Gene/GBM_LGG_ULK_conbined_gene_score_by_RW_pvalue_FC_GBM_DepMap_filter_n20.tsv"
#nodes_to="/home/wch23/Project/LifeArc/ITGB5/result_2/PAAD/Sig.Genes/PAAD_ITGB5_gene_score_by_RW_pvalue_0.007_FC_PAAD_DepMap_filter.tsv"

disease_name="LUSC"
#disease_name="CHOL"
#disease_name="PAAD"
#target_name="ITGB5"
target_name="SOX2"

#output_prefix="/home/wch23/Project/LifeArc/ULK/result/GBM_LGG_v2/Drug/round2/${disease_name}_DepMap_drug_proximity_ULK_mp"
output_prefix="/home/wch23/Project/LifeArc/SOX2/result/SOX2_STAT1/Drug/round2/${disease_name}_DepMap_drug_proximity_SOX2_mp"
#output_prefix="/home/wch23/Project/LifeArc/ITGB5/result_2/CHC/Drug/round4/${disease_name}_DepMap_drug_proximity_${target_name}_mp"
#output_prefix="/home/wch23/Project/LifeArc/ITGB5/result_2/PAAD/Drug/round1/${disease_name}_DepMap_drug_proximity_${target_name}_mp"


python p01_calculate_proximity_from_nature.py -e ${network_file} -s ${nodes_from_prefix} -t $nodes_to -d $disease_name -o ${output_prefix}_db06176.txt &
