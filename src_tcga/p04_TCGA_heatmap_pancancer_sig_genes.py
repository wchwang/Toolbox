# Created by woochanghwang at 21/10/2019


# Created by woochanghwang at 15/10/2018
'''
Modified by woochanghwang at 22/10/2018
# -- to add anova test, and select significant genes
# -- plot heatmap with signicant genes
-get_TCGA_FC_matrix_for_all_tumourTypes()
'''



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import toolbox.data_handlers as dh
import csv
import sys
import statsmodels.api as sm
from statsmodels.formula.api import ols

sys.setrecursionlimit(100000)


def get_TCGA_cancer_type_matrix(tumour_tcga_addr,normal_tcga_addr, cancer_type,gene_type):

    all_tumour_table = pd.read_table(tumour_tcga_addr, sep='\t')
    all_normal_table = pd.read_table(normal_tcga_addr, sep='\t')


    all_tumour_table = all_tumour_table.replace('\n', '', regex=True)
    all_normal_table = all_normal_table.replace('\n', '', regex=True)

    # print(all_tumour_table.head())
    # print(all_normal_table["Entrez",cancer_type].head())

    filter_tumour_col = [col for col in all_tumour_table if col.startswith(cancer_type)]
    filter_normal_col = [col for col in all_normal_table if col.startswith(cancer_type)]

    cancer_type_tumour_df = all_tumour_table[filter_tumour_col]
    cancer_type_normal_df = all_normal_table[filter_normal_col]

    tcga_symbol_entrez = all_tumour_table[['Symbol','Entrez','new_Symbol']]


    print(cancer_type_tumour_df.head())
    print(cancer_type_normal_df.head())
    print(tcga_symbol_entrez.head())

    all_tumour_table_with_symbol_df = pd.concat([tcga_symbol_entrez,cancer_type_tumour_df],axis = 1,sort=False)
    all_normal_table_with_symbol_df = pd.concat([tcga_symbol_entrez,cancer_type_normal_df],axis = 1,sort=False)

    # print(all_normal_table_with_symbol_df)

    data_index = 2  # new symbol
    all_normal_table_with_symbol_df.to_csv("/Users/woochanghwang/PycharmProjects/TCGA/data/RNAseqV2/processed.gene.norm/TCGA_"+cancer_type+"_normal"+gene_type+".tsv",sep='\t',index=data_index)
    all_tumour_table_with_symbol_df.to_csv(
        "/Users/woochanghwang/PycharmProjects/TCGA/data/RNAseqV2/processed.gene.norm/TCGA_" + cancer_type + "_tumour"+gene_type+".tsv",
        sep='\t', index=data_index)

    return

def get_TCGA_cancer_types_in_both_tumour_normal(tcga_tumour_matrix_addr, tcga_normal_matrix_addr, gene_type):
    all_tumour_columns = pd.read_table(tcga_tumour_matrix_addr, sep='\t',nrows = 1).columns
    all_normal_columns = pd.read_table(tcga_normal_matrix_addr, sep='\t',nrows = 1).columns

    all_tumour_types =list(set([x.split('.')[0] for x in list(all_tumour_columns)[3:]]))
    all_normal_types = list(set([x.split('.')[0] for x in list(all_normal_columns)[3:]]))

    print("tumour:", len(all_tumour_types))
    print("normal:",len(all_normal_types))
    print("both", len(set(all_tumour_types) & set(all_normal_types)))

    tumor_types = list(set(all_tumour_types) & set(all_normal_types))

    return tumor_types

def get_TCGA_FC_matrix_for_all_tumourTypes(tcga_tumour_matrix_addr, tcga_normal_matrix_addr, tumour_types):
    all_tumour_table = pd.read_table(tcga_tumour_matrix_addr, sep='\t')
    all_normal_table = pd.read_table(tcga_normal_matrix_addr, sep='\t')

    all_tumour_table = all_tumour_table.replace('\n', '', regex=True)
    all_normal_table = all_normal_table.replace('\n', '', regex=True)

    tcga_symbol_entrez = all_tumour_table[['Symbol', 'Entrez', 'new_Symbol']]
    tcga_symbol_entrez['Entrez'] = tcga_symbol_entrez['Entrez'].astype(str)

    tcga_symbol_entrez['Symbol|Entrez'] = tcga_symbol_entrez[['new_Symbol', 'Entrez']].apply(lambda x: '|'.join(x),
                                                                                             axis=1)

    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_entrez['Symbol|Entrez']

    # print(tcga_symbol_and_entrez_all_tumours_df.head())
    for tumour_type in tumour_types:
        filter_tumour_col = [col for col in all_tumour_table if col.startswith(tumour_type)]
        filter_normal_col = [col for col in all_normal_table if col.startswith(tumour_type)]

        tcga_tumour_a_type_df = all_tumour_table[filter_tumour_col]
        tcga_normal_a_type_df = all_normal_table[filter_normal_col]

        tcga_normal_a_type_df['Normal_mean'] = tcga_normal_a_type_df.mean(axis=1)

        tcga_tumour_a_type_diff_df = tcga_tumour_a_type_df.div(tcga_normal_a_type_df.mean(axis=1), axis="index")
        tcga_tumour_a_type_diff_df = np.log2(tcga_tumour_a_type_diff_df)

        tcga_symbol_and_entrez_all_tumours_df = pd.concat(
            (tcga_symbol_and_entrez_all_tumours_df, tcga_tumour_a_type_diff_df), axis=1)

    print(tcga_symbol_and_entrez_all_tumours_df.head())

    # # column_names = pd.Series(tcga_symbol_and_entrez_all_tumours_df.columns,name='Cancer_Type')
    column_names = pd.Series(tcga_symbol_and_entrez_all_tumours_df.columns)
    # # column_names = column_names[1:]
    print(column_names.head())
    column_names = column_names.str.replace('\.\d+', '')  # remove digital in colume names for groupby
    column_names[0] = 'Cancer_Type'
    column_names_df = pd.DataFrame(column_names).T
    column_names_df.columns = tcga_symbol_and_entrez_all_tumours_df.columns

    print(column_names_df.head())
    #
    column_names_df.reset_index()
    tcga_symbol_and_entrez_all_tumours_df.reset_index()

    tcga_symbol_and_entrez_all_tumours_df = pd.concat((tcga_symbol_and_entrez_all_tumours_df, column_names_df),
                                                      ignore_index=True, axis=0)

    tcga_symbol_and_entrez_all_tumours_df.set_index('Symbol|Entrez', inplace=True)
    ######################
    ## To remove inf, NA
    ######################
    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_and_entrez_all_tumours_df.replace([np.inf, -np.inf], np.nan)
    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_and_entrez_all_tumours_df.dropna()

    print(tcga_symbol_and_entrez_all_tumours_df)

    tcga_symbol_and_entrez_all_tumours_df.to_csv('../result/TCGA_log_FC.tsv', sep='\t', quoting=csv.QUOTE_NONE)
    return tcga_symbol_and_entrez_all_tumours_df

def get_TCGA_FC_matrix_for_all_tumourTypes_pre(tcga_tumour_matrix_addr, tcga_normal_matrix_addr, tumour_types):
    all_tumour_table = pd.read_table(tcga_tumour_matrix_addr, sep='\t')
    all_normal_table = pd.read_table(tcga_normal_matrix_addr, sep='\t')

    all_tumour_table = all_tumour_table.replace('\n', '', regex=True)
    all_normal_table = all_normal_table.replace('\n', '', regex=True)


    tcga_symbol_entrez = all_tumour_table[['Symbol', 'Entrez', 'new_Symbol']]
    tcga_symbol_entrez['Entrez'] = tcga_symbol_entrez['Entrez'].astype(str)

    tcga_symbol_entrez['Symbol|Entrez'] = tcga_symbol_entrez[['new_Symbol','Entrez']].apply(lambda x: '|'.join(x), axis = 1)

    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_entrez[['Symbol|Entrez']]

    for tumour_type in tumour_types:
        # tumour_type = tumour_types[0]

        filter_tumour_col = [col for col in all_tumour_table if col.startswith(tumour_type)]
        filter_normal_col = [col for col in all_normal_table if col.startswith(tumour_type)]

        tcga_tumour_a_type_df = all_tumour_table[filter_tumour_col]
        tcga_normal_a_type_df = all_normal_table[filter_normal_col]

        tcga_tumour_a_type_df['Tumour_mean'] = tcga_tumour_a_type_df.mean(axis = 1)
        tcga_normal_a_type_df['Normal_mean'] = tcga_normal_a_type_df.mean(axis = 1)

        # print(tcga_tumour_a_type_df.head())
        # print("row_count_tcga_tumour", tcga_normal_a_type_df.count())
        a_type_tcga_tumour_normal_means_df = pd.concat([tcga_tumour_a_type_df['Tumour_mean'],tcga_normal_a_type_df['Normal_mean']],axis=1)
        tcga_symbol_and_entrez_all_tumours_df[tumour_type] = np.log2(a_type_tcga_tumour_normal_means_df['Tumour_mean']/a_type_tcga_tumour_normal_means_df['Normal_mean'])

    tcga_symbol_and_entrez_all_tumours_df.set_index('Symbol|Entrez',inplace = True)
    # print(tcga_symbol_and_entrez_all_tumours_df.head())

    # print("after combine:", a_type_tcga_tumour_normal_means_df.count())

    ## To remove inf, NA
    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_and_entrez_all_tumours_df.replace([np.inf, -np.inf], np.nan)
    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_and_entrez_all_tumours_df.dropna()

    tcga_symbol_and_entrez_all_tumours_df.to_csv('../result/TCGA_log_FC.tsv', sep='\t', quoting=csv.QUOTE_NONE)
    return tcga_symbol_and_entrez_all_tumours_df

def get_TCGA_Diff_matrix_for_all_tumourTypes(tcga_tumour_matrix_addr, tcga_normal_matrix_addr, tumour_types):
    all_tumour_table = pd.read_table(tcga_tumour_matrix_addr, sep='\t')
    all_normal_table = pd.read_table(tcga_normal_matrix_addr, sep='\t')

    all_tumour_table = all_tumour_table.replace('\n', '', regex=True)
    all_normal_table = all_normal_table.replace('\n', '', regex=True)


    tcga_symbol_entrez = all_tumour_table[['Symbol', 'Entrez', 'new_Symbol']]
    tcga_symbol_entrez['Entrez'] = tcga_symbol_entrez['Entrez'].astype(str)

    tcga_symbol_entrez['Symbol|Entrez'] = tcga_symbol_entrez[['new_Symbol','Entrez']].apply(lambda x: '|'.join(x), axis = 1)

    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_entrez[['Symbol|Entrez']]

    for tumour_type in tumour_types:
        # tumour_type = tumour_types[0]

        filter_tumour_col = [col for col in all_tumour_table if col.startswith(tumour_type)]
        filter_normal_col = [col for col in all_normal_table if col.startswith(tumour_type)]

        tcga_tumour_a_type_df = all_tumour_table[filter_tumour_col]
        tcga_normal_a_type_df = all_normal_table[filter_normal_col]

        # tcga_tumour_a_type_mean_df['Tumour_mean'] = tcga_tumour_a_type_df.mean(axis = 1)
        tcga_normal_a_type_df['Normal_mean'] = tcga_normal_a_type_df.mean(axis = 1)

        tcga_tumour_a_type_diff_df = tcga_tumour_a_type_df.subtract(tcga_normal_a_type_df.mean(axis = 1),axis = "index")

        # print(tcga_tumour_a_type_df.head())
        # print(tcga_normal_a_type_df.head())
        # print(tcga_tumour_a_type_diff_df.head())


        # tcga_symbol_and_entrez_all_tumours_df.merge(tcga_tumour_a_type_diff_df)
        # print(tcga_symbol_and_entrez_all_tumours_df.head())
        tcga_symbol_and_entrez_all_tumours_df = pd.concat((tcga_symbol_and_entrez_all_tumours_df,tcga_tumour_a_type_diff_df),axis = 1)

        # a_type_tcga_tumour_normal_means_df = pd.concat([tcga_tumour_a_type_df['Tumour_mean'],tcga_normal_a_type_df['Normal_mean']],axis=1)
        # tcga_symbol_and_entrez_all_tumours_df[tumour_type] = np.log2(a_type_tcga_tumour_normal_means_df['Tumour_mean']/a_type_tcga_tumour_normal_means_df['Normal_mean'])

        ################################
        # # column_names = pd.Series(tcga_symbol_and_entrez_all_tumours_df.columns,name='Cancer_Type')
    column_names = pd.Series(tcga_symbol_and_entrez_all_tumours_df.columns)
    # # column_names = column_names[1:]
    print(column_names.head())
    column_names = column_names.str.replace('\.\d+', '')  # remove digital in colume names for groupby
    column_names[0] = 'Cancer_Type'
    column_names_df = pd.DataFrame(column_names).T
    column_names_df.columns = tcga_symbol_and_entrez_all_tumours_df.columns

    print(column_names_df.head())
    #
    column_names_df.reset_index()
    tcga_symbol_and_entrez_all_tumours_df.reset_index()

    tcga_symbol_and_entrez_all_tumours_df = pd.concat((tcga_symbol_and_entrez_all_tumours_df, column_names_df),
                                                      ignore_index=True, axis=0)

    tcga_symbol_and_entrez_all_tumours_df.set_index('Symbol|Entrez', inplace=True)

    # tcga_symbol_and_entrez_all_tumours_df.set_index('Symbol|Entrez',inplace = True)
    print(tcga_symbol_and_entrez_all_tumours_df.head())

    # print("after combine:", a_type_tcga_tumour_normal_means_df.count())

    ## To remove inf, NA
    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_and_entrez_all_tumours_df.replace([np.inf, -np.inf], np.nan)
    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_and_entrez_all_tumours_df.dropna()

    tcga_symbol_and_entrez_all_tumours_df.to_csv('../result/tmp_TCGA_diff.tsv', sep='\t', quoting=csv.QUOTE_NONE)
    return tcga_symbol_and_entrez_all_tumours_df


def get_TCGA_Diff_matrix_for_all_tumourTypes_pre(tcga_tumour_matrix_addr, tcga_normal_matrix_addr, tumour_types):
    all_tumour_table = pd.read_table(tcga_tumour_matrix_addr, sep='\t')
    all_normal_table = pd.read_table(tcga_normal_matrix_addr, sep='\t')

    all_tumour_table = all_tumour_table.replace('\n', '', regex=True)
    all_normal_table = all_normal_table.replace('\n', '', regex=True)


    tcga_symbol_entrez = all_tumour_table[['Symbol', 'Entrez', 'new_Symbol']]
    tcga_symbol_entrez['Entrez'] = tcga_symbol_entrez['Entrez'].astype(str)

    tcga_symbol_entrez['Symbol|Entrez'] = tcga_symbol_entrez[['new_Symbol','Entrez']].apply(lambda x: '|'.join(x), axis = 1)

    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_entrez[['Symbol|Entrez']]

    for tumour_type in tumour_types:
        # tumour_type = tumour_types[0]

        filter_tumour_col = [col for col in all_tumour_table if col.startswith(tumour_type)]
        filter_normal_col = [col for col in all_normal_table if col.startswith(tumour_type)]

        tcga_tumour_a_type_df = all_tumour_table[filter_tumour_col]
        tcga_normal_a_type_df = all_normal_table[filter_normal_col]

        tcga_tumour_a_type_df['Tumour_mean'] = tcga_tumour_a_type_df.mean(axis = 1)
        tcga_normal_a_type_df['Normal_mean'] = tcga_normal_a_type_df.mean(axis = 1)

        # a_type_tcga_tumour_normal_means_df = pd.concat([tcga_tumour_a_type_df['Tumour_mean'],tcga_normal_a_type_df['Normal_mean']],axis=1)
        # tcga_symbol_and_entrez_all_tumours_df[tumour_type] = np.log2(a_type_tcga_tumour_normal_means_df['Tumour_mean']/a_type_tcga_tumour_normal_means_df['Normal_mean'])

    tcga_symbol_and_entrez_all_tumours_df.set_index('Symbol|Entrez',inplace = True)
    # print(tcga_symbol_and_entrez_all_tumours_df.head())

    # print("after combine:", a_type_tcga_tumour_normal_means_df.count())

    ## To remove inf, NA
    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_and_entrez_all_tumours_df.replace([np.inf, -np.inf], np.nan)
    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_and_entrez_all_tumours_df.dropna()

    tcga_symbol_and_entrez_all_tumours_df.to_csv('../result/TCGA_log_FC.tsv', sep='\t', quoting=csv.QUOTE_NONE)
    return tcga_symbol_and_entrez_all_tumours_df


def get_selected_genes_using_oneway_ANOVA(TCGA_df, p_value_threshold, cancer_type):
    anova_res = pd.DataFrame(columns=['Gene','F-score','P-value'])

    TCGA_df_T = TCGA_df.T

    significant_gene = []

    for gene in TCGA_df_T.columns:
        TCGA_df_T[gene]
        a_gene_cancertype_df = pd.concat((TCGA_df_T[gene],cancer_type), axis=1)
        a_gene_cancertype_df['Cancer_Type'].astype("category")

        anova_model = 'Q(\"{}\") ~ C(Cancer_Type)'.format(gene)
        # print(anova_model)
        # print(a_gene_cancertype_df.head())
        mod = ols(anova_model, data = a_gene_cancertype_df).fit()
        aov_table = sm.stats.anova_lm(mod, typ=2)
        # print (aov_table['F'][0], aov_table['PR(>F)'][0])
        # print(aov_table)
        pvalue_anova = aov_table['PR(>F)'][0]
        if pvalue_anova < p_value_threshold:
            # print ("thr, pvalue:", p_value_threshold, pvalue_anova)
            significant_gene.append(gene)

    print("sig gene:", len(significant_gene))
    print("whole", TCGA_df_T.columns)

    return significant_gene

def main_before():
    gene_type = "_new_genesymbol"
    gene_value_mode = "Diff"  # FC , Diff
    tcga_tumour_matrix_addr = "/Users/woochanghwang/PycharmProjects/TCGA/data/RNAseqV2/processed.gene.norm/TCGA_all_tumour" + gene_type + ".tsv"
    tcga_normal_matrix_addr = "/Users/woochanghwang/PycharmProjects/TCGA/data/RNAseqV2/processed.gene.norm/TCGA_all_normal" + gene_type + ".tsv"

    tumour_types_in_TCGA = get_TCGA_cancer_types_in_both_tumour_normal(tcga_tumour_matrix_addr, tcga_normal_matrix_addr,
                                                                       gene_type)


    if gene_value_mode == "FC":
        TCGA_df = get_TCGA_FC_matrix_for_all_tumourTypes(tcga_tumour_matrix_addr, tcga_normal_matrix_addr,
                                                         tumour_types_in_TCGA)
    elif gene_value_mode == "Diff":
        TCGA_df = get_TCGA_Diff_matrix_for_all_tumourTypes(tcga_tumour_matrix_addr, tcga_normal_matrix_addr,
                                                           tumour_types_in_TCGA)
    # print(TCGA_FC_df.head())

    TCGA_df, cancer_type = TCGA_df.drop(TCGA_df.tail(1).index), TCGA_df.tail(1)
    TCGA_df = TCGA_df.astype('float64')

    cancer_type = cancer_type.iloc[0]

    ######################
    ## ANOVA
    ######################

    # p_value_threshold = 0.000000000000000000000000000000000001
    p_value_threshold = 0.1e-300    # -350 = significant gene(0) , -320 = significant(8300)
    selected_genes_df = get_selected_genes_using_oneway_ANOVA(TCGA_df, p_value_threshold, cancer_type)

    ######################
    ## set multiple colors
    ######################
    # rgb_colors = sns.color_palette("Set2",
    #                                len(cancer_type.unique()))  # http://seaborn.pydata.org/tutorial/color_palettes.html
    #
    # cancer_type_color = dict(zip(cancer_type.unique(), rgb_colors))
    #
    # # print (cancer_type_color)
    # col_colors = cancer_type.map(cancer_type_color)
    # #
    #
    # print(TCGA_df.info())
    #
    # g = sns.clustermap(TCGA_df, metric="correlation", cmap="RdBu_r", robust=True, method="average",
    #                    col_colors=col_colors)  # Average is best
    #
    # g.savefig('../result/tmp_TCGA_diff_heatmap_single_robust_colColor.pdf')
    #
    # plt.show()

def main():
    # organism = '9606'
    # string_node_addr = "/home/wch23/Project/LifeArc/General/data/STRING/{}.protein.links.v11.0.400.nodes.txt".format(
    #     organism)
    # with open(string_node_addr) as string_node_f:
    #     string_node = [x.strip() for x in string_node_f.readlines()]

    ####
    # get sig genes from anova reslut
    ##########

    target_cancer = 'PAAD'
    # target_cancer = "LUSC"
    gene_value_mode = 'Diff'
    # anova_result_addr = "/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_{}_vs_other_{}_anova_result_with_ensembl.tsv".format(target_cancer,gene_value_mode)
    anova_result_addr = "/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_result.csv".format(target_cancer, gene_value_mode)
    anova_result_df = pd.read_csv(anova_result_addr)

    p_value_threshold = 1.0E-150

    anova_sig_result_df = anova_result_df.loc[anova_result_df['pvalue']<p_value_threshold]
    print(anova_sig_result_df.shape)
    print(anova_sig_result_df.head())
    significant_genes = anova_sig_result_df['gene'].tolist()

    print(significant_genes)

    anova_sig_result_df.to_csv("/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_sig_genes_pv_{}.tsv".format(target_cancer, gene_value_mode,p_value_threshold),sep='\t')

    with open ("/home/wch23/Project/LifeArc/SOX2/result/Sig.Genes/{}_anova_sig_genes.txt".format(target_cancer),'w') as sig_f:
        sig_f.write('\n'.join(significant_genes))
    #####
    # draw heatmap with sig genes
    #####
    tcga_cancer_diff_df = dh.load_obj("/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original")

    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')

    tcga_cancer_diff_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'] != 'None']
    tcga_cancer_diff_df = tcga_cancer_diff_df.drop(columns=['ensembl_id', 'ensembl_gene'])

    tcga_cancer_original_order = list(tcga_cancer_diff_df)

    tcga_cancer_diff_sig_gene_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'].isin(significant_genes)]

    print(tcga_cancer_diff_sig_gene_df.head())
    print(cancer_id_df.head())
    cancer_id_df = cancer_id_df.set_index('fullcode')
    tcga_cancer_diff_sig_gene_df = tcga_cancer_diff_sig_gene_df.set_index('gene_symbol')
    tcga_cancer_diff_sig_gene_cancer_type_df = pd.concat([tcga_cancer_diff_sig_gene_df.T, cancer_id_df],axis=1, join='inner' )
    tcga_cancer_diff_sig_gene_cancer_type_df = tcga_cancer_diff_sig_gene_cancer_type_df.T
    print(tcga_cancer_diff_sig_gene_cancer_type_df)
    cancer_type = tcga_cancer_diff_sig_gene_cancer_type_df.loc['Abbreviation']
    print(cancer_type)

    rgb_colors = sns.color_palette("hls",  len(cancer_type.unique()))

    cancer_type_color = dict(zip(cancer_type.unique(), rgb_colors))

    print (cancer_type_color)
    col_colors = cancer_type.map(cancer_type_color)
    #

    print(tcga_cancer_diff_sig_gene_df.info())

    g = sns.clustermap(tcga_cancer_diff_sig_gene_df, metric="correlation", cmap="RdBu_r", robust=True, method="average",z_score=0,
                       col_colors=col_colors,xticklabels=False)  # Average is best
    # g = sns.clustermap(tcga_cancer_diff_sig_gene_df, metric="correlation", cmap="RdBu_r", robust=True, method="average",
    #                    col_colors=col_colors, xticklabels=False)  # Average is best

    g.savefig('/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_diff_heatmap_single_robust_colColor_norm.png'.format(target_cancer))

    print(g.dendrogram_col.reordered_ind)
    print(g.dendrogram_row.reordered_ind)

    clustred_col = g.dendrogram_col.reordered_ind

    cancer_type_df = tcga_cancer_diff_sig_gene_cancer_type_df.T[['Abbreviation']]
    cancer_type_df.index.name = 'fullcode'
    cancer_type_df = cancer_type_df.reset_index()
    print(cancer_type_df)

    clusterd_cancer_type_df = cancer_type_df.reindex(clustred_col)
    print(clusterd_cancer_type_df)

    clusterd_cancer_type_df.to_csv('/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_diff_heatmap_tcga_clustred_result.csv'.format(target_cancer))



def TCGA_heatmap_pancacer_sig_genes(target_cancer, p_value_threshold, anova_result_addr):
    # organism = '9606'
    # string_node_addr = "/home/wch23/Project/LifeArc/General/data/STRING/{}.protein.links.v11.0.400.nodes.txt".format(
    #     organism)
    # with open(string_node_addr) as string_node_f:
    #     string_node = [x.strip() for x in string_node_f.readlines()]

    ####
    # get sig genes from anova reslut
    ##########

    # target_cancer = 'PAAD'
    # target_cancer = "LUSC"
    gene_value_mode = 'Diff'
    # anova_result_addr = "/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_{}_vs_other_{}_anova_result_with_ensembl.tsv".format(target_cancer,gene_value_mode)
    # anova_result_addr = "/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_result.csv".format(target_cancer, gene_value_mode)
    anova_result_df = pd.read_csv(anova_result_addr)

    # p_value_threshold = 1.0E-150

    anova_sig_result_df = anova_result_df.loc[anova_result_df['pvalue']<p_value_threshold]
    print(anova_sig_result_df.shape)
    print(anova_sig_result_df.head())
    significant_genes = anova_sig_result_df['gene'].tolist()

    print(significant_genes)

    anova_sig_result_df.to_csv("/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_sig_genes_pv_{}.tsv".format(target_cancer, gene_value_mode,p_value_threshold),sep='\t')

    # with open ("/home/wch23/Project/LifeArc/SOX2/result/Sig.Genes/{}_anova_sig_genes.txt".format(target_cancer),'w') as sig_f:
    #     sig_f.write('\n'.join(significant_genes))
    #####
    # draw heatmap with sig genes
    #####
    tcga_cancer_diff_df = dh.load_obj("/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original")

    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')

    tcga_cancer_diff_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'] != 'None']
    tcga_cancer_diff_df = tcga_cancer_diff_df.drop(columns=['ensembl_id', 'ensembl_gene'])

    tcga_cancer_original_order = list(tcga_cancer_diff_df)

    tcga_cancer_diff_sig_gene_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'].isin(significant_genes)]

    print(tcga_cancer_diff_sig_gene_df.head())
    print(cancer_id_df.head())
    cancer_id_df = cancer_id_df.set_index('fullcode')
    tcga_cancer_diff_sig_gene_df = tcga_cancer_diff_sig_gene_df.set_index('gene_symbol')
    tcga_cancer_diff_sig_gene_cancer_type_df = pd.concat([tcga_cancer_diff_sig_gene_df.T, cancer_id_df],axis=1, join='inner' )
    tcga_cancer_diff_sig_gene_cancer_type_df = tcga_cancer_diff_sig_gene_cancer_type_df.T
    print(tcga_cancer_diff_sig_gene_cancer_type_df)
    cancer_type = tcga_cancer_diff_sig_gene_cancer_type_df.loc['Abbreviation']
    print(cancer_type)

    rgb_colors = sns.color_palette("hls",  len(cancer_type.unique()))

    cancer_type_color = dict(zip(cancer_type.unique(), rgb_colors))

    print (cancer_type_color)
    col_colors = cancer_type.map(cancer_type_color)
    #

    print(tcga_cancer_diff_sig_gene_df.info())

    file_dir = '/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}'.format(target_cancer)
    from pathlib import Path
    Path(file_dir).mkdir(parents=True, exist_ok=True)

    tcga_cancer_diff_sig_gene_df.to_csv('/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_diff_heatmap_matrix.tsv'.format(target_cancer),sep='\t')

    g = sns.clustermap(tcga_cancer_diff_sig_gene_df, metric="correlation", cmap="RdBu_r", robust=True, method="average",z_score=0,
                       col_colors=col_colors,xticklabels=False)  # Average is best
    # # g = sns.clustermap(tcga_cancer_diff_sig_gene_df, metric="correlation", cmap="RdBu_r", robust=True, method="average",
    # #                    col_colors=col_colors, xticklabels=False)  # Average is best
    #
    g.savefig('/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_diff_heatmap_single_robust_colColor_norm.png'.format(target_cancer))

    print(g.dendrogram_col.reordered_ind)
    print(g.dendrogram_row.reordered_ind)

    clustred_col = g.dendrogram_col.reordered_ind

    cancer_type_df = tcga_cancer_diff_sig_gene_cancer_type_df.T[['Abbreviation']]
    cancer_type_df.index.name = 'fullcode'
    cancer_type_df = cancer_type_df.reset_index()
    print(cancer_type_df)

    clusterd_cancer_type_df = cancer_type_df.reindex(clustred_col)
    print(clusterd_cancer_type_df)

    clusterd_cancer_type_df.to_csv('/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_diff_heatmap_tcga_clustred_result.csv'.format(target_cancer))


def TCGA_heatmap_pancancer_exist_sigGenes(target_cancer, sigGene_addr):
    # # organism = '9606'
    # # string_node_addr = "/home/wch23/Project/LifeArc/General/data/STRING/{}.protein.links.v11.0.400.nodes.txt".format(
    # #     organism)
    # # with open(string_node_addr) as string_node_f:
    # #     string_node = [x.strip() for x in string_node_f.readlines()]
    #
    # ####
    # # get sig genes from anova reslut
    # ##########
    #
    # # target_cancer = 'PAAD'
    # # target_cancer = "LUSC"
    # gene_value_mode = 'Diff'
    # # anova_result_addr = "/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_{}_vs_other_{}_anova_result_with_ensembl.tsv".format(target_cancer,gene_value_mode)
    # anova_result_addr = "/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_result.csv".format(target_cancer, gene_value_mode)
    # anova_result_df = pd.read_csv(anova_result_addr)
    #
    # # p_value_threshold = 1.0E-150
    #
    # anova_sig_result_df = anova_result_df.loc[anova_result_df['pvalue']<p_value_threshold]
    # print(anova_sig_result_df.shape)
    # print(anova_sig_result_df.head())
    with open(sigGene_addr) as sigGene_f:
        significant_genes = [x.strip() for x in sigGene_f.readlines()]
    # significant_genes = anova_sig_result_df['gene'].tolist()

    print(significant_genes)

    # anova_sig_result_df.to_csv("/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_sig_genes_pv_{}.tsv".format(target_cancer, gene_value_mode,p_value_threshold),sep='\t')

    # with open ("/home/wch23/Project/LifeArc/SOX2/result/Sig.Genes/{}_anova_sig_genes.txt".format(target_cancer),'w') as sig_f:
    #     sig_f.write('\n'.join(significant_genes))
    #####
    # draw heatmap with sig genes
    #####
    tcga_cancer_diff_df = dh.load_obj("/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original")

    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')

    tcga_cancer_diff_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'] != 'None']
    tcga_cancer_diff_df = tcga_cancer_diff_df.drop(columns=['ensembl_id', 'ensembl_gene'])

    tcga_cancer_original_order = list(tcga_cancer_diff_df)

    tcga_cancer_diff_sig_gene_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'].isin(significant_genes)]

    print(tcga_cancer_diff_sig_gene_df.head())
    print(cancer_id_df.head())
    cancer_id_df = cancer_id_df.set_index('fullcode')
    tcga_cancer_diff_sig_gene_df = tcga_cancer_diff_sig_gene_df.set_index('gene_symbol')
    tcga_cancer_diff_sig_gene_cancer_type_df = pd.concat([tcga_cancer_diff_sig_gene_df.T, cancer_id_df],axis=1, join='inner' )
    tcga_cancer_diff_sig_gene_cancer_type_df = tcga_cancer_diff_sig_gene_cancer_type_df.T
    print(tcga_cancer_diff_sig_gene_cancer_type_df)
    cancer_type = tcga_cancer_diff_sig_gene_cancer_type_df.loc['Abbreviation']
    print(cancer_type)

    rgb_colors = sns.color_palette("hls",  len(cancer_type.unique()))

    cancer_type_color = dict(zip(cancer_type.unique(), rgb_colors))

    print (cancer_type_color)
    col_colors = cancer_type.map(cancer_type_color)
    #

    print(tcga_cancer_diff_sig_gene_df.info())

    tcga_cancer_diff_sig_gene_df.to_csv('/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_diff_heatmap_matrix.tsv'.format(target_cancer),sep='\t')

    g = sns.clustermap(tcga_cancer_diff_sig_gene_df, metric="correlation", cmap="RdBu_r", robust=True, method="average",z_score=0,
                       col_colors=col_colors,xticklabels=False)  # Average is best
    # g = sns.clustermap(tcga_cancer_diff_sig_gene_df, metric="correlation", cmap="RdBu_r", robust=True, method="average",
    #                    col_colors=col_colors, xticklabels=False)  # Average is best

    g.savefig('/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_diff_heatmap_single_robust_colColor_norm.png'.format(target_cancer))

    print(g.dendrogram_col.reordered_ind)
    print(g.dendrogram_row.reordered_ind)

    clustred_col = g.dendrogram_col.reordered_ind

    cancer_type_df = tcga_cancer_diff_sig_gene_cancer_type_df.T[['Abbreviation']]
    cancer_type_df.index.name = 'fullcode'
    cancer_type_df = cancer_type_df.reset_index()
    print(cancer_type_df)

    clusterd_cancer_type_df = cancer_type_df.reindex(clustred_col)
    print(clusterd_cancer_type_df)

    clusterd_cancer_type_df.to_csv('/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_diff_heatmap_tcga_clustred_result.csv'.format(target_cancer))

def TCGA_heatmap_pancancer_exist_sigGenes_private(target_cancer, sigGene_addr,tcga_cancer_diff_sig_gene_addr, tcga_heatmap_pancancer_addr, tcga_pancancer_cluster_addr):

    with open(sigGene_addr) as sigGene_f:
        significant_genes = [x.strip() for x in sigGene_f.readlines()]
    # significant_genes = anova_sig_result_df['gene'].tolist()

    print(significant_genes)

    # anova_sig_result_df.to_csv("/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_sig_genes_pv_{}.tsv".format(target_cancer, gene_value_mode,p_value_threshold),sep='\t')

    # with open ("/home/wch23/Project/LifeArc/SOX2/result/Sig.Genes/{}_anova_sig_genes.txt".format(target_cancer),'w') as sig_f:
    #     sig_f.write('\n'.join(significant_genes))
    #####
    # draw heatmap with sig genes
    #####
    tcga_cancer_diff_df = dh.load_obj("/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original")

    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')

    tcga_cancer_diff_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'] != 'None']
    tcga_cancer_diff_df = tcga_cancer_diff_df.drop(columns=['ensembl_id', 'ensembl_gene'])

    tcga_cancer_original_order = list(tcga_cancer_diff_df)

    tcga_cancer_diff_sig_gene_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'].isin(significant_genes)]

    print(tcga_cancer_diff_sig_gene_df.head())
    print(cancer_id_df.head())
    cancer_id_df = cancer_id_df.set_index('fullcode')
    tcga_cancer_diff_sig_gene_df = tcga_cancer_diff_sig_gene_df.set_index('gene_symbol')
    tcga_cancer_diff_sig_gene_cancer_type_df = pd.concat([tcga_cancer_diff_sig_gene_df.T, cancer_id_df],axis=1, join='inner' )
    tcga_cancer_diff_sig_gene_cancer_type_df = tcga_cancer_diff_sig_gene_cancer_type_df.T
    print(tcga_cancer_diff_sig_gene_cancer_type_df)
    cancer_type = tcga_cancer_diff_sig_gene_cancer_type_df.loc['Abbreviation']
    print(cancer_type)

    rgb_colors = sns.color_palette("hls",  len(cancer_type.unique()))

    cancer_type_color = dict(zip(cancer_type.unique(), rgb_colors))

    print (cancer_type_color)
    col_colors = cancer_type.map(cancer_type_color)
    #

    print(tcga_cancer_diff_sig_gene_df.info())

    tcga_cancer_diff_sig_gene_df.to_csv(tcga_cancer_diff_sig_gene_addr,sep='\t')

    # g = sns.clustermap(tcga_cancer_diff_sig_gene_df, metric="correlation", cmap="RdBu_r", robust=True, method="average",z_score=0,
    #                    col_colors=col_colors,xticklabels=False)  # Average is best
    # # g = sns.clustermap(tcga_cancer_diff_sig_gene_df, metric="correlation", cmap="RdBu_r", robust=True, method="average",
    # #                    col_colors=col_colors, xticklabels=False)  # Average is best
    #
    # g.savefig(tcga_heatmap_pancancer_addr)
    #
    # print(g.dendrogram_col.reordered_ind)
    # print(g.dendrogram_row.reordered_ind)
    #
    # clustred_col = g.dendrogram_col.reordered_ind
    #
    # cancer_type_df = tcga_cancer_diff_sig_gene_cancer_type_df.T[['Abbreviation']]
    # cancer_type_df.index.name = 'fullcode'
    # cancer_type_df = cancer_type_df.reset_index()
    # print(cancer_type_df)
    #
    # clusterd_cancer_type_df = cancer_type_df.reindex(clustred_col)
    # print(clusterd_cancer_type_df)
    #
    # clusterd_cancer_type_df.to_csv(tcga_pancancer_cluster_addr)


def change_col_name(matrix_addr, col_name_addr, matrix_new_addr,index_file_addr, column_file_addr):

    matrix_df= pd.read_csv(matrix_addr, sep='\t',index_col=0)

    col_name_df = pd.read_csv(col_name_addr, index_col=0)

    print(list(col_name_df))
    col_name_df['type_code'] = col_name_df.Abbreviation+'-'+col_name_df.fullcode

    # print(col_name_df)
    col_name_change_dict = col_name_df.set_index('fullcode').to_dict()['type_code']
    print(col_name_change_dict)

    new_matrix_df = matrix_df.rename(columns=col_name_change_dict)
    matrix_index = list(new_matrix_df.index.to_list())
    matrix_col = list(new_matrix_df)

    print(len(matrix_index))
    print(len(matrix_col))
    print(matrix_index)
    new_matrix_df.to_csv(matrix_new_addr, sep='\t',index=False, header=False)
    with open(index_file_addr, 'w' ) as index_f:
        index_f.write('\n'.join(matrix_index))
    with open(column_file_addr,'w') as col_f:
        col_f.write('\n'.join(matrix_col))


def main_PAAD():
    target_cancer = 'PAAD'
    p_value_threshold = 1.0E-40
    TCGA_heatmap_pancacer_sig_genes(target_cancer, p_value_threshold)

def main_change_col(cancer_type):

    target_cancer = cancer_type
    matrix_addr = "/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_diff_heatmap_matrix.tsv".format(target_cancer)
    col_name_addr = "/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_diff_heatmap_tcga_clustred_result.csv".format(target_cancer)
    matrix_new_addr=  "/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_{}_diff_heatmap_matrix.tsv".format(target_cancer,target_cancer)
    index_file_addr = "/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_{}_diff_heatmap_matrix_index.tsv".format(target_cancer,target_cancer)
    column_file_addr = "/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_{}_diff_heatmap_matrix_column.tsv".format(
        target_cancer, target_cancer)
    # matrix_addr = "/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_Frank_diff_heatmap_matrix.tsv".format(target_cancer)
    # col_name_addr = "/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_Frank_diff_heatmap_tcga_clustred_result.csv".format(target_cancer)
    # matrix_new_addr=  "/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_Frank_{}_diff_heatmap_matrix_with_cancertype.tsv".format(target_cancer,target_cancer)
    # index_file_addr = "/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_Frank_{}_diff_heatmap_matrix_index.tsv".format(target_cancer,target_cancer)
    # column_file_addr = "/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_Frank_{}_diff_heatmap_matrix_column.tsv".format(
    #     target_cancer, target_cancer)
    change_col_name(matrix_addr, col_name_addr,matrix_new_addr, index_file_addr, column_file_addr)

def main_PAAD_ex():
    target_cancer = 'PAAD'
    # p_value_threshold = 1.0E-40
    sigGenes_addr = "/home/wch23/Project/LifeArc/TCGA/Result/significant_genes_only_original_Diff_from_PAAD_anova_1e-50.csv"
    TCGA_heatmap_pancancer_exist_sigGenes(target_cancer, sigGenes_addr)


def main_LIHC():
    target_cancer = 'LIHC'
    sigGenes_addr = "/home/wch23/Project/LifeArc/ITGB5/result_2/CHC/Sig.Genes/tcga_anova_FDR_result_significant_{}_Diff.csv".format(target_cancer)
    TCGA_heatmap_pancancer_exist_sigGenes(target_cancer, sigGenes_addr)

def main_LUSC():
    target_cancer = 'LUSC'
    sigGenes_addr = "/home/wch23/Project/LifeArc/TCGA/Result/significant_genes_only_original_Diff_from_LUSC_anova_1e-150.csv"
    TCGA_heatmap_pancancer_exist_sigGenes(target_cancer, sigGenes_addr)

def main_CHOL_ex():
    target_cancer = 'CHOL'
    sigGenes_addr = "/home/wch23/Project/LifeArc/ITGB5/result_2/CHC/Sig.Genes/tcga_anova_FDR_result_significant_CHOL_Diff.csv"
    TCGA_heatmap_pancancer_exist_sigGenes(target_cancer, sigGenes_addr)

def main_CHOL():
    target_cancer = 'CHOL'
    p_value_threshold = 1.0E-20
    TCGA_heatmap_pancacer_sig_genes(target_cancer, p_value_threshold)

def main_GBM():
    target_cancer = 'GBM'
    # p_value_threshold = 1.0E-100
    p_value_threshold = 1.0E-80
    gene_value_mode = 'Diff'
    anova_result_addr = "/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_result_filter_min_max.csv".format(
        target_cancer, gene_value_mode)
    TCGA_heatmap_pancacer_sig_genes(target_cancer, p_value_threshold, anova_result_addr)

def main_TCGA():

    cancer_pvalue_dict= {
        # 'GBM' : 1.0E-80,
        # 'LGG' : 1.0E-200,
        # 'LAML': 1.0E-300,
        # 'DLBC': 1.0E-24,
        # 'OV'  : 1.0E-150,
        # 'SARC': 1.0E-70,
        # 'PRAD': 1.0E-50
        # 'KICH' : 1.0E-20,
        # 'KIRC' : 1.0E-90,
        # 'KIRP' : 1.0E-35,
        # 'THCA' : 1.0E-100,
        "BRCA" : 1.0E-50     # to add ERBB2 modified, 1.0E-90


    }
    # target_cancer = 'GBM'
    # p_value_threshold = 1.0E-100
    # p_value_threshold = 1.0E-80
    for target_cancer, p_value_threshold in cancer_pvalue_dict.items():
        gene_value_mode = 'Diff'
        anova_result_addr = "/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_result_filter_min_max.csv".format(
            target_cancer, gene_value_mode)
        TCGA_heatmap_pancacer_sig_genes(target_cancer, p_value_threshold, anova_result_addr)

def main_LUSC_frank():
    target_cancer = "LUSC"

    sigGenes_addr = "/home/wch23/Project/LifeArc/SOX2/result/SOX2_STAT1/Sig.Genes/SOX2_STAT1/Frank_D12_D16_logfc__fdr_0.05.txt"
    tcga_cancer_diff_sig_gene_addr = '/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_Frank_diff_heatmap_matrix.tsv'.format(target_cancer)
    tcga_heatmap_pancancer_addr = '/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_Frank_diff_heatmap_single_robust_colColor_norm.png'.format(target_cancer)
    tcga_pancancer_cluster_addr = '/home/wch23/Project/LifeArc/TCGA/Result/heatmap/{}/TCGA_Frank_diff_heatmap_tcga_clustred_result.csv'.format(
        target_cancer)
    TCGA_heatmap_pancancer_exist_sigGenes_private(target_cancer, sigGenes_addr,tcga_cancer_diff_sig_gene_addr, tcga_heatmap_pancancer_addr, tcga_pancancer_cluster_addr)

if __name__ == '__main__':
    # main()
    # main_PAAD()
    # main_temp()
    # main_LIHC()
    # main_PAAD_ex()
    # main_LUSC()
    # main_CHOL_ex()
    # target_cancer_type = ["GBM","LGG","LAML",'DLBC','OV','SARC','PRAD']
    # target_cancer_type = ["KICH","KIRC","KIRP","THCA"]
    target_cancer_type = ["BRCA"]
    for cancer_type in target_cancer_type:
        main_change_col(cancer_type)
    # main_CHOL()
    # main_GBM()
    # main_TCGA()
    # main_LUSC_frank()
