# Created by woochanghwang at 21/10/2019

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import csv
import toolbox.data_handlers as dh


def get_TCGA_aGene_Diff_matrix_for_all_tumourTypes(tcga_tumour_matrix_addr, tcga_normal_matrix_addr, tumour_types,gene_symbol):
    all_tumour_table = pd.read_table(tcga_tumour_matrix_addr, sep='\t')
    all_normal_table = pd.read_table(tcga_normal_matrix_addr, sep='\t')

    all_tumour_table = all_tumour_table.replace('\n', '', regex=True)
    all_normal_table = all_normal_table.replace('\n', '', regex=True)


    tcga_symbol_entrez = all_tumour_table[['Symbol', 'Entrez', 'new_Symbol']]
    tcga_symbol_entrez['Entrez'] = tcga_symbol_entrez['Entrez'].astype(str)

    # tcga_symbol_entrez['Symbol|Entrez'] = tcga_symbol_entrez[['new_Symbol','Entrez']].apply(lambda x: '|'.join(x), axis = 1)

    tcga_symbol_and_entrez_all_tumours_df = tcga_symbol_entrez[['new_Symbol']]

    selected_gene_index = tcga_symbol_entrez[tcga_symbol_entrez['new_Symbol']== gene_symbol].index.item()
    print("Gene Index:",selected_gene_index)

    for tumour_type in tumour_types:

        filter_tumour_col = [col for col in all_tumour_table if col.startswith(tumour_type)]
        filter_normal_col = [col for col in all_normal_table if col.startswith(tumour_type)]

        tcga_tumour_a_type_df = all_tumour_table[filter_tumour_col]
        tcga_normal_a_type_df = all_normal_table[filter_normal_col]

        tcga_normal_a_type_df['Normal_mean'] = tcga_normal_a_type_df.mean(axis = 1)

        tcga_tumour_a_type_diff_df = tcga_tumour_a_type_df.subtract(tcga_normal_a_type_df.mean(axis = 1),axis = "index")

        tcga_symbol_and_entrez_all_tumours_df = pd.concat((tcga_symbol_and_entrez_all_tumours_df,tcga_tumour_a_type_diff_df),axis = 1)


    # print(tcga_symbol_and_entrez_all_tumours_df.head())
    column_names = pd.Series(tcga_symbol_and_entrez_all_tumours_df.columns[1:],name='Cancer_Type')
    # column_names = column_names[1:]
    # print(column_names)
    column_names = column_names.str.replace('\.\d+', '')    #remove digital in colume names for groupby
    column_names_df = pd.DataFrame(column_names)

    # column_names_df= column_names_df.drop(column_names_df.index[0])
    # print(column_names_df.head())
    # print("row count")
    # print(column_names_df.shape)

    tcga_symbol_and_entrez_all_tumours_df=tcga_symbol_and_entrez_all_tumours_df.set_index('new_Symbol')

    selected_gene_series = tcga_symbol_and_entrez_all_tumours_df.iloc[selected_gene_index]
    selected_gene_df = pd.DataFrame(selected_gene_series)
    # print(selected_gene_df)
    # print("row count")
    # print(selected_gene_df.shape)

    column_names_df.reset_index(inplace=True)
    selected_gene_df.reset_index(inplace=True)

    selected_gene_df = pd.concat([column_names_df['Cancer_Type'], selected_gene_df[gene_symbol]], axis=1)
    # print (selected_gene_df)


    return selected_gene_df

# def draw_TCGA_boxplot(gene_set):
#     gene_type = "_new_genesymbol"
#
#     tcga_tumour_matrix_addr = "/Users/woochanghwang/PycharmProjects/TCGA/data/RNAseqV2/processed.gene.norm/TCGA_all_tumour" + gene_type + ".tsv"
#     tcga_normal_matrix_addr = "/Users/woochanghwang/PycharmProjects/TCGA/data/RNAseqV2/processed.gene.norm/TCGA_all_normal" + gene_type + ".tsv"
#
#     # tcga_tumour_matrix_addr = "../data/tmp_tcga_tumour_data.tsv"
#     # tcga_normal_matrix_addr = "../data/tmp_tcga_normal_data.tsv"
#
#     tumour_types_in_TCGA = get_TCGA_cancer_types_in_both_tumour_normal(tcga_tumour_matrix_addr, tcga_normal_matrix_addr,
#                                                                        gene_type)
#
#     tumour_types_in_TCGA = sorted(tumour_types_in_TCGA)
#     # gene_symbol_list = ['SOX2']  # ['SRC', 'FYN', 'PRKCA', 'VTN', 'MYL12A','EPHA2' ,'SDC1','YES1']
#     gene_symbol_list = gene_set
#     for gene_symbol in gene_symbol_list:
#         TCGA_theGene_df = get_TCGA_aGene_Diff_matrix_for_all_tumourTypes(tcga_tumour_matrix_addr,
#                                                                          tcga_normal_matrix_addr,
#                                                                          tumour_types_in_TCGA, gene_symbol)
#
#         print("TCGA the Gene")
#         print(TCGA_theGene_df.head())
#
#         print("Tumour types:", tumour_types_in_TCGA)
#         TCGA_theGene_df['Cancer_Type'] = TCGA_theGene_df['Cancer_Type'].astype("category")
#         TCGA_theGene_df['Cancer_Type'].cat.set_categories(tumour_types_in_TCGA, inplace=True)
#
#         print(TCGA_theGene_df.head())
#
#         ax = sns.boxplot(x='Cancer_Type', y=gene_symbol, data=TCGA_theGene_df, order=tumour_types_in_TCGA,
#                          showfliers=False);
#         ax.set_title(gene_symbol)
#         ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
#         ax.get_figure().savefig('../result/boxplot/TCGA_boxplot_' + gene_symbol + '.pdf')
#
#         plt.show()

def main():
    tcga_cancer_diff_df = dh.load_obj("/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original")
    with open('/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_types_has_normal.tsv') as cancer_type_f:
        tcga_cancer_types = [x.strip() for x in cancer_type_f.readlines()]
    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')
    # print(cancer_id_df)
    # target_cancer_list=["PAAD"] # pancreas, ITGB5

    # cancer_id_df['Case_Ctrl'] = ['Case' if x in target_cancer_list else 'Ctrl' for x in cancer_id_df['Abbreviation']]
    # print(cancer_id_df)
    #
    tcga_cancer_diff_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'] != 'None']
    tcga_cancer_diff_df = tcga_cancer_diff_df.drop(columns=['ensembl_id', 'ensembl_gene'])

    # print(cancer_id_df)
    # gene_symbol_list = ["AC128709.3", "AL035258.1", "TCF4-AS1", "AC012498.2", "GSTA8P", "AC128709.2", "LINC01932",
    #                     "POU6F2-AS2", "SCGB3A2", "KRT74", "AC022031.1", "LINC01206", "GBP6", "AC134043.2", "SERPINB13",
    #                     "AC012498.1", "ADH7", "SFTPA1",
    #                     "SFTPB"]  # ['SRC', 'FYN', 'PRKCA', 'VTN', 'MYL12A','EPHA2' ,'SDC1','YES1']
    # gene_symbol_list = ["CYTIP","B2M","AZI2","LINC00487","TNFSF12-TNFSF13","AL365203.2","HSPA1B"]
    # gene_symbol_list = ["ITGA1","ITGA2","ITGA2B","ITGA3","ITGA4","ITGA5","ITGA6","ITGA7","ITGA8","ITGA9","ITGA10","ITGA11","ITGAD","ITGAE","ITGAL","ITGAM",
    #                     "ITGAV","ITGAX","ITGB1","ITGB2","ITGB3","ITGB4","ITGB5","ITGB6","ITGB7","ITGB8","ITGBL1"]
    # gene_symbol_list = ["SMG8","CDH12P1",'AC131392.1',]
    # gene_symbol_list = ["DLGAP1","DTNB","BHLHE40","PXN","CYLD","RYBP","GSC","LEFTY2","BMP7","NRP2","BMP4","OTX2","PRDM1","DENND2A","RAD51C","JARID2","KLF9","SOCS3","TGIF1","KANK1","MKRN1","FOS","DNMT3A","STAB2","RABIF","SNCG","ZIC3","FOXD3"]
    # gene_symbol_list= ["SOX2"]
    gene_symbol_list = ["FGFR1","FGFR2", "CDK1", "ABL1", "VDR", "MTOR", "PARP1", "KLF5"]
    cancer_type = 'ESCA'
    cancer_id_df = cancer_id_df.set_index('fullcode')

    anova_result_df = pd.read_csv("/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_Diff_anova_result_filter_min_max.csv".format(cancer_type),index_col=0)

    for gene_symbol in gene_symbol_list:

        aGene_TCGA_df = tcga_cancer_diff_df.loc[tcga_cancer_diff_df['gene_symbol']==gene_symbol]
        # aGene_TCGA_df = aGene_TCGA_df
        aGene_TCGA_df = aGene_TCGA_df.set_index('gene_symbol')
        aGene_TCGA_df = aGene_TCGA_df.T
        # aGene_TCGA_df = aGene_TCGA_df.reset_index()

        print("Gene:", gene_symbol )
        print(aGene_TCGA_df.head())

        aGene_TCGA_cancer_id_df = pd.concat([aGene_TCGA_df, cancer_id_df], join='inner', axis=1)
        # aGene_TCGA_cancer_id_df = pd.merge(aGene_TCGA_df,cancer_id_df,left_on='index', right_on='fullcode')
        print(aGene_TCGA_cancer_id_df)

        anova_p_value = anova_result_df.loc[gene_symbol].values[0]
        print(anova_p_value)

        aGene_TCGA_cancer_id_df['Abbreviation'] = aGene_TCGA_cancer_id_df['Abbreviation'].astype('category')
        cancer_types = aGene_TCGA_cancer_id_df['Abbreviation'].tolist()
        cancer_types = list(set(cancer_types))
        cancer_types = sorted(cancer_types)


        pkmn_type_colors = ['#78C850',  # Grass
                            '#F08030',  # Fire
                            '#6890F0',  # Water
                            '#A8B820',  # Bug
                            '#A8A878',  # Normal
                            '#A040A0',  # Poison
                            '#F8D030',  # Electric
                            '#E0C068',  # Ground
                            '#EE99AC',  # Fairy
                            '#C03028',  # Fighting
                            '#F85888',  # Psychic
                            '#B8A038',  # Rock
                            '#705898',  # Ghost
                            '#98D8D8',  # Ice
                            '#7038F8',  # Dragon
                            ]

        # ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types,
        #                  showfliers=False);
        # ax.set_title(gene_symbol)
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        # ax.get_figure().savefig('/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}/TCGA_boxplot_{}_{}_original.pdf'.format(cancer_type,gene_symbol, anova_p_value))

        ## for other type boxplot

        plt.figure(figsize=(10,5))
        # ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types, palette=pkmn_type_colors,
        #                  showfliers=False);
        ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types,
                         # palette=pkmn_type_colors,
                         color='w',
                         linewidth=1.5,
                         showfliers=False
                         );
        # ax = sns.swarmplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types, color=".25")
        ax.set_title(gene_symbol)
        ax.set_xlabel('TCGA')
        ax.set_ylabel('Tumour vs Normal')
        ax.xaxis.grid(True)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

        #####################
        ## Set xlabe color
        ####################
        # colors = ['k','k','k','r','k','b','k','b','r','k',
        #           'k','k','k','k','r','k','r','r','k','k',
        #           'k','k','b','k','k','b','b','k','k','k','k']
        # for xtick, color in zip(ax.get_xticklabels(), colors):
        #     xtick.set_color(color)


        ################
        # this is for making color brighter
        #################
        # for patch in ax.artists:
        #     r, g, b, a = patch.get_facecolor()
        #     patch.set_facecolor((r, g, b, .3))
        #########################

        plt.setp(ax.artists, edgecolor='k', facecolor='w')
        plt.setp(ax.lines, color='k')
        file_dir = '/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}'.format(cancer_type)
        from pathlib import Path
        Path(file_dir).mkdir(parents=True, exist_ok=True)

        ax.get_figure().savefig(
            '/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}/TCGA_boxplot_{}_{}_original_new.pdf'.format(cancer_type,
                                                                                                            gene_symbol,
                                                                                                            anova_p_value))

        plt.show()


def main_for_combine_genes():
    tcga_cancer_diff_df = dh.load_obj("/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original")
    with open('/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_types_has_normal.tsv') as cancer_type_f:
        tcga_cancer_types = [x.strip() for x in cancer_type_f.readlines()]
    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')
    # print(cancer_id_df)
    # target_cancer_list=["PAAD"] # pancreas, ITGB5

    # cancer_id_df['Case_Ctrl'] = ['Case' if x in target_cancer_list else 'Ctrl' for x in cancer_id_df['Abbreviation']]
    # print(cancer_id_df)
    #
    tcga_cancer_diff_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'] != 'None']
    tcga_cancer_diff_df = tcga_cancer_diff_df.drop(columns=['ensembl_id', 'ensembl_gene'])

    # print(cancer_id_df)
    # gene_symbol_list = ["AC128709.3", "AL035258.1", "TCF4-AS1", "AC012498.2", "GSTA8P", "AC128709.2", "LINC01932",
    #                     "POU6F2-AS2", "SCGB3A2", "KRT74", "AC022031.1", "LINC01206", "GBP6", "AC134043.2", "SERPINB13",
    #                     "AC012498.1", "ADH7", "SFTPA1",
    #                     "SFTPB"]  # ['SRC', 'FYN', 'PRKCA', 'VTN', 'MYL12A','EPHA2' ,'SDC1','YES1']
    # gene_symbol_list = ["CYTIP","B2M","AZI2","LINC00487","TNFSF12-TNFSF13","AL365203.2","HSPA1B"]
    # gene_symbol_list = ["ITGA1","ITGA2","ITGA2B","ITGA3","ITGA4","ITGA5","ITGA6","ITGA7","ITGA8","ITGA9","ITGA10","ITGA11","ITGAD","ITGAE","ITGAL","ITGAM",
    #                     "ITGAV","ITGAX","ITGB1","ITGB2","ITGB3","ITGB4","ITGB5","ITGB6","ITGB7","ITGB8","ITGBL1"]
    # gene_symbol_list = ["TMEM52","EPB41L4B","AC011754.1","RBPJL","AC096633.1","PNLIP","CELP","LHFPL5","AC092535.1","TMED6"]
    # gene_symbol_list = ["DLGAP1","DTNB","BHLHE40","PXN","CYLD","RYBP","GSC","LEFTY2","BMP7","NRP2","BMP4","OTX2","PRDM1","DENND2A","RAD51C","JARID2","KLF9","SOCS3","TGIF1","KANK1","MKRN1","FOS","DNMT3A","STAB2","RABIF","SNCG","ZIC3","FOXD3"]
    gene_symbol_list= ["SOX2","STAT1"]
    gene_symbol = '_'.join(gene_symbol_list)
    cancer_type = 'LUSC'
    cancer_id_df = cancer_id_df.set_index('fullcode')

    anova_result_df = pd.read_csv("/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_Diff_anova_result.csv".format(cancer_type),index_col=0)

    aGene_TCGA_df = tcga_cancer_diff_df.loc[tcga_cancer_diff_df['gene_symbol'].isin(gene_symbol_list)]
    # aGene_TCGA_df = aGene_TCGA_df
    aGene_TCGA_df = aGene_TCGA_df.set_index('gene_symbol')
    aGene_TCGA_df = aGene_TCGA_df.T
    # aGene_TCGA_df[gene_symbol] = aGene_TCGA_df.sum(axis=1)
    aGene_TCGA_df[gene_symbol] = aGene_TCGA_df.mean(axis=1)

    # aGene_TCGA_df = aGene_TCGA_df.reset_index()
    aGene_TCGA_df = aGene_TCGA_df.drop(columns=gene_symbol_list)
    print("Gene:", gene_symbol)
    print(aGene_TCGA_df.head())

    aGene_TCGA_cancer_id_df = pd.concat([aGene_TCGA_df, cancer_id_df], join='inner', axis=1)
    # aGene_TCGA_cancer_id_df = pd.merge(aGene_TCGA_df,cancer_id_df,left_on='index', right_on='fullcode')
    print(aGene_TCGA_cancer_id_df)

    # anova_p_value = anova_result_df.loc[gene_symbol].values[0]
    # print(anova_p_value)

    aGene_TCGA_cancer_id_df['Abbreviation'] = aGene_TCGA_cancer_id_df['Abbreviation'].astype('category')
    cancer_types = aGene_TCGA_cancer_id_df['Abbreviation'].tolist()
    cancer_types = list(set(cancer_types))
    cancer_types = sorted(cancer_types)

    pkmn_type_colors = ['#78C850',  # Grass
                        '#F08030',  # Fire
                        '#6890F0',  # Water
                        '#A8B820',  # Bug
                        '#A8A878',  # Normal
                        '#A040A0',  # Poison
                        '#F8D030',  # Electric
                        '#E0C068',  # Ground
                        '#EE99AC',  # Fairy
                        '#C03028',  # Fighting
                        '#F85888',  # Psychic
                        '#B8A038',  # Rock
                        '#705898',  # Ghost
                        '#98D8D8',  # Ice
                        '#7038F8',  # Dragon
                        ]

    # ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types,
    #                  showfliers=False);
    # ax.set_title(gene_symbol)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # ax.get_figure().savefig('/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}/TCGA_boxplot_{}_{}_original.pdf'.format(cancer_type,gene_symbol, anova_p_value))

    ## for other type boxplot

    plt.figure(figsize=(10, 5))
    # ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types, palette=pkmn_type_colors,
    #                  showfliers=False);
    ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types,
                     # palette=pkmn_type_colors,
                     color='w',
                     linewidth=1.5,
                     showfliers=False
                     );
    # ax = sns.swarmplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types, color=".25")
    ax.set_title(gene_symbol)
    ax.set_xlabel('TCGA')
    ax.set_ylabel('Tumour vs Normal')
    ax.xaxis.grid(True)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

    ################
    # this is for making color brighter
    #################
    # for patch in ax.artists:
    #     r, g, b, a = patch.get_facecolor()
    #     patch.set_facecolor((r, g, b, .3))
    #########################

    plt.setp(ax.artists, edgecolor='k', facecolor='w')
    plt.setp(ax.lines, color='k')
    ax.get_figure().savefig(
        '/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}/TCGA_boxplot_{}_original_new_v2.pdf'.format(cancer_type,
                                                                                                           gene_symbol
                                                                                                           # anova_p_value
                                                                                                           ))

    plt.show()


def main_for_cancer(miRNA_list):
    # cancer_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_xena_to_original_100.csv"
    # normal_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_gepia_normal_xena_to_original.csv"

    # tcga_cancer_df = pd.read_csv(cancer_matrix_addr)

    tcga_cancer_df = dh.load_obj("/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_xena_to_original")
    print(tcga_cancer_df.head())

    print(tcga_cancer_df)

    with open('/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_types_has_normal.tsv') as cancer_type_f:
        tcga_cancer_types = [x.strip() for x in cancer_type_f.readlines()]
    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')
    # print(cancer_id_df)
    # target_cancer_list=["PAAD"] # pancreas, ITGB5

    # cancer_id_df['Case_Ctrl'] = ['Case' if x in target_cancer_list else 'Ctrl' for x in cancer_id_df['Abbreviation']]
    # print(cancer_id_df)
    #
    tcga_cancer_df = tcga_cancer_df[tcga_cancer_df['gene_symbol'] != 'None']
    tcga_cancer_df = tcga_cancer_df.drop(columns=['ensembl_id', 'ensembl_gene'])
    # tcga_cancer_df = tcga_cancer_df.drop(columns=['gene_symbol', 'ensembl_gene'])

    # print(cancer_id_df)

    # gene_symbol_list = ["SBNO2", "AL031587.5", "SMC6", "DUSP12", "CBWD6", "SLC39A2", "IVD", "GDNF-AS1", "RUSC1-AS1",
    #                     "PABPC4L", "SCARB2", "FAT3", "AC127024.5", "EXOSC8"]
    gene_symbol_list = miRNA_list
    cancer_type = "Cancer"
    cancer_id_df = cancer_id_df.set_index('fullcode')


    gene_list_from_TCGA = tcga_cancer_df['gene_symbol'].to_list()
    # gene_symbols_in_TCGA = list(set(gene_list_from_TCGA)&set(gene_symbol_list))
    gene_symbols_in_TCGA = ["MIR3648-2", "AL513534.2", "MIR3648-1", "MIR6753", "AC099677.4", "EXOSC8", "AC127024.5", "FAT3",
                        "SCARB2", "PABPC4L", "RUSC1-AS1", "GDNF-AS1", "IVD", "SLC39A2", "CBWD6", "DUSP12", "SMC6",
                        "AL031587.5", "SBNO2"]

    # print(len(gene_symbol_list), len(gene_symbols_in_TCGA))
    for gene_symbol in gene_symbols_in_TCGA:

        aGene_TCGA_df = tcga_cancer_df.loc[tcga_cancer_df['gene_symbol'] == gene_symbol]
        # aGene_TCGA_df = aGene_TCGA_df
        aGene_TCGA_df = aGene_TCGA_df.set_index('gene_symbol')
        # aGene_TCGA_df = tcga_cancer_df.loc[tcga_cancer_df['ensembl_id'] == gene_symbol]
        # # aGene_TCGA_df = aGene_TCGA_df
        # aGene_TCGA_df = aGene_TCGA_df.set_index('ensembl_id')

        aGene_TCGA_df = aGene_TCGA_df.T
        # aGene_TCGA_df = aGene_TCGA_df.reset_index()

        print("Gene:", gene_symbol)
        print(aGene_TCGA_df.head())

        aGene_TCGA_cancer_id_df = pd.concat([aGene_TCGA_df, cancer_id_df], join='inner', axis=1)
        # aGene_TCGA_cancer_id_df = pd.merge(aGene_TCGA_df,cancer_id_df,left_on='index', right_on='fullcode')
        print(aGene_TCGA_cancer_id_df)

        aGene_TCGA_cancer_id_df['Abbreviation'] = aGene_TCGA_cancer_id_df['Abbreviation'].astype('category')
        cancer_types = aGene_TCGA_cancer_id_df['Abbreviation'].tolist()
        cancer_types = list(set(cancer_types))
        cancer_types = sorted(cancer_types)

        ## for other type boxplot

        plt.figure(figsize=(10, 5))
        # ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types, palette=pkmn_type_colors,
        #                  showfliers=False);
        ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types,
                         # palette=pkmn_type_colors,
                         color='w',
                         linewidth=1.5,
                         showfliers=False
                         );
        # ax = sns.swarmplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types, color=".25")
        ax.set_title(gene_symbol)
        ax.set_xlabel('TCGA')
        ax.set_ylabel('Tumour FC')
        ax.xaxis.grid(True)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

        plt.setp(ax.artists, edgecolor='k', facecolor='w')
        plt.setp(ax.lines, color='k')
        file_dir = '/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}'.format(cancer_type)
        from pathlib import Path
        Path(file_dir).mkdir(parents=True, exist_ok=True)

        ax.get_figure().savefig(
            '/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}/expressed/TCGA_boxplot_{}_original_for_cancer.pdf'.format(
                cancer_type,
                gene_symbol))

    # plt.show()



def main_for_normal(miRNA_list):
    # cancer_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_xena_to_original_100.csv"
    # normal_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_gepia_normal_xena_to_original.csv"

    # tcga_cancer_df = pd.read_csv(cancer_matrix_addr)

    tcga_cancer_df = dh.load_obj("/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_gepia_normal_xena_to_original")
    print(tcga_cancer_df.head())

    print(tcga_cancer_df)

    with open('/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_types_has_normal.tsv') as cancer_type_f:
        tcga_cancer_types = [x.strip() for x in cancer_type_f.readlines()]
    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/normal_id_info.tsv"
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')
    # print(cancer_id_df)
    # target_cancer_list=["PAAD"] # pancreas, ITGB5

    # cancer_id_df['Case_Ctrl'] = ['Case' if x in target_cancer_list else 'Ctrl' for x in cancer_id_df['Abbreviation']]
    # print(cancer_id_df)
    #
    tcga_cancer_df = tcga_cancer_df[tcga_cancer_df['gene_symbol'] != 'None']
    tcga_cancer_df = tcga_cancer_df.drop(columns=['ensembl_id', 'ensembl_gene'])
    # tcga_cancer_df = tcga_cancer_df.drop(columns=['gene_symbol', 'ensembl_gene'])
    # print(cancer_id_df)

    # gene_symbol_list = ["SBNO2", "AL031587.5", "SMC6", "DUSP12", "CBWD6", "SLC39A2", "IVD", "GDNF-AS1", "RUSC1-AS1",
    #                     "PABPC4L", "SCARB2", "FAT3", "AC127024.5", "EXOSC8"]
    gene_symbol_list = miRNA_list
    cancer_type = "Normal"

    cancer_id_df = cancer_id_df.set_index('fullcode')

    print(len(list(tcga_cancer_df)), len(list(cancer_id_df)))

    gene_list_from_TCGA = tcga_cancer_df['gene_symbol'].to_list()
    gene_symbols_in_TCGA = list(set(gene_list_from_TCGA)&set(gene_symbol_list))

    for gene_symbol in gene_symbols_in_TCGA:

        aGene_TCGA_df = tcga_cancer_df.loc[tcga_cancer_df['gene_symbol'] == gene_symbol]
        # aGene_TCGA_df = aGene_TCGA_df
        aGene_TCGA_df = aGene_TCGA_df.set_index('gene_symbol')
        # aGene_TCGA_df = tcga_cancer_df.loc[tcga_cancer_df['ensembl_id'] == gene_symbol]
        # # aGene_TCGA_df = aGene_TCGA_df
        # aGene_TCGA_df = aGene_TCGA_df.set_index('ensembl_id')

        aGene_TCGA_df = aGene_TCGA_df.T
        # aGene_TCGA_df = aGene_TCGA_df.reset_index()

        print("Gene:", gene_symbol)
        print(aGene_TCGA_df.head())

        aGene_TCGA_cancer_id_df = pd.concat([aGene_TCGA_df, cancer_id_df], join='inner', axis=1)
        # aGene_TCGA_cancer_id_df = pd.merge(aGene_TCGA_df,cancer_id_df,left_on='index', right_on='fullcode')
        print(aGene_TCGA_cancer_id_df)

        aGene_TCGA_cancer_id_df['Abbreviation'] = aGene_TCGA_cancer_id_df['Abbreviation'].astype('category')
        cancer_types = aGene_TCGA_cancer_id_df['Abbreviation'].tolist()
        cancer_types = list(set(cancer_types))
        cancer_types = sorted(cancer_types)

        ## for other type boxplot

        plt.figure(figsize=(10, 5))
        # ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types, palette=pkmn_type_colors,
        #                  showfliers=False);
        ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types,
                         # palette=pkmn_type_colors,
                         color='w',
                         linewidth=1.5,
                         showfliers=False
                         );
        # ax = sns.swarmplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types, color=".25")
        ax.set_title(gene_symbol)
        ax.set_xlabel('TCGA')
        ax.set_ylabel('Tumour FC')
        ax.xaxis.grid(True)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

        plt.setp(ax.artists, edgecolor='k', facecolor='w')
        plt.setp(ax.lines, color='k')
        file_dir = '/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}'.format(cancer_type)
        from pathlib import Path
        Path(file_dir).mkdir(parents=True, exist_ok=True)

        ax.get_figure().savefig(
            '/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}/TCGA_boxplot_{}_original_for_cancer.pdf'.format(
                cancer_type,
                gene_symbol))

    # plt.show()

# def main_pre():
#     gene_type = "_new_genesymbol"
#
#     tcga_tumour_matrix_addr = "/Users/woochanghwang/PycharmProjects/TCGA/data/RNAseqV2/processed.gene.norm/TCGA_all_tumour" + gene_type + ".tsv"
#     tcga_normal_matrix_addr = "/Users/woochanghwang/PycharmProjects/TCGA/data/RNAseqV2/processed.gene.norm/TCGA_all_normal" + gene_type + ".tsv"
#
#     # tcga_tumour_matrix_addr = "../data/tmp_tcga_tumour_data.tsv"
#     # tcga_normal_matrix_addr = "../data/tmp_tcga_normal_data.tsv"
#
#     tumour_types_in_TCGA = get_TCGA_cancer_types_in_both_tumour_normal(tcga_tumour_matrix_addr, tcga_normal_matrix_addr,
#                                                                        gene_type)
#
#     tumour_types_in_TCGA = sorted(tumour_types_in_TCGA)
#     cancer_type = 'UCSC'
#     gene_symbol_list = ["AC128709.3","AL035258.1","TCF4-AS1","AC012498.2","GSTA8P","AC128709.2","LINC01932","POU6F2-AS2","SCGB3A2","KRT74","AC022031.1","LINC01206","GBP6","AC134043.2","SERPINB13","AC012498.1","ADH7","SFTPA1","SFTPB"] #['SRC', 'FYN', 'PRKCA', 'VTN', 'MYL12A','EPHA2' ,'SDC1','YES1']
#     for gene_symbol in gene_symbol_list:
#         TCGA_theGene_df = get_TCGA_aGene_Diff_matrix_for_all_tumourTypes(tcga_tumour_matrix_addr, tcga_normal_matrix_addr,
#                                                            tumour_types_in_TCGA,gene_symbol)
#
#         print("TCGA the Gene")
#         print(TCGA_theGene_df.head())
#
#         print("Tumour types:", tumour_types_in_TCGA)
#         TCGA_theGene_df['Cancer_Type'] = TCGA_theGene_df['Cancer_Type'].astype("category")
#         TCGA_theGene_df['Cancer_Type'].cat.set_categories(tumour_types_in_TCGA,inplace = True)
#
#         print(TCGA_theGene_df.head())
#
#         ax = sns.boxplot(x='Cancer_Type', y=gene_symbol, data=TCGA_theGene_df, order = tumour_types_in_TCGA, showfliers=False);
#         ax.set_title(gene_symbol)
#         ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
#         ax.get_figure().savefig('../result/boxplot/{}/TCGA_boxplot_{}.pdf'.format(cancer_type,gene_symbol))
#
#         plt.show()


def read_miRNA_gene_set():
    miRNA_gene_addr= "/mnt/raid0_data/MTI_DATA/miRNA/miRNA_enslOD_mart_export_20200304.txt"
    miRNA_gene_df = pd.read_csv(miRNA_gene_addr, sep='\t')
    print(miRNA_gene_df)
    print(list(miRNA_gene_df))
    return miRNA_gene_df['Gene name'].to_list()

def main_for_nonType():
    tcga_cancer_diff_df = dh.load_obj("/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original")
    with open('/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_types_has_normal.tsv') as cancer_type_f:
        tcga_cancer_types = [x.strip() for x in cancer_type_f.readlines()]
    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')
    # print(cancer_id_df)
    # target_cancer_list=["PAAD"] # pancreas, ITGB5

    # cancer_id_df['Case_Ctrl'] = ['Case' if x in target_cancer_list else 'Ctrl' for x in cancer_id_df['Abbreviation']]
    # print(cancer_id_df)
    #
    tcga_cancer_diff_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol'] != 'None']
    tcga_cancer_diff_df = tcga_cancer_diff_df.drop(columns=['ensembl_id', 'ensembl_gene'])


    # gene_symbol_list = ["MIR3648-2","AL513534.2","MIR3648-1","MIR6753","AC099677.4","EXOSC8","AC127024.5","FAT3","SCARB2","PABPC4L","RUSC1-AS1","GDNF-AS1","IVD","SLC39A2","CBWD6","DUSP12","SMC6","AL031587.5","SBNO2"]
    gene_symbol_list = ["ERBB2"]
    cancer_type = 'All'
    cancer_id_df = cancer_id_df.set_index('fullcode')


    for gene_symbol in gene_symbol_list:

        aGene_TCGA_df = tcga_cancer_diff_df.loc[tcga_cancer_diff_df['gene_symbol']==gene_symbol]
        # aGene_TCGA_df = aGene_TCGA_df
        aGene_TCGA_df = aGene_TCGA_df.set_index('gene_symbol')
        aGene_TCGA_df = aGene_TCGA_df.T
        # aGene_TCGA_df = aGene_TCGA_df.reset_index()

        print("Gene:", gene_symbol )
        print(aGene_TCGA_df.head())

        aGene_TCGA_cancer_id_df = pd.concat([aGene_TCGA_df, cancer_id_df], join='inner', axis=1)
        # aGene_TCGA_cancer_id_df = pd.merge(aGene_TCGA_df,cancer_id_df,left_on='index', right_on='fullcode')
        print(aGene_TCGA_cancer_id_df)

        aGene_TCGA_cancer_id_df['Abbreviation'] = aGene_TCGA_cancer_id_df['Abbreviation'].astype('category')
        cancer_types = aGene_TCGA_cancer_id_df['Abbreviation'].tolist()
        cancer_types = list(set(cancer_types))
        cancer_types = sorted(cancer_types)




        ## for other type boxplot

        plt.figure(figsize=(10,5))
        # ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types, palette=pkmn_type_colors,
        #                  showfliers=False);
        ax = sns.boxplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types,
                         # palette=pkmn_type_colors,
                         color='w',
                         linewidth=1.5,
                         showfliers=False
                         );
        # ax = sns.swarmplot(x='Abbreviation', y=gene_symbol, data=aGene_TCGA_cancer_id_df, order=cancer_types, color=".25")
        ax.set_title(gene_symbol)
        ax.set_xlabel('TCGA')
        ax.set_ylabel('Tumour vs Normal')
        ax.xaxis.grid(True)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)



        plt.setp(ax.artists, edgecolor='k', facecolor='w')
        plt.setp(ax.lines, color='k')
        file_dir = '/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}'.format(cancer_type)
        from pathlib import Path
        Path(file_dir).mkdir(parents=True, exist_ok=True)

        ax.get_figure().savefig(
            '/home/wch23/Project/LifeArc/TCGA/Result/boxplot/{}/TCGA_boxplot_{}_original_new_v2.pdf'.format(cancer_type,
                                                                                                            gene_symbol))

        plt.show()

if __name__ == '__main__':
    main()
    # main_for_nonType()
    # main_for_combine_genes()
    # miRNA_list = read_miRNA_gene_set()
    # main_for_cancer(miRNA_list)
    # main_for_normal(miRNA_list)