# Created by woochanghwang at 2019-06-18

# Created by woochanghwang at 2019-06-13

import pandas as pd
import toolbox.data_handler as dah

def broad_depmap_handle():
    '''
    Col: Gene Symbol(Entrez)
    Row: Cell line
    :return:
    '''
    sample_info_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/DepMap/broad_ccle_depmap_info.csv"
    depmap_result_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/DepMap/broad_LOF_score.csv"

    sample_info_df = pd.read_csv(sample_info_addr)
    depmap_result_df = pd.read_csv(depmap_result_addr)

    print(list(sample_info_df))
    # print(list(depmap_result_df))

    # print(sample_info_df.head())

    sample_id_name_df = sample_info_df[['DepMap_ID','stripped_cell_line_name']]

    # print(sample_id_name_df)
    depmap_result_df = depmap_result_df.set_index('ID')
    sample_id_name_df = sample_id_name_df.set_index('DepMap_ID')

    depmap_sample_df = pd.concat([depmap_result_df,sample_id_name_df],axis=1, join='inner')

    # print(depmap_sample_df)

    depmap_sample_df = depmap_sample_df.reset_index()

    # rearrange_column
    depmap_sample_df_col = depmap_sample_df.columns.tolist()
    depmap_sample_df_col = [depmap_sample_df_col[0]]+[depmap_sample_df_col[-1]]+depmap_sample_df_col[1:-1]

    depmap_sample_df = depmap_sample_df[depmap_sample_df_col]

    depmap_sample_df.rename(columns = lambda x:x.split(' ')[0].strip(), inplace=True)

    depmap_sample_T_df = depmap_sample_df.T

    # set cellline name as columne name

    depmap_sample_T_df.rename(columns = depmap_sample_T_df.iloc[1], inplace=True)

    depmap_sample_T_df = depmap_sample_T_df.drop(['stripped_cell_line_name','index'])

    print(depmap_sample_T_df)

    depmap_sample_T_df.to_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/DepMap/broad_ccle_LOF_score.csv")

def sanger_depmap_handle():
    '''
    Col: Cell line
    Row: Gene
    :return:
    Cell line rename same as Broad cellline
    '''
    depmap_result_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/DepMap/sanger_LOF_score.tsv"
    sanger_depmap_df = pd.read_csv(depmap_result_addr,sep='\t')

    #rename columns

    sanger_cols = sanger_depmap_df.columns.tolist()
    # sanger_cols = [x.replace('-','').upper() for x in sanger_cols]

    sanger_depmap_df.rename(columns = lambda x:x.replace('-','').upper(), inplace=True)

    ########
    ## conver value to -(value)
    ## #######
    sanger_depmap_df = sanger_depmap_df.set_index('GENE')
    sanger_depmap_df = sanger_depmap_df.apply(lambda x:x*-1)

    sanger_depmap_df = sanger_depmap_df.reset_index()
    sanger_depmap_df.to_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/DepMap/sanger_ccle_LOF_score.csv")

def sanger_depmap_selected_celllines(selected_celllines):

    sanger_depmap_df = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/DepMap/sanger_ccle_LOF_score.csv")

    sanger_depmap_df = sanger_depmap_df.set_index('GENE')

    sanger_depmap_selected_celllines_df = sanger_depmap_df[sanger_depmap_df.columns & selected_celllines]

    # print(sanger_depmap_selected_celllines_df)

    return sanger_depmap_selected_celllines_df

def broad_depmap_selected_celllines(selected_celllines=None):
    broad_depmap_df = pd.read_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/DepMap/broad_ccle_LOF_score.csv", index_col=0)
    # broad_depmap_df = broad_depmap_df.set_index('ID')

    if selected_celllines == None:
        return broad_depmap_df
    else:
        broad_depmap_selected_celllines_df = broad_depmap_df[broad_depmap_df.columns & selected_celllines]

    # print(broad_depmap_selected_celllines_df)

    return broad_depmap_selected_celllines_df

def calc_significance_of_geneSet_A_cellline(geneSet, cellline, depmap_df,boot_num):
    '''

    :param geneSet:
    :param cellline:
    :param depmap_df:
    :return: significance of selected geneset (LOF)
    '''
    import random
    gene_length = len(geneSet)

    whole_gene = list(depmap_df.index)

    # cols = depmap_df.columns
    # depmap_df = depmap_df[cols].astype('float')

    guide_depmap_df = depmap_df.loc[geneSet]

    guide_depmap_df = guide_depmap_df.dropna()
    # cols = guide_depmap_df.columns
    # guide_depmap_df = guide_depmap_df[cols].astype('float')
    # print(guide_depmap_df[cellline])
    guide_depmap_significat_df = guide_depmap_df[guide_depmap_df[cellline] >= 0]


    effect_score = len(guide_depmap_significat_df)/gene_length
    print("effect score:",effect_score)
    bootstrap = []
    for i in range(boot_num):
        random_genes = random.choices(whole_gene,k=gene_length)
        # print(random_genes)
        random_depmap_df = depmap_df.loc[random_genes]
        random_depmap_significat_df = random_depmap_df[random_depmap_df[cellline]>=0]
        # print(random_depmap_significat_df)
        effect_geens_length = len(random_depmap_significat_df)
        bootstrap.append(effect_geens_length/gene_length)

    significant_list = [i for i in bootstrap if i >= effect_score]

    print(significant_list)
    print("p-value", 1-len(significant_list)/boot_num)

def get_cellline_subtype(selected_celllines):

    broad_ccle_info_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/DepMap/broad_ccle_depmap_info.csv"
    broad_ccle_info_df = pd.read_csv(broad_ccle_info_addr)

    broad_ccle_info_df = broad_ccle_info_df.set_index('stripped_cell_line_name')

    print(list(broad_ccle_info_df))
    broad_ccle_info_selected_df = broad_ccle_info_df.loc[selected_celllines]

    broad_ccle_info_selected_df = broad_ccle_info_selected_df.reset_index()
    broad_cellline_subtype_df = broad_ccle_info_selected_df[['stripped_cell_line_name','disease_sub_subtype']]
    broad_cellline_subtype_df = broad_cellline_subtype_df.set_index('stripped_cell_line_name')

    ## Fill NA
    broad_cellline_subtype_df = dah.recode_empty_cells(broad_cellline_subtype_df,'NA')

    print(broad_cellline_subtype_df)

    return broad_cellline_subtype_df
def broad_ccle_celllines_selected_disease(disease_types_dict):

    broad_ccle_info_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/DepMap/broad_ccle_depmap_info.csv"
    broad_ccle_info_df = pd.read_csv(broad_ccle_info_addr)

    ## add filter
    for key, value in disease_types_dict.items():
        broad_ccle_info_df = broad_ccle_info_df[broad_ccle_info_df[key]==value]

    # print(broad_ccle_info_df)
    # print(len(broad_ccle_info_df))

    return list(broad_ccle_info_df['stripped_cell_line_name'])

def make_binary_up_down(lof):
    if lof < 0 : return 0
    else:   return 1

def make_binary_up_down_broad(lof):
    if lof < -0.5 : return 0
    else:   return 1

def broad_sanger_effective_genes_on_selected_celllines(selected_celllines):
    '''
    binary value dataframe
    row : genes
    binary:
        0 : effective
        1 : non-effective
    :param selected_celllines:
    :return:
    '''
    sanger_depmap_df = sanger_depmap_selected_celllines(selected_celllines)
    broad_depmap_df = broad_depmap_selected_celllines(selected_celllines)

    sanger_depmap_df['mean'] = sanger_depmap_df.mean(axis=1)
    broad_depmap_df['mean'] = broad_depmap_df.mean(axis=1)

    sanger_depmap_df['binary_sanger'] = sanger_depmap_df['mean'].apply(make_binary_up_down)
    broad_depmap_df['binary_broad'] = broad_depmap_df['mean'].apply(make_binary_up_down_broad)

    print(sanger_depmap_df)
    print(broad_depmap_df)

    sanger_broad_depmap_binary_df = pd.concat([sanger_depmap_df['binary_sanger'], broad_depmap_df['binary_broad']],
                                              axis=1, join='inner')
    print(sanger_broad_depmap_binary_df)

    print(len(sanger_depmap_df[sanger_depmap_df['binary_sanger'] == 0]))
    print(len(broad_depmap_df[broad_depmap_df['binary_broad'] == 0]))

    return sanger_broad_depmap_binary_df


def get_broad_depmap_selected_disease_genes(disease_type_dict, genes):
    selected_celllines = broad_ccle_celllines_selected_disease(disease_type_dict)

    # print(selected_celllines)
    # cellines_subtypes = get_cellline_subtype(selected_celllines)
    broad_depmap_df = broad_depmap_selected_celllines(selected_celllines)
    selected_celllines_genes_df = broad_depmap_df.loc[genes]

    return selected_celllines_genes_df

def get_broad_depmap_selectec_genes(genes):
    broad_depmap_df = broad_depmap_selected_celllines()
    selected_celllines_genes_df = broad_depmap_df.loc[genes]

    return selected_celllines_genes_df

def main():

    # broad_depmap_handle()

    # sanger_depmap_handle()

    # ####################
    # ## selected_ Disease
    # ####################
    #
    disease_type_dict = {
        'disease': "central_nervous_system",
        'disease_sutype': 'glioma'
    }



    novel_genes_df = pd.read_csv(
        '/Users/woochanghwang/PycharmProjects/LifeArc/ULK/result_2/GBM_LGG/Sig_Gene/GBM_LGG_OT/GBM_LGG_ULK_novel_genes.tsv',
        sep='\t')

    novel_genes = novel_genes_df['Gene']

    selected_celllines_genes_df = get_broad_depmap_selected_disease_genes(disease_type_dict,novel_genes)
    #
    # print(selected_celllines_genes_df)

    #############################################

    # selected_celllines = broad_ccle_celllines_selected_disease(disease_type_dict)
    # # selected_celllines = ["U251MG","KNS42","SNU201"]
    # # selected_celllines = ["KNS42"]

    ####
    ## get celline subtype info
    ###

    # cellines_subtypes = get_cellline_subtype(selected_celllines)



    # sanger_depmap_df = sanger_depmap_selected_celllines(selected_celllines)
    # broad_depmap_df = broad_depmap_selected_celllines(selected_celllines)
    #
    #
    # ##############
    # ## boosting , signiricant
    # ##############
    # guide_result = "/home/wch23/Project/LifeArc/ULK/result/GBM_LGG/Sig_Gene/GBM_LGG_OT/GBM_LGG_ULK_conbined_gene_score_by_RW_pvalue_FC.tsv"
    # guide_result_df = pd.read_csv(guide_result, sep='\t')
    # selected_genes = list(guide_result_df['Gene'])
    # boot_num = 1000
    # calc_significance_of_geneSet_A_cellline(selected_genes, selected_celllines[0], broad_depmap_df, boot_num)

    #################
    ## DepMap, disease selection
    ################

    # sanger_cellines = list(sanger_depmap_df)
    # broad_celllines = list(broad_depmap_df)
    #
    # print(len(sanger_cellines), sanger_cellines)
    # print(len(broad_celllines), broad_celllines)
    # print(set(sanger_cellines)&set(broad_celllines))

    #
    # sanger_depmap_df['mean'] = sanger_depmap_df.mean( axis=1)
    # broad_depmap_df['mean'] = broad_depmap_df.mean( axis =1)
    #
    # sanger_depmap_df['binary_sanger'] = sanger_depmap_df['mean'].apply(make_binary_up_down)
    # broad_depmap_df['binary_broad'] = broad_depmap_df['mean'].apply(make_binary_up_down_broad)
    #
    # print(sanger_depmap_df)
    # print(broad_depmap_df)
    #
    # sanger_broad_depmap_binary_df = pd.concat([sanger_depmap_df['binary_sanger'],broad_depmap_df['binary_broad']],axis=1,join='inner')
    # print(sanger_broad_depmap_binary_df)
    #
    # print(len(sanger_depmap_df[sanger_depmap_df['binary_sanger']==0]))
    # print(len(broad_depmap_df[broad_depmap_df['binary_broad']==0]))





if __name__ == '__main__':
    main()