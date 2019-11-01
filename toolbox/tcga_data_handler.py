# Created by woochanghwang at 2019-05-29
import pandas as pd
import csv

def get_TCGA_Tumour_mean_value_for_selected_geneSet( tumour_types, selected_gene_list, cancer_type, result_addr=None):
    '''

    :param tumour_types: TCGA Tumour code ex)CHOL
    :param selected_gene_list:
    :param cancer_type: basically same as tumour type, but if you want multi tumour types, you need to create one name( CHC : CHOL+LICH)
    :param result_addr: for writing result, result addr
    :return: Dataframe(Row: seletect genes, Col: Mean
    '''
    tcga_tumour_addr = "/Users/woochanghwang/PycharmProjects/TCGA/data/RNAseqV2/processed.gene.norm/TCGA_{}_tumour_new_genesymbol.tsv".format(cancer_type)

    all_tumour_table = pd.read_table(tcga_tumour_addr, sep='\t')

    all_tumour_table = all_tumour_table.replace('\n', '', regex=True)

    tcga_symbol_entrez = all_tumour_table[['Symbol', 'Entrez', 'new_Symbol']]
    # tcga_symbol_entrez['Entrez'] = tcga_symbol_entrez['Entrez'].astype(str)
    #
    # tcga_symbol_entrez['Symbol|Entrez'] = tcga_symbol_entrez[['new_Symbol', 'Entrez']].apply(lambda x: '|'.join(x),
    #                                                                                          axis=1)

    tcga_symbol_selected_tumours_df = tcga_symbol_entrez['new_Symbol']

    print(tcga_symbol_selected_tumours_df.head())

    filter_tumour_col = [col for col in all_tumour_table if col.startswith(tumour_types)]

    tcga_tumour_a_type_df = all_tumour_table[filter_tumour_col]

    # tcga_normal_a_type_df['normal_mean'] = tcga_normal_a_type_df.mean(axis=1)
    # print("tumor_type:",tcga_normal_a_type_df.mean(axis=1))

    tcga_tumour_a_type_df = tcga_tumour_a_type_df.mean(axis=1)


    tcga_symbol_selected_tumours_df = pd.concat(
        (tcga_symbol_selected_tumours_df,tcga_tumour_a_type_df), axis=1)


    print(tcga_symbol_selected_tumours_df.head())

    tcga_symbol_selected_tumours_df = tcga_symbol_selected_tumours_df[tcga_symbol_selected_tumours_df.new_Symbol !='?']
    # print(tcga_symbol_selected_tumours_df['new_Symbol'].isin(selected_gene_list))

    tcga_symbol_selected_tumours_df = tcga_symbol_selected_tumours_df[tcga_symbol_selected_tumours_df['new_Symbol'].isin(selected_gene_list)]


    tcga_symbol_selected_tumours_df = tcga_symbol_selected_tumours_df.rename(columns={'new_Symbol':"Symbol",0:cancer_type})
    print(tcga_symbol_selected_tumours_df)

    if result_addr != None:
        tcga_symbol_selected_tumours_df.to_csv(result_addr, sep='\t', quoting=csv.QUOTE_NONE, index=False)
    return tcga_symbol_selected_tumours_df

def main():
    pass


if __name__ == '__main__':
    main()