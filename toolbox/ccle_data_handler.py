# Created by woochanghwang at 2019-05-29

from cmapPy.pandasGEXpress.parse import parse
import toolbox.tcga_data_handler as tcga_dh
import pandas as pd
import mygene
import toolbox.data_handler as data_h
import csv

CCLE_RNA_TPM_SYMBOL_ADDR = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/CCLE/CCLE_RNAseq_rsem_genes_tpm_20180929_symbol"

def find_gene_rid(ccle_GCToo, target_gene):
    try:
        return ccle_GCToo.row_metadata_df[ccle_GCToo.row_metadata_df['Description']==target_gene].index.item()
    except:
        pass


def get_CCLE_for_selecte_geneSet(selected_gene_list,result_addr=None):
    '''

    :param selected_gene_list: Gene Symbol list
    :param result_addr: For writing result, result address
    :return: DataFrame(Row : Genes, Col : Cell lines
    '''
    ccle_GCToo = parse(
        "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/CCLE/CCLE_RNAseq_genes_rpkm_20180929.gct")
    print(ccle_GCToo.row_metadata_df[:10])
    print(ccle_GCToo.col_metadata_df[:10])
    print(ccle_GCToo.data_df.head())

    selected_gene_rid_list = []
    for target_gene in selected_gene_list:
        target_gene_rid = find_gene_rid(ccle_GCToo, target_gene)
        selected_gene_rid_list.append(target_gene_rid)

    selected_gene_expression_df = ccle_GCToo.data_df.loc[selected_gene_rid_list, :]
    selected_gene_expression_df_T = selected_gene_expression_df.T
    print(selected_gene_expression_df_T)

    gene_map_dict = dict()
    for gene, rid in zip(selected_gene_list, selected_gene_rid_list):
        gene_map_dict[rid] = gene

    # rename
    selected_gene_expression_df_T = selected_gene_expression_df_T.rename(index=str, columns=gene_map_dict)
    selected_gene_expression_df = selected_gene_expression_df_T.T

    if result_addr != None:
        selected_gene_expression_df.to_csv(result_addr,sep='\t')

    return selected_gene_expression_df


def get_CCLE_Exp_gct_from_selected_genes(selected_gene_list,result_addr=None):
    ccle_GCToo = parse(
        "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/CCLE/CCLE_RNAseq_genes_rpkm_20180929.gct")
    # print(ccle_GCToo.row_metadata_df[:10])
    # print(ccle_GCToo.col_metadata_df[:10])
    # print(ccle_GCToo.data_df.head())

    selected_gene_rid_list = []
    for target_gene in selected_gene_list:
        target_gene_rid = find_gene_rid(ccle_GCToo, target_gene)
        selected_gene_rid_list.append(target_gene_rid)

    selected_gene_expression_df = ccle_GCToo.data_df.loc[selected_gene_rid_list, :]
    selected_gene_expression_df_T = selected_gene_expression_df.T
    print(selected_gene_expression_df_T)

    gene_map_dict = dict()
    for gene, rid in zip(selected_gene_list, selected_gene_rid_list):
        gene_map_dict[rid] = gene

    # rename
    selected_gene_expression_df_T = selected_gene_expression_df_T.rename(index=str, columns=gene_map_dict)
    selected_gene_expression_df = selected_gene_expression_df_T.T

    if result_addr != None:
        selected_gene_expression_df.to_csv(result_addr, sep='\t', quoting=csv.QUOTE_NONE, index=False)
    return selected_gene_expression_df

def get_CCLE_RNA_tpm_for_selected_genes(selected_genelist,result_addr=None):


    ccle_rna_tpm_symbol_df = data_h.load_obj(CCLE_RNA_TPM_SYMBOL_ADDR)

    ccle_rna_tpm_symbol_df = ccle_rna_tpm_symbol_df.set_index('symbol')

    ccle_rna_tmp_for_selected_genes_df = ccle_rna_tpm_symbol_df.loc[selected_genelist]

    print(ccle_rna_tmp_for_selected_genes_df)

    if result_addr != None:
        ccle_rna_tmp_for_selected_genes_df =  ccle_rna_tmp_for_selected_genes_df.reset_index()
        ccle_rna_tmp_for_selected_genes_df.to_csv(result_addr,sep='\t', index=False)

    return ccle_rna_tmp_for_selected_genes_df

def get_CCLE_RNA_tpm_for_specific_histology(histology,result_addr=None):
    '''

    :param histology: ex)glioma
    :return: rows: symbol, cols = ccle_ids
    '''
    ccle_info_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/CCLE/Cell_lines_annotations_20181226.txt"

    ccle_info_df = pd.read_csv(ccle_info_addr, sep='\t')

    # print(list(ccle_info_df))
    ccle_specific_histology = ccle_info_df[ccle_info_df['Histology'] == histology]

    ccle_ids_for_histology = list(ccle_specific_histology['CCLE_ID'])

    ccle_rna_tpm_symbol_df = data_h.load_obj(CCLE_RNA_TPM_SYMBOL_ADDR)

    ccle_rna_tpm_symbol_df = ccle_rna_tpm_symbol_df.set_index('symbol')

    ccle_samples = list(ccle_rna_tpm_symbol_df)
    # print(ccle_ids_for_histology)

    ccle_ids_for_histology_in = set(ccle_ids_for_histology) & set(ccle_samples)

    ccle_rna_tmp_for_histology_df = ccle_rna_tpm_symbol_df[ccle_ids_for_histology_in]

    print(ccle_rna_tmp_for_histology_df)

    if result_addr != None:
        ccle_rna_tmp_for_histology_df = ccle_rna_tmp_for_histology_df.reset_index()
        ccle_rna_tmp_for_histology_df.to_csv(result_addr, index=False)

    return ccle_rna_tmp_for_histology_df

def get_CCLE_RNA_tmp_for_specific_histology_and_genes(histology, genes, result_addr=None):
    ccle_for_specific_histology_df = get_CCLE_RNA_tpm_for_specific_histology(histology)

    ccle_for_specific_histology_genes_df = ccle_for_specific_histology_df.loc[genes]

    print(ccle_for_specific_histology_genes_df)

    if result_addr !=None:
        ccle_for_specific_histology_genes_df = ccle_for_specific_histology_genes_df.reset_index()
        ccle_for_specific_histology_genes_df.to_csv(result_addr, index=False)

    return ccle_for_specific_histology_genes_df

def make_CCLE_RNA_TPM_with_GeneSymbol():
    ccle_tpm_esemble_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/CCLE/CCLE_RNAseq_rsem_genes_tpm_20180929.txt"

    ccle_tpm_esemble_df = pd.read_csv(ccle_tpm_esemble_addr,sep='\t')
    # print(ccle_tpm_esemble_df)

    ccle_genes = ccle_tpm_esemble_df['gene_id'].to_list()
    ccle_genes = [x.split('.')[0] for x in ccle_genes]
    new_gene_id= pd.Series(ccle_genes)
    ccle_tpm_esemble_df.insert(loc=0, column='gene_new_id',value=new_gene_id)
    ccle_tpm_esemble_df = ccle_tpm_esemble_df.drop(columns=['gene_id','transcript_ids'])

    # print(ccle_tpm_esemble_df)
    ccle_tpm_esemble_df = ccle_tpm_esemble_df.set_index('gene_new_id')
    ### I have alread made gene mapping file
    # mg = mygene.MyGeneInfo()
    #
    # gene_ids_df = mg.querymany(ccle_genes,scopes="ensembl.gene,symbol",species='human',as_dataframe=True)
    #
    # gene_ids_df.to_csv('/Users/woochanghwang/PycharmProjects/LifeArc/General/data/CCLE/ccle_ensembl_to_symbol.csv')

    ccle_ensembl_to_symbol_df = pd.read_csv('/Users/woochanghwang/PycharmProjects/LifeArc/General/data/CCLE/ccle_ensembl_to_symbol.csv', index_col=0)
    ccle_ensembl_to_symbol_df = ccle_ensembl_to_symbol_df.reset_index().drop_duplicates(subset='query',keep='first').set_index('query')
    ccle_symbol_series = ccle_ensembl_to_symbol_df['symbol']
    # print(ccle_ensembl_to_symbol_df)
    # print(ccle_tpm_esemble_df)


    # ccle_tpm_esemble_df = pd.concat([ccle_tpm_esemble_df,ccle_symbol_series],axis=1)
    ccle_tpm_esemble_df.insert(loc=0, column='symbol',value=ccle_symbol_series)
    # print(ccle_tpm_esemble_df)

    ccle_tpm_symbol_df = ccle_tpm_esemble_df.reset_index().drop(columns=['gene_new_id'])
    ccle_tpm_symbol_df.to_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/CCLE/CCLE_RNAseq_rsem_genes_tpm_20180929_symbol.csv",index=False)

    data_h.save_obj(ccle_tpm_symbol_df, CCLE_RNA_TPM_SYMBOL_ADDR)



def main():
    # make_CCLE_RNA_TPM_with_GeneSymbol()
    # get_CCLE_Exp_tpm_for_specific_Histology('glioma')

    novel_genes_df = pd.read_csv('/Users/woochanghwang/PycharmProjects/LifeArc/ULK/result_2/GBM_LGG/Sig_Gene/GBM_LGG_OT/GBM_LGG_ULK_novel_genes.tsv',sep='\t')

    novel_genes = novel_genes_df['Gene']
    novel_genes_ccle = get_CCLE_RNA_tpm_for_selected_genes(novel_genes)

if __name__ == '__main__':
    main()