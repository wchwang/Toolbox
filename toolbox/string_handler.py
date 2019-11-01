# Created by woochanghwang at 21/10/2019
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import toolbox.visual_utilities as vu
import toolbox.gprofiler_utilities as gp_util
import toolbox.data_handler as dh
import mygene

def make_sig_ppi_for_mouse(string_addr, string_sig_addr):
    sig_val = 400

    string_all = pd.read_csv(string_addr, sep=' ')
    string_sig = string_all[string_all['combined_score']>=sig_val]

    # protein_a_ensembl = string_sig['protein1'].str.split('.', n =1, expand=True)
    # protein_b_ensembl = string_sig['protein2'].str.split('.', n =1, expand=True)
    #
    # string_sig['protein1_ensembl'] = protein_a_ensembl[1]
    # string_sig['protein2_ensembl'] = protein_b_ensembl[1]
    #
    # string_sig.drop(columns=['protein1','protein2'], inplace=True)
    #
    # string_sig = string_sig[['protein1_ensembl','protein2_ensembl','combined_score']]
    print(string_sig.head())
    string_sig.to_csv(string_sig_addr,index=False,sep='\t',quoting=None)

def make_sig_ppi_with_symbol_for_mouse(string_sig_addr,string_info_addr, string_sig_symbol_addr):

    string_info_df = pd.read_csv(string_info_addr, sep='\t')
    print(string_info_df.columns)

    string_info_df = string_info_df.set_index('protein_external_id')

    print(string_info_df.head())
    string_id_symbol_dict = string_info_df.to_dict()['preferred_name']

    string_sig_df = pd.read_csv(string_sig_addr, sep='\t')
    string_sig_df['protein1_symbol'] = string_sig_df['protein1'].map(string_id_symbol_dict)
    string_sig_df['protein2_symbol'] = string_sig_df['protein2'].map(string_id_symbol_dict)

    # print(string_sig_df.head())
    string_sig_symbol_df = string_sig_df[["protein1_symbol","protein2_symbol"]]
    print(string_sig_symbol_df.head())

    string_sig_symbol_df.to_csv(string_sig_symbol_addr,index=False,sep='\t',quoting=None)

def make_sig_ppi_with_ensembl_for_mouse(string_sig_addr,string_info_addr, string_sig_ensembl_addr):
    string_info_df = pd.read_csv(string_info_addr, sep='\t')
    print(string_info_df.columns)

    string_info_df = string_info_df.set_index('protein_external_id')

    print(string_info_df.head())
    string_id_ensembl_dict = string_info_df.to_dict()['converted']

    string_sig_df = pd.read_csv(string_sig_addr, sep='\t')
    string_sig_df['protein1_ensembl'] = string_sig_df['protein1'].map(string_id_ensembl_dict)
    string_sig_df['protein2_ensembl'] = string_sig_df['protein2'].map(string_id_ensembl_dict)

    # print(string_sig_df.head())
    string_sig_ensembl_df = string_sig_df[["protein1_ensembl", "protein2_ensembl"]]
    string_sig_ensembl_df = string_sig_ensembl_df.dropna()
    print(string_sig_ensembl_df.head())

    string_sig_ensembl_df.to_csv(string_sig_ensembl_addr, index=False, sep='\t', quoting=None)

def add_ensembl_gene_into_string_info(string_info_addr):
    string_info_df = pd.read_csv(string_info_addr, sep='\t')
    protein_tax_ensembl = string_info_df['protein_external_id'].str.split(".",n=1,expand=True)
    string_info_df['protein_ensembl.protein'] = protein_tax_ensembl[1]

    protein_ensembl = protein_tax_ensembl[1].tolist()
    print(protein_ensembl)
    # string_info_df = string_info_df.iloc[:-1]

    # mg = mygene.MyGeneInfo()
    # ensembl_protein_to_gene_df = mg.querymany(protein_ensembl, scopes='ensembl.protein', fields='ensembl.gene',
    #                                           species=10090,returnall=False, as_dataframe=True)
    #
    # # ensembl_protein_to_gene_df = mg.getgenes(protein_ensembl, fields='ensembl.gene',
    # #                                           species=10090,as_dataframe=True)
    #
    # print(ensembl_protein_to_gene_df)
    #
    #
    # ensembl_protein_to_gene_df.to_csv("/Users/woochanghwang/PycharmProjects/CIMR/Data/STRING/ensembl.id.tsv",sep='\t')

    from gprofiler import GProfiler
    gp = GProfiler(return_dataframe=True)

    ensembl_protein_to_gene_df = gp.orth(organism='mmusculus',
                       target='ENSG',
            query=protein_ensembl)


    ensembl_protein_to_gene_df = ensembl_protein_to_gene_df.reset_index()
    ensembl_protein_to_gene_df = ensembl_protein_to_gene_df.set_index('incoming')
    # ensembl_protein_to_gene_df = ensembl_protein_to_gene_df.drop('index')
    # ensembl_protein_to_gene_df.to_csv("/Users/woochanghwang/PycharmProjects/CIMR/Data/STRING/string_ensembl_protein_to_gene_gp.tsv",sep='\t',index=False)
    # string_info_df = string_info_df.drop_duplicates()
    # ensembl_protein_to_gene_df = ensembl_protein_to_gene_df.drop_duplicates()
    string_info_df = string_info_df.set_index('protein_ensembl.protein')

    print(string_info_df.head())
    print(ensembl_protein_to_gene_df.head())
    string_info_ensembl_df = pd.concat([string_info_df,ensembl_protein_to_gene_df], axis=1, sort=False)

    # string_info_ensembl_df = pd.merge(string_info_df, ensembl_protein_to_gene_df)
    string_info_ensembl_df.to_csv("/Users/woochanghwang/PycharmProjects/CIMR/Data/STRING/10090.protein.info.v11.0.ensembl.txt",sep='\t')


def main():
    string_addr = "/Users/woochanghwang/PycharmProjects/CIMR/Data/STRING/10090.protein.links.v11.0.txt"
    string_sig_addr = "/Users/woochanghwang/PycharmProjects/CIMR/Data/STRING/10090.protein.links.v11.0.400.txt"
    string_info_addr = "/Users/woochanghwang/PycharmProjects/CIMR/Data/STRING/10090.protein.info.v11.0.ensembl.txt"
    make_sig_ppi_for_mouse(string_addr, string_sig_addr)

    string_sig_symbol_addr = "/Users/woochanghwang/PycharmProjects/CIMR/Data/STRING/10090.protein.links.v11.0.400.symbol.txt"
    string_sig_ensembl_addr = "/Users/woochanghwang/PycharmProjects/CIMR/Data/STRING/10090.protein.links.v11.0.400.ensembl.txt"

    # make_sig_ppi_with_symbol_for_mouse(string_sig_addr,string_info_addr, string_sig_symbol_addr)
    make_sig_ppi_with_ensembl_for_mouse(string_sig_addr, string_info_addr, string_sig_ensembl_addr)
    # add_ensembl_gene_into_string_info(string_info_addr)


if __name__ == '__main__':
    main()