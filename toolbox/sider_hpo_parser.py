# Created by woochanghwang at 12/09/2019
'''
SIDER and HPO data parser
to make Drug - phenotype - gene
'''

import pandas as pd

def make_hpo_groupby_pheno(hpo_phono_to_genes_addr,hpo_pheno_to_genes_groupby_addr ):
    '''
    HPO_ID  HPO_Name    Gene_ID Gene_Name
    :param hpo_phono_to_genes_addr:
    :return:
    '''

    hpo_df = pd.read_csv(hpo_phono_to_genes_addr, sep='\t')

    # hpo_groupby_df = hpo_df.groupby(['HPO_ID','HPO_Name'],as_index=False)['Gene_ID','Gene_Name'].agg(lambda x:list(x))
    hpo_groupby_df = hpo_df.groupby(['HPO_ID', 'HPO_Name'], as_index=False)['Gene_ID', 'Gene_Name'].agg(lambda x:'|'.join(map(str,x)))
    print(hpo_groupby_df.head())

    hpo_groupby_df.to_csv(hpo_pheno_to_genes_groupby_addr ,sep='\t', index=False)


def make_join_sider_hpo_by_phenotype(sider_addr, hpo_addr, sider_hpo_addr, sider_comp_name_addr):
    # ['STITCH_Comp_ID', 'UMLS_ID', 'Method_of_detection', 'Concept_name', 'Concept_type', 'UMLS_ID_MedDRA','MedDRA_concept_name']
    sider_df = pd.read_csv(sider_addr,sep='\t')

    sider_pt_df = sider_df[sider_df['Concept_type']=='PT']

    # print(sider_pt_df.head())
    hpo_df = pd.read_csv(hpo_addr,sep='\t')

    # sider_pt_df = sider_pt_df.set_index('MedDRA_concept_name')

    # hpo_df = hpo_df.rename(columns={'HPO_Name':'MedDRA_concept_name'})
    # hpo_df = hpo_df.set_index('MedDRA_concept_name')

    sider_hpo_df = sider_pt_df.merge(hpo_df, left_on='MedDRA_concept_name', right_on='HPO_Name')
    # sider_hpo_df = pd.concat([sider_pt_df,hpo_df],axis=1, join='inner', ignore_index=False)
    # #
    print(sider_hpo_df.head())

    print(sider_pt_df.head())
    print(hpo_df.head())

    # sider_hpo_df = sider_hpo_df[['STITCH_Comp_ID','Method_of_detection','UMLS_ID_MedDRA','HPO_ID', 'HPO_Name','Gene_ID', 'Gene_Name']]
    sider_hpo_df = sider_hpo_df[
        ['STITCH_Comp_ID', 'UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'Gene_ID', 'Gene_Name']]

    sider_comp_name_df = pd.read_csv(sider_comp_name_addr, sep='\t')
    sider_hpo_drugname_df = pd.merge(sider_hpo_df,sider_comp_name_df,on='STITCH_Comp_ID')

    # print(sider_hpo_drugname_df.head())
    sider_hpo_drugname_df = sider_hpo_drugname_df[
        ['STITCH_Comp_ID', 'Comp_Name','UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'Gene_ID', 'Gene_Name']]
    sider_hpo_drugname_df.to_csv(sider_hpo_addr, index=False, sep='\t')


def main():
    hpo_pheno_to_genes_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/HPO/HPO_phenotype_to_genes.txt"
    sider_compound_to_indication_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/meddra_all_indications.tsv"
    sider_compound_to_se_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/meddra_all_se.tsv"
    sider_compound_name_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/drug_names.tsv"

    ##HPO data groupby HPO_ID
    hpo_pheno_to_genes_groupby_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/HPO/HPO_phenotype_to_genes_groupby.txt"
    # make_hpo_groupby_pheno(hpo_pheno_to_genes_addr, hpo_pheno_to_genes_groupby_addr  )

    ##Join Sider to HPO by phenotype
    sider_indication_hpo_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_indications_hpo_mapped.tsv"
    make_join_sider_hpo_by_phenotype(sider_compound_to_indication_addr, hpo_pheno_to_genes_groupby_addr, sider_indication_hpo_addr, sider_compound_name_addr)

    ##Join Sider to HPO by phenotype
    sider_se_hpo_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_se_hpo_mapped.tsv"
    # make_join_sider_hpo_by_phenotype(sider_compound_to_se_addr, hpo_pheno_to_genes_groupby_addr,
    #                                  sider_se_hpo_addr, sider_compound_name_addr)


def sider_hpo_grouby_genes():
    hpo_pheno_to_genes_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/HPO/HPO_phenotype_to_genes.txt"
    sider_compound_to_indication_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/meddra_all_indications.tsv"
    sider_compound_to_se_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/meddra_all_se.tsv"
    sider_compound_name_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/drug_names.tsv"


    ##Join Sider to HPO by phenotype
    # sider_indication_hpo_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_indications_mapped_hpo_all.tsv"
    sider_se_hpo_addr =  "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_se_hpo_mapped.tsv"
    sider_df = pd.read_csv(sider_compound_to_se_addr, sep='\t')

    sider_pt_df = sider_df[sider_df['Concept_type'] == 'PT']
    hpo_df = pd.read_csv(hpo_pheno_to_genes_addr, sep='\t')

    sider_hpo_df = sider_pt_df.merge(hpo_df, left_on='MedDRA_concept_name', right_on='HPO_Name')

    print(sider_hpo_df.head())
    print(list(sider_hpo_df))

    sider_hpo_selected_df = sider_hpo_df[['Gene_ID','Gene_Name','MedDRA_concept_name']]

    sider_hpo_selected_df = sider_hpo_selected_df.drop_duplicates()

    sider_hpo_selected_df.to_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_hpo_all.tsv",sep='\t',index=False)

    # sider_groupby_df = sider_hpo_selected_df.groupby(['Gene_ID', 'Gene_Name'], as_index=False)['MedDRA_concept_name'].agg(
    #     lambda x: '|'.join(map(str, x)))
    #
    # sider_groupby_df = sider_hpo_selected_df.groupby(['Gene_ID', 'Gene_Name'], as_index=False)[
    #     'MedDRA_concept_name'].agg(
    #     lambda x:list(x))
    #
    # print(sider_groupby_df.head())
    #
    # sider_groupby_df['concept_length'] = sider_groupby_df['MedDRA_concept_name'].apply(lambda x:len(x) )
    #
    #
    # sider_hpo_gene_phenotype_len_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_se_hpo_mapped_concept_length.tsv"
    # sider_groupby_df.to_csv(sider_hpo_gene_phenotype_len_addr, sep='\t', index=False)

def find_gene_count_on_se_indi():
    hpo_pheno_to_genes_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/HPO/HPO_phenotype_to_genes.txt"
    sider_compound_to_indication_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/meddra_all_indications.tsv"
    sider_compound_to_se_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/meddra_all_se.tsv"
    sider_compound_name_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/drug_names.tsv"

    ##Join Sider to HPO by phenotype
    # sider_indication_hpo_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_indications_mapped_hpo_all.tsv"
    sider_se_hpo_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_se_hpo_mapped.tsv"
    sider_df = pd.read_csv(sider_compound_to_se_addr, sep='\t')

    sider_pt_df = sider_df[sider_df['Concept_type'] == 'PT']
    hpo_df = pd.read_csv(hpo_pheno_to_genes_addr, sep='\t')

    sider_hpo_df = sider_pt_df.merge(hpo_df, left_on='MedDRA_concept_name', right_on='HPO_Name')

    print(sider_hpo_df.head())
    print(list(sider_hpo_df))

    sider_hpo_selected_df = sider_hpo_df[['STITCH_Comp_ID','Gene_ID', 'Gene_Name', 'MedDRA_concept_name']]

    sider_hpo_selected_df = sider_hpo_selected_df.drop_duplicates()

    print(sider_hpo_selected_df.head())


def add_atc_code_on_sider_hpo_mapped():
    '''
    Compound- one ATC Code - one phenotype
    :return:
    '''
    atc_code_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drugbank/drugbank_atc_code_hierachi.tsv"
    sider_drug_atc_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/drug_atc.tsv"
    sider_indication_hpo_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_indications_hpo_mapped.tsv"
    sider_se_hpo_addr="/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_se_hpo_mapped.tsv"

    sider_drug_atc_df = pd.read_csv(sider_drug_atc_addr,sep='\t',names=['STITCH_Comp_ID','ATC_CODE'])
    # sider_indi_hpo_df = pd.read_csv(sider_indication_hpo_addr,sep='\t')
    sider_indi_hpo_df = pd.read_csv(sider_se_hpo_addr, sep='\t')


    print(list(sider_drug_atc_df))
    print(sider_drug_atc_df.head())

    print(list(sider_indi_hpo_df))
    print(sider_indi_hpo_df.size)
    # ['STITCH_Comp_ID', 'Comp_Name', 'UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'Gene_ID', 'Gene_Name']

    sider_indi_hpo_atc_df = sider_indi_hpo_df.merge(sider_drug_atc_df, on='STITCH_Comp_ID')
    sider_indi_hpo_atc_df = sider_indi_hpo_atc_df[['STITCH_Comp_ID', 'Comp_Name','ATC_CODE', 'UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'Gene_ID', 'Gene_Name']]
    sider_indi_hpo_atc_df = sider_indi_hpo_atc_df.drop_duplicates()
    print(sider_indi_hpo_atc_df.size)

    sider_indi_hpo_atc_df['ATC_Code_parent'] = sider_indi_hpo_atc_df.apply(lambda x:x['ATC_CODE'][:-2],axis=1)
    sider_indi_hpo_atc_df = sider_indi_hpo_atc_df.drop_duplicates()
    print(list(sider_indi_hpo_atc_df))
    print(sider_indi_hpo_atc_df.head())
    print(sider_indi_hpo_atc_df.size)

    drug_bank_atc_code_df = pd.read_csv(atc_code_addr, sep='\t', names=["ATC_Code","Description"])

    sider_indi_hpo_atc_drugbank_df = sider_indi_hpo_atc_df.merge(drug_bank_atc_code_df,how='left',left_on="ATC_Code_parent",right_on="ATC_Code")

    sider_indi_hpo_atc_drugbank_df = sider_indi_hpo_atc_drugbank_df[['STITCH_Comp_ID', 'Comp_Name', 'ATC_CODE','Description', 'UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'Gene_ID', 'Gene_Name']]
    sider_indi_hpo_atc_drugbank_df = sider_indi_hpo_atc_drugbank_df.drop_duplicates()
    print(list(sider_indi_hpo_atc_drugbank_df))
    print(sider_indi_hpo_atc_drugbank_df.head())
    print(sider_indi_hpo_atc_drugbank_df.size)

    # sider_indi_hpo_atc_drugbank_df.to_csv("/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_indications_hpo_mapped_atc_code.tsv",sep='\t',index=False)
    sider_indi_hpo_atc_drugbank_df.to_csv(
        "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_se_hpo_mapped_atc_code.tsv",
        sep='\t', index=False)


def make_unique_genes(genes):
    genes = genes.split('|')
    genes = list(set(genes))
    genes = '|'.join(genes)
    return genes

def make_unique_compound_ATC_codes_phenotypes():
    '''
    Compound- multi ATC Code - multi phenotype
    :return:
    '''
    atc_code_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drugbank/drugbank_atc_code_hierachi.tsv"
    sider_drug_atc_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/drug_atc.tsv"
    sider_indication_hpo_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_indications_hpo_mapped.tsv"
    sider_se_hpo_addr="/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_se_hpo_mapped.tsv"

    sider_indication_hpo_final_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_indications_multi_hpo_multi_atc_mapped.tsv"
    sider_se_hpo_final_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_ses_multi_hpo_multi_atc_mapped.tsv"

    sider_drug_atc_df = pd.read_csv(sider_drug_atc_addr,sep='\t',names=['STITCH_Comp_ID','ATC_CODE'])

    #############
    # Indication
    ##############
    # sider_hpo_addr = sider_indication_hpo_addr  # indication or side effect
    # sider_hpo_final_addr = sider_indication_hpo_final_addr
    #############
    # SE
    #############
    sider_hpo_addr = sider_se_hpo_addr  # indication or side effect
    sider_hpo_final_addr = sider_se_hpo_final_addr
    #######
    sider_hpo_df = pd.read_csv(sider_hpo_addr,sep='\t')

    # hpo_df.groupby(['HPO_ID', 'HPO_Name'], as_index=False)['Gene_ID', 'Gene_Name'].agg(lambda x: '|'.join(map(str, x)))

    ######
    # Sider , atc code, groupy compound
    ######
    sider_drug_atcs_df = sider_drug_atc_df.groupby('STITCH_Comp_ID',as_index=False).agg(lambda x:'|'.join(map(str,x)))

    #####
    # sider-hpo groupby compound
    #######
    print(list(sider_hpo_df))

    # ['STITCH_Comp_ID', 'Comp_Name', 'UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'Gene_ID', 'Gene_Name']

    sider_hpo_groupby_comp_df = sider_hpo_df.groupby(['STITCH_Comp_ID','Comp_Name'], as_index=False)['UMLS_ID_MedDRA','HPO_ID','HPO_Name','Gene_ID','Gene_Name'].agg(lambda x:'|'.join(map(str,x)))

    print(sider_hpo_groupby_comp_df.head())

    sider_hpo_groupby_comp_df['unique_Gene_ID'] = sider_hpo_groupby_comp_df['Gene_ID'].apply(make_unique_genes)
    sider_hpo_groupby_comp_df['unique_Gene_Name'] = sider_hpo_groupby_comp_df['Gene_Name'].apply(make_unique_genes)

    sider_hpo_groupby_comp_df = sider_hpo_groupby_comp_df[['STITCH_Comp_ID', 'Comp_Name', 'UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'unique_Gene_ID', 'unique_Gene_Name']]
    sider_hpo_groupby_comp_df = sider_hpo_groupby_comp_df.rename(columns={'unique_Gene_ID':'Gene_ID','unique_Gene_Name':'Gene_Name'})
    print(sider_hpo_groupby_comp_df.head())

    sider_hpo_atc_df = sider_hpo_groupby_comp_df.merge(sider_drug_atcs_df,on='STITCH_Comp_ID')

    print(list(sider_hpo_atc_df))
    print(sider_hpo_atc_df.head())

    sider_hpo_atc_df = sider_hpo_atc_df[['STITCH_Comp_ID', 'Comp_Name', 'ATC_CODE','UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'Gene_ID', 'Gene_Name']]

    sider_hpo_atc_df.to_csv(sider_hpo_final_addr,sep='\t',index=False)

def make_unique_compound_phenotypes():
    '''
    Compound- multi ATC Code - multi phenotype
    :return:
    '''
    atc_code_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drugbank/drugbank_atc_code_hierachi.tsv"
    sider_drug_atc_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/drug_atc.tsv"
    sider_indication_hpo_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_indications_hpo_mapped.tsv"
    sider_se_hpo_addr="/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_se_hpo_mapped.tsv"

    sider_indication_hpo_final_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_indications_multi_hpo_multi.tsv"
    sider_se_hpo_final_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_ses_multi_hpo_multi.tsv"

    sider_drug_atc_df = pd.read_csv(sider_drug_atc_addr,sep='\t',names=['STITCH_Comp_ID','ATC_CODE'])

    #############
    # Indication
    ##############
    sider_hpo_addr = sider_indication_hpo_addr  # indication or side effect
    sider_hpo_final_addr = sider_indication_hpo_final_addr
    #############
    # SE
    #############
    # sider_hpo_addr = sider_se_hpo_addr  # indication or side effect
    # sider_hpo_final_addr = sider_se_hpo_final_addr
    #######
    sider_hpo_df = pd.read_csv(sider_hpo_addr,sep='\t')


    #####
    # sider-hpo groupby compound
    #######
    print(list(sider_hpo_df))

    # ['STITCH_Comp_ID', 'Comp_Name', 'UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'Gene_ID', 'Gene_Name']

    sider_hpo_groupby_comp_df = sider_hpo_df.groupby(['STITCH_Comp_ID','Comp_Name'], as_index=False)['UMLS_ID_MedDRA','HPO_ID','HPO_Name','Gene_ID','Gene_Name'].agg(lambda x:'|'.join(map(str,x)))

    print(sider_hpo_groupby_comp_df.head())

    sider_hpo_groupby_comp_df['unique_Gene_ID'] = sider_hpo_groupby_comp_df['Gene_ID'].apply(make_unique_genes)
    sider_hpo_groupby_comp_df['unique_Gene_Name'] = sider_hpo_groupby_comp_df['Gene_Name'].apply(make_unique_genes)

    sider_hpo_groupby_comp_df = sider_hpo_groupby_comp_df[['STITCH_Comp_ID', 'Comp_Name', 'UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'unique_Gene_ID', 'unique_Gene_Name']]
    sider_hpo_groupby_comp_df = sider_hpo_groupby_comp_df.rename(columns={'unique_Gene_ID':'Gene_ID','unique_Gene_Name':'Gene_Name'})
    print(sider_hpo_groupby_comp_df.head())

    # sider_hpo_atc_df = sider_hpo_groupby_comp_df.merge(sider_drug_atcs_df,on='STITCH_Comp_ID')

    print(list(sider_hpo_groupby_comp_df))
    print(sider_hpo_groupby_comp_df.head())

    sider_hpo_groupby_comp_df = sider_hpo_groupby_comp_df[['STITCH_Comp_ID', 'Comp_Name','UMLS_ID_MedDRA', 'HPO_ID', 'HPO_Name', 'Gene_ID', 'Gene_Name']]

    sider_hpo_groupby_comp_df.to_csv(sider_hpo_final_addr,sep='\t',index=False)

def select_se_for_specific_drugs(selected_drugs):

    # sider_indication_hpo_final_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_indications_multi_hpo_multi.tsv"
    sider_se_hpo_final_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/data/SIDER/sider_ses_multi_hpo_multi.tsv"

    sider_se_df = pd.read_csv(sider_se_hpo_final_addr,sep='\t')

    sider_se_selected_drug_df = sider_se_df[sider_se_df['Comp_Name'].isin(selected_drugs)]

    print(sider_se_selected_drug_df)

    sider_se_hpo = sider_se_selected_drug_df['HPO_Name']
    hpos = []
    for i in range(len(sider_se_hpo)):
        hpos.append(sider_se_hpo.iloc[i].split('|'))
    print(hpos)

    for hpo in hpos:
        print(len(hpo))
    print(len(set(hpo[0])&set(hpo[1])&set(hpo[2])&set(hpo[3])))
    # result_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/General/src_side_effect/Result/sider_se_selected_drug.tsv"
    # sider_se_selected_drug_df.to_csv(result_addr, sep='\t', index=False)

if __name__ == '__main__':
    # main()
    # sider_hpo_grouby_genes()
    # find_gene_count_on_se_indi()
    # add_atc_code_on_sider_hpo_mapped()

    # make_unique_compound_ATC_codes_phenotypes()

    # make_unique_compound_phenotypes()
    selected_drugs = ['gefitinib','erlotinib','cediranib','vandetanib','canertinib','bosutinib']
    select_se_for_specific_drugs(selected_drugs)



    # a