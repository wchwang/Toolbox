# Created by woochanghwang at 15/10/2019

'''
Make TCGA + GTEx (Gepia)
'''
import pandas as pd
from gprofiler  import GProfiler

def make_TCGA_xena_samples_verbose():
    tcga_verbose_addr = "/home/wch23/Project/LifeArc/TCGA/data/TGCA_xena_samples_verbose.csv"
    tcga_abbreviation_addr = "/home/wch23/Project/LifeArc/TCGA/data/TCGA_Study_Abbreviation.tsv"

    tcga_verbose_df = pd.read_csv(tcga_verbose_addr)
    print(tcga_verbose_df.head())

    tcga_abbreviation_df = pd.read_csv(tcga_abbreviation_addr, sep='\t', names=['Abbreviation','Study Name'])
    print(tcga_abbreviation_df.head())

    tcga_verbose_df['Abbreviation'] = tcga_verbose_df[['Study Name']].merge(tcga_abbreviation_df,how='left').Abbreviation

    print(tcga_verbose_df.head())

    tcga_verbose_df.to_csv("/home/wch23/Project/LifeArc/TCGA/data/TCGA_xena_samples_verbose_with_abbreviation.csv", index=False)

def make_GTEx_xena_samples_verbose():
    # gtex_verbose_addr = "/home/wch23/Project/LifeArc/TCGA/data/GTEX_xena_samples_verbose.csv"
    gtex_verbose_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/GTEx_v7_Annotations_SampleAttributesDS.csv"
    gtex_tcga_mapping_addr = "/home/wch23/Project/LifeArc/TCGA/data/Pair_match.txt"

    gtex_verbose_df = pd.read_csv(gtex_verbose_addr)
    gtex_tcga_mapping_df = pd.read_csv(gtex_tcga_mapping_addr, sep='\t', names=['Abbreviation','SMTSD'])

    # gtex_verbose_df['Abbreviation'] = gtex_verbose_df[['SMTSD']].merge(gtex_tcga_mapping_df,how='left').Abbreviation
    gtex_verbose_df = pd.merge(gtex_verbose_df,gtex_tcga_mapping_df, on='SMTSD', how='right')

    print(gtex_verbose_df[['SMTSD','Abbreviation']].head())
    print(gtex_tcga_mapping_df)
    gtex_verbose_df = gtex_verbose_df.fillna('NA')
    gtex_verbose_df.to_csv("/home/wch23/Project/LifeArc/TCGA/data/GTEX_v7_xena_samples_verbose_with_abbreviation.csv")

def make_tcga_gtex_id_mapping_file(tcga_gtex_id_df, tcga_gtex_id_addr):

    # print(tcga_gtex_id_df)
    ensembl_id = tcga_gtex_id_df['sample'].str.split(".", n=1, expand=True)
    tcga_gtex_id_df['ensembl_gene'] = ensembl_id[0]

    # print(tcga_gtex_id_df)
    gp = GProfiler(return_dataframe = True)
    ensembl_2_symbol = gp.convert(organism='hsapiens',
                                  query = tcga_gtex_id_df['ensembl_gene'].tolist(),
                                  target_namespace='ENSG')
    # print(ensembl_2_symbol[['incoming','name']])

    tcga_gtex_id_df['gene_symbol'] = tcga_gtex_id_df[['ensembl_gene']].merge(ensembl_2_symbol,how='left',right_on='incoming',left_on='ensembl_gene').name

    # print(tcga_gtex_id_df)

    tcga_gtex_id_df.to_csv(tcga_gtex_id_addr,sep='\t', index=False)




def make_TCGA_GTEx_matrix_each(tcga_gtex_raw_matrix_addr, tcga_verbose_addr, gtex_verbose_addr,
                               cancer_matrix_addr, normal_matrix_addr, tcga_gtex_id_addr, cancer_id_addr, normal_id_addr):
    tcga_gtex_raw_df = pd.read_csv(tcga_gtex_raw_matrix_addr)
    tcga_verbose_df = "/home/wch23/Project/LifeArc/TCGA/data/TGCA_xena_samples_verbose_with_abbreviation.csv"
    # print(list(tcga_gtex_raw_df))
    # print(tcga_gtex_raw_df.head())

    ################
    # Make ID mapping file
    #################
    # tcga_gtex_id_df = tcga_gtex_raw_df[['sample']]
    # make_tcga_gtex_id_mapping_file(tcga_gtex_id_df, tcga_gtex_id_addr)

    ########################
    # Make TCGA Cancer Matrix
    # Remove CNTL ( Controls)
    ########################
    tcga_gtex_raw_df = tcga_gtex_raw_df.set_index('sample')
    # print(list(tcga_gtex_raw_df))
    tcga_col = [col for col in tcga_gtex_raw_df if col.startswith('TCGA')]
    gtex_col = [col for col in tcga_gtex_raw_df if col.startswith('GTEX')]
    gtex_col += [col for col in tcga_gtex_raw_df if col.startswith('K')]

    tcga_raw_df = tcga_gtex_raw_df[tcga_col]
    gtex_raw_df = tcga_gtex_raw_df[gtex_col]

    # print(tcga_raw_df.head())
    # print(gtex_raw_df.head())

    ###########################################
    # Remove control or nomapped columns
    ###########################################
    ## make cancer matrix, Remove CNTL Samples , TCGA
    tcga_verbose_df = pd.read_csv(tcga_verbose_addr)
    tcga_samples = tcga_verbose_df.loc[(tcga_verbose_df['Abbreviation'] == 'CNTL')|(tcga_verbose_df['Short.Letter.Code'] =='NT')]
    cntl_samples = tcga_samples['fullcode']
    # tcga_cancer_df = tcga_raw_df.drop(columns=cntl_samples)
    # tcga_cancer_df = tcga_cancer_df.reset_index()
    # tcga_cancer_df.to_csv(cancer_matrix_addr, index=False)
    # print(tcga_cancer_df.head())

    # ## Remove GTEx unmapped data (to TCGA)
    ## make normal matrix
    gtex_verbose_df = pd.read_csv(gtex_verbose_addr)
    gtex_verbose_df = gtex_verbose_df.fillna('NA')
    print(set(gtex_verbose_df['Abbreviation']))
    gtex_samples_mapping_tcga = gtex_verbose_df.loc[gtex_verbose_df['Abbreviation']!='NA' ]

    tcga_mapped_gtex_samples_all = gtex_samples_mapping_tcga['SAMPID'].tolist()
    tcga_mapped_gtex_samples = [col for col in gtex_raw_df.columns if col in tcga_mapped_gtex_samples_all]

    gtex_mapped_tcga_df = gtex_raw_df[tcga_mapped_gtex_samples]

    tcga_normal_verbose_df = tcga_verbose_df.loc[tcga_verbose_df['Short.Letter.Code'] == 'NT']
    tcga_normal_samples = tcga_normal_verbose_df['fullcode'].tolist()
    tcga_normal_df = tcga_raw_df[tcga_normal_samples]
    # print(tcga_normal_df)

    tcga_gtex_normal_df = pd.concat([tcga_normal_df,gtex_mapped_tcga_df], axis=1)
    # tcga_gtex_normal_df = gtex_mapped_tcga_df.append(tcga_normal_df, ignore_index=True)
    print(tcga_gtex_normal_df)
    tcga_gtex_normal_df = tcga_gtex_normal_df.reset_index()
    tcga_gtex_normal_df.to_csv(normal_matrix_addr,index=False)
    print(tcga_gtex_normal_df.head())
    #
    # ################################
    # # make cancer normal sample_ids
    # #############################
    tcga_cancer_samples = tcga_verbose_df.loc[
        (tcga_verbose_df['Abbreviation'] != 'CNTL') & (tcga_verbose_df['Short.Letter.Code'] != 'NT')]

    tcga_cancer_samples_info_df = tcga_cancer_samples[['fullcode', 'Abbreviation']]
    print(tcga_cancer_samples_info_df)
    tcga_cancer_samples_info_df.to_csv(cancer_id_addr,sep='\t', index=False)

    gtex_normal_samples = gtex_verbose_df[gtex_verbose_df['SAMPID'].isin(tcga_mapped_gtex_samples) ]
    gtex_normal_samples_info_df = gtex_normal_samples[['SAMPID','Abbreviation']]
    gtex_normal_samples_info_df = gtex_normal_samples_info_df.rename(columns={'SAMPID': 'fullcode'})

    tcga_normal_samples_info_df = tcga_normal_verbose_df[['fullcode','Abbreviation']]
    tcga_gtex_normal_samples_info_df = pd.concat([gtex_normal_samples_info_df,tcga_normal_samples_info_df], ignore_index=True)
    print(tcga_gtex_normal_samples_info_df)
    tcga_gtex_normal_samples_info_df.to_csv(normal_id_addr, sep='\t', index=False)

def make_cancer_normal_sample_ids(tcga_verbose_addr, gtex_verbose_addr, cancer_id_addr, normal_id_addr):
    tcga_verbose_df = pd.read_csv(tcga_verbose_addr)
    tcga_cancer_samples = tcga_verbose_df.loc[
        (tcga_verbose_df['Abbreviation'] != 'CNTL') | (tcga_verbose_df['Short.Letter.Code'] != 'NT')]

    tcga_cancer_samples_info_df = tcga_cancer_samples[['fullcode','Abbreviation']]
    print(tcga_cancer_samples_info_df)

    gtex_verbose_df = pd.read_csv(gtex_verbose_addr)
    gtex_verbose_df = gtex_verbose_df.fillna('NA')
    # print(set(gtex_verbose_df['Abbreviation']))
    gtex_samples_mapping_tcga = gtex_verbose_df.loc[gtex_verbose_df['Abbreviation'] != 'NA']

    gtex_normal_samples_info_df = gtex_samples_mapping_tcga[['SAMPID','Abbreviation']]
    gtex_normal_samples_info_df = gtex_normal_samples_info_df.rename(columns={'SAMPID':'fullcode'})

    print(gtex_normal_samples_info_df)

def main():
    #########
    # Step 1 : Add TCGA Abbreviation on TCGA_xena_samples_verbose
    #          GTEx tissue --> TCGA Abbreviation
    # I have already done.
    #########
    # make_TCGA_xena_samples_verbose()
    # make_GTEx_xena_samples_verbose()
    ##########
    # Step 2 : Change Head of TCGA + GTEx matrix --> Seperate to Cancer , Normal(GTEx)
    ##########
    tcga_gtex_raw_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/master_gepis_xena_gtex_tcga.csv"
    tcga_verbose_addr = "/home/wch23/Project/LifeArc/TCGA/data/TCGA_xena_samples_verbose_with_abbreviation.csv"
    gtex_verbose_addr = "/home/wch23/Project/LifeArc/TCGA/data/GTEX_v7_xena_samples_verbose_with_abbreviation.csv"
    cancer_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_xena.csv"
    normal_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_gepia_normal_xena.csv"
    tcga_gtex_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_id_mapping.csv"
    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    normal_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/normal_id_info.tsv"
    # make_TCGA_GTEx_matrix_each(tcga_gtex_raw_matrix_addr, tcga_verbose_addr, gtex_verbose_addr, cancer_matrix_addr, normal_matrix_addr,
    #                            tcga_gtex_id_addr, cancer_id_addr, normal_id_addr)

    ############
    # Step 3 : log value to expression value
    ############

    cancer_matrix_original_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_xena_to_original.csv"
    cancer_matrix_df = pd.read_csv(cancer_matrix_addr,index_col=0)
    print(cancer_matrix_df.head())
    print(cancer_matrix_df.columns)

    cancer_matrix_original_df = cancer_matrix_df.rpow(2)
    cancer_matrix_original_df = cancer_matrix_original_df.replace(0.001,0)
    cancer_matrix_original_df = cancer_matrix_original_df.reset_index()
    cancer_matrix_original_df.to_csv(cancer_matrix_original_addr, index = False)
    print("Cancer Done")
    normal_matrix_original_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_gepia_normal_xena_to_original.csv"
    normal_matrix_df = pd.read_csv(normal_matrix_addr, index_col=0)
    print(normal_matrix_df.head())
    normal_matrix_original_df = normal_matrix_df.rpow(2)
    normal_matrix_original_df = normal_matrix_original_df.replace(0.001, 0)
    normal_matrix_original_df = normal_matrix_original_df.reset_index()
    normal_matrix_original_df.to_csv(normal_matrix_original_addr, index=False)
    print("Normal Done")

    # make_cancer_normal_sample_ids(tcga_verbose_addr, gtex_verbose_addr, cancer_id_addr, normal_id_addr)

def test():
    import numpy as np

    df1 = pd.DataFrame({"A": [14, 4, 5, 4, 1],
                        "B": [5, 2, 54, 3, 0.001],
                        "C": [20, 20, 7, 3, 0.8],
                        "D": [14, 3, 6, 2, 0.6]})

    print(df1)
    df2 = pd.DataFrame({"A": [1, 5, 3, 4, 2]})

    df3 = df1.applymap(np.log2)
    df4 = df3.rpow(2)
    print(df3)
    print(df4)
    df5 = df4.replace(0.001,0)
    print(df5)


if __name__ == '__main__':
    main()
    # test()