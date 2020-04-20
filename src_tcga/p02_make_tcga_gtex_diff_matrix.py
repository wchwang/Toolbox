# Created by woochanghwang at 17/10/2019

import pandas as pd


def make_TCGA_Cancer_Diff_matrix_for_all_cancer(cancer_matrix_addr, normal_matrix_addr, cancer_id_addr, normal_id_addr, tcga_gtex_id_addr, tcga_cancer_diff_addr):

    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')
    normal_id_df = pd.read_csv(normal_id_addr, sep='\t')
    cancer_types = list(set(cancer_id_df['Abbreviation'].tolist()))

    # print(normal_id_df)
    tcga_cancer_diff_df = pd.read_csv(tcga_gtex_id_addr, sep='\t', index_col=0)

    cancer_matrix_df = pd.read_csv(cancer_matrix_addr, index_col=0)
    normal_matrix_df = pd.read_csv(normal_matrix_addr, index_col=0)


    # print(cancer_matrix_df.head())
    print("Start making cancer matrix with diff")
    cancer_with_normal = []

    for cancer_type in cancer_types:
        a_tcga_cancer_samples_df = cancer_id_df.loc[cancer_id_df['Abbreviation']==cancer_type]
        a_tcga_cancer_samples = a_tcga_cancer_samples_df['fullcode'].tolist()
        print('cancer: ', cancer_type, len(a_tcga_cancer_samples))

        a_tcga_gtex_normal_samples_df = normal_id_df.loc[normal_id_df['Abbreviation']==cancer_type]
        a_tcga_gtex_normal_samples = a_tcga_gtex_normal_samples_df['fullcode'].tolist()
        print('normal: ', cancer_type, len(a_tcga_gtex_normal_samples) )
        if len(a_tcga_gtex_normal_samples) == 0 : continue  # skip if there is no normal samples
        print("Cancer Type: {}".format(cancer_type))
        cancer_with_normal.append(cancer_type)

        a_tcga_cancer_matrix_df = cancer_matrix_df[a_tcga_cancer_samples]
        a_tcga_gtex_normal_matrix_df = normal_matrix_df[a_tcga_gtex_normal_samples]

        # print(a_tcga_cancer_matrix_df.head())
        # print(a_tcga_gtex_normal_matrix_df.head())
        a_tcga_cancer_normal_diff_df = a_tcga_cancer_matrix_df.subtract(a_tcga_gtex_normal_matrix_df.mean(axis=1), axis="index" )

        print(a_tcga_cancer_normal_diff_df.head())

        tcga_cancer_diff_df = pd.concat([tcga_cancer_diff_df,a_tcga_cancer_normal_diff_df], axis=1)


    with open('/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_types_has_normal.tsv','w') as cancer_w_n_f:
        cancer_w_n_f.write('\n'.join(cancer_with_normal))

    tcga_cancer_diff_df.index.name = 'ensembl_id'
    tcga_cancer_diff_df = tcga_cancer_diff_df.reset_index()
    print(tcga_cancer_diff_df.iloc[:,:7])
    print("All Done")
    tcga_cancer_diff_df.to_csv(tcga_cancer_diff_addr, sep='\t', index=False)
    import toolbox.data_handlers as dh
    dh.save_obj(tcga_cancer_diff_df, "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original")


def make_TCGA_Cancer_Diff_matrix_over_basemean_for_all_cancer(cancer_matrix_addr, normal_matrix_addr, cancer_id_addr, normal_id_addr, tcga_gtex_id_addr, tcga_cancer_diff_addr,basemean_th):

    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')
    normal_id_df = pd.read_csv(normal_id_addr, sep='\t')
    cancer_types = list(set(cancer_id_df['Abbreviation'].tolist()))

    # print(normal_id_df)
    tcga_cancer_diff_df = pd.read_csv(tcga_gtex_id_addr, sep='\t', index_col=0)

    cancer_matrix_df = pd.read_csv(cancer_matrix_addr, index_col=0)
    normal_matrix_df = pd.read_csv(normal_matrix_addr, index_col=0)


    # print(cancer_matrix_df.head())
    print("Start making cancer matrix with diff")
    cancer_with_normal = []

    for cancer_type in cancer_types[:2]:
        a_tcga_cancer_samples_df = cancer_id_df.loc[cancer_id_df['Abbreviation']==cancer_type]
        a_tcga_cancer_samples = a_tcga_cancer_samples_df['fullcode'].tolist()
        print('cancer: ', cancer_type, len(a_tcga_cancer_samples))

        a_tcga_gtex_normal_samples_df = normal_id_df.loc[normal_id_df['Abbreviation']==cancer_type]
        a_tcga_gtex_normal_samples = a_tcga_gtex_normal_samples_df['fullcode'].tolist()
        print('normal: ', cancer_type, len(a_tcga_gtex_normal_samples) )
        if len(a_tcga_gtex_normal_samples) == 0 : continue  # skip if there is no normal samples
        print("Cancer Type: {}".format(cancer_type))
        cancer_with_normal.append(cancer_type)

        a_tcga_cancer_matrix_df = cancer_matrix_df[a_tcga_cancer_samples]
        a_tcga_gtex_normal_matrix_df = normal_matrix_df[a_tcga_gtex_normal_samples]

        # print(a_tcga_cancer_matrix_df.head())
        # print(a_tcga_gtex_normal_matrix_df.head())
        a_tcga_cancer_normal_diff_df = a_tcga_cancer_matrix_df.subtract(a_tcga_gtex_normal_matrix_df.mean(axis=1), axis="index" )

        print(a_tcga_cancer_normal_diff_df.head())

        tcga_cancer_diff_df = pd.concat([tcga_cancer_diff_df,a_tcga_cancer_normal_diff_df], axis=1)


    with open('/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_types_has_normal.tsv','w') as cancer_w_n_f:
        cancer_w_n_f.write('\n'.join(cancer_with_normal))

    tcga_cancer_diff_df.index.name = 'ensembl_id'
    tcga_cancer_diff_df = tcga_cancer_diff_df.reset_index()
    print(tcga_cancer_diff_df.iloc[:,:7])
    print("All Done")
    tcga_cancer_diff_df.to_csv(tcga_cancer_diff_addr, sep='\t', index=False)
    import toolbox.data_handlers as dh
    dh.save_obj(tcga_cancer_diff_df, "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original")


def make_TCGA_Cancer_matrix_for_a_cancer(cancer_matrix_addr, cancer_id_addr, tcga_gtex_id_addr, cancer_type, tcga_a_cancer_addr):
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')

    cancer_types = list(set(cancer_id_df['Abbreviation'].tolist()))

    # print(normal_id_df)
    tcga_a_cancer_df = pd.read_csv(tcga_gtex_id_addr, sep='\t', index_col=0)

    cancer_matrix_df = pd.read_csv(cancer_matrix_addr, index_col=0)

    # print(cancer_matrix_df.head())
    print("Start making cancer matrix")

    a_tcga_cancer_samples_df = cancer_id_df.loc[cancer_id_df['Abbreviation'] == cancer_type]
    a_tcga_cancer_samples = a_tcga_cancer_samples_df['fullcode'].tolist()
    print('cancer: ', cancer_type, len(a_tcga_cancer_samples))

    print("Cancer Type: {}".format(cancer_type))

    a_tcga_cancer_matrix_df = cancer_matrix_df[a_tcga_cancer_samples]
    tcga_a_cancer_df = pd.concat([tcga_a_cancer_df, a_tcga_cancer_matrix_df], axis=1)

    print(tcga_a_cancer_df.head())

    tcga_a_cancer_df.to_csv(tcga_a_cancer_addr, sep='\t', index=False)


def make_TCGA_Cancer_for_cancer_normal_matrix(cancer_matrix_addr, normal_matrix_addr, cancer_id_addr, normal_id_addr, tcga_gtex_id_addr, tcga_cancer_addr, tcga_gtex_normal_addr):
    import toolbox.data_handlers as dh

    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')
    normal_id_df = pd.read_csv(normal_id_addr, sep='\t')
    cancer_types = list(set(cancer_id_df['Abbreviation'].tolist()))

    print(normal_id_df)
    tcga_cancer_id_df = pd.read_csv(tcga_gtex_id_addr, sep='\t', index_col=0)

    # cancer_matrix_df = pd.read_csv(cancer_matrix_addr, index_col=0)
    normal_matrix_df = pd.read_csv(normal_matrix_addr, index_col=0)

    # tcga_cancer_df = pd.concat([tcga_cancer_id_df,cancer_matrix_df], axis=1)
    tcga_gtex_normal_df = pd.concat([tcga_cancer_id_df, normal_matrix_df], axis=1)

    # tcga_cancer_df.index.name = 'ensembl_id'
    # tcga_cancer_df = tcga_cancer_df.reset_index()
    # print(tcga_cancer_df.iloc[:,:7])
    # print("Cancer Done")
    # tcga_cancer_df.to_csv(tcga_cancer_addr, sep='\t', index=False)

    # dh.save_obj(tcga_cancer_df, tcga_cancer_addr)

    # tcga_gtex_normal_df.index.name = 'ensembl_id'
    # tcga_gtex_normal_df = tcga_gtex_normal_df.reset_index()
    # normal_id_df = pd.read_csv(normal_id_addr,sep='\t')
    # print(list(normal_id_df))
    normal_samples = normal_id_df['fullcode'].to_list()
    print(normal_samples)
    tcga_gtex_normal_df = tcga_gtex_normal_df[normal_samples]
    print(tcga_gtex_normal_df.iloc[:, :7])
    print("Normal Done")
    # tcga_cancer_df.to_csv(tcga_cancer_addr, sep='\t', index=False)

    dh.save_obj(tcga_gtex_normal_df, tcga_gtex_normal_addr)



def main():
    # cancer_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_xena.csv"
    # normal_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_gepia_normal_xena.csv"
    cancer_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_xena_to_original.csv"
    normal_matrix_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_gepia_normal_xena_to_original.csv"
    tcga_gtex_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_id_mapping.csv"
    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    normal_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/normal_id_info.tsv"
    tcga_cancer_diff_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original.tsv"

    cancer_type = 'CHOL'
    tcga_a_cancer_addr =  "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_{}_matrix.tsv".format(cancer_type)
    # make_TCGA_Cancer_matrix_for_a_cancer(cancer_matrix_addr, cancer_id_addr, tcga_gtex_id_addr, cancer_type, tcga_a_cancer_addr)

    # make_TCGA_Cancer_Diff_matrix_for_all_cancer(cancer_matrix_addr, normal_matrix_addr, cancer_id_addr, normal_id_addr, tcga_gtex_id_addr, tcga_cancer_diff_addr)

    tcga_cancer_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_xena_to_original"
    tcga_gtex_normal_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_gepia_normal_xena_to_original"

    make_TCGA_Cancer_for_cancer_normal_matrix(cancer_matrix_addr, normal_matrix_addr, cancer_id_addr, normal_id_addr, tcga_gtex_id_addr, tcga_cancer_addr, tcga_gtex_normal_addr)




if __name__ == '__main__':
    main()