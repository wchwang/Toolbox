# Created by woochanghwang at 17/10/2019

import sys
sys.path.append('/home/wch23/Project/LifeArc/General')

import toolbox.data_handlers as dh
import statsmodels.api as sm
from statsmodels.formula.api import ols
import pandas as pd


def get_selected_genes_using_oneway_ANOVA(TCGA_df, p_value_threshold, cancer_id_df, target_cancer, gene_value_mode, anova_result_addr):
    # anova_res = pd.DataFrame(columns=['Gene','F-score','P-value'])
    TCGA_df = TCGA_df.set_index('gene_symbol')
    TCGA_df_T = TCGA_df.T
    print(TCGA_df_T.head())

    cancer_id_df = cancer_id_df.set_index('fullcode')
    cancer_id_df = cancer_id_df[['Case_Ctrl']]
    print(cancer_id_df.head())
    significant_gene = []

    anova_res_list = []
    for gene in TCGA_df_T.columns:
        try:
            a_gene_cancertype_df = pd.concat((TCGA_df_T[gene], cancer_id_df), axis=1)
            a_gene_cancertype_df['Case_Ctrl'].astype("category")

            ##############
            #model_name = ols('outcome_variable ~ C(group_variable)', data=your_data).fit()
            ##############
            anova_model = 'Q(\"{}\") ~ C(Case_Ctrl)'.format(gene)
            # print("anova_model:", anova_model)
            # print(a_gene_cancertype_df.head())
            mod = ols(anova_model, data = a_gene_cancertype_df).fit()
            aov_table = sm.stats.anova_lm(mod, typ=2)
            # print (aov_table['F'][0], aov_table['PR(>F)'][0])
            # print(aov_table)
            pvalue_anova = aov_table['PR(>F)'][0]
            anova_res_list.append([gene,pvalue_anova])
            if pvalue_anova < p_value_threshold:
                # print ("thr, pvalue:", p_value_threshold, pvalue_anova)
                # significant_gene.append(gene)
                significant_gene.append([gene,pvalue_anova])
        except:
            print('error:' , gene)

    print("sig gene:", len(significant_gene))
    print("whole", TCGA_df_T.columns)

    anova_res_df = pd.DataFrame(anova_res_list,columns=["gene",'pvalue'])


    anova_res_df.to_csv(anova_result_addr,index=False)


    return significant_gene

def get_selected_genes_using_oneway_ANOVA(TCGA_df, p_value_threshold, cancer_id_df, target_cancer, gene_value_mode, anova_result_addr):
    # anova_res = pd.DataFrame(columns=['Gene','F-score','P-value'])
    TCGA_df = TCGA_df.set_index('gene_symbol')
    TCGA_df_T = TCGA_df.T
    print(TCGA_df_T.head())

    cancer_id_df = cancer_id_df.set_index('fullcode')
    cancer_id_df = cancer_id_df[['Case_Ctrl']]
    print(cancer_id_df.head())
    significant_gene = []

    anova_res_list = []
    for gene in TCGA_df_T.columns:
        try:
            a_gene_cancertype_df = pd.concat((TCGA_df_T[gene], cancer_id_df), axis=1)
            a_gene_cancertype_df['Case_Ctrl'].astype("category")
            # print(TCGA_df_T[gene])
            # print(cancer_id_df)
            # print(a_gene_cancertype_df)
            ##############
            #model_name = ols('outcome_variable ~ C(group_variable)', data=your_data).fit()
            ##############
            anova_model = 'Q(\"{}\") ~ C(Case_Ctrl)'.format(gene)
            # print("anova_model:", anova_model)
            # print(a_gene_cancertype_df.head())
            mod = ols(anova_model, data = a_gene_cancertype_df).fit()
            aov_table = sm.stats.anova_lm(mod, typ=2)
            # print (aov_table['F'][0], aov_table['PR(>F)'][0])
            # print(aov_table)
            pvalue_anova = aov_table['PR(>F)'][0]
            anova_res_list.append([gene,pvalue_anova])
            if pvalue_anova < p_value_threshold:
                # print ("thr, pvalue:", p_value_threshold, pvalue_anova)
                # significant_gene.append(gene)
                significant_gene.append([gene,pvalue_anova])
        except:
            print('error:' , gene)

    print("sig gene:", len(significant_gene))
    print("whole", TCGA_df_T.columns)
    # print(anova_res_list)

    anova_res_df = pd.DataFrame(anova_res_list,columns=["gene",'pvalue'])

    anova_res_df.to_csv(anova_result_addr,index=False)


    return significant_gene

def make_anova_result_with_symbol_ensembl(anova_result_addr, tcga_id_mapping_addr, anova_result_with_ensembl_addr):
    tcga_id_mapping_df = pd.read_csv(tcga_id_mapping_addr, sep='\t')

    anova_result_df = pd.read_csv(anova_result_addr)

    anova_id_result= pd.concat([anova_result_df, tcga_id_mapping_df])
    anova_result_df['ensembl_gene'] = anova_result_df[['gene']].merge(tcga_id_mapping_df,left_on='gene', right_on='gene_symbol', how='left').ensembl_gene
    # tcga_verbose_df['Abbreviation'] = tcga_verbose_df[['Study Name']].merge(tcga_abbreviation_df,how='left').Abbreviation

    anova_result_df = anova_result_df[['gene','ensembl_gene','pvalue']]
    anova_result_df.to_csv(anova_result_with_ensembl_addr,sep='\t',index=False)
    # print(anova_result_df.head())
    # print(tcga_id_mapping_df.head())
    # print(anova_id_result.head())

def get_sig_genes_by_anova(target_cancers, p_value_th):
    tcga_cancer_diff_df = dh.load_obj("/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_diff_matrix_original")
    with open('/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_cancer_types_has_normal.tsv') as cancer_type_f:
        tcga_cancer_types = [x.strip() for x in cancer_type_f.readlines()]
    cancer_id_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/cancer_id_info.tsv"
    cancer_id_df = pd.read_csv(cancer_id_addr, sep='\t')
    # print(cancer_id_df)
    # target_cancer_list=["PAAD"] # pancreas, ITGB5
    target_cancer_list=target_cancers
    # target_cancer_list = ["LUSC"]  # Lung squamous cell carcinoma , Sox2
    target_cancer = '_'.join(target_cancer_list)
    cancer_id_df['Case_Ctrl'] = ['Case' if x in target_cancer_list else 'Ctrl' for x in cancer_id_df['Abbreviation']]
    # print(cancer_id_df)
    #
    tcga_cancer_diff_df = tcga_cancer_diff_df[tcga_cancer_diff_df['gene_symbol']!='None']
    tcga_cancer_diff_df = tcga_cancer_diff_df.drop(columns=['ensembl_id','ensembl_gene'])


    # # # print(tcga_cancer_diff_df.head())
    # #
    gene_value_mode = 'Diff'
    # p_value_threshold = 1.0e-40  # -350 = significant gene(0) , -320 = significant(8300)
    p_value_threshold = p_value_th

    base_mean_threshold = 0
    anova_result_addr = "/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_result.csv".format(target_cancer, gene_value_mode)
    selected_genes_l = get_selected_genes_using_oneway_ANOVA(tcga_cancer_diff_df, p_value_threshold, cancer_id_df,
                                                             target_cancer, gene_value_mode, anova_result_addr)
    selected_genes_addr = "/home/wch23/Project/LifeArc/TCGA/Result/significant_genes_original_{}_from_{}_anova_{}.csv".format(gene_value_mode, target_cancer,
                                                                                       str(p_value_threshold))


    sig_genes_df = pd.DataFrame(selected_genes_l, columns=["gene", 'pvalue'])
    sig_genes_df.to_csv(selected_genes_addr, index=False)

    tcga_id_mapping_addr = "/mnt/raid0_data/MTI_DATA/TCGA/GepiaTCGA/tcga_gtex_id_mapping.csv"
    anova_result_with_ensembl_addr = "/home/wch23/Project/LifeArc/TCGA/Result/TCGA_GTEx_original_{}_vs_other_{}_anova_result_with_ensembl.csv".format(target_cancer,gene_value_mode)
    make_anova_result_with_symbol_ensembl(anova_result_addr, tcga_id_mapping_addr, anova_result_with_ensembl_addr)

    ################
    # save significant genes
    ################
    significant_genes_addr = "/home/wch23/Project/LifeArc/TCGA/Result/significant_genes_only_original_{}_from_{}_anova_{}.csv".format(
        gene_value_mode, target_cancer, str(p_value_threshold))
    anova_result_df = pd.read_csv(anova_result_addr)
    sig_genes_df = anova_result_df[anova_result_df['pvalue']<=p_value_threshold]
    # sig_genes = sig_genes_df['gene'].tolist()
    #
    # with open(significant_genes_addr, 'w') as sig_genes_f:
    #     sig_genes_f.write('\n'.join(sig_genes))
    sig_genes_df.to_csv(significant_genes_addr, sep='\t')


def main_PAAD():
    target_cancer_list = ["PAAD"]
    p_value_threshold = 1.0e-40
    get_sig_genes_by_anova(target_cancer_list, p_value_threshold)

def main_LUSC():
    target_cancer_list = ["LUSC"]
    p_value_threshold = 1.0e-100
    get_sig_genes_by_anova(target_cancer_list, p_value_threshold)

def main_LIHC():
    target_cancer_list = ['LIHC']
    p_value_threshold = 1.0e-40
    get_sig_genes_by_anova(target_cancer_list, p_value_threshold)

def main_CHOL():
    target_cancer_list = ['CHOL']
    p_value_threshold = 1.0e-40
    get_sig_genes_by_anova(target_cancer_list, p_value_threshold)

def main_CHC():
    target_cancer_list = ['CHOL','LIHC']
    p_value_threshold = 1.0e-40
    get_sig_genes_by_anova(target_cancer_list, p_value_threshold)

def main_GBM():
    target_cancer_list = ['GBM']
    p_value_threshold = 1.0e-100
    get_sig_genes_by_anova(target_cancer_list, p_value_threshold)

def main_LGG():
    target_cancer_list = ['LGG']
    p_value_threshold = 1.0e-100
    get_sig_genes_by_anova(target_cancer_list, p_value_threshold)

def main():
    tcga_cancers = ["KICH","DLBC","PCPG"]
    for cancer_type in tcga_cancers:
        print(cancer_type)
        target_cancer_list = [cancer_type]
        p_value_threshold = 1.0e-50
        get_sig_genes_by_anova(target_cancer_list, p_value_threshold)

if __name__ == '__main__':
    # main_PAAD()
    # main_LIHC()
    # main_CHOL()
    # main_CHC()
    # main_GBM()
    # main_LGG()
    main()