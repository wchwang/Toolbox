# Created by woochanghwang at 19/09/2019

# Created by woochanghwang at 13/11/2018

# Created by woochanghwang at 12/11/2018
# Kaplan - Meier plot

import xenaPython as xena
import statistics as stat
import pandas as pd
import numpy as np

from lifelines.statistics import multivariate_logrank_test

def get_codes(host, dataset, fields, data):
    "get codes for enumerations"
    codes = xena.field_codes(host, dataset, fields)
    codes_idx = dict([(x['name'], x['code'].split('\t')) for x in codes if x['code'] is not None])
    for i in range(len(fields)):
        if fields[i] in codes_idx:
            data[i] = [None if v == 'NaN' else codes_idx[fields[i]][int(v)] for v in data[i]]
    return data

def get_fields(host, dataset, samples, fields):
    "get field values"
    data = xena.dataset_fetch(host, dataset, samples, fields)
    return data

def get_fields_and_codes(host, dataset, samples, fields):
    "get fields and resolve codes"
    return get_codes( host, dataset, fields, get_fields( host, dataset, samples, fields))

def  get_TCGA_rna_expressions(host,samples, cohort_symbol,mode, genes=['all']):

    rna_seq_dataset = "TCGA.{}.sampleMap/HiSeqV2_PANCAN".format(cohort_symbol)
    if mode == 'all':
        genes = xena.dataset_field(host, rna_seq_dataset)
    elif mode == 'selected':
        genes = genes

    rna_expresion = get_fields(host, rna_seq_dataset, samples, genes)
    rna_expression_by_genes = dict(zip(genes, rna_expresion))

    # print("expression:", rna_expresion)
    return rna_expression_by_genes

def get_TCGA_OS_phenotypes(host, samples, cohort_symbol):
    os_dataset = "TCGA.{}.sampleMap/{}_clinicalMatrix".format(cohort_symbol, cohort_symbol)

    fields = ['OS', 'OS.time','sample_type']
    values = get_fields_and_codes(host, os_dataset, samples, fields)

    phenotypes = dict(zip(fields, values))

    return phenotypes

def get_TCGA_os_groups(expression, os_group_type):
    first_quartile = np.percentile(list(expression),75)
    median_value = stat.median(list(expression))
    last_quartile = np.percentile(list(expression),25)

    l_expression = list(expression)

    # print(first_quartile)
    # print(median_value)
    # print(last_quartile)

    if os_group_type == 'median':

        os_groups = []
        for sample in l_expression:
            if sample >= median_value:
                os_groups.append('High')
            else:
                os_groups.append('Low')
    elif os_group_type == 'quartiles':
        os_groups = []
        for sample in l_expression:
            if sample >= first_quartile:
                os_groups.append('First')
            elif sample > last_quartile and sample < first_quartile:
                os_groups.append('Middle')
            elif sample <= last_quartile:
                os_groups.append('Last')
        # print("first",first_quartile)

    return os_groups


def tutorial_main():
    ## TCGA Hub  Data

    host = xena.PUBLIC_HUBS['tcgaHub']
    excludeTypes = ['probeMap', 'probemap', 'genePredExt']

    cohort_list = xena.all_cohorts(host, excludeTypes)

    cohort_symbol_list = [x.split('(')[1][:-1] for x in cohort_list]
    # print(cohort_symbol_list)

    cohort= 'TCGA Bile Duct Cancer (CHOL)'
    cohort_symbol = 'CHOL'

    # cohort = "GDC TCGA Bile Duct Cancer (CHOL)"
    # cohort_symbol = 'CHOL'

    samples = xena.cohort_samples(host, cohort, None)

    ## mode = 'all' , 'selected'
    mode = 'selected'
    genes = ['ITGB5']

    rna_expressions_by_genes = get_TCGA_rna_expressions(host, samples, cohort_symbol,mode,genes)

    print(rna_expressions_by_genes[genes[0]])

    # Overall Survival data
    os_phenotypes = get_TCGA_OS_phenotypes(host, samples, cohort_symbol)
    os_result = clean_samples_for_tumour(rna_expressions_by_genes[genes[0]], os_phenotypes)

    os_time = pd.to_numeric(os_result['os_time'],errors='coerce')
    os_event = pd.to_numeric(os_result['os_event'], errors='coerce')

    rna_expression = list(os_result['expression'])

    os_groups = get_TCGA_os_groups(os_result['expression'], os_group_type ='quartiles')


    print(os_groups)
    os_df = pd.DataFrame(
        {
            'os_time' : os_time,
            'os_group': os_groups,
            'os_event': os_event,
        }
    )

    os_df = make_df_for_kaplan(os_result, 'quartiles')
    # os_df = make_df_for_kaplan(os_result, 'median')

    print(os_df.head())
    results = multivariate_logrank_test(os_df['os_time'],os_df['os_group'],os_df['os_event'])
    print (results)
    print(results.p_value)
    print(len(os_df))
    #
    # from lifelines import KaplanMeierFitter
    # from matplotlib import pyplot as plt
    #
    # ax = plt.subplot(111)
    #
    # kmf = KaplanMeierFitter()
    # for group in os_df['os_group'].unique():
    #     data = os_result.get_group(group)
    #     kmf.fit(data["T"], data["E"], label=group)
    #     kmf.plot(ax=ax)

def clean_samples_for_tumour(rna_expression_by_a_gene, os_phenotypes):
    tcga_sample_df = pd.DataFrame(
        {
            'os_time' : os_phenotypes['OS.time'],
            'os_event' : os_phenotypes['OS'],
            'expression' : rna_expression_by_a_gene,
            'sample_type' : os_phenotypes['sample_type']
        }
    )

    sample_type_normal = ['Solid Tissue Normal']


    # print(tcga_sample_df.info())

    #fill up None --> NaN, 'NaN' to NaN for using dropna
    tcga_sample_df.fillna(value=pd.np.nan, inplace=True)
    tcga_sample_df['expression'] = tcga_sample_df['expression'].replace('NaN',np.nan)
    tcga_sample_df = tcga_sample_df[~tcga_sample_df.sample_type.isin(sample_type_normal)]
    # print(list(tcga_sample_df['sample_type']))
    tcga_sample_df.dropna(inplace=True)

    # print(tcga_sample_df['sample_type'].unique())
    # print(list(tcga_sample_df['sample_type']))
    # print(tcga_sample_df.info())

    return tcga_sample_df

def make_df_for_kaplan(os_result, grouping_mode):
    # rna_expression = list(os_result['expression'])

    os_time = pd.to_numeric(os_result['os_time'], errors='coerce')
    os_event = pd.to_numeric(os_result['os_event'], errors='coerce')

    if grouping_mode == 'median':
        os_groups = get_TCGA_os_groups(os_result['expression'], os_group_type=grouping_mode)

        os_df = pd.DataFrame(
            {
                'os_time': os_time,
                'os_group': os_groups,
                'os_event': os_event,
            }
        )
    elif grouping_mode == 'quartiles':
        os_groups = get_TCGA_os_groups(os_result['expression'], os_group_type=grouping_mode)

        # print(os_groups)
        os_df = pd.DataFrame(
            {
                'os_time': os_time,
                'os_group': os_groups,
                'os_event': os_event,
            }
        )
        # select 1/4 , 4/4 samples
        os_df = os_df[~os_df.os_group.isin(['Middle'])]

    return os_df

def kaplan_meier_a_gene_all_TCGA(genes,grouping_mode, result_dir):
    ## TCGA Hub  Data

    host = xena.PUBLIC_HUBS['tcgaHub']
    excludeTypes = ['probeMap', 'probemap', 'genePredExt']

    cohort_list = xena.all_cohorts(host, excludeTypes)

    cohort_list.remove('TCGA Pan-Cancer (PANCAN)')
    cohort_list.remove('TCGA Formalin Fixed Paraffin-Embedded Pilot Phase II (FPPP)')
    cohort_symbol_list = [x.split('(')[1][:-1] for x in cohort_list]

    # cohort = 'TCGA Bile Duct Cancer (CHOL)'
    # cohort_symbol = 'CHOL'

    tcga_os_kaplan_result_dict = dict()
    #################
    ## Temp
    #################
    # cohort_list = ['TCGA Bile Duct Cancer (CHOL)']
    # cohort_symbol_list =['CHOL']


    for cohort in cohort_list:
        cohort_symbol = cohort.split('(')[1][:-1]
        samples = xena.cohort_samples(host, cohort, None)

        # print(cohort, cohort_symbol)

        ## mode = 'all' , 'selected'
        mode = 'selected'

        rna_expressions_by_genes = get_TCGA_rna_expressions(host, samples, cohort_symbol, mode, genes)

        # Overall Survival data
        os_phenotypes = get_TCGA_OS_phenotypes(host, samples, cohort_symbol)

        # remove NaN samples, and select only tumour samples
        os_result = clean_samples_for_tumour(rna_expressions_by_genes[genes[0]], os_phenotypes)

        os_df = make_df_for_kaplan(os_result, grouping_mode)

        # print("sample:",len(os_df))
        results = multivariate_logrank_test(os_df['os_time'], os_df['os_group'], os_df['os_event'])
        # print("p-value:",results.p_value)
        tcga_os_kaplan_result_dict[cohort] = results.p_value

    # for key,value in tcga_os_kaplan_result_dict.items():
    #     print(key, value)

    tcga_os_kaplan_result_df = pd.DataFrame.from_dict(tcga_os_kaplan_result_dict, orient='index')

    tcga_os_kaplan_result_df.to_csv("{}/{}_kaplan_meier_all_TCGA.tsv".format(result_dir, genes[0]),
                                        sep='\t')

def kaplan_meier_all_gene_a_TCGA(cohort, grouping_mode, result_dir):
    ## TCGA Hub  Data

    host = xena.PUBLIC_HUBS['tcgaHub']
    excludeTypes = ['probeMap', 'probemap', 'genePredExt']

    # cohort_list = xena.all_cohorts(host, excludeTypes)
    cohort_symbol = cohort.split('(')[1][:-1]

    rna_seq_dataset = "TCGA.{}.sampleMap/HiSeqV2_PANCAN".format(cohort_symbol)

    all_genes  = xena.dataset_field(host, rna_seq_dataset)  ## all gene sets



    # for a_cohort in cohort_list:
    #     if cohort_symbol in a_cohort:
    #         cohort = a_cohort

    samples = xena.cohort_samples(host, cohort, None)

    mode = 'all'
    rna_expressions_by_genes = get_TCGA_rna_expressions(host, samples, cohort_symbol, mode)

    tcga_os_kaplan_result_dict = dict()

    for gene in  all_genes:
        print("cohort , gene", cohort, gene)
        ## mode = 'all' , 'selected'
        mode = 'selected'

        # Overall Survival data
        os_phenotypes = get_TCGA_OS_phenotypes(host, samples, cohort_symbol)

        # remove NaN samples, and select only tumour samples
        os_result = clean_samples_for_tumour(rna_expressions_by_genes[gene], os_phenotypes)

        os_df = make_df_for_kaplan(os_result, grouping_mode)

        # print("sample:", len(os_df))
        results = multivariate_logrank_test(os_df['os_time'], os_df['os_group'], os_df['os_event'])
        # print("p-value:", results.p_value)
        tcga_os_kaplan_result_dict[gene] = results.p_value


    for key, value in tcga_os_kaplan_result_dict.items():
        print(key, value)

    tcga_os_kaplan_result_df = pd.DataFrame.from_dict(tcga_os_kaplan_result_dict, orient='index')

    tcga_os_kaplan_result_df.to_csv("{}/{}_kaplan_meier_all_genes.tsv".format(result_dir, cohort_symbol),sep='\t')


def main():


    grouping_mode = "quartiles"   # median, quartiles (top25%, down25%)
    # genes = ['SCN5A']  #EIF4A2, PDCD4
    # kaplan_meier_a_gene_all_TCGA(genes, grouping_mode)

    cohort = 'TCGA Glioblastoma (GBM)'
    kaplan_meier_all_gene_a_TCGA(cohort, grouping_mode)


if __name__ == '__main__':
    # tutorial_main()
    main()