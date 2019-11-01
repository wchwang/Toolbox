# Created by woochanghwang at 19/09/2019
'''
gene set functional annotation
'''

from gprofiler import GProfiler

def Functional_profiling(gene_list,
                         organism='hsapiens',
                         sources=["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP"],
                         user_threshold=0.05):
    gp = GProfiler(return_dataframe=True)

    gp_result_df = gp.profile(query=gene_list,organism=organism, user_threshold=user_threshold,no_iea=True, sources=sources)

    return gp_result_df

def main():
    gene_list_addr = "/Users/woochanghwang/PycharmProjects/LifeArc/MALAT/Data/RBP.TableS6.tsv"

    with open(gene_list_addr) as gene_list_f:
        gene_list = [x.strip() for x in gene_list_f.readlines()]

    sources = ['KEGG']
    gp_result_df = Functional_profiling(gene_list[:50],sources=sources)

    print(gp_result_df.head())

if __name__ == '__main__':
    main()