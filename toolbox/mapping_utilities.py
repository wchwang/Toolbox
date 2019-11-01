# Created by woochanghwang at 2019-03-29

import mygene

def main():
    pass

def mapping_genes_id_to(gene_list,id_from, id_to, species=9606, as_dataframe = False):
    mg = mygene.MyGeneInfo()
    return mg.querymany(gene_list, scopes=id_from ,fields= id_to ,species=species,as_dataframe=as_dataframe)

def get_gene_info(gene):
    #entrezgene, ensembl.gene,symbol
    mg = mygene.MyGeneInfo()
    return mg.getgene(gene,fields='name,symbol,entrezgene,taxid')

def convert_ensembl_to_entrez(source_ensembl_list):
    '''
    convert id (ENSEMBL --> Gene Entrez)
    :param source_ensembl_list: ENSEMBL LIST [ 'ENSXXX','ENSXXX',...]
    :return: dictionary[ENSEMBL] = Gene Entrez
    '''
    mg = mygene.MyGeneInfo()

    # tmp_gene = ['ENSP00000000233', 'ENSP00000263431', 'ENSP00000353863', 'ENSP00000342026', 'ENSP00000240874']
    # tmp_gene2 = ['ENSG00000148795', 'ENSG00000165359', 'ENSG00000150676']

    source_geneEntrez = mg.querymany(source_ensembl_list, scopes='ensembl.protein', fields='entrezgene', species='human')
    ensembl_to_symbol_dict = dict()

    na_genes = []
    for g in source_geneEntrez:

        s_query = g['query']
        s_symbol = g.get('entrezgene','NA')
        if s_symbol == 'NA': na_genes.append(s_query)
        ensembl_to_symbol_dict[s_query] = s_symbol

    print("NA:", len(na_genes), na_genes)
    return ensembl_to_symbol_dict

if __name__ == '__main__':
    main()