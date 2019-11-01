##############################################################################
# DrugBank XML parser to parse targets of drugs
#
# Created by woochanghwang at 2019-05-29
##############################################################################

import re, io, gzip
import xml.etree.ElementTree as ET
import collections
import requests
import pandas
import json

def parse( drugbank_file ):

    with open(drugbank_file) as xml_file:
        tree = ET.parse(xml_file)
    root = tree.getroot()

    ns = '{http://www.drugbank.ca}'
    inchikey_template = "{ns}calculated-properties/{ns}property[{ns}kind='InChIKey']/{ns}value"
    inchi_template = "{ns}calculated-properties/{ns}property[{ns}kind='InChI']/{ns}value"
    smiles_template = "{ns}calculated-properties/{ns}property[{ns}kind='SMILES']/{ns}value"
    class_template = "{ns}classification/{ns}class"
    target_template = "{ns}targets/{ns}target/{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='GenAtlas']/{ns}identifier"

    rows = list()
    for i, drug in enumerate(root):
        row = collections.OrderedDict()
        assert drug.tag == ns + 'drug'
        row['type'] = drug.get('type')
        row['drugbank_id'] = drug.findtext(ns + "drugbank-id[@primary='true']")
        row['name'] = drug.findtext(ns + "name")
        row['description'] = drug.findtext(ns + "description")
        row['groups'] = [group.text for group in
                         drug.findall("{ns}groups/{ns}group".format(ns=ns))]
        row['atc_codes'] = [code.get('code') for code in
                            drug.findall("{ns}atc-codes/{ns}atc-code".format(ns=ns))]
        row['categories'] = [x.findtext(ns + 'category') for x in
                             drug.findall("{ns}categories/{ns}category".format(ns=ns))]
        row['inchi'] = drug.findtext(inchi_template.format(ns=ns))
        row['inchikey'] = drug.findtext(inchikey_template.format(ns=ns))
        row['smiles'] = drug.findtext(smiles_template.format(ns=ns))
        row['indication'] = drug.findtext(ns + "indication")
        row['class'] = drug.findtext(class_template.format(ns=ns))
        row['target'] = [target.text for target in
                         drug.findall(target_template.format(ns=ns))]
        # Add drug aliases
        aliases = {
            elem.text for elem in
            drug.findall("{ns}international-brands/{ns}international-brand".format(ns=ns)) +
            drug.findall("{ns}synonyms/{ns}synonym[@language='English']".format(ns=ns)) +
            drug.findall("{ns}international-brands/{ns}international-brand".format(ns=ns)) +
            drug.findall("{ns}products/{ns}product/{ns}name".format(ns=ns))

        }
        aliases.add(row['name'])
        row['aliases'] = sorted(aliases)

        rows.append(row)

    alias_dict = {row['drugbank_id']: row['aliases'] for row in rows}
    with open('../data/Drugbank/aliases.json', 'w') as fp:
        json.dump(alias_dict, fp, indent=2, sort_keys=True)

    rows = list(map(collapse_list_values, rows))

    columns = ['drugbank_id', 'name', 'type', 'groups','class', 'atc_codes','target', 'categories', 'inchikey', 'inchi','smiles', 'description','indication']
    drugbank_df = pandas.DataFrame.from_dict(rows)[columns]

    drugbank_df.to_csv("../data/DrugBank/drugbank.tsv",sep='\t',index=False)
    print(drugbank_df.head())

    drugbank_slim_df = drugbank_df[
        drugbank_df.groups.map(lambda x: 'approved' in x) &
        drugbank_df.inchi.map(lambda x: x is not None) &
        drugbank_df.type.map(lambda x: x == 'small molecule')
        ]
    drugbank_slim_df.to_csv("../data/DrugBank/drugbank_slim.tsv", sep='\t', index=False)
    print(drugbank_df.head())

def collapse_list_values(row):
    for key, value in row.items():
        if isinstance(value, list):
            row[key] = '|'.join(value)
    return row


def parse_protein_info(drugbank_file):
    with open(drugbank_file) as xml_file:
        tree = ET.parse(xml_file)
    root = tree.getroot()

    ns = '{http://www.drugbank.ca}'

    protein_rows = list()
    for i, drug in enumerate(root):
        drugbank_id = drug.findtext(ns + "drugbank-id[@primary='true']")
        for category in ['target', 'enzyme', 'carrier', 'transporter']:
            proteins = drug.findall('{ns}{cat}s/{ns}{cat}'.format(ns=ns, cat=category))
            for protein in proteins:
                row = {'drugbank_id': drugbank_id, 'category': category}
                row['organism'] = protein.findtext('{}organism'.format(ns))
                row['known_action'] = protein.findtext('{}known-action'.format(ns))
                actions = protein.findall('{ns}actions/{ns}action'.format(ns=ns))
                row['actions'] = '|'.join(action.text for action in actions)
                uniprot_ids = [polypep.text for polypep in protein.findall(
                    "{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier".format(
                        ns=ns))]
                if len(uniprot_ids) != 1:
                    continue

                row['uniprot_id'] = uniprot_ids[0]
                gene_symbols = [polypep.text for polypep in protein.findall(
                    "{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='GenAtlas']/{ns}identifier".format(
                        ns=ns))]
                if len(gene_symbols) >= 1:
                    row['gene_symbol'] = gene_symbols[0]
                else:
                    row['gene_symbol'] = ''
                # ref_text = protein.findtext("{ns}references[@format='textile']".format(ns=ns))
                # pmids = re.findall(r'pubmed/([0-9]+)', str(ref_text))
                pmids = [pmid.text for pmid in protein.findall("{ns}references/{ns}articles/{ns}article/{ns}pubmed-id".format(ns=ns))]


                if len(pmids) >= 1 and None not in pmids:
                    # print(pmids)

                    row['pubmed_ids'] = '|'.join(pmids)
                else:
                    row['pubmed_ids'] = ''
                protein_rows.append(row)


    protein_df = pandas.DataFrame.from_dict(protein_rows)

    # Read our uniprot to entrez_gene mapping
    response = requests.get('http://git.dhimmel.com/uniprot/data/map/GeneID.tsv.gz', stream=True)
    text = io.TextIOWrapper(gzip.GzipFile(fileobj=response.raw))
    uniprot_df = pandas.read_table(text, engine='python')
    # print(list(uniprot_df))
    uniprot_df.rename(columns={'uniprot': 'uniprot_id', 'GeneID': 'entrez_gene_id'}, inplace=True)

    # merge uniprot mapping with protein_df
    entrez_df = protein_df.merge(uniprot_df, how='inner')

    columns = ['drugbank_id', 'category', 'uniprot_id', 'entrez_gene_id','gene_symbol', 'organism',
               'known_action', 'actions', 'pubmed_ids']
    entrez_df = entrez_df[columns]


    entrez_df.to_csv("../data/DrugBank/drugbank_protein.tsv", sep='\t', index=False)

def drugbank_to_pubchem_mapping():
    import pubchempy

    # Read DrugBank compounds
    drugbank_df = pandas.read_table('/Users/woochanghwang/PycharmProjects/LifeArc/General/data/Drugbank/drugbank.tsv')
    drugbank_df = drugbank_df[-drugbank_df.inchi.isnull()]

    print(drugbank_df.head())

    # map DrugBank compounds to pubchem using InChI
    rows = list()
    for i, row in drugbank_df.iterrows():
        try:
            compounds = pubchempy.get_compounds(row.inchi, namespace='inchi')
        except pubchempy.BadRequestError:
            print('BadRequestError', row)
            continue
        try:
            compound, = compounds
        except ValueError:
            print(row, compounds)
            continue
        row['pubchem_cid'] = compound.cid
        rows.append(row)

    # Create a DataFrame of the mapping
    mapped_df = pandas.DataFrame(rows)
    mapping_df = mapped_df[['drugbank_id', 'pubchem_cid']].dropna()
    mapping_df['pubchem_cid'] = mapping_df['pubchem_cid'].astype(int)
    print(mapping_df.head())

    # Save mapping
    mapping_df.to_csv('/Users/woochanghwang/PycharmProjects/LifeArc/General/data/drugbank_pubchem-mapping.tsv', index=False, sep='\t')

def parse_atc_code(drugbank_file):
    with open(drugbank_file) as xml_file:
        tree = ET.parse(xml_file)
    root = tree.getroot()

    ns = '{http://www.drugbank.ca}'
    # inchikey_template = "{ns}calculated-properties/{ns}property[{ns}kind='InChIKey']/{ns}value"
    # inchi_template = "{ns}calculated-properties/{ns}property[{ns}kind='InChI']/{ns}value"
    # smiles_template = "{ns}calculated-properties/{ns}property[{ns}kind='SMILES']/{ns}value"
    # class_template = "{ns}classification/{ns}class"
    # target_template = "{ns}targets/{ns}target/{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='GenAtlas']/{ns}identifier"

    rows = list()
    for i, drug in enumerate(root):
        row = collections.OrderedDict()
        # assert drug.tag == ns + 'drug'
        # row['type'] = drug.get('type')
        row['drugbank_id'] = drug.findtext(ns + "drugbank-id[@primary='true']")
        row['name'] = drug.findtext(ns + "name")
        # row['description'] = drug.findtext(ns + "description")
        row['groups'] = [group.text for group in
                         drug.findall("{ns}groups/{ns}group".format(ns=ns))]
        row['atc_codes'] = [code.get('code') for code in
                            drug.findall("{ns}atc-codes/{ns}atc-code".format(ns=ns))]
        row['levels'] = [code.get('code')+':'+code.text for code in
                            drug.findall("{ns}atc-codes/{ns}atc-code/{ns}level".format(ns=ns))]
        # row['categories'] = [x.findtext(ns + 'category') for x in
        #                      drug.findall("{ns}categories/{ns}category".format(ns=ns))]
        # row['inchi'] = drug.findtext(inchi_template.format(ns=ns))
        # row['inchikey'] = drug.findtext(inchikey_template.format(ns=ns))
        # row['smiles'] = drug.findtext(smiles_template.format(ns=ns))
        # row['indication'] = drug.findtext(ns + "indication")
        # row['class'] = drug.findtext(class_template.format(ns=ns))
        # row['target'] = [target.text for target in
        #                  drug.findall(target_template.format(ns=ns))]

        # print(row)
        rows.append(row)



    rows = list(map(collapse_list_values, rows))

    columns = ['drugbank_id', 'name', 'groups', 'atc_codes','levels']
    drugbank_df = pandas.DataFrame.from_dict(rows)[columns]

    # drugbank_df.to_csv("../data/DrugBank/drugbank_atc_code.tsv",sep='\t',index=False)
    # print(drugbank_df.head())

    drugbank_slim_df = drugbank_df[
        drugbank_df.groups.map(lambda x: 'approved' in x) &
        drugbank_df.atc_codes.map(lambda x: x is not None)
        ]
    drugbank_slim_df.to_csv("../data/DrugBank/drugbank_atc_code.tsv", sep='\t', index=False)
    print(drugbank_df.head())

def make_atc_code_catagory():
    drug_bank_atc_code_addr = "../data/DrugBank/drugbank_atc_code.tsv"

    with open(drug_bank_atc_code_addr) as atc_code_f:
        drugbank_drugs_atc_code = [x.strip().split('\t') for x in atc_code_f.readlines()]

    atc_codes = []

    for drug in drugbank_drugs_atc_code[1:]:
        if len(drug)<5: continue
        drug[-1] = drug[-1].replace(':','\t')
        atc_levels_for_drug = drug[-1].split('|')
        atc_codes += atc_levels_for_drug
        atc_code_for_drug = drug[3].split('|')
        atc_code_for_drug = [x+'\t'+drug[1] for x in atc_code_for_drug]
        atc_codes += atc_code_for_drug

    atc_codes = list(set(atc_codes))

    atc_codes.sort()

    with open("../data/DrugBank/drugbank_atc_code_hierachi.tsv",'w') as atc_code_f:
        atc_code_f.write('\n'.join(atc_codes))
def main():
    drugbank_file  = "../data/Drugbank/full database.xml"  # test.xml
    # parse(drugbank_file)
    # parse_protein_info(drugbank_file)

    # drugbank_to_pubchem_mapping()

    parse_atc_code(drugbank_file)
    make_atc_code_catagory()
    return




if __name__ == "__main__":
    main()
