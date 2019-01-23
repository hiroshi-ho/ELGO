import json
import re
import pickle

def Taxid_remover(one_gene_id_in_Pubtator):
    new_id = re.sub('\(Tax:.*','',one_gene_id_in_Pubtator)
    return new_id

def many_entity_candidate_splitter(one_raw_gene_and_candidate_gene_id_set):
    raw_gene_name_in_Pubtator_txt = one_raw_gene_and_candidate_gene_id_set[0]
    one_gene_or_many_gene = one_raw_gene_and_candidate_gene_id_set[1].split(',')
    if len(one_gene_or_many_gene) == 1:
        return one_raw_gene_and_candidate_gene_id_set
    else:
        new_returned_set =list()
        for one_gene in one_gene_or_many_gene:
            new_returned_set +=[raw_gene_name_in_Pubtator_txt,one_gene]

        return new_returned_set

def from_Pubtatorfile_Gene_extractor(pubtator_filepath,ontology_json_filepath):
    Pubtator_geneset_without_abstract = list()

    with open(pubtator_filepath,'r') as pf:
        for oneline in pf:
            if ('|t|' in oneline) or ('|a|' in oneline) or( oneline.strip() == ''):
                continue
            else:
                if oneline.split('\t')[4] == 'Gene':
                    Pubtator_geneset_without_abstract.append([oneline.split('\t')[3],oneline.split('\t')[5].strip()])

    new_raw_gene_in_Pubtator_and_GOid_set_list = list()

    with open(ontology_json_filepath,'r') as ojf:
        ontology_json = json.load(ojf)

    for oneset in Pubtator_geneset_without_abstract:
        new_oneset = list()
        gene_in_text = oneset[0]
        id_of_Entrezgene = Taxid_remover(oneset[1])

        ontology_ids = ontology_json.keys()
        if id_of_Entrezgene not in ontology_ids: #if no correct ans, just ignore
            continue

        new_oneset.append(gene_in_text)
        new_oneset.append(id_of_Entrezgene)

        if not ',' in id_of_Entrezgene:
            new_raw_gene_in_Pubtator_and_GOid_set_list.append(new_oneset)
        else:
            candidate_gene = id_of_Entrezgene.split(',')
            for candidate in candidate_gene:
                new_raw_gene_in_Pubtator_and_GOid_set_list.append([gene_in_text,candidate])

    return new_raw_gene_in_Pubtator_and_GOid_set_list

def Entrez_gene_and_gene_name_dictionary_maker(Entrez_gene_Allinfo_filepath):
    counter = 0
    Entrez_gene_id_to_gene_name_dictionary = {}

    with open(Entrez_gene_Allinfo_filepath,'r') as EgA:
        for line in EgA:
            counter += 1
            if counter == 1:
                continue
            oneline_list = line.strip().split('\t')
            Entrez_gene_id_to_gene_name_dictionary[oneline_list[1]] = oneline_list[2] # Entrez_gene_id : gene_name
    return Entrez_gene_id_to_gene_name_dictionary

if __name__ =='__main__':

    Entrez_gene_id_to_gene_dump_json_path = './dataset_dir/All_Data.gene_info.json'
    Entrez_gene_id_to_gene_name_dictionary = Entrez_gene_and_gene_name_dictionary_maker(Entrez_gene_Allinfo_filepath='./dataset_dir/All_Data.gene_info')
    Pubtator_format_trainfile_filepath = './dataset_dir/GNormPlusCorpus/BC2GNtrain.PubTator.txt'
    Pubtator_format_test_file_filepath = './dataset_dir/GNormPlusCorpus/BC2GNtest.PubTator.txt'

    with open(Entrez_gene_id_to_gene_dump_json_path,'w') as Eg:
        json.dump(Entrez_gene_id_to_gene_name_dictionary, Eg,ensure_ascii=False, indent=4)


    train_data = from_Pubtatorfile_Gene_extractor(pubtator_filepath=Pubtator_format_trainfile_filepath,ontology_json_filepath=Entrez_gene_id_to_gene_dump_json_path)
    test_data = from_Pubtatorfile_Gene_extractor(pubtator_filepath=Pubtator_format_test_file_filepath,ontology_json_filepath=Entrez_gene_id_to_gene_dump_json_path)


    with open('./dataset_dir/BC2GNtrain_gene.pkl','wb') as BC2tr:
        pickle.dump(train_data,BC2tr)

    with open('./dataset_dir/BC2GNtest_gene.pkl','wb') as BC2te:
        pickle.dump(train_data,BC2te)