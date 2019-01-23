import json,pickle,scipy, numpy
from collections import defaultdict
from Entrezgene_char_ngram_dataloader import ngram_char_bow_returner
import torch

def make_ngram_from_Entrez_gene_ontology(max_feature,ngram_maxnum,Entrez_gene_ontology_json_filepath):
    with open(Entrez_gene_ontology_json_filepath,'r') as Eg:
        Gene_id_and_Gene_json = json.load(Eg)

    feature_ngram_id = {}

    for one_gene in Gene_id_and_Gene_json.values():
        ngram_feature_list = ngram_char_bow_returner(one_Entrez_gene_raw_name=one_gene,ngrammax_num=ngram_maxnum)
        for one_feature in ngram_feature_list:
            if one_feature in feature_ngram_id:
                feature_ngram_id[one_feature] +=1
            else:
                feature_ngram_id[one_feature] = 1

    sorted_feature_ngram_id = sorted(feature_ngram_id.items(), key=lambda x: -x[1])[0:max_feature]
    return sorted_feature_ngram_id

def make_feature_id_dict_from_sorted_feature_ngram_id(sorted_feature_id):
    feature_id_default_dict =  defaultdict(lambda: len(feature_id_default_dict))
    for one_sorted_feature in sorted_feature_id:
        feature_id_default_dict[one_sorted_feature[0]]

    return feature_id_default_dict

def make_one_feature_vector_from_one_gene_and_feature_id(one_gene_raw_text,ngram_maxnum,feature_id_default_dict):
    ngram_dict = defaultdict(lambda:0)
    for one_ngram in ngram_char_bow_returner(one_Entrez_gene_raw_name=one_gene_raw_text,ngrammax_num=ngram_maxnum):
        if one_ngram in feature_id_default_dict.kyes():
            ngram_dict[feature_id_default_dict[one_ngram]] +=1

    # n-gram dict
    #  {409: 1, 1: 2, 78: 1, 11: 1, 81: 1, 210: 1, ...} ただしdefaultdictのインスタンスとして
    I = numpy.zeros(len(ngram_dict.keys()))
    J = numpy.array(list(ngram_dict.keys()))
    V = numpy.array(list(ngram_dict.values()))

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
    vector = scipy.sparse.coo_matrix((V, (I, J)), shape=(1, len(feature_id_default_dict)))
    # http://techeten.xyz/1004
    return vector

def train_dataset_maker(feature_id_dict,ngram_maxnum,train_dataset_pkl,Entrez_gene_ontology_json_path):
    vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list = [] # [[vec_from_text,vec_from_Entrez_gene_id],[...]...]

    with open(Entrez_gene_ontology_json_path,'r') as Egoj:
        Entrez_gene_ontology_json = json.load(Egoj)

    with open(train_dataset_pkl,'rb') as tdp:
        train_gene_in_pubtator_and_correct_set_list = pickle.load(tdp)
        for one_set in train_gene_in_pubtator_and_correct_set_list:
            gene_text_in_Pubtator = one_set[0]
            gene_text_in_Entrez_gene = Entrez_gene_ontology_json[one_set[1]]

            vector_in_text_in_Pubtator = make_one_feature_vector_from_one_gene_and_feature_id(one_gene_raw_text=gene_text_in_Pubtator,
                                                                                              ngram_maxnum=ngram_maxnum,
                                                                                              feature_id_default_dict=feature_id_dict)

            vector_in_Entrez_gene = make_one_feature_vector_from_one_gene_and_feature_id(one_gene_raw_text=gene_text_in_Entrez_gene,
                                                                                        ngram_maxnum=ngram_maxnum,
                                                                                        feature_id_default_dict=feature_id_dict)

            vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list.append([torch.tensor(vector_in_text_in_Pubtator),torch.tensor(vector_in_Entrez_gene)])

    return vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list

def Entrez_gene_text_picklizer(ngram_maxnum,feature_id_dict,Entrez_gene_id_json_path):
    return

def ontology_dataset_loader_for_test_evaluator(test_dataset_which_contains_list_of_list_which_has_one_gene_text_and_correct_gene_id,
                                               gene_in_Entrez_ontology_json_path,
                                               pickled_Entrez_gene_vector_path):
    return




if __name__ == '__main__':
    '''
    
    MAX_FEATURE = 300000
    NGRAM_MAXNUM = 5
    Entrez_gene_ontology_json_filepath = './dataset_dir/All_Data.gene_info.json'

    Sorted_feature_ngram_id = make_ngram_from_Entrez_gene_ontology(max_feature=MAX_FEATURE,
                                                                   ngram_maxnum=NGRAM_MAXNUM,
                                                                   Entrez_gene_ontology_json_filepath=Entrez_gene_ontology_json_filepath)
    with open('./dataset_dir/feature_from_ontology_feature_300000.pkl','wb') as ffop:
        pickle.dump(Sorted_feature_ngram_id,ffop)
    
    '''
    FEATURE_ID_SORTED_DICT_PATH = './dataset_dir/feature_from_ontology_feature_300000.pkl'
    with open(FEATURE_ID_SORTED_DICT_PATH,'rb') as FISD:
        sorted_feature_id_set_list = pickle.load(FISD)

    feature_id_default_dict = make_feature_id_dict_from_sorted_feature_ngram_id(sorted_feature_id=sorted_feature_id_set_list)

