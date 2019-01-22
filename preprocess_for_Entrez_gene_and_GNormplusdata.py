import os
import dill
import numpy
import json
import sklearn
import regex
import scipy
import re
from collections import defaultdict
from sklearn.externals import joblib
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.svm import SVC
import MeCab
import subprocess
import pickle
import codecs
from sklearn.metrics import accuracy_score, precision_score, recall_score, make_scorer, f1_score
from sklearn.model_selection import cross_val_score
from sklearn.metrics import classification_report
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.linear_model import LogisticRegression
from sklearn.multiclass import OneVsRestClassifier

def Taxid_remover(one_gene_id_in_Pubtator):
    new_id = re.sub(pattern=r'Tax:[0-9]+',repl='',string=one_gene_id_in_Pubtator)
    new_id = new_id.replace('()','')
    print(new_id)
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
            print(new_returned_set)
        return new_returned_set

def from_Pubtatorfile_Gene_extractor(pubtator_filepath):
    Pubtator_geneset_without_abstract = list()

    with open(pubtator_filepath,'r') as pf:
        for oneline in pf:
            if ('|t|' in oneline) or ('|a|' in oneline) or( oneline.strip() == ''):
                continue
            else:
                if oneline.split('\t')[4] == 'Gene':
                    Pubtator_geneset_without_abstract.append([oneline.split('\t')[3],oneline.split('\t')[5].strip()])

    new_raw_gene_in_Pubtator_and_GOid_set_list = list()

    for oneset in Pubtator_geneset_without_abstract:
        new_oneset = list()
        new_oneset.append(oneset[0])
        new_oneset.append(Taxid_remover(oneset[1]))
        if not ',' in oneset[1]:
            new_raw_gene_in_Pubtator_and_GOid_set_list.append(new_oneset)
        else:
            candidate_gene = oneset[1].split(',')
            for candidate in candidate_gene:
                new_raw_gene_in_Pubtator_and_GOid_set_list.append([oneset[0],candidate])

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
    train_data = from_Pubtatorfile_Gene_extractor(pubtator_filepath=Pubtator_format_trainfile_filepath)
    test_data = from_Pubtatorfile_Gene_extractor(pubtator_filepath=Pubtator_format_test_file_filepath)

    with open(Entrez_gene_id_to_gene_dump_json_path,'w') as Eg:
        json.dump(Entrez_gene_id_to_gene_name_dictionary, Eg,ensure_ascii=False, indent=4)

    with open('./dataset_dir/BC2GNtrain_gene.pkl','wb') as BC2tr:
        pickle.dump(train_data,BC2tr)

    with open('./dataset_dir/BC2GNtest_gene.pkl','wb') as BC2te:
        pickle.dump(train_data,BC2te)