import os
import dill
import numpy
import json
import sklearn
import regex
import scipy
import re
from collections import defaultdict
import pickle
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import time
import numpy as np
from torch.autograd import Variable

def Entrez_gene_id2raw_gene_converter(one_Entrez_gene_id,Entrez_gene_id2rawgene_json_filepath):
    with open(Entrez_gene_id2rawgene_json_filepath,'r') as Egj:
        json_loaded = json.load(Egj)

    return json_loaded[one_Entrez_gene_id]

def ngram_char_bow_returner(one_Entrez_gene_raw_name,ngrammax_num):

    to_be_returned_list = list()

    for i in range(1,ngrammax_num+1):
        try :
            for x in range(len(one_Entrez_gene_raw_name)):
                if len(one_Entrez_gene_raw_name) + 1 >= x:
                    n = one_Entrez_gene_raw_name[x:x + i]
                    if len(n) == i:
                        to_be_returned_list.append(n)
                else:
                    continue
        except:
            pass

    return to_be_returned_list

def make_feature_id_dict(ngram_list_list):
    feature_id = defaultdict(lambda: len(feature_id))
    for ngram_list in ngram_list_list:
        for ngram in ngram_list:
            feature_id[ngram]
    return feature_id

def make_feature_ngram_list_list_from_raw_data(train_data_filepath,Entrez_gene_id2rawgene_json_filepath):
    ngram_list_list = list()

    with open(train_data_filepath,'rb') as td:
        train_data = pickle.load(td)

    for one_gene_in_PubTator_and_correct_Entrez_gene_id_set in train_data:
        ngram_list_list.append([one_gene_in_PubTator_and_correct_Entrez_gene_id_set[0],
                                Entrez_gene_id2raw_gene_converter(one_gene_in_PubTator_and_correct_Entrez_gene_id_set[1],
                                                                  Entrez_gene_id2rawgene_json_filepath=Entrez_gene_id2rawgene_json_filepath)])

    return ngram_list_list

def make_vector_from_ngram_list(feature_id, ngram_list):
    # input は n-gram list
    #ex.['派手', 'に', '降っ', 'て', 'き', 'まし', 'た', '伊万里', '市', 'で', 'は', ...2-gram,...3-gram...]

    # defaultdict で語彙を入れると自動的に対応int付けしてくれる
    ngram_dict = defaultdict(lambda:0)
    for ngram in ngram_list:
        if ngram in feature_id.keys():
            #ngram→BoWでベクトル（よくある書き方）
            ngram_dict[feature_id[ngram]] += 1

    #n-gram dict
    #  {409: 1, 1: 2, 78: 1, 11: 1, 81: 1, 210: 1, ...} ただしdefaultdictのインスタンスとして
    I = numpy.zeros(len(ngram_dict.keys()))
    J = numpy.array(list(ngram_dict.keys()))
    V = numpy.array(list(ngram_dict.values()))

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
    vector = scipy.sparse.coo_matrix((V,(I,J)),shape=(1,len(feature_id)))
    # http://techeten.xyz/1004
    return vector

def batch_loader_for_train_or_test(feature_id,train_or_test_file_path,batch_num):

    return # batched_input_tensor_from_raw_gene_id, batched_answe_tensor_from_correct_Entrez_gene_id


