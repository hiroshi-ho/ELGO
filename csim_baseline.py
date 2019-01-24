from Autoencoder import AffineLearner
import json,pickle,scipy.sparse, numpy, random
import scipy as sp
from collections import defaultdict
from Entrezgene_char_ngram_dataloader import ngram_char_bow_returner
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import time
import numpy as np
from torch.autograd import Variable

from Autoencoder import *

def All_gene_extractor_for_tfidf_vector_maked(Allgene_filepath):
    '''
    :param Allgene_filepath: Entrez_gene_ontology_json_filepath
    :return: [gene,gene,...]
    '''
    All_gene_list = []

    with open(Allgene_filepath,'r') as Af:
        for geneid, gene in json.load(Af).items():
            All_gene_list.append(gene)

    return All_gene_list

def Countdict_test(All_gene_list):
    Genecountdict = {}

    for one_gene in All_gene_list:
        if one_gene in Genecountdict:
            Genecountdict[one_gene] +=1
        else:
            Genecountdict[one_gene] = 1

    Sorted_Genecountdict = sorted(Genecountdict.items(), key=lambda x: x[1], reverse=True)
    return Sorted_Genecountdict

if __name__ == '__main__':
    ENTREZ_GENE_ID_JSON = './dataset_dir/All_Data.gene_info.json'
    Dumped_SORTED_DICT_PATH = './dataset_dir/sorted_dict.pkl'

    # sorted_dict_tuple_list = Countdict_test(All_gene_list=All_gene_extractor_for_tfidf_vector_maked(Allgene_filepath=ENTREZ_GENE_ID_JSON))
    #
    # with open(Dumped_SORTED_DICT_PATH,'wb') as DSD:
    #     pickle.dump(sorted_dict_tuple_list,DSD)

    with open(Dumped_SORTED_DICT_PATH,'rb') as DSDload:
        sorted_dict_tuple_list = pickle.load(DSDload)


    print(sorted_dict_tuple_list[0:500])