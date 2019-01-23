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

def make_ngram_from_Entrez_gene_ontology(max_feature,ngram_minnum,ngram_maxnum,Entrez_gene_ontology_json_filepath):
    with open(Entrez_gene_ontology_json_filepath,'r') as Eg:
        Gene_id_and_Gene_json = json.load(Eg)

    feature_ngram_id = {}

    for one_gene in Gene_id_and_Gene_json.values():
        ngram_feature_list = ngram_char_bow_returner(one_Entrez_gene_raw_name=one_gene,ngram_minnum=ngram_minnum,ngrammax_num=ngram_maxnum)
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

def make_one_feature_vector_from_one_gene_and_feature_id(one_gene_raw_text,ngram_minnum,ngram_maxnum,feature_id_default_dict):
    ngram_dict = defaultdict(lambda:0)
    for one_ngram in ngram_char_bow_returner(one_Entrez_gene_raw_name=one_gene_raw_text,ngram_minnum=ngram_minnum,ngrammax_num=ngram_maxnum):
        if one_ngram in feature_id_default_dict.keys():
            ngram_dict[feature_id_default_dict[one_ngram]] +=1

    # n-gram dict
    #  {409: 1, 1: 2, 78: 1, 11: 1, 81: 1, 210: 1, ...} ただしdefaultdictのインスタンスとして
    I = numpy.zeros(len(ngram_dict.keys()))
    J = numpy.array(list(ngram_dict.keys()))
    V = numpy.array(list(ngram_dict.values()))

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html
    vector = scipy.sparse.coo_matrix((V, (I, J)), shape=(1, len(feature_id_default_dict)), dtype='float64')
    # http://techeten.xyz/1004

    # scipy to numpy
    # https://stackoverflow.com/questions/26576524/how-do-i-transform-a-scipy-sparse-matrix-to-a-numpy-matrix
    return vector.todense()

def tensor_set_maker(feature_id_dict,ngram_minnum,ngram_maxnum,train_or_test_dataset_pkl_path,Entrez_gene_ontology_json_path):
    vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list = [] # [[vec_from_text,vec_from_Entrez_gene_id],[...]...]

    with open(Entrez_gene_ontology_json_path,'r') as Egoj:
        Entrez_gene_ontology_json = json.load(Egoj)

    with open(train_or_test_dataset_pkl_path,'rb') as tdp:
        train_gene_in_pubtator_and_correct_set_list = pickle.load(tdp)
        for one_set in train_gene_in_pubtator_and_correct_set_list:
            gene_text_in_Pubtator = one_set[0]
            gene_text_in_Entrez_gene = Entrez_gene_ontology_json[one_set[1]]

            vector_in_text_in_Pubtator = make_one_feature_vector_from_one_gene_and_feature_id(one_gene_raw_text=gene_text_in_Pubtator,
                                                                                              ngram_minnum=ngram_minnum,
                                                                                              ngram_maxnum=ngram_maxnum,
                                                                                              feature_id_default_dict=feature_id_dict)

            vector_in_Entrez_gene = make_one_feature_vector_from_one_gene_and_feature_id(one_gene_raw_text=gene_text_in_Entrez_gene,
                                                                                         ngram_minnum=ngram_minnum,
                                                                                        ngram_maxnum=ngram_maxnum,
                                                                                        feature_id_default_dict=feature_id_dict)

            vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list.append([torch.LongTensor(vector_in_text_in_Pubtator).to(device),torch.LongTensor(vector_in_Entrez_gene).to(device)])

    return vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list

def Entrez_gene_text_id_and_tensor_set(ngram_minnum,ngram_maxnum,feature_id_dict,Entrez_gene_id_json_path):
    '''
    This func cant be conducted directly, cause of memory exprosion at tensor_list_of_gene
    :param ngram_minnum:
    :param ngram_maxnum:
    :param feature_id_dict:
    :param Entrez_gene_id_json_path:
    :return:
    '''

    tensor_list_of_gene = list()
    id_list = list()

    with open(Entrez_gene_id_json_path, 'r') as EGij:
        Entrez_gene_id_json = json.load(EGij)

    for Entrez_gene_id, gene_itself in Entrez_gene_id_json.items():
        id_list.append(Entrez_gene_id)
        tensor_list_of_gene.append(torch.LongTensor(make_one_feature_vector_from_one_gene_and_feature_id(one_gene_raw_text=gene_itself,
                                                                                                     ngram_minnum=ngram_minnum,
                                                                                ngram_maxnum=ngram_maxnum,
                                                                                feature_id_default_dict=feature_id_dict
                                                                                )).to(device))

    return id_list, tensor_list_of_gene

### batch loader for train and text
def batch_index_list_list_loader(train_or_test_one_by_one_vector_set,batch_num):
    batched_index_list_list = list()
    maxindex = len(train_or_test_one_by_one_vector_set)
    index_list = [i for i in range(maxindex)]
    shuffle_index_list = random.sample(index_list,len(index_list))

    sub_list = list()
    for i in shuffle_index_list:
        sub_list.append(i)
        if len(sub_list) % batch_num == 0 :
            batched_index_list_list.append(sub_list)
            sub_list = list()

    return batched_index_list_list

def one_batch_loader(one_indexes_of_batch,one_tensor_from_text_and_tensor_from_correct_gene_set):
    batched_tensor_from_text = list()
    batched_tensor_from_correct_gene = list()

    for idx in one_indexes_of_batch:
        batched_tensor_from_text.append(one_tensor_from_text_and_tensor_from_correct_gene_set[idx][0])
        batched_tensor_from_correct_gene.append(one_tensor_from_text_and_tensor_from_correct_gene_set[idx][1])


    print(batched_tensor_from_text,batched_tensor_from_correct_gene)
    print(one_indexes_of_batch)
    print(len(one_tensor_from_text_and_tensor_from_correct_gene_set))
    return torch.LongTensor(batched_tensor_from_text).to(device), torch.LongTensor(batched_tensor_from_correct_gene).to(device)
### batch loader end

### Linear Projection model
class AffineLearner(nn.Module):

    def __init__(self,vec_dim_init_batched,vec_dim_projected_batched,batch_size):
        super(AffineLearner,self).__init__()
        self.batch_size = batch_size
        self.fc = nn.Linear(vec_dim_init_batched,vec_dim_projected_batched)
        self.feature_vec_dim = vec_dim_init_batched

    def forward(self, to_be_projected_vector):
        projected = self.fc(to_be_projected_vector.float()).float()
        return projected

    def loss_custom(self,projected_tensor_from_Pubtator_Text_batched,correct_text_tensor_batched):
        correct_text_tensor_batched.view(self.batch_size,-1,1)
        ce_loss = - torch.bmm(projected_tensor_from_Pubtator_Text_batched,correct_text_tensor_batched)

        return ce_loss

def test_loss_evaluator(model,test_one_by_one_vector_set,batch_num):
    model.eval()
    test_batch_index_list_of_list = batch_index_list_list_loader(train_or_test_one_by_one_vector_set=test_one_by_one_vector_set,
                                                                 batch_num=batch_num)
    test_loss_sum = list()

    for idx, one_test_batch_index_list in enumerate(test_batch_index_list_of_list):
        wrapped_test_tensor_set_from_pupator_text, wrapped_test_tensor_set_from_Entrez_gene = one_batch_loader(
            one_indexes_of_batch=one_test_batch_index_list,
            one_tensor_from_text_and_tensor_from_correct_gene_set=test_one_by_one_vector_set)

        test_batched_projected_vectors = model(wrapped_test_tensor_set_from_pupator_text)
        loss_for_test_batch = model.loss_custom(projected_tensor_from_Pubtator_Text_batched=test_batched_projected_vectors,
                                                correct_text_tensor_batched=wrapped_test_tensor_set_from_Entrez_gene)
        test_loss_sum.append(loss_for_test_batch.data[0].cpu().numpy())

    return np.mean(test_loss_sum)

def logger(logger_path,epoch,one_epoch_time,train_loss,test_loss):
    with open(logger_path,'a') as logf:
        written_str = "epoch " + str(epoch) + " train loss: " + str(train_loss) + ' test loss: ' + str(test_loss) + ' epoch time: ' + str(round(one_epoch_time,2))  + '\n'
        logf.write(written_str)


if __name__ == '__main__':
    '''
    Entrez_gene_ontology_json_filepath = './dataset_dir/All_Data.gene_info.json'
    Sorted_feature_ngram_id = make_ngram_from_Entrez_gene_ontology(max_feature=MAX_FEATURE,
                                                                   ngram_maxnum=NGRAM_MAXNUM,
                                                                   Entrez_gene_ontology_json_filepath=Entrez_gene_ontology_json_filepath)
    with open('./dataset_dir/feature_from_ontology_feature_300000.pkl','wb') as ffop:
        pickle.dump(Sorted_feature_ngram_id,ffop)
    '''

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # filepath
    NGRAM_MINNUM = 2
    NGRAM_MAXNUM = 3
    MAX_FEATURE = 1000
    TRAIN_DATASET_PKL = './dataset_dir/BC2GNtrain_gene.pkl'
    TEST_DATASET_PKL = './dataset_dir/BC2GNtest_gene.pkl'
    FEATURE_ID_SORTED_DICT_PATH = './dataset_dir/feature_from_ontology_feature_300000.pkl'
    ENTREZ_GENE_ID_JSON = './dataset_dir/All_Data.gene_info.json'

    TRAIN_TENSOR_DATASET = './model_data/train_tensor_dataset.pkl'
    TEST_TENSOR_DATASET = './model_data/test_tensor_dataset.pkl'
    FEATURE_ID_DEFAULT_DICT = './model_data/feature_id_default_dict.pkl'
    ONTOLOGY_ID_LIST = './model_data/ontology_id_list.pkl'
    ONTOLOGY_TENSOR_LIST_OF_GENE = './model_data/onotology_tensor_list_of_gene.pkl'

    LOG_FILE = './model_data/logs.log'
    MODEL_FILEPATH = './model_data/model.model'

    ### model params
    BATCH_SIZE = 10
    LR = 0.5
    EPOCH_NUM = 100
    ### model params end

    with open(FEATURE_ID_SORTED_DICT_PATH,'rb') as FISD:
        sorted_feature_id_set_list = pickle.load(FISD)

    feature_id_default_dict = make_feature_id_dict_from_sorted_feature_ngram_id(sorted_feature_id=sorted_feature_id_set_list)
    # print(len(feature_id_default_dict))

    # train tensor dataset
    train_vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list = tensor_set_maker(feature_id_dict=feature_id_default_dict,
                                                                                                       ngram_minnum=NGRAM_MINNUM,
                                                                                                    ngram_maxnum=NGRAM_MAXNUM,
                                                                                                    train_or_test_dataset_pkl_path=TRAIN_DATASET_PKL,
                                                                                                    Entrez_gene_ontology_json_path=ENTREZ_GENE_ID_JSON)

    print("TRAIN  VEC TENSORIZED")
    # test tensor dataset
    test_vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list = tensor_set_maker(feature_id_dict=feature_id_default_dict,
                                                                                                    ngram_minnum=NGRAM_MINNUM,
                                                                                                    ngram_maxnum=NGRAM_MAXNUM,
                                                                                                    train_or_test_dataset_pkl_path=TEST_DATASET_PKL,
                                                                                                    Entrez_gene_ontology_json_path=ENTREZ_GENE_ID_JSON)

    print("TEST  VEC TENSORIZED")
    # loaded for test validation
    # ontology makes memory explosion, so timely this is supended.
    # ontology_id_list, ontology_tensor_list_of_gene = Entrez_gene_text_id_and_tensor_set(ngram_minnum=NGRAM_MINNUM,ngram_maxnum=NGRAM_MAXNUM,
    #                                                                                     feature_id_dict=feature_id_default_dict,
    #                                                                                     Entrez_gene_id_json_path=ENTREZ_GENE_ID_JSON)

    # dump
    print("Dump start")
    # torch.save(train_vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list,TRAIN_TENSOR_DATASET)
    # torch.save(test_vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list,TEST_TENSOR_DATASET)
    print("Dump skipped")

    # model
    model = AffineLearner(vec_dim_init_batched=MAX_FEATURE,
                          vec_dim_projected_batched=MAX_FEATURE,
                          batch_size=BATCH_SIZE)

    model = model.to(device=device)
    # loss_function = nn.MSELoss() # L1Loss?
    optimizer = optim.SGD(model.parameters(),lr=LR)

    epoch_list = []
    test_loss_list = []
    train_loss_list = []

    for epoch in range(EPOCH_NUM):

        t1 = time.time()
        batch_index_list_of_list = batch_index_list_list_loader(train_or_test_one_by_one_vector_set=train_vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list,
                                                                batch_num=BATCH_SIZE)
        loss_sum = 0
        batch_tot = 0

        for idx, one_index_list in enumerate(batch_index_list_of_list):
            model.train()
            optimizer.zero_grad()


            wrapped_train_tensor_set_from_pupator_text, wrapped_train_tensor_set_from_Entrez_gene = one_batch_loader(one_indexes_of_batch=one_index_list,
                                                                                                                     one_tensor_from_text_and_tensor_from_correct_gene_set=train_vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list)
            projected_by_model_tensor_batched = model(wrapped_train_tensor_set_from_pupator_text)
            loss = model.loss_custom(projected_tensor_from_Pubtator_Text_batched=projected_by_model_tensor_batched,
                                     correct_text_tensor_batched=wrapped_train_tensor_set_from_Entrez_gene)
            loss.backword()
            optimizer.step()
            loss_sum += loss.data
            batch_tot += BATCH_SIZE
            if batch_tot % 200 == 0:
                print(batch_tot,'trained')

        epoch_list.append(epoch+1)
        train_loss_list.append(loss_sum / len(batch_index_list_of_list) * BATCH_SIZE)
        print("epoch",epoch + 1)
        print("train_loss", loss_sum / (len(batch_index_list_of_list) * BATCH_SIZE))

        ## test_loss_evaluation
        model.eval()
        test_loss = test_loss_evaluator(model=model,
                                        test_one_by_one_vector_set=test_vector_from_in_Pubtator_text_and_correct_vector_from_Entrez_gene_set_list,
                                        batch_num=BATCH_SIZE)
        test_loss_list.append(test_loss)
        t2 = time.time()

        # logger
        time_for_one_epoch = t2 - t1
        one_epoch_train_loss = loss_sum / len(batch_index_list_of_list) * BATCH_SIZE
        one_epoch_test_loss = test_loss
        logger(logger_path=LOG_FILE,one_epoch_time=time_for_one_epoch,train_loss=loss_sum / (len(batch_index_list_of_list) * BATCH_SIZE),
               test_loss=test_loss)
    torch.save(model.state_dict(),MODEL_FILEPATH)
