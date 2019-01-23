Train_data_path = 'BC2GNtrain_gene.pkl'
Test_data_path = 'BC2GNtest_gene.pkl'

import re, pickle

with open(Train_data_path,'rb') as Tra:
    train_data = pickle.load(Tra)
    counter = 0
    for one_set in train_data:
        counter += 1
        if 'Tax' in one_set[1]:
            print(one_set[1])
            print(re.sub('\(Tax:.*','',one_set[1]))
    print(counter)

with open(Test_data_path,'rb') as Tra:
    train_data = pickle.load(Tra)
    counter = 0
    for one_set in train_data:
        counter += 1
        if 'Tax' in one_set[1]:
            print(one_set[1])
            print(re.sub('\(Tax:.*','',one_set[1]))
    print(counter)