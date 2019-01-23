import pickle

with open('feature_from_ontology.pkl','rb') as ff:
    feature_list = pickle.load(ff)
    print(feature_list[0:10000])