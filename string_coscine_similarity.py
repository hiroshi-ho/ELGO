from featureid_generator_from_all_data_gene_info_json import AffineLearner
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

