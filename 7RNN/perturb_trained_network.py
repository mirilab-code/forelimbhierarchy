# .venv\scripts\activate

import sys
from cmath import inf
import matplotlib.pyplot as plt
import numpy as np
import psychrnn
import seaborn as sb
import tensorflow as tf
from matplotlib.colors import Normalize
from psychrnn.backend.models.basic import Basic
from psychrnn.backend.simulation import BasicSimulator
from sklearn.model_selection import train_test_split

from reachingtask import *
from utils import *

################ code ###################################################

# load params
new_network_params = {'N_batch': 32,
                  'N_in': 3,
                  'N_out': 3,
                  'dt': 1,
                  'tau': 1,
                  'T': 151,
                  'N_steps': 151,
                  'N_rec': 1000,
                  'name': 'perturb',
                  'load_weights_path': './saved_model_params/sym_bad.npz'}

trained_weights = np.load('./saved_model_params/sym_good.npz')
input_batch = np.load('./saved_model_params/input_set.npy')
conditions = np.load('./saved_model_params/input_set_conditions.npy')

trained_keys = list(trained_weights.keys())
for key in trained_keys:
    new_network_params[key] = trained_weights[key]


model = Basic(new_network_params)
simulator = BasicSimulator(rnn_model = model)

print('loaded!')

# run the model on just normal inputs that it's been trained on.
model_output, model_state = model.test(input_batch)
output,states = simulator.run_trials(input_batch)      # I don't know what the difference between these are. I think the simulation does it entirely in numpy(?)