# .venv\scripts\activate

import matplotlib.pyplot as plt
import numpy as np
import psychrnn
import seaborn as sb
import tensorflow as tf
from matplotlib.colors import Normalize
from psychrnn.backend.models.basic import Basic
from sklearn.model_selection import train_test_split

from reachingtask import *
from utils import *

################ code ###################################################

W_mask = np.load('data/w.npy')
np.fill_diagonal(W_mask,0)
# muscle_control = np.load('data/muscle_control.npy').flatten()
# muscle_inactCFA = np.load('data/muscle_inactCFA.npy').flatten()
# muscle_inactRFA = np.load('data/muscle_inactRFA.npy').flatten()
# CFA_nostim = np.load('data/CFA_nostim.npy').flatten()
# CFA_inact = np.load('data/CFA_inact.npy').flatten()
# RFA_nostim = np.load('data/RFA_nostim.npy').flatten()
# RFA_inact = np.load('data/RFA_inact.npy').flatten()
# laser_on = np.load('data/laser_on.npy').flatten()
# laser_off = laser_on * 0

Nneurons = np.shape(W_mask)[0]

# make a mask with 1s of the weights we don't want to change during training. in this case, any weight that is 0 we want to stay 0.
W_rec_fixed = W_mask*0
W_rec_fixed[W_mask==0] = 1
# Now fix the intERregional weights, just to see what the model does when constraining that way.
# W_rec_fixed[:500,501:] = 1
# W_rec_fixed[501:,:500] = 1

# Dale principle, specify which recurrent neurons are exc or inh so the training doesn't change them.
sign = np.mean(W_mask,0)
sign[sign>0] = 1
sign[sign<0] = -1
dale_rec = np.diag(sign)

dr = directedreach(dt = 1, tau = 1, T = 151, N_batch = 32)
network_params = dr.get_task_params() # get the params passed in and defined in dr
network_params['name'] = 'model'
network_params['N_rec'] = Nneurons
N_in = network_params['N_in']
N_rec = network_params['N_rec']
N_out = network_params['N_out']

randweights = np.random.rand(Nneurons,Nneurons) * 0.1
W = np.multiply(randweights,W_mask)


network_params['autapses'] = False
# set the network's recurrent connectivity to this new architecture
network_params['W_rec'] = W
network_params['rec_connectivity'] = np.abs(W_mask)
# to satisfy dale's principle, make is so that the exc and inh neurons remain of the same type and don't change during training.
network_params['Dale_rec'] = dale_rec
network_params['dale_ratio'] = 0.8

##### INPUTS #######
# make two inputs go to either cfa or rfa
input_connectivity = np.zeros((N_rec, N_in))
input_connectivity[400:500,0] = 1   # first input only goes to CFA interneurons
input_connectivity[900:1000,1] = 1   # second input only goes to RFA interneurons
input_connectivity[:,2] = 1   # go cue input goes to all the neurons
network_params['input_connectivity'] = input_connectivity
W_in = input_connectivity
network_params['W_in'] = W_in

##### OUTPUTS #######
# output connectivity
# and the three outputs are either from CFA exc neurons, RFA exc neurons, or (CFA&RFA) exc neurons for muscle EMG
output_connectivity = np.zeros((N_out, N_rec))
output_connectivity[0,:100] = 1     # CFA neurons to CFA output
output_connectivity[1,500:600] = 1  # RFA neurons to RFA output
output_connectivity[2,100:200] = 1
output_connectivity[2,600:700] = 1  # 100 neurons from CFA and 100 from RFA are the muscle outputs
network_params['output_connectivity'] = output_connectivity
# output weights
W_out = output_connectivity
network_params['W_out'] = W_out



# fix the weights
W_in_fixed = np.ones((N_rec,N_in)) # fix these since we want the laser to only have positive weights (ie only activate in inh neurons)
W_out_fixed = np.ones((N_out, N_rec))
W_out_fixed[2,100:200] = 0    # don't fix the weights that creates the muscle outputs
W_out_fixed[2,600:700] = 0    # "

# define the custom loss function
network_params["loss_function"] = "mylossfunc"
network_params[network_params["loss_function"]] = mylossfunc


train_params = {}
train_params['fixed_weights'] = {
    'W_in': W_in_fixed,
    'W_out': W_out_fixed,
    'W_rec': W_rec_fixed
}
#### stuff from the documentation, not sure what should be used or not
train_params['save_weights_path'] =  'saved_models' # Where to save the model after training. Default: None
train_params['training_iters'] = 400000 # number of iterations to train for Default: 50000
train_params['learning_rate'] = .006 # Sets learning rate if use default optimizer Default: .001
# train_params['loss_epoch'] = 10 # Compute and record loss every 'loss_epoch' epochs. Default: 10
# train_params['verbosity'] = False # If true, prints information as training progresses. Default: True
# train_params['save_training_weights_epoch'] = 100 # save training weights every 'save_training_weights_epoch' epochs. Default: 100
# train_params['training_weights_path'] = None # where to save training weights as training progresses. Default: None
# train_params['optimizer'] = tf.compat.v1.train.AdamOptimizer(learning_rate=train_params['learning_rate']) # What optimizer to use to compute gradients. Default: tf.train.AdamOptimizer(learning_rate=train_params['learning_rate'])
# train_params['clip_grads'] = True # If true, clip gradients by norm 1. Default: True


model = Basic(network_params)
w_pretrain = model.get_weights()


# ---------------------- Train a basic model ---------------------------
model.train(dr, train_params) # train model to perform pd task
w_posttrain = model.get_weights()

# ---------------------- Test the trained model ---------------------------
x,target_output,mask, trial_params = dr.get_trial_batch() # get dr task inputs and outputs
model_output, model_state = model.test(x) # run the model on input x
conds = [_['condition'] for _ in trial_params]

# ---------------------- check some things --------------------------
Q = w_posttrain['W_rec'].copy()
bQ = Q.copy()
bQ[bQ>0] = 1
bQ[bQ<0] = -1

# check to see if the signs of the weights stayed the same during training.
check = W_mask - bQ
print('sum of the differences between W_mask and signs of new weights is: {}'.format(np.sum(check)))

# ---------------------- save stuff ---------------------------------
random_id = str(np.random.randint(10000))
savename = 'newW_' + random_id + '.npy'
np.save(savename,Q)
np.save('oldW.npy',w_pretrain['W_rec'])


# ---------------------- Plot the results ---------------------------
trial_num = 0
print(conds[trial_num])
plot_model_vs_data(trial_params,model_output)
plot_activity(model_state[trial_num,:,:])




# ---------------------- Teardown the model -------------------------
model.destruct()
