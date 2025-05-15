import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
import tensorflow as tf
from matplotlib.colors import Normalize


################ utility functions ######################################
def plot_weights(weights, title=""):
    h = sb.heatmap(weights, center=0)
    plt.title(title)
    plt.show()

def addd(x,y):
    return x+y

def performance_measure(trial_batch, trial_y, output_mask, output, epoch, losses, verbosity):
    # return pds[len(pds)-1].accuracy_function(trial_y, output, output_mask)

    if(len(losses)==0):
        last_loss = [0]
        perf = -100
        min_loss = -100
    else:
        last_loss = losses[-1]
        perf = -1*last_loss
        min_loss = np.min(losses)


    print(losses[-5:],min_loss)


    return perf

def sym_reg(model,params):
    w = model.get_effective_W_rec()
    AtoB = w[500:,:500]
    BtoA = w[:500,500:]

    # if you just want to do the sum do this
    AtoB_sum = tf.math.reduce_sum(AtoB)
    BtoA_sum = tf.math.reduce_sum(BtoA)
    dif_sum = tf.math.abs(AtoB_sum - BtoA_sum)
    reg = dif_sum

    # Let's try to do the abs difference of the cumulative sum. This... does not work.
    # AtoB_cs = tf.math.cumsum(AtoB)
    # BtoA_cs = tf.math.cumsum(BtoA)
    # abs_dif_cs = tf.math.abs(AtoB_cs - BtoA_cs)
    # reg = abs_dif_cs

    return reg


def mylossfunc(predictions, y, output_mask):
    """ error .

    Args:
        predictions (*tf.Tensor(dtype=float, shape =(*:attr:`N_batch`, :attr:`N_steps`, :attr:`N_out` *))*): Network output.
        y (*tf.Tensor(dtype=float, shape =(*?, :attr:`N_steps`, :attr:`N_out` *))*): Target output.
        output_mask (*tf.Tensor(dtype=float, shape =(*?, :attr:`N_steps`, :attr:`N_out` *))*): Output mask for :attr:`N_batch` trials. True when the network should aim to match the target output, False when the target output can be ignored.

    Returns:
        tf.Tensor(dtype=float): error of experimental data to model output.

    """
    return tf.reduce_mean(input_tensor=tf.square(output_mask * (predictions - y)))      # this is what's built in so I know it works
    # return tf.reduce_sum(input_tensor=tf.square(output_mask * (predictions - y)))

def plot_trial(trial_num,trial_params,model_output):
    plt.subplot(3,1,1)
    plt.plot(trial_params[trial_num]['cfa_output'])
    plt.plot(model_output[trial_num,:,0])
    plt.title('cfa')

    plt.subplot(3,1,2)
    plt.plot(trial_params[trial_num]['rfa_output'])
    plt.plot(model_output[trial_num,:,1])
    plt.title('rfa')

    plt.subplot(3,1,3)
    plt.plot(trial_params[trial_num]['muscle_output'])
    plt.plot(model_output[trial_num,:,2])
    plt.title('EMG')

    plt.suptitle(trial_params[trial_num]['condition'])
    
    plt.show()

def plot_activity(state):
    plt.subplot(2,1,1)
    plt.plot(state[:,:400], '-b', alpha=0.1)
    plt.plot(state[:,400:500], '-c', alpha=0.1)
    plt.title('CFA')

    plt.subplot(2,1,2)
    plt.plot(state[:,500:900], '-r', alpha=0.1)
    plt.plot(state[:,900:1000], '-m', alpha=0.1)
    plt.title('RFA')

    plt.show()


def plot_model_vs_data(trial_params,model_output,sym_constr):
    conds = [_['condition'] for _ in trial_params]
    conds = np.array(conds)
    no_stim_idx = np.where(conds=='no_stim')[0][0]
    measurecfa_inactrfa_idx = np.where(conds=='measurecfa_inactrfa')[0][0]
    measurerfa_inactcfa_idx = np.where(conds=='measurerfa_inactcfa')[0][0]

    data_nostim = trial_params[no_stim_idx]
    data_measurecfa_inactrfa = trial_params[measurecfa_inactrfa_idx]
    data_measurerfa_inactcfa = trial_params[measurerfa_inactcfa_idx]

    # muscles
    plt.subplot(3,3,1)
    plt.plot(data_nostim['muscle_output'], '-k')
    plt.plot(model_output[no_stim_idx,:,2], '-c')
    muscle_nostim_perf = np.sum(abs(data_nostim['muscle_output'] - model_output[no_stim_idx,:,2]))
    plt.ylabel('EMG')
    plt.title('no stim')

    plt.subplot(3,3,2)
    plt.plot(data_measurecfa_inactrfa['muscle_output'], '-k')
    plt.plot(model_output[measurecfa_inactrfa_idx,:,2], '-c')
    muscle_inactrfa_perf = np.sum(abs(data_measurecfa_inactrfa['muscle_output'] - model_output[measurecfa_inactrfa_idx,:,2]))
    plt.title('measure cfa inactivate rfa')

 
    plt.subplot(3,3,3)
    plt.plot(data_measurerfa_inactcfa['muscle_output'], '-k')
    plt.plot(model_output[measurerfa_inactcfa_idx,:,2], '-c')
    muscle_inactcfa_perf = np.sum(abs(data_measurerfa_inactcfa['muscle_output'] - model_output[measurerfa_inactcfa_idx,:,2]))
    plt.title('measure rfa inactivate cfa')


    # cfa
    plt.subplot(3,3,4)
    plt.plot(data_nostim['cfa_output'], '-k')
    plt.plot(model_output[no_stim_idx,:,0], '-c')
    cfa_nostim_perf = np.sum(abs(data_nostim['cfa_output'] - model_output[no_stim_idx,:,0]))
    plt.ylabel('cfa')

    plt.subplot(3,3,5)
    plt.plot(data_measurecfa_inactrfa['cfa_output'], '-k')
    plt.plot(model_output[measurecfa_inactrfa_idx,:,0], '-c')
    cfa_inactrfa_perf = np.sum(abs(data_measurecfa_inactrfa['cfa_output'] - model_output[measurecfa_inactrfa_idx,:,0]))

    
    plt.subplot(3,3,6)
    plt.plot(model_output[measurerfa_inactcfa_idx,:,0], '-c')
    
    # rfa
    plt.subplot(3,3,7)
    plt.plot(data_nostim['rfa_output'], '-k')
    plt.plot(model_output[no_stim_idx,:,1], '-c')
    rfa_nostim_perf = np.sum(abs(data_nostim['rfa_output'] - model_output[no_stim_idx,:,1]))
    plt.ylabel('rfa')


    plt.subplot(3,3,8)
    plt.plot(model_output[measurecfa_inactrfa_idx,:,1], '-c')
    
    plt.subplot(3,3,9)
    plt.plot(data_measurerfa_inactcfa['rfa_output'], '-k')
    plt.plot(model_output[measurerfa_inactcfa_idx,:,1], '-c')
    rfa_inactcfa_perf = np.sum(abs(data_measurerfa_inactcfa['rfa_output'] - model_output[measurerfa_inactcfa_idx,:,1]))




    performance = muscle_nostim_perf + muscle_inactrfa_perf + muscle_inactcfa_perf + cfa_nostim_perf + cfa_inactrfa_perf + rfa_nostim_perf + rfa_inactcfa_perf
    performance = np.round(performance,3)

    fname = 'performance_{}.png'.format(performance)
    plt.show(block=False)
    if(sym_constr == True):
        plt.savefig('output_symmetric/' + fname)
    else:
        plt.savefig('output_unconstrained/' + fname)

    return(performance)


def save_activity(model_state,trial_params):

    conds = [_['condition'] for _ in trial_params]
    ns = np.where(np.array(conds)=='no_stim')[0][0]
    mCiR = np.where(np.array(conds)=='measurecfa_inactrfa')[0][0]
    mRiC = np.where(np.array(conds)=='measurerfa_inactcfa')[0][0]

    no_stim = model_state[ns,:,:]
    measurecfa_inactrfa = model_state[mCiR,:,:]
    measurerfa_inactcfa = model_state[mRiC,:,:]

    act = np.dstack((no_stim,measurecfa_inactrfa,measurerfa_inactcfa))

    return(act)
    