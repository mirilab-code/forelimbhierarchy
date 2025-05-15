from psychrnn.tasks.task import Task
import numpy as np

# https://psychrnn.readthedocs.io/en/latest/notebooks/NewTask.html
# https://psychrnn.readthedocs.io/en/latest/apidoc/backend.html#module-psychrnn.backend.initializations


class directedreach(Task):
    def __init__(self, dt, tau, T, N_batch):
        super(directedreach, self).__init__(3, 3, dt, tau, T, N_batch)
        # First param: since there are three inputs: laser to cfa and laser to rfa and then the input signal (go cue)
        # Second param: and three outputs: output of cfa, output of rfa, output of both to muscles

    def generate_trial_params(self, batch, trial):
        """"Define parameters for each trial.

        Using a combination of randomness, presets, and task attributes, define the necessary trial parameters.

        Args:
            batch (int): The batch number that this trial is part of.
            trial (int): The trial number of the trial within the batch.

        Returns:
            dict: Dictionary of trial parameters.

        """

        # ----------------------------------
        # Define parameters of a trial
        # ----------------------------------
        params = dict()
        c = np.random.choice([0,1,2])
        laser_amp = 6 # idk maybe the laser should be stronger than the go cue? 6x is from Andrew's 2017 paper
        laser_on = np.load('data/laser_on.npy').flatten()   
        if(c==0):
            # measure from cfa, inactivate rfa. So the cfa output is what we measure when inactivating rfa.
            params['condition'] = 'measurecfa_inactrfa'
            params['cfa_input'] = laser_on * 0
            params['rfa_input'] = laser_on * laser_amp
            params['cfa_output'] = np.load('data/CFA_inactrfa.npy').flatten()

            params['rfa_output'] = np.load('data/RFA_nostim.npy').flatten() * (1-(laser_on*0.9))    # the 0.9 because we don't want cortical activity to go to 0, just like 10% of the max
            # params['rfa_output'] = np.load('data/RFA_nostim.npy').flatten()
            
            params['muscle_output'] = np.load('data/muscle_inactRFA.npy').flatten()
        if(c==1):
            # measure from rfa, inactivate cfa
            params['condition'] = 'measurerfa_inactcfa'
            params['rfa_input'] = laser_on * 0
            params['cfa_input'] = laser_on * laser_amp
            
            params['cfa_output'] = np.load('data/CFA_nostim.npy').flatten() * (1-(laser_on*0.9))
            # params['cfa_output'] = np.load('data/CFA_nostim.npy').flatten()


            params['rfa_output'] = np.load('data/RFA_inactcfa.npy').flatten()
            params['muscle_output'] = np.load('data/muscle_inactCFA.npy').flatten()
        if(c==2):
            # no laser
            params['condition'] = 'no_stim'
            params['cfa_input'] = laser_on * 0
            params['rfa_input'] = laser_on * 0
            params['cfa_output'] = np.load('data/CFA_nostim.npy').flatten()
            params['rfa_output'] = np.load('data/RFA_nostim.npy').flatten()
            params['muscle_output'] = np.load('data/muscle_control.npy').flatten()

        params['go_cue'] = np.load('data/gocue.npy').flatten()

        return params

    def trial_function(self, time, params):
        """ Compute the trial properties at the given time.

        Based on the params compute the trial stimulus (x_t), correct output (y_t), and mask (mask_t) at the given time.

        Args:
            time (int): The time within the trial (0 <= time < T).
            params (dict): The trial params produced generate_trial_params()

        Returns:
            tuple:

            x_t (ndarray(dtype=float, shape=(N_in,))): Trial input at time given params.
            y_t (ndarray(dtype=float, shape=(N_out,))): Correct trial output at time given params.
            mask_t (ndarray(dtype=bool, shape=(N_out,))): True if the network should train to match the y_t, False if the network should ignore y_t when training.

        """
        # stim_noise = 0.1
        # onset = self.T/4.0
        # stim_dur = self.T/2.0

        # ----------------------------------
        # Initialize 
        # ----------------------------------
        x_t = np.zeros(self.N_in)
        y_t = np.zeros(self.N_out)
        mask_t = np.ones(self.N_out)

        noise = np.random.normal(0,0.05,(3,151)) 
        

        # ----------------------------------
        # Retrieve parameters
        # ----------------------------------
        condition = params['condition']
        cfa_input = params['cfa_input']
        rfa_input = params['rfa_input']
        cfa_output = params['cfa_output'] + noise[0,:]
        rfa_output = params['rfa_output'] + noise[1,:]
        muscle_output = params['muscle_output'] + noise[2,:]
        gocue = params['go_cue']

        # ----------------------------------
        # Compute values
        # ----------------------------------

        # time>100 when we care about what the inactivated region is doing DURING inactivation (implementing the cfa_nostim * (1-laser_on) above),
        # time>50 when we only care about what the inactivated region is doing BEFORE inactivation
        if time > 100:
            if(condition == 'measurerfa_inactcfa'):
                mask_t[0] = 0
            if(condition == 'measurecfa_inactrfa'):
                mask_t[1] = 0
            if(condition == 'no_stim'):
                q=1             # we DO care about both

        # print("time={}".format(time))
        x_t[0] = cfa_input[time]
        x_t[1] = rfa_input[time]
        x_t[2] = gocue[time]
        y_t[0] = cfa_output[time]
        y_t[1] = rfa_output[time]
        y_t[2] = muscle_output[time]


        return x_t, y_t, mask_t

############### the built in simple PD task ##################################################

class SimplePD(Task):
    def __init__(self, dt, tau, T, N_batch):
        super(SimplePD, self).__init__(2, 2, dt, tau, T, N_batch)

    def generate_trial_params(self, batch, trial):
        """"Define parameters for each trial.

        Using a combination of randomness, presets, and task attributes, define the necessary trial parameters.

        Args:
            batch (int): The batch number that this trial is part of.
            trial (int): The trial number of the trial within the batch.

        Returns:
            dict: Dictionary of trial parameters.

        """

        # ----------------------------------
        # Define parameters of a trial
        # ----------------------------------
        params = dict()
        params['coherence'] = np.random.exponential(scale=1/5)
        params['direction'] = np.random.choice([0, 1])

        return params

    def trial_function(self, time, params):
        """ Compute the trial properties at the given time.

        Based on the params compute the trial stimulus (x_t), correct output (y_t), and mask (mask_t) at the given time.

        Args:
            time (int): The time within the trial (0 <= time < T).
            params (dict): The trial params produced generate_trial_params()

        Returns:
            tuple:

            x_t (ndarray(dtype=float, shape=(N_in,))): Trial input at time given params.
            y_t (ndarray(dtype=float, shape=(N_out,))): Correct trial output at time given params.
            mask_t (ndarray(dtype=bool, shape=(N_out,))): True if the network should train to match the y_t, False if the network should ignore y_t when training.

        """
        stim_noise = 0.1
        onset = self.T/4.0
        stim_dur = self.T/2.0

        # ----------------------------------
        # Initialize with noise
        # ----------------------------------
        x_t = np.sqrt(2*self.alpha*stim_noise*stim_noise)*np.random.randn(self.N_in)
        y_t = np.zeros(self.N_out)
        mask_t = np.ones(self.N_out)

        # ----------------------------------
        # Retrieve parameters
        # ----------------------------------
        coh = params['coherence']
        direction = params['direction']

        # ----------------------------------
        # Compute values
        # ----------------------------------
        if onset < time < onset + stim_dur:
            x_t[direction] += 1 + coh
            x_t[(direction + 1) % 2] += 1

        if time > onset + stim_dur + 20:
            y_t[direction] = 1.

        if time < onset + stim_dur:
            mask_t = np.zeros(self.N_out)

        return x_t, y_t, mask_t