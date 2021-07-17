###########################################################################
### Execution Flags
###########################################################################
load_flg  = 0       # Loading model     0: no, 1: saved model, 2: ckpt
train_flg = 1       # Training          0: no, 1: yes
test_flg  = 1       # Evaluating        0: no, 1: yes
plot_flg  = 2       # Plotting          0: no, 1: after training, 2: also during training
# video_flg = 0       # Video history     0: no, 1: yes
save_flg  = 1       # Saving model      0: no, 1: yes


###########################################################################
### ODE Properties
###########################################################################
name  = 'Shock1D_OpNN'

x0, x_end = 0., 0.3         # Space/time limits


###########################################################################
### NN Architecture
###########################################################################
net = {
    'architecture': 'FFNN',

    'args': {

        'dim': {
            "inp": [2, 1],     # Branch-Trunk input dimensionality
            "out": 2          # Output dimensionality
        },

        'layer_size': {
            "branch": [dim['inp'][0], 32, 32, 16],
            "trunk":  [dim['inp'][1], 32, 32, 16]
        },

        'activation': {
            "branch": 'elu',
            "trunk":  'tanh'
        },

        'regularization': {
            'name':       'l1+l2',
            'l1':         1.e-5,
            'l2':         1.e-5
        },

        'batch_norm': {
            "branch": None,
            "trunk":  None
        },

        'dropout_rate': 0,

        'kernel_initializer': "GlorotNormal",

        'block': "FFBlock"
    }
}

# I/O Transformation ------------------------------------------------------
input_fn  = None
output_fn = None


###########################################################################
### Plotting Quantities
###########################################################################
# Variables Name ----------------------------------------------------------
var_names = [ r'$O_2^{({%s})}$' % (i+1) for i in range(dim['out']-1) ]
var_names.append(r'$O$')

# Axis --------------------------------------------------------------------
labels   = [r'$t\quad[s]$', r'$Y$']         # x, y labels
scales   = ['log', 'linear']                # x, y scales
log_plot = True
limit_y  = False


###########################################################################
### Data
###########################################################################
data = {
    'type': 'DataSet',

    'args': {
    
        'file_train': None,

        'file_valid': None,

        'file_test': None,

        'valid_split': 0.2,

        'col_x': [[0,2],[2,3]],

        'col_y': [3,5],

        'batch_size': 0
    }
}


###########################################################################
### Training
###########################################################################
epochs       = 10
optimizer    = 'adam'
lr           = 1.e-4

# Learning rate scheduler -------------------------------------------------
decay        = ["exponential", 10000, 0.98]

# Losses ------------------------------------------------------------------
# Uppercase names are custom losses, while lowercase ones are tf.keras.losses
loss         = 'MALPE'
metrics      = ['mse']

# Callbacks ---------------------------------------------------------------
callbacks = {

    'early_stopping': {
        'monitor':              'val_tot_loss',
        'min_delta':            1.e-6,
        'patience':             20000,
        'restore_best_weights': True,
        'verbose':              1
    },

    'lr_tracker': {
        'verbose': 1
    },

    'model_ckpt': {
        'monitor':           'val_tot_loss',
        'save_best_only':    True, 
        'save_weights_only': True
    },

    'plotter': {
        'var_names': var_names,
        'labels':    labels,
        'scales':    scales,
        'log_flg':   log_plot,
        'limit_y':   limit_y,
        'period':    10,
    },

    'tensorboard': {
        'histogram_freq': 0,
        'write_graph':    True,
        'write_images':   True,
        'profile_batch':  0
    },

    'weighted_loss': {
        'name':  'EmpiricalWeightsAdapter',
        'alpha': 0.9,
        'freq':  100
    }

}
