# Input file
# >> Chemistry surrogate model
###############################################################################

import prode as pro

from tensorflow.keras import layers
from prode.utilities.model import layers as pro_layers

# Execution Flags
# -----------------------------------------------------------------------------
load_flg  = 1   # Loading model 0: no, 1: saved model, 2: ckpt
train_flg = 0   # Training      0: no, 1: yes
test_flg  = 1   # Evaluating    0: no, 1: yes
plot_flg  = 2   # Plotting      0: no, 1: after training, 2: `1` + during training
save_flg  = 1   # Saving model  0: no, 1: yes

verbose = 0   # Saving model  0: no, 1: yes

# NN Architecture
# -----------------------------------------------------------------------------
net = {
    'architecture': 'FFNN',

    'args': {

        'dim': {
            "inp": 1,
            "out": 4
        },

        'layer_size': [1, 16, 16, 4],

        'activation': 'sigmoid',

        'regularization': {
            'name':       'l1+l2',
            'l1':         1.e-3,
            'l2':         1.e-3
        },

        'kernel_initializer': "GlorotNormal",

        'block': "FFBlock"
    }
}

# >> I/O Transformation
input_fn = pro_layers.Scaling(init_scale=1.e7, name='input_fn', trainable=False)
output_fn = layers.Softmax(axis=-1, name='output_fn')

# Execution Flags
# -----------------------------------------------------------------------------
name = 'Shock1D'
path = './chemSurrogateModel/'

# Plotting
# -----------------------------------------------------------------------------
labels   = [r'$x\quad[m]$', r'$Y$']         # x, y labels
scales   = ['log', 'linear']                # x, y scales
log_plot = True
limit_y  = False

# Data
# -----------------------------------------------------------------------------
path_to_data = './dataGenerator/output_shock/shock.dat'
columns = ['x', 'X_O', 'X_O2_1', 'X_O2_2', 'X_O2_3', 'u', 'T', 'rho', 'p', 'nd', 'H', 'Mf']
data = {
    'type': 'DataSet',

    'args': {

        'source': 'array',

        'valid_split': 0.2,

        'batch_size': 64
    }
}

# Training
# -----------------------------------------------------------------------------
epochs = 30000
optimizer = 'adam'
lr = 1.e-3

# >> Learning rate scheduler
decay = ["exponential", 10000, 0.98]

# >> Losses
loss = {'identifier': 'malpe', 'axis': -1}
metrics = ['mape']

# >> Callbacks
callbacks = {

    'early_stopping': {
        'monitor':              'val_loss',
        'min_delta':            1.e-6,
        'patience':             20000,
        'restore_best_weights': True,
        'verbose':              1
    },

    'lr_tracker': {
        'verbose': 1
    },

    'model_ckpt': {
        'monitor':           'val_loss',
        'save_best_only':    True, 
        'save_weights_only': True
    },

    'plotter': {
        'labels':    labels,
        'scales':    scales,
        'log_flg':   log_plot,
        'limit_y':   limit_y,
        'period':    500,
    },

    'tensorboard': {
        'histogram_freq': 0,
        'write_graph':    True,
        'write_images':   True,
        'profile_batch':  0
    }

}
