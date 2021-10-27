# *************************************************************************** #
#                                    PrODE                                    #
# --------------------------------------------------------------------------- #
#                     Machine-Learning-Based Approximators                    #
#                  for Ordinary Differential Equations (ODEs)                 #
#                                                                             #
# *************************************************************************** #

# --------------------------- CONFIGURATION FILE ---------------------------- #


# HyperNet --------------------------------------------------------------------
headers = {
    'main':    '[HyperNet]: ',
    'warning': '[HyperNet]:   WARNING! ',
    'val_err': 'from HyperNet\n>>> '
}

values = {
    'SMALL':    1.e-15,
    'BIG':      1.e+15
}

# TensorFlow ------------------------------------------------------------------
tf_setup = {
    'EPSILON':              1.e-15,
    'DTYPE':                'float64',
    'NUM_THREADS':          16,
    'TF_CPP_MIN_LOG_LEVEL': '3'
}
