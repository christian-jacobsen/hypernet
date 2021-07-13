### Configuration file

# PrODE
headers = {
    'main':    '\n[PrODE]: ',
    'submain': '[PrODE]:   ',
    'warning': '[PrODE]:   WARNING! ',
    'val_err': 'from PrODE\n'
}

# TensorFlow setup
tf_setup = {
    'EPSILON':              1.e-15,
    'DTYPE':                'float64',
    'NUM_THREADS':          8,
    'TF_CPP_MIN_LOG_LEVEL': '3'
}

# Multiprocessing setup
mp_setup = {
    'NUM_CPU':  8
}
