# *************************************************************************** #
#                                   HyperNet                                  #
#                                                                             #
# *************************************************************************** #

# -------------------------------- INPUT FILE ------------------------------- #
# Description:
# >> Defining global input parameters.
# --------------------------------------------------------------------------- #


# Thermo
# =============================================================================
ambient = {
    'p': 5000.e0,
    'T': 10000.e0,
}

mixture = {
    'var': 'X',
    'name': 'MultiComponent',
    'composition': {
        'O2':   [0.95, 1.e-13, 1.e-13],
        'O':    0.05
    }
}

thermo = 'CoupledEnergyModes'

EOS = 'PerfectGas'

# *** Species
specie = {
    'O2': {
        'rovib': {
            'system': 'O3',
            'PES': 'UMN',
            'grouping': 'CBM3'
        }
    },
    'O': None
}

# Chemistry
# =============================================================================
chemistry = {
    'processFlags': {      # Reactive processes
        'excit':  0,
        'diss':   1
    },
    'model': 'Standard',
    'solver': 'Standard',
    'reactionsList': None,
    'constPV': 'V'
}

# ODE parameters
# =============================================================================
ode = {
    'start': 0.,
    'end': 1.e-3,
    'delta': {
        'min': 1.e-12,
        'max': 1.e-4,
        'str': 1.015
    }
}
