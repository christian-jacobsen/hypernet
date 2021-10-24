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
    'p': 1000.e0,
    'T': 10000.e0,
}

mixture = {
    'var': 'X',
    'name': 'MultiComponent',
    'species': {
        'O2':   [0.95, 0., 0.],
        'O':    0.05
    }
}

thermo = 'CoupledEnergyModes'

# *** Species
species = {
    'O2': {
        'rovib': {
            'system': 'O3',
            'PES': 'UMN',
            'grouping': 'RVE3'
        }
    },
    'O': None
}

# Chemistry
# =============================================================================
processes_mask = {      # Reactive processes
    'excit':  1,
    'diss':   1,
    'recomb': 1
}

# ODE parameters
# =============================================================================
ode = {
    'start': 0.,
    'end': 1.,
    'delta': {
        'min': 1.e-14,
        'max': 1.e-4,
        'str': 1.015
    }
}
