# *************************************************************************** #
#                                   HyperNet                                  #
# --------------------------------------------------------------------------- #
#                 Machine Learning-Based library for modeling                 #
#           multi-component non-equilibrium thermochemical processes          #
#                                                                             #
# *************************************************************************** #

# -------------------------------- INPUT FILE ------------------------------- #
# Description:
# >> Defining global input parameters for the `fitRates` app
# --------------------------------------------------------------------------- #


# Thermochemistry
# =============================================================================
method = 'CBM' # 'RVE' 'CBM' 'DPM'
n_bins = 3
grouping = method + str(n_bins)

# *** Temperature
T = [
    1500.0e0,
    2500.0e0,
    5000.0e0,
    6000.0e0,
    8000.0e0,
    10000.0e0,
    12000.0e0,
    14000.0e0,
    15000.0e0,
    20000.0e0
]

# *** Species
species = {
    'O2': {
        'rovib': {
            'system': 'O3',
            'PES': 'UMN',
            'grouping': grouping
        }
    },
    'O': None
}

# *** Reactions
reacReader = {
    'path': '/O3_UMN/' + grouping + '/',
    'files': {
        'diss': 'Diss_Corrected.dat',
        'exch': 'Exch_Type1.dat',
        'inel': 'Inel.dat'
    },
    'file_version': 'old'
}

reacWriter = {
    'indices': {
        'diss': 2,
        'exch': 5,
        'inel': 6
    },
    'rate': 'Arrhenius',
    'type': 'MicroReversible'
}
