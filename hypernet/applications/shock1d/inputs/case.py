# Input file
# >> Case set up
###############################################################################

# Paths
# -----------------------------------------------------------------------------
hypernet = '/home/zanardi/WORKSPACE/CFD/HyperNet/hypernet/'
databases = hypernet + 'hypernet/databases/'
postprocess = './postprocessing/'

# Thermo
# -----------------------------------------------------------------------------
freestream = {
    'p': 500.e0,
    'T': 300.e0,
    'u': 5000.e0
}

mixture = {
    'var': 'X',
    'name': 'MultiComponent',
    'species': {
        'O':    0.05,
        'O2':   [0.95, 0., 0.]
    }
}

thermo = 'CoupledEnergyModes'

# Species
# -----------------------------------------------------------------------------
grouping = 'CBM3'
specie = {
    'thermo_path': databases + 'air/thermo/',
    'grouping': {
        'O2': {
            'name': 'CBM3',
            'path': {
                'rovib_levels': databases + 'grouping/O2/levels/O3_UMN/O2/FromUMN_Sorted.inp',
                'lev_to_bin': databases + 'grouping/O2/grouping/O3_UMN/O2/LevelsMap_' + grouping + '.csv'
            }
        }
    }
}

# Root finding algorithm
# -----------------------------------------------------------------------------
algorithm = {

    'method': 'hybr',

    'x': {
        'start': 0.,
        'end': 0.3
    },

    'dx': {
        'min': 1.e-7,
        'max': 1.e-1,
        'str': 5.e0
    }
}
