
hypernet = '/home/zanardi/WORKSPACE/CFD/HyperNet/hypernet/'
databases = hypernet + 'hypernet/databases/'

data = './dataGen/output_shock/shock.dat'
columns = ['x', 'X_O', 'X_O2_1', 'X_O2_2', 'X_O2_3', 'u', 'T', 'rho', 'p', 'nd', 'H', 'Mf']

#---------------------------
#Free-strem conditions (pressure, temperature, velocity and mass fractions)
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

#---------------------------
# Species database
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

#---------------------------
#Space-grid parameters
algorithm = {

    'method':    'hybr',

    'x': {
        'start': 0.,
        'end':   0.3
    },

    'dx': {
        'min':   1.e-14,
        'max':   5.e-3,
        'str':   1.01e0
    }
}
