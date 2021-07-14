
hypernet = '/home/zanardi/WORKSPACE/CFD/HyperNet/hypernet/'
databases = hypernet + 'hypernet/databases/'

data = './shock.dat'
columns = ['x', 'X_O', 'X_O2', 'u', 'T', 'rho', 'p', 'nd', 'H', 'Mf']

#---------------------------
#Free-strem conditions (pressure, temperature, velocity and mass fractions)
freestream = {
    'p': 50.e0,
    'T': 300.e0,
    'u': 6000.e0
}

mixture = {
    'var': 'X',
    'name': 'MultiComponent',
    'species': {
        'O':    0.,
        'O2':   1.
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
            'name': None,
            'path': {
                'rovib_levels': databases + 'O2_grouping/levels/O3_UMN/O2/FromUMN_Sorted.inp',
                'lev_to_bin': databases + 'O2_grouping/grouping/O3_UMN/O2/LevelsMap_' + grouping + '.csv'
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