# Input file
###############################################################################

# Paths
# -----------------------------------------------------------------------------
hypernet = '/home/zanardi/WORKSPACE/CFD/HyperNet/hypernet/'
databases = hypernet + 'hypernet/databases/'

# Thermo
# -----------------------------------------------------------------------------
mixture = ['O2', 'O']
T = [1500.0e0,2500.0e0,5000.0e0,6000.0e0,8000.0e0,10000.0e0,12000.0e0,14000.0e0,15000.0e0,20000.0e0]

# Species
# -----------------------------------------------------------------------------
grouping_name = 'CBM3'
specie = {
    'thermo_path': databases + 'air/thermo/',
    'grouping': {
        'O2': {
            'name': grouping_name,
            'path': {
                'rovib_levels': databases + 'O2_grouping/levels/O3_UMN/O2/FromUMN_Sorted.inp',
                'lev_to_bin': databases + 'O2_grouping/grouping/O3_UMN/O2/LevelsMap_' + grouping_name + '.csv'
            }
        }
    }
}

# Chemistry
# -----------------------------------------------------------------------------
system = 'O3'
chemistry = {
    'rates_path': databases + 'O2_grouping/kinetics/O3_UMN_' + grouping_name + '/',
    'reac_files': {
        'diss': 'Diss_Corrected.dat',
        'exch': 'Exch_Type1.dat',
        'inel': 'Inel.dat'
    },
    'reac_idx': {
        'diss': 0,
        'exc': 5
    }
}

# Writing
# -----------------------------------------------------------------------------
path = '../'
