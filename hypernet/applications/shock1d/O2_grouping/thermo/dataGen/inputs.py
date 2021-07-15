# Input file
###############################################################################

# Paths
# -----------------------------------------------------------------------------
hypernet = '/home/zanardi/WORKSPACE/CFD/HyperNet/hypernet/'
databases = hypernet + 'hypernet/databases/'

# Thermo
# -----------------------------------------------------------------------------
species = ['O2', 'O']
thermo = 'CoupledEnergyModes'
T = [300., 30000.]

# Species
# -----------------------------------------------------------------------------
grouping_name = {
    'O2': 'CBM3'
}
specie = {
    'thermo_path': databases + 'air/thermo/',
    'grouping': {
        'O2': {
            'name': grouping_name['O2'],
            'path': {
                'rovib_levels': databases + 'O2_grouping/levels/O3_UMN/O2/FromUMN_Sorted.inp',
                'lev_to_bin': databases + 'O2_grouping/grouping/O3_UMN/O2/LevelsMap_' + grouping_name['O2'] + '.csv'
            }
        }
    }
}

# Writing
# -----------------------------------------------------------------------------
path = '../'
tab_name = {
    'temp': ['Tg', 'Tt'],
    'groups': ['eint', 'qint', 'cvint']
}
file_type = 'bin'