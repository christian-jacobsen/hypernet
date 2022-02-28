# *************************************************************************** #
#                                   HyperNet                                  #
# --------------------------------------------------------------------------- #
#                 Machine-Learning-Based library for modeling                 #
#           multi-component non-equilibrium thermochemical processes          #
#                                                                             #
# *************************************************************************** #

# -------------------------------- INPUT FILE ------------------------------- #
# Description:
# >> Defining global input variables for the `box` app
# --------------------------------------------------------------------------- #


# Thermo
# =============================================================================
ambient = {
    'p': 5000.e0,                       # Pressure [Pa]
    'T': 10000.e0,                      # Translational Temperature [K]
}

mixture = {
    'var': 'X',                         # 'X' or 'Y' (molar/mass fraction)
    'name': 'MultiComponent',           # DO NOT CHANGE!
    'composition': {                    # Mass/Molar fractions of each species
        'O2':   [0.95, 1.e-13, 1.e-13], # Num. values = Num. groups
        'O':    0.05
    }
}

# *** Species
specie = {
    'O2': {                             # If molecule, define rovibrational
                                        #   parameters (used to create the
                                        #   path to the proper grouping
                                        #   folder in 'database')
        'rovib': {
            'system': 'O3',             # O3 system <=> O2+O reactive mixture
            'PES': 'UMN',               # PES from Minnesota
            'grouping': 'CBM3'          # Grouping strategy:
                                        #  'CBM3' = name (CBM) + groups (3)
        }
    },
    'O': None                           # If atom, put `None`
}

thermo = 'CoupledEnergyModes'           # DO NOT CHANGE!
EOS = 'PerfectGas'                      # DO NOT CHANGE!
constVP = 'V'                           # DO NOT CHANGE!

# Chemistry
# =============================================================================
chemistry = {
    'processFlags': {                   # Reactive processes:
        'excit':  1,                    #   Thermal excitation/relaxation proc.
        'diss':   1                     #   Dissociation/Recombination proc.
    },
    'model': 'Standard',                # DO NOT CHANGE!
    'solver': 'Standard',               # DO NOT CHANGE!
    'reactionsList': None,              # <path-to-file> or `None`
                                        # If `None`, it will recover the rates
                                        #   from the grouping information given
                                        #   in `specie` above.
    'heatBath': 'adiabatic'             # 'adiabatic' or 'isothermal'
}

# ODE parameters
# =============================================================================
setup = {
    'start': 0.,
    'end': 1.e-1,
    'delta': {
        'min': 1.e-14,
        'max': 1.e-4,
        'str': 1.015
    }
}
