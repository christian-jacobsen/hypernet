import os
import numpy as np
import pandas as pd

from hypernet.src.general import const
from hypernet.src.general import utils

import hypernet.database as db
thermo_db = os.path.dirname(db.__file__) + '/air/thermo/'
grouping_db = os.path.dirname(db.__file__) + '/air/grouping/'


class Specie(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        name,
        rovib=None,
        *args,
        **kwargs
    ):
        # Read specie properties ----------------------------------------------
        self.properties = {
            'NAME': [str, 'name'],              # Specie name
            'MOLAR_MASS': [float, 'M'],         # Molecular mass [kg/mol]
            'NB_ATOMS': [int, 'n_at'],          # Number of atoms
            'ELEM_NAME': [str, 'elem_name'],    # Chemical elements (symbols)
            'ELEM_QUANT': [int, 'elem_num'],    # Chemical elements (quantities)
            'EF': [float, 'Ef'],                # Formation energy [J/mol]
            'DE': [float, 'D'],                 # Dissociation energy [eV]
            'CHARGE': [int, 'C'],               # Electric charge
            'NUCLEAR_CHARGE': [int, 'C_nucl'],  # Nuclear charge (or atomic number)
            'IP': [float, 'Eion'],              # Ionization potential [eV]
            'LIN': [int, 'lin'],                # Linearity
            'SYM': [float, 'sym'],              # Symmetry factor
            'LJ_DIAM': [float, 'lj_diam'],      # Lennard-Jones parameter: diameter [A]
            'LJ_EPS': [float, 'lj_eps'],        # Lennard-Jones parameter: potential well-depth [K]
            'POLAR': [float, 'pol'],            # Polarizability [A^3]s
            'MODEL': [str, 'model'],            # Thermal model for internal energy (e.g. RR, HO, EL, RR_HO,...):
                                                #   RR = Rigid-rotor
                                                #   HO = Harmomic oscillator
                                                #   EL = Electronic energy computed based on a Boltzmann distribution of electronic levels 
                                                #   NONE = pure State-to-State approach for internal energy levels
            'THETA_ROT': [float, 'theta_r'],    # Characteristic rotational temperature [K]
            'THETA_VIB': [float, 'theta_v']     # Characteristic vibrational temperature(s) [K]
        }
        self.read_properties(thermo_db + name)

        # Read specie electronic levels ---------------------------------------
        self.el_lev, self.g_e = self.read_elec_levels(thermo_db + name)

        # Read specie ro-vibrational levels -----------------------------------
        self.rovib = rovib
        if self.n_at > 1 and self.rovib is not None:
            path = grouping_db + self.rovib['system'] + '_' \
                + self.rovib['PES'] + '/' + name + '/'
            # Retrieve ro-vib levels
            self.rv_lev = self.read_rovib_levels(path + '/levels.inp')
            # Retrieve ro-vib level to bin mapping
            map_file = path + '/LevelsMap_' + rovib['grouping'] + '.csv' \
                if self.rovib['grouping'] else None
            self.lev_to_bin, self.n_bins = self.read_grouping(
                map_file, self.rovib['grouping']
            )

    # Properties
    ###########################################################################
    @property
    def Y(self):
        return self._Y
    @Y.setter
    def Y(self, value):
        self._Y = value

    @property
    def X(self):
        return self._X
    @X.setter
    def X(self, value):
        self._X = value

    @property
    def n(self):
        return self._n
    @n.setter
    def n(self, value):
        self._n = value

    @property
    def rho(self):
        return self._rho
    @rho.setter
    def rho(self, value):
        self._rho = value

    # Methods
    ###########################################################################
    # Global properties -------------------------------------------------------
    def read_properties(self, file):
        for name, (data_type, var) in self.properties.items():
            val = self.read_property(file, name, data_type)
            setattr(self, var, val)

        # Specific gas constant R [J/(kg K)]
        self.R = const.URG / self.M

        # Mass [kg]
        self.m = self.M / const.UNA

        # Conversions
        if self.D:
            # Dissociation energy: [eV] -> [J]
            self.D = self.D * const.EV_to_J
        if self.Eion:
            # Ionization potential: [eV] -> [J]
            self.Eion = self.Eion * const.EV_to_J

    def read_property(self, file, name, data_type):
        result = utils.search_string(file, name)
        if result:
            num, line = result
            k, v = line.strip().split(" = ")
            if " " in v:
                v = v.split(" ")
            if isinstance(v, (list,tuple)):
                return [
                    data_type(utils.check_format(data_type, i)) for i in v
                ]
            else:
                return data_type(utils.check_format(data_type, v))
        else:
            return None

    # Electronic levels -------------------------------------------------------
    def read_elec_levels(self, file):
        '''Retrieve all the electronic levels.'''
        # Initilize input_file
        db = {}
        # Number of electronic levels
        db['num'] = self.read_property(file, 'NB_ELEC_LEVELS', int)
        # Electronic levels degeneracies 'g' and energies 'E' [J]
        num, line = utils.search_string(file, 'NB_ELEC_LEVELS')
        deg, en = [], []
        with open(file, 'r') as f:
            for l in f.readlines()[num+1:]:
                g, E = l.strip().split("  ")
                deg.append(float(g))
                en.append(float(E) * const.EH_to_J)
        db['g'] = np.array(deg)
        db['E'] = np.array(en) - en[0]
        g_e = deg[0]
        return db, g_e

    # Molecular ro-vibrational levels -----------------------------------------
    def read_rovib_levels(self, file):
        '''Retrieve all the ro-vib levels.'''
        # Initilize input_file
        db = {}
        # Initilize input_file
        data = pd.read_csv(
            file, header=None, skiprows=15, delim_whitespace=True
        )
        db['vqn'] = np.array(data[0])                     # Vibrational Q.N.
        db['jqn'] = np.array(data[1])                     # Rotational  Q.N.
        db['num'] = len(db['vqn'])
        db['g'] = self.g_e/2.0*(2.0*db['jqn']+1.0)        # Degeneracies
        E = np.array(data[2]) * const.EH_to_J       # Energy [J]
        db['E'] = E - E[0]
        return db

    # Grouping ----------------------------------------------------------------
    def read_grouping(self, file, name):
        '''Retrieve mapping of each rot-vib level to the respective bin.'''
        if name is not None:
            data = pd.read_csv(file, header=None, skiprows=1)
            lev_to_bin = np.array(data[1]) - 1
            n_bins = max(lev_to_bin) + 1
        else:
            lev_to_bin = np.arange(self.rv_lev['num'], dtype=np.int32)
            n_bins = self.rv_lev['num']
        return lev_to_bin, n_bins
