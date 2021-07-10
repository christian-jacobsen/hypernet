import hypernet.src.general.const as const

class specie(object):

    # Initialization
    ###########################################################################
    def __init__(self, path_to_database, mixture):

        self.path_to_database = path_to_database

        self.properties = {
            'NAME': str,            # Specie name
            'MOLAR_MASS': float,    # Molecular mass [kg/mol]
            'NB_ATOMS': int,        # Number of atoms
            'ELEM_NAME': str,       # Chemical elements (symbols)
            'ELEM_QUANT': int,      # Chemical elements (quantities)
            'EF': float,            # Formation energy [J/mol]
            'DE': float,            # Dissociation energy [eV]
            'CHARGE': int,          # Electric charge
            'NUCLEAR_CHARGE': int,  # Nuclear charge (or atomic number)
            'IP': float,            # Ionization potential [eV]
            'LIN': int,             # Linearity
            'SYM': float,           # Symmetry factor
            'LJ_DIAM': float,       # Lennard-Jones parameter: diameter [A]
            'LJ_EPS': float,        # Lennard-Jones parameter: potential well-depth [K]
            'POLAR': float,         # Polarizability [A^3]
            'MODEL': str,           # Thermal model for internal energy (e.g. RR, HO, EL, RR_HO,...):
                                    #   RR = Rigid-rotor
                                    #   HO = Harmomic oscillator
                                    #   EL = Electronic energy computed based on a Boltzmann distribution of electronic levels 
                                    #   NONE = pure State-to-State approach for internal energy levels
            'THETA_ROT': float,     # Characteristic rotational temperature [K]
            'THETA_VIB': float,     # Characteristic vibrational temperature(s) [K]
            'NB_ELEC_LEVELS': int,  # Number of electronic levels
        }

        self.database = {}
        for specie in mixture:
            self.read_properties(specie)


    # Methods
    ###########################################################################
    # Reading -----------------------------------------------------------------
    def read_properties(self, specie):
        file = self.path_to_database + specie

        if specie not in self.database:

            # Fill the database
            self.database[specie] = {
                prop: self.read_property(file, prop, data_type) \
                    for prop, data_type in self.properties.items()
            }

            # Define specie molar mass [kg/mol]
            m = self.database[specie]['MOLAR_MASS']

            # Additions -------------------------------------------------------
            # Electronic degeneracies 'g' and energies 'E' [J] levels
            if 'NB_ELEC_LEVELS' in self.database[specie]:
                self.database[specie]['ELEC_LEVELS'] = \
                    self.read_elec_levels(file)

            # Specific gas constant R [J/(kg K)]
            self.database[specie]['R'] = const.URG / m

            # Conversions -------------------------------------------------------
            if self.database[specie]['EF']:
                # Formation energy: [J/mol] -> [J/kg]
                self.database[specie]['EF'] = self.database[specie]['EF'] * m
            if self.database[specie]['DE']:
                # Dissociation energy: [eV] -> [J]
                self.database[specie]['DE'] = self.database[specie]['DE'] * const.EV_to_J
            if self.database[specie]['IP']:
                # Ionization potential: [eV] -> [J]
                self.database[specie]['IP'] = self.database[specie]['IP'] * const.EV_to_J

    def read_property(self, file, name, data_type):
        result = self.search_string(file, name)
        if result:
            num, line = result
            k, v = line.strip().split(" = ")
            if " " in v:
                v = v.split(" ")
            if isinstance(v, (list,tuple)):
                return [ data_type(self.check_format(data_type, i)) for i in v ]
            else:
                return data_type(self.check_format(data_type, v))
        else:
            return None

    def read_elec_levels(self, file):
        num, line = self.search_string(file, 'NB_ELEC_LEVELS')
        degeneracies, energies = [], []
        with open(file, 'r') as f:
            for l in f.readlines()[num+1:]:
                g, e = l.strip().split("  ")
                degeneracies.append(float(g))
                energies.append(float(e) * const.EV_to_J)
            return {
                'g': degeneracies, 'E': energies
            }

    # Utils -------------------------------------------------------------------
    def search_string(self, file, string):
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                if string+' ' in line:
                    return i, line

    def check_format(self, data_type, string):
        if data_type == float and 'd' in string:
            return string.replace("d", "e")
        else:
            return string
