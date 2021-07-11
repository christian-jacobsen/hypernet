import hypernet.src.general.const as const

class specie(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        database,
        name,
        Y
    ):
        # Read specie properties ==============================================
        self.properties = {
            'NAME': [str, 'name'],              # Specie name
            'MOLAR_MASS': [float, 'm'],         # Molecular mass [kg/mol]
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
            'POLAR': [float, 'pol'],            # Polarizability [A^3]
            'MODEL': [str, 'model'],            # Thermal model for internal energy (e.g. RR, HO, EL, RR_HO,...):
                                                #   RR = Rigid-rotor
                                                #   HO = Harmomic oscillator
                                                #   EL = Electronic energy computed based on a Boltzmann distribution of electronic levels 
                                                #   NONE = pure State-to-State approach for internal energy levels
            'THETA_ROT': [float, 'theta_r'],    # Characteristic rotational temperature [K]
            'THETA_VIB': [float, 'theta_v']     # Characteristic vibrational temperature(s) [K]
        }
        self.read_properties(database.path_to_properties + name)

        # Read specie electronic levels =======================================
        self.el_lev, g_e = self.read_elec_levels(
            database.path_to_properties + name
        )

        # Read specie ro-vibrational levels ===================================
        if self.n_at > 1:
            self.rv_lev = self.read_rovib_levels(
                database.path_to_inter_levels, g_e
            )
            self.lev_to_bin, self.n_bins = self.read_grouping(
                database.grouping_file, database.grouping
            )

        # Initilize specie mass fraction ======================================
        self.Y_(Y)


    # Methods
    ###########################################################################
    # Mass Fraction ===========================================================
    def Y_(self, Y):
        self.Y = Y

    # Reading methods
    ###########################################################################
    # Global properties =======================================================
    def read_properties(self, file):
        for name, (data_type, var) in self.properties.items():
            val = self.read_property(file, name, data_type)
            setattr(self, var, val)

        # Specific gas constant R [J/(kg K)]
        self.R = const.URG / self.m

        # Conversions ---------------------------------------------------------
        if self.Ef:
            # Formation energy: [J/mol] -> [J/kg]
            self.Ef = self.Ef / self.m
        if self.D:
            # Dissociation energy: [eV] -> [J]
            self.D = self.D * const.EV_to_J
        if self.Eion:
            # Ionization potential: [eV] -> [J]
            self.Eion = self.Eion * const.EV_to_J

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

    # Electronic levels =======================================================
    def read_elec_levels(self, file):
        '''Retrieve all the electronic levels.'''
        # Initilize database
        db = {}
        # Number of electronic levels
        db['num'] = self.read_property(file, 'NB_ELEC_LEVELS', int)
        # Electronic levels degeneracies 'g' and energies 'E' [J]
        num, line = self.search_string(file, 'NB_ELEC_LEVELS')
        deg, en = [], []
        with open(file, 'r') as f:
            for l in f.readlines()[num+1:]:
                g, E = l.strip().split("  ")
                deg.append(float(g))
                en.append(float(E) * const.EH_to_J)
        db['g'] = deg
        db['E'] = en

        return db, deg[0]

    # Molecular ro-vibrational levels =========================================
    def read_rovib_levels(self, file, g_e):
        '''Retrieve all the ro-vib levels.
        g_e: electronic ground state degeneracy

        '''
        # Initilize database
        db = {}
        # Initilize database
        data = pd.read_csv(file, header=None, skiprows=15, delim_whitespace=True).values
        vqn = np.array(data[0])                     # Vibrational Q.N.
        jqn = np.array(data[1])                     # Rotational  Q.N.
        db['num'] = len(vqn)
        db['g'] = g_e/2.0*(2.0*jqn+1.0)       # Degeneracies
        E = np.array(data[2]) * const.EH_to_J       # Energy [J]
        db['E'] = E - E[0]
        return db

    # Grouping ================================================================
    def read_grouping(self, file, name):
        '''Retrieve mapping of each rot-vib level to the respective bin.'''
        if name is not None:
            data = pd.read_csv(file, header=None, skiprows=1).values
            lev_to_bin = np.array(data[1]) - 1
            n_bins = max(bins) + 1
        else:
            lev_to_bin = np.zeros(self.rv_lev['num'], dtype=np.int32)
            n_bins = 1
        return lev_to_bin, n_bins


    # Util methods
    ###########################################################################
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
