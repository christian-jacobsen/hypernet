import numpy as np
import hypernet.src.general.const as const

class energyModes(object):

    # Initialization
    ###########################################################################
    def __init__(self, database, mixture):

        self.database = database    # Species database
        self.mixture = mixture      # Initial mixture (dict = {name: Y})

        self.R = 0.
        # Updata initial specifc gas constant of the mixture
        self.R(mixture)

        self.levels = {}
        for specie in mixture:
            self.read_internal_levels(specie)


    # Methods
    ###########################################################################
    # Enthalpy ----------------------------------------------------------------
    def cp(self,R,T):
        # [J/(kg K)]
        return self.cv(R,T) + R

    def h(self,R,T):
        # [J/kg]
        return self.e(R,T) + R*T

    # Internal Energy ---------------------------------------------------------
    def cv(self,R,T):
        # [J/(kg K)]
        return self.cv_tr(R) + self.cv_rv(R,T)

    def e(self,R,T):
        # [J/kg]
        return self.e_tr(R) + self.e_rv(R,T) + self.e_f(R,T)

    # Translational Internal Energy -------------------------------------------
    def cv_tr(self,R):
        # [J/(kg K)]
        return 3./2.*R

    def e_tr(self,R,T):
        # [J/kg]
        return 3./2.*R*T

    # Ro-Vibrational Internal Energy ------------------------------------------
    def cv_rv(self,R):
        # [J/(kg K)]
        return 3./2.*R

    def e_rv(self,R,T):
        # [J/kmol]
        return 3./2.*R*T


    # Partition functions
    def levels_partition_fn(self, T):
        return self.g * np.exp( - self.E / (T * const.UKB) )

    def groups_partition_fn(self, T):
        '''Compute bins partition function.'''
        q = self.levels_partition_fn(T)
        Q_bins = np.zeros(self.num_bins)
        for bin_i in range(self.num_bins):
            Q_bins[bin_i] = np.sum(q[self.lev_to_bin == bin_i])
        return Q_bins

    # Reading
    def read_internal_levels(self):
        '''Retrieve all the ro-vib levels.'''
        data = pd.read_csv(self.levels_file, header=None, skiprows=15, delim_whitespace=True)
        data = data.apply(pd.to_numeric, errors='coerce')
        vqn  = np.array(data[0])                    # Vibrational Q.N.
        jqn  = np.array(data[1])                    # Rotational  Q.N.
        g    = 1.5*(2.0*jqn+1.0)                    # Degeneracies (electronic included)
        E    = np.array(data[2]) * const.EH_to_J    # Energy in Joule
        E    = E - E[0]
        num  = len(vqn)
        return g, E, num