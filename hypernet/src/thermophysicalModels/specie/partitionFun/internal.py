import numpy as np

from hypernet.src.general import const
from hypernet.src.thermophysicalModels.specie.partitionFun import Basic


class Internal(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        super(Internal, self).__init__(specie)

    # Methods
    ###########################################################################
    # Update ------------------------------------------------------------------
    def update(self, T):
        self.z = self.z_(T)
        self.dzdT = self.dzdT_(T)
        self.d2zdT2 = self.d2zdT2_(T)
        self.q = self.q_()
        self.dqdT = self.dqdT_()
        self.d2qdT2 = self.d2qdT2_()
        self.Q = self.Q_()
        self.dQdT = self.dQdT_()
        self.d2QdT2 = self.d2QdT2_()

    # Internal partition functions --------------------------------------------
    def z_(self, T):
        if self.specie.n_at > 1:
            return - self.specie.rv_lev['E'] / (T * const.UKB)
        else:
            return np.ones(1)

    def dzdT_(self, T):
        if self.specie.n_at > 1:
            return - self.z / T
        else:
            return np.zeros(1)

    def d2zdT2_(self, T):
        if self.specie.n_at > 1:
            return - 2 * self.dzdT / T
        else:
            return np.zeros(1)

    def q_(self):
        if self.specie.n_at > 1:
            return self.specie.rv_lev['g'] * np.exp(self.z)
        else:
            return np.zeros(1)

    def dqdT_(self):
        if self.specie.n_at > 1:
            return self.q * self.dzdT 
        else:
            return np.zeros(1)

    def d2qdT2_(self):
        if self.specie.n_at > 1:
            return self.dqdT * self.dzdT + self.q * self.d2zdT2
        else:
            return np.zeros(1)

    def Q_(self):
        if self.specie.n_at > 1:
            Q = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                Q[bin_i] = np.sum(
                    self.q[self.specie.lev_to_bin == bin_i]
                )
            return Q
        else:
            return np.ones(1)*self.specie.g_e

    def dQdT_(self):
        if self.specie.n_at > 1:
            dQdT = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                dQdT[bin_i] = np.sum(
                    self.dqdT[self.specie.lev_to_bin == bin_i]
                )
            return dQdT
        else:
            return np.zeros(1)

    def d2QdT2_(self):
        if self.specie.n_at > 1:
            d2QdT2 = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                d2QdT2[bin_i] = np.sum(
                    self.d2qdT2[self.specie.lev_to_bin == bin_i]
                )
            return d2QdT2
        else:
            return np.zeros(1)

    # Levels distribution -----------------------------------------------------
    def levels_distribution(self, T=0.):
        # [J/mol]
        n, g, E = [], [], []
        if self.specie.n_at > 1:
            g_sp = self.specie.rv_lev['g']
            E_sp = self.specie.rv_lev['E'] / const.EV_to_J
            n_over_Q = self.specie.n / self.Q
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                n.append(n_over_Q[bin_i] * self.q[mask])
                g.append(g_sp[mask])
                E.append(E_sp[mask])
        return n, g, E
