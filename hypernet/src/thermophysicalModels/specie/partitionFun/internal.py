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

    # Properties
    ###########################################################################
    @property
    def q(self):
        return self._q
    @q.setter
    def q(self, value):
        self._q = value

    @property
    def dqdT(self):
        return self._dqdT
    @dqdT.setter
    def dqdT(self, value):
        self._dqdT = value

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

    # Partition functions -----------------------------------------------------
    def z_(self, T):
        if self.specie.n_at > 1:
            return - self.specie.rv_lev['E'] / (T * const.UKB)
        else:
            return np.zeros(1)

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
            _Q = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                _Q[bin_i] = np.sum(
                    self.q[self.specie.lev_to_bin == bin_i]
                )
            return _Q
        else:
            return np.ones(1)*self.specie.g_e

    def dQdT_(self):
        if self.specie.n_at > 1:
            _dQdT = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                _dQdT[bin_i] = np.sum(
                    self.dqdT[self.specie.lev_to_bin == bin_i]
                )
            return _dQdT
        else:
            return np.zeros(1)

    def d2QdT2_(self):
        if self.specie.n_at > 1:
            _d2QdT2 = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                _d2QdT2[bin_i] = np.sum(
                    self.d2qdT2[self.specie.lev_to_bin == bin_i]
                )
            return _d2QdT2
        else:
            return np.zeros(1)
