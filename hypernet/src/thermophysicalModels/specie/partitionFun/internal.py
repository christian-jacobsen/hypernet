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
        super(Internal, self).update(T)
        self.q = self.q_(T)
        self.dqdT = self.dqdT_(T)

    # Partition functions -----------------------------------------------------
    def z_(self, T):
        return - self.specie.rv_lev['E'] / (T * const.UKB)

    def dzdT_(self, T):
        return - self.z_(T) / T

    def q_(self, T):
        if self.specie.n_at > 1:
            return self.specie.rv_lev['g'] * np.exp(self.z_(T))
        else:
            return np.zeros(1)

    def dqdT_(self, T):
        if self.specie.n_at > 1:
            return self.q_(T) * self.dzdT_(T)
        else:
            return np.zeros(1)

    def Q_(self, T):
        if self.specie.n_at > 1:
            _q = self.q_(T)
            _Q = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                _Q[bin_i] = np.sum(_q[self.specie.lev_to_bin == bin_i])
            return _Q
        else:
            return np.ones(1)*self.specie.g_e

    def dQdT_(self, T):
        if self.specie.n_at > 1:
            _dqdT = self.dqdT_(T)
            _dQdT = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                _dQdT[bin_i] = np.sum(_dqdT[self.specie.lev_to_bin == bin_i])
            return _dQdT
        else:
            return np.zeros(1)
