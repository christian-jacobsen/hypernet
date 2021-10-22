import numpy as np

from hypernet.src.general import const
from hypernet.src.specie.partitionFun.basicPF import PartitionFun


class Internal(PartitionFun):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        super(Internal, self).__init__(specie)

    # Internal partition functions --------------------------------------------
    def z(self, T):
        return - self.specie.rv_lev['E'] / (T * const.UKB)

    def dz_dT(self, T):
        return - self.z(T) / T

    def q(self, T):
        return self.specie.rv_lev['g'] * np.exp(self.z(T))

    def dq_dT(self, T):
        return self.q(T) * self.dz_dT(T)

    def Q_(self, T):
        if self.specie.n_at > 1:
            _q = self.q(T)
            _Q = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                _Q[bin_i] = np.sum(_q[self.specie.lev_to_bin == bin_i])
            return _Q
        else:
            return np.ones(1)*self.specie.g_e

    def dQ_dT_(self, T):
        if self.specie.n_at > 1:
            _dq_dT = self.dq_dT(T)
            _dQ_dT = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                _dQ_dT[bin_i] = np.sum(_dq_dT[self.specie.lev_to_bin == bin_i])
            return _dQ_dT
        else:
            return np.zeros(1)
