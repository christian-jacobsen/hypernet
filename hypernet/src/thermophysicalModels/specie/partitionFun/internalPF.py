# TODO: add electronic partition function calculations

import numpy as np

from hypernet.src.general import const
from hypernet.src.specie.partitionFun.basicPF import PartitionFun


class Internal(PartitionFun):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        electronic=False, # TODO
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

    def Q(self, T):
        if self.specie.n_at > 1:
            q_ = self.q(T)
            Q_ = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                Q_[bin_i] = np.sum(q_[self.specie.lev_to_bin == bin_i])
            return Q_
        else:
            return self.specie.g_e

    def dQ_dT(self, T):
        if self.specie.n_at > 1:
            dq_dT_ = self.dq_dT(T)
            dQ_dT_ = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                dQ_dT_[bin_i] = np.sum(dq_dT_[self.specie.lev_to_bin == bin_i])
            return dQ_dT_
        else:
            return 0.
