import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils


class CoupledEnergyModes(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        # Specie Properties ===================================================
        self.specie = specie

    # Methods
    ###########################################################################
    # Enthalpy ================================================================
    def cp(self, T):
        # [J/(mol K)]
        return self.cv(T) + const.URG

    def h(self, T):
        # [J/mol]
        return self.e(T) + const.URG*T

    # Internal Energy =========================================================
    def cv(self, T):
        # [J/(mol K)]
        return self.cv_tr() + self.cv_rv(T)

    def e(self, T):
        # [J/mol]
        return self.e_tr(T) + self.e_rv(T) + self.e_f()

    # Energy of formation -----------------------------------------------------
    def e_f(self):
        # [J/mol]
        return self.specie.Ef

    # Translational Internal Energy -------------------------------------------
    def cv_tr(self):
        # [J/(mol K)]
        return 3./2.*const.URG

    def e_tr(self, T):
        # [J/mol]
        return 3./2.*const.URG*T

    # Ro-Vibrational Internal Energy ------------------------------------------
    def cv_rv(self, T):
        # [J/(mol K)]
        if self.specie.n_at > 1:
            q_, Q_ = self.q(T), self.Q(T)
            dq_dT_, dQ_dT_ = self.dq_dT(T), self.dQ_dT(T)
            cv_rv_ = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                f1 = self.specie.rv_lev['E'][mask] / Q_[bin_i]
                f2 = dq_dT_[mask] - q_[mask] * dQ_dT_[bin_i] / Q_[bin_i]
                cv_rv_[bin_i] = np.sum(f1 * f2)
            return cv_rv_ * const.UNA
        else:
            return 0.

    def e_rv(self, T):
        # [J/mol]
        if self.specie.n_at > 1:
            q_, Q_ = self.q(T), self.Q(T)
            e_rv_ = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                e_rv_[bin_i] = np.sum(self.specie.rv_lev['E'][mask] * q_[mask] / Q_[bin_i])
            return e_rv_ * const.UNA
        else:
            return 0.

    # Partition functions
    def z(self, T):
        return - self.specie.rv_lev['E'] / (T * const.UKB)

    def dz_dT(self, T):
        return - self.z(T) / T

    def q(self, T):
        return self.specie.rv_lev['g'] * np.exp(self.z(T))

    def dq_dT(self, T):
        return self.q(T) * self.dz_dT(T)

    def Q(self, T):
        '''Compute bins partition function.'''
        q_ = self.q(T)
        Q_ = np.zeros(self.specie.n_bins)
        for bin_i in range(self.specie.n_bins):
            Q_[bin_i] = np.sum(q_[self.specie.lev_to_bin == bin_i])
        return Q_

    def dQ_dT(self, T):
        '''Compute bins partition function.'''
        dq_dT_ = self.dq_dT(T)
        dQ_dT_ = np.zeros(self.specie.n_bins)
        for bin_i in range(self.specie.n_bins):
            dQ_dT_[bin_i] = np.sum(dq_dT_[self.specie.lev_to_bin == bin_i])
        return dQ_dT_
