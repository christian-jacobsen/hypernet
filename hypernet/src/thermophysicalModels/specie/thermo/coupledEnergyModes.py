import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils
from hypernet.src.thermophysicalModels.specie.thermo import Basic


class CoupledEnergyModes(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        internalPF,
        *args,
        **kwargs
    ):
        super(CoupledEnergyModes, self).__init__(specie)
        self.intPF = internalPF

    # Methods
    ###########################################################################
    # Energy ------------------------------------------------------------------
    def cv(self, T):
        # [J/(mol K)]
        return self.cv_tr() + self.cv_int(T)

    def e(self, T):
        # [J/mol]
        return self.e_tr(T) + self.e_int(T) + self.e_f()

    # Translational Energy
    def cv_tr(self):
        # [J/(mol K)]
        return 3./2.*const.URG

    def e_tr(self, T):
        # [J/mol]
        return 3./2.*const.URG*T

    # Ro-vibrational Energy
    def cv_int(self, T):
        # [J/(mol K)]
        if self.specie.n_at > 1:
            q, Q = self.intPF.q, self.intPF.Q
            dqdT, dQdT = self.intPF.dqdT, self.intPF.dQdT
            cv_int_ = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                f1 = self.specie.rv_lev['E'][mask] / Q[bin_i]
                f2 = dqdT[mask] - q[mask] * dQdT[bin_i] / Q[bin_i]
                cv_int_[bin_i] = np.sum(f1 * f2)
            return cv_int_ * const.UNA
        else:
            return 0.

    def e_int(self, T):
        # [J/mol]
        if self.specie.n_at > 1:
            q, Q = self.intPF.q, self.intPF.Q
            e_int_ = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                e_int_[bin_i] = np.sum(
                    self.specie.rv_lev['E'][mask] * q[mask] / Q[bin_i]
                )
            return e_int_ * const.UNA
        else:
            return 0.

    # Levels populations ------------------------------------------------------
    def n_g_E_lev(self, T):
        # [J/mol]
        if self.specie.n_at > 1:
            q, Q = self.intPF.q_(T), self.intPF.Q_(T)
            n_lev, g_lev, E_lev = [], [], []
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                n_lev.append(self.specie.n[bin_i] * q[mask] / Q[bin_i])
                g_lev.append(self.specie.rv_lev['g'][mask])
                E_lev.append(self.specie.rv_lev['E'][mask]*1./const.EH_to_J)
            return n_lev, g_lev, E_lev
