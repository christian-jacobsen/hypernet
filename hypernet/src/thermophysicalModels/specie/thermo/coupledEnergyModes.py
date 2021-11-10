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
    # Enthalpy ================================================================
    def dcpdT(self, T):
        # [J/(mol K)]
        return self.dcvdT(T)

    def cp(self, T):
        # [J/(mol K)]
        return self.cv(T) + const.URG

    def h(self, T):
        # [J/mol]
        return self.e(T) + const.URG*T

    # Energy ------------------------------------------------------------------
    def dcvdT(self, T=0.):
        # [J/(mol K)]
        return self.dcv_intdT(T)

    def cv(self, T=0.):
        # [J/(mol K)]
        return self.cv_tr() + self.cv_int(T)

    def e(self, T=0.):
        # [J/mol]
        return self.e_tr(T) + self.e_int(T) + self.e_f()

    # Translational Energy
    def cv_tr(self):
        # [J/(mol K)]
        return 3./2.*const.URG

    def e_tr(self, T=0.):
        # [J/mol]
        return 3./2.*const.URG*T

    # Ro-vibrational Energy
    def dcv_intdT(self, T=0.):
        # [J/(mol K^2)]
        if self.specie.n_at > 1:
            e_int_ = self.e_int(T) / const.UNA
            cv_int_ = self.cv_int(T) / const.UNA
            Q, dQdT = self.intPF.Q, self.intPF.dQdT
            d2qdT2, d2QdT2 = self.intPF.d2qdT2, self.intPF.d2QdT2
            dcv_intdT_ = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                f1 = np.sum(self.specie.rv_lev['E'][mask] * d2qdT2[mask])
                f2 = - 2 * cv_int_[bin_i] * dQdT[bin_i]
                f3 = - e_int_[bin_i] * d2QdT2[bin_i]
                dcv_intdT_[bin_i] = (f1 + f2 + f3) / Q[bin_i]
            return dcv_intdT_ * const.UNA
        else:
            return 0.

    def cv_int(self, T=0.):
        # [J/(mol K)]
        if self.specie.n_at > 1:
            e_int_ = self.e_int(T) / const.UNA
            q, Q = self.intPF.q, self.intPF.Q
            dqdT, dQdT = self.intPF.dqdT, self.intPF.dQdT
            cv_int_ = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                f1 = np.sum(self.specie.rv_lev['E'][mask] * dqdT[mask])
                f2 = - e_int_[bin_i] * dQdT[bin_i]
                cv_int_[bin_i] = (f1 + f2) / Q[bin_i]
            return cv_int_ * const.UNA
        else:
            return 0.

    def e_int(self, T=0.):
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
    def n_g_E_lev(self, T=0.):
        # [J/mol]
        if self.specie.n_at > 1:
            q, Q = self.intPF.q, self.intPF.Q
            n_lev, g_lev, E_lev = [], [], []
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                n_lev.append(self.specie.n[bin_i] * q[mask] / Q[bin_i])
                g_lev.append(self.specie.rv_lev['g'][mask])
                E_lev.append(self.specie.rv_lev['E'][mask]*1./const.EV_to_J)
            return n_lev, g_lev, E_lev
