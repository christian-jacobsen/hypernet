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
    # Update ------------------------------------------------------------------
    def update(self, T):
        # Update energy modes
        self.e_tr = self.e_tr_(T)
        self.cv_tr = self.cv_tr_()
        self.e_int = self.e_int_()
        self.cv_int = self.cv_int_()
        super(CoupledEnergyModes, self).update(T)

    # Energy ------------------------------------------------------------------
    def e_(self, T):
        # [J/kg]
        return self.e_tr + self.e_int + self.e_f_()

    def cv_(self, T):
        # [J/(kg K)]
        return self.cv_tr + self.cv_int

    def dcvdT_(self, T):
        # [J/(kg K^2)]
        return self.dcv_intdT_()

    # Translational Energy
    def e_tr_(self, T):
        # [J/kg]
        return 3./2. * self.specie.R * T

    def cv_tr_(self):
        # [J/(kg K)]
        return 3./2. * self.specie.R

    # Ro-vibrational Energy
    def e_int_(self):
        # [J/kg]
        if self.specie.n_at > 1:
            # Get partition functions
            q, Q = self.intPF.q, self.intPF.Q
            # Evaluate internal energy [J]
            e_int = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                e_int[bin_i] = np.sum(
                    self.specie.rv_lev['E'][mask] * q[mask]
                )
            e_int = e_int / Q
            # Return the specific internal energy per unit mass [J/kg]
            return e_int / self.specie.m
        else:
            return 0.

    def cv_int_(self):
        # [J/(kg K)]
        if self.specie.n_at > 1:
            # Get internal energy value [J]
            e_int = self.e_int * self.specie.m
            # Get partition functions
            Q = self.intPF.Q
            dqdT, dQdT = self.intPF.dqdT, self.intPF.dQdT
            # Evaluate internal heat capacities [J/K]
            cv_int = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                cv_int[bin_i] = np.sum(
                    self.specie.rv_lev['E'][mask] * dqdT[mask]
                )
            cv_int = (cv_int - e_int * dQdT) / Q
            # Return the specific internal heat capacities
            # per unit mass [J/(kg K)]
            return cv_int / self.specie.m
        else:
            return 0.

    def dcv_intdT_(self):
        # [J/(kg K^2)]
        if self.specie.n_at > 1:
            # Get internal energy value [J]
            e_int = self.e_int * self.specie.m
            # Get internal heat capacities [J/K]
            cv_int = self.cv_int * self.specie.m
            # Get internal partition functions (and derivatives)
            Q, dQdT = self.intPF.Q, self.intPF.dQdT
            d2qdT2, d2QdT2 = self.intPF.d2qdT2, self.intPF.d2QdT2
            # Evaluate derivative of the internal heat 
            # capacities derivatives [J/K^2]
            dcv_intdT = np.zeros(self.specie.n_bins)
            for bin_i in range(self.specie.n_bins):
                mask = self.specie.lev_to_bin == bin_i
                dcv_intdT = np.sum(
                    self.specie.rv_lev['E'][mask] * d2qdT2[mask]
                )
            dcv_intdT = (dcv_intdT - 2*cv_int*dQdT - e_int*d2QdT2) / Q
            return dcv_intdT / self.specie.m
        else:
            return 0.
