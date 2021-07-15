import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils

from scipy import interpolate


class SplineSurrogate(Surrogate):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        data_x,
        data_Y,
        *args,
        **kwargs
    ):
        # Specie Properties ===================================================
        self.x = data_x
        self.Y = data_Y
        self.n_species = data_Y.shape[1]

    # Methods
    ###########################################################################
    def fit(self):
        self.params = []
        for i in range(self.n_species-1):
            self.params.append(interpolate.splrep(self.x, self.Y[:,i]))

    def update(self, x):
        # [J/mol]
        Y = []
        for i in range(self.n_species-1):
            Y.append(interpolate.splev(x, self.params[i]))

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
