import numpy as np

from hypernet.src.general import const
from hypernet.src.general import utils


class ArrheniusReactionRate(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        A,
        beta,
        Ta
    ):
        self.A = A
        self.beta = beta
        self.Ta = Ta

    # Master Equation matrices
    ###########################################################################
    def k(self, T=0.):
        k_ = self.A
        if T > 0.:
            k_ *= np.power(T, self.beta)
            k_ *= np.exp(-self.Ta / T)
        return k_

    def dkdT(self, T=0.):
        dkdT_ = self.k(T)
        if T > 0.:
            dkdT_ *= self.beta / T
            dkdT_ *= self.Ta / T**2
        return dkdT_
