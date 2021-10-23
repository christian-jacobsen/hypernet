import numpy as np

from hypernet.src.thermophysicalModels.specie.equationOfState import Basic


class PerfectGas(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        super(PerfectGas, self).__init__(specie)

    # Methods
    ###########################################################################
    def p(self, T):
        return self.specie.rho * self.specie.R * T

    def psi(self, T):
        return 1. / (self.specie.R * T)
