import numpy as np

from hypernet.src.specie.equationOfState.basicEOS import EOS


class PerfectGas(EOS):

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
        return self.rho*self.specie.R*T

    def psi(self, T):
        return 1./(self.specie.R*T)
