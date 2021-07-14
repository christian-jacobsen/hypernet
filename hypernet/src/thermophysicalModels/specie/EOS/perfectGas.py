class PerfectGas(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        R,
        *args,
        **kwargs
    ):
        self.R = R

    # Methods
    ###########################################################################
    def rho(self, p, T):
        return p/(self.R*T)

    def p(self, rho, T):
        return rho*self.R*T

    def psi(self, T):
        return 1./(self.R*T)
