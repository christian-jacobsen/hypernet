class PerfectGas(object):

    # Initialization
    ###########################################################################
    def __init__(
        self
    )

    # Methods
    ###########################################################################
    def rho(self, p, T):
        return p/(R*T)

    def p(self, rho, T):
        return rho*R*T

    def psi(self, T):
        return 1./(R*T)
