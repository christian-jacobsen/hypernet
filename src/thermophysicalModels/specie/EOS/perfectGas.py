class perfectGas(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie
    ):
        self.sp = specie

    # Methods
    ###########################################################################
    def rho(self, p, T):
        return p/(self.sp.R*T)

    def p(self, rho, T):
        return rho*self.sp.R*T

    def psi(self, T):
        return 1./(self.sp.R*T)
