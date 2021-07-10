class perfectGas(object):

    # Initialization
    ###########################################################################
    def __init__(self, database, mixture):

        self.database = database    # Species database
        self.mixture = mixture      # Initial mixture (dict = {name: Y})

        self.R = 0.
        # Updata initial specifc gas constant of the mixture
        self.R(mixture)


    # Methods
    ###########################################################################
    def R(self, Y):
        for specie in mixture:
            self.R += self.database[specie]['R'] * Y[specie]

    def rho(self, p, T):
        return p/(self.R*T)

    def p(self, rho, T):
        return rho*self.R*T

    def psi(self, T):
        return 1./(self.R*T)
