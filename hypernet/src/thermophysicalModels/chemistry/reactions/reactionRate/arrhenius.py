import numpy as np

from hypernet.src.thermophysicalModels.chemistry.reactions.reactionRate import Basic


class Arrhenius(Basic):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        reactionsDatabase,
        *args,
        **kwargs
    ):
        super(ArrheniusReactionRate, self).__init__(
            reactionsDatabase,
            *args,
            **kwargs
        )
        self.A = self.reacDB['A']
        self.beta = self.reacDB['beta']
        self.Ta = self.reacDB['Ta']

    # Methods
    ###########################################################################
    # Forward reaction rates --------------------------------------------------
    def k_(self, T=0.):
        _k = self.A
        if T > 0.:
            _k *= np.power(T, self.beta)
            _k *= np.exp(-self.Ta / T)
        return _k

    def dkdT_(self, T=0.):
        _dkdT = self.k(T)
        if T > 0.:
            _dkdT *= self.beta / T
            _dkdT *= self.Ta / T**2
        return _dkdT
