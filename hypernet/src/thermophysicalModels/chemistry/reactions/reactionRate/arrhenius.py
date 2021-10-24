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
        super(Arrhenius, self).__init__(
            reactionsDatabase,
            *args,
            **kwargs
        )
        self.A = self.reacDB['A'].to_numpy()
        self.beta = self.reacDB['beta'].to_numpy()
        self.Ta = self.reacDB['Ta'].to_numpy()

    # Methods
    ###########################################################################
    # Forward reaction rates --------------------------------------------------
    def k_(self, T):
        _k = self.A
        _k *= np.power(T, self.beta)
        _k *= np.exp(-self.Ta / T)
        return _k

    def dkdT_(self, T):
        _dkdT = self.beta / T
        _dkdT += self.Ta / T**2
        _dkdT *= self.k
        return _dkdT
