import numpy as np

from hypernet.src.thermophysicalModels.chemistryModel.reactions.reactionRate import BasicReactionRate


class ArrheniusReactionRate(BasicReactionRate):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        reactionsDataFrame,
        *args,
        **kwargs
    ):
        super(ArrheniusReactionRate, self).__init__(
            reactionsDataFrame,
            *args,
            **kwargs
        )
        self.A = self.reacDB['A']
        self.beta = self.reacDB['beta']
        self.Ta = self.reacDB['Ta']

    # Master Equation matrices
    ###########################################################################
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
