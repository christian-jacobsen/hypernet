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
        return self.A * np.power(T, self.beta) * np.exp(-self.Ta / T)

    def dkdT_(self, T):
        return (self.beta + self.Ta / T) * self.k  / T
