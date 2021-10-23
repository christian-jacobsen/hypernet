import abc


class BasicReactionRate(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        reactionsDataFrame,
        *args,
        **kwargs
    ):
        self.reacDB = reactionsDataFrame

    # Update method -----------------------------------------------------------
    def update(self, T):
        self.k = self.k_(T)
        self.dkdT = self.dkdT_(T)

    # Forward reaction rates --------------------------------------------------
    @abc.abstractmethod
    def k_(self, T=0.):
        pass

    @abc.abstractmethod
    def dkdT_(self, T=0.):
        pass

    # Properties --------------------------------------------------------------
    @property
    def k(self):
        return self._k
    @k.setter
    def k(self, value):
        self._k = value

    @property
    def dkdT(self):
        return self._dkdT
    @dkdT.setter
    def dkdT(self, value):
        self._dkdT = value
