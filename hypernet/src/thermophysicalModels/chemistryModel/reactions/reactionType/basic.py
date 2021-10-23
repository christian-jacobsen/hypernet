import abc


class BasicReactionType(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specieThermos,
        reactionRate,
        processIndeces,
        *args,
        **kwargs
    ):
        self.spTh = specieThermos
        self.reacRate = reactionRate
        self.processIndeces = processIndeces

    # Update method -----------------------------------------------------------
    def update(self, T):
        for name, spTh_ in self.spTh.items():
            spTh_.intPF.update(T)
            spTh_.transPF.update(T)
        self.reacRate.update(T)
        self.kf = self.reacRate.k
        self.dkfdT = self.reacRate.dkdT

    # Properties --------------------------------------------------------------
    @property
    def kf(self):
        return self._kf
    @kf.setter
    def kf(self, value):
        self._kf = value

    @property
    def dkfdT(self):
        return self._dkfdT
    @dkfdT.setter
    def dkfdT(self, value):
        self._dkfdT = value

    # Forward reaction rates --------------------------------------------------
    def kf_(self, T):
        return self.reactionRate.k(T)

    def dkfdT_(self, T):
        return self.reactionRate.dkdT(T)

    # Reverse reaction rates --------------------------------------------------
    @abc.abstractmethod
    def kr(self, kf, *args, **kwargs):
        pass

    @abc.abstractmethod
    def dkrdT(self, kr, dkfdT, *args, **kwargs):
        pass
