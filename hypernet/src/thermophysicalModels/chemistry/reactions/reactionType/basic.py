import abc


class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specieThermos,
        reactionRate,
        processIndices,
        *args,
        **kwargs
    ):
        self.spTh = specieThermos
        self.reacRate = reactionRate
        self.processIndices = processIndices

    # Properties
    ###########################################################################
    # Forward reaction rates --------------------------------------------------
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

    # Methods
    ###########################################################################
    # Update method -----------------------------------------------------------
    def update(self, T):
        # Update partition functions
        for name, spTh_ in self.spTh.items():
            spTh_.intPF.update(T)
            spTh_.transPF.update(T)
        # Update reactions
        self.reacRate.update(T)
        # Get forward reactions rates
        self.kf = self.reacRate.k
        self.dkfdT = self.reacRate.dkdT

    # Reverse reaction rates --------------------------------------------------
    @abc.abstractmethod
    def kr_(self, kf, *args, **kwargs):
        pass

    @abc.abstractmethod
    def dkrdT_(self, kr, dkfdT, *args, **kwargs):
        pass
