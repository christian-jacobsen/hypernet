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
        # Thermodynamic specie properties
        self.spTh = specieThermos
        for name, spTh in self.spTh.items():
            if spTh.specie.n_at > 1:
                self.molecule = name
            else:
                self.atom = name

        # Reactive processes
        self.processIndices = processIndices

        # Reaction rate
        self.reacRate = reactionRate

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
