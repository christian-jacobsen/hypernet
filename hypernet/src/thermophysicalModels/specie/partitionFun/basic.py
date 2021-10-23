import abc

from hypernet.src.general import const


class Basic(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        # Specie Properties
        self.specie = specie

        self.Q = 1.
        self.dQ_dT = 0.

    # Properties
    ###########################################################################
    @property
    def Q(self):
        return self._Q
    @Q.setter
    def Q(self, value):
        self._Q = value

    @property
    def dQ_dT(self):
        return self._dQ_dT
    @dQ_dT.setter
    def dQ_dT(self, value):
        self._dQ_dT = value

    # Methods
    ###########################################################################
    # Update ------------------------------------------------------------------
    def update(self, T):
        self.Q = self.Q_(T)
        self.dQ_dT = self.dQ_dT_(T)

    # Partition functions -----------------------------------------------------
    @abc.abstractmethod
    def Q_(self, T):
        pass

    @abc.abstractmethod
    def dQ_dT_(self, T):
        pass