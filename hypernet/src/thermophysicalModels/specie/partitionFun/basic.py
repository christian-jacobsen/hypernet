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
        self.dQdT = 0.

    # Properties
    ###########################################################################
    @property
    def Q(self):
        return self._Q
    @Q.setter
    def Q(self, value):
        self._Q = value

    @property
    def dQdT(self):
        return self._dQdT
    @dQdT.setter
    def dQdT(self, value):
        self._dQdT = value

    # Methods
    ###########################################################################
    # Update ------------------------------------------------------------------
    def update(self, T):
        self.Q = self.Q_(T)
        self.dQdT = self.dQdT_(T)

    # Partition functions -----------------------------------------------------
    @abc.abstractmethod
    def Q_(self, T):
        pass

    @abc.abstractmethod
    def dQdT_(self, T):
        pass
