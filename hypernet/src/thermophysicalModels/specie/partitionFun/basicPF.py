import abc


class PartitionFun(object):

    # Initialization
    ###########################################################################
    def __init__(
        self,
        specie,
        *args,
        **kwargs
    ):
        # Specie Properties ===================================================
        self.specie = specie

    # Partition functions -----------------------------------------------------
    @abc.abstractmethod
    def Q_(self, T):
        pass

    @abc.abstractmethod
    def dQ_dT_(self, T):
        pass

    # Properties --------------------------------------------------------------
    @property
    def Q(self):
        return self._Q
    @Q.setter
    def Q(self, T):
        self._Q = self.Q_(T)

    @property
    def dQ_dT(self):
        return self._dQ_dT
    @dQ_dT.setter
    def dQ_dT(self, T):
        self._dQ_dT = self.dQ_dT_(T)